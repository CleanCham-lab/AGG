#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AGG - ANNOVAR Gene Gleaner
==========================
Gene-Centric Variant Extraction from ANNOVAR

A GATK -> ANNOVAR (table_annovar.pl) run leaves one row per variant for the
whole exome/genome. AGG keeps only the rows whose gene column
(default ``Gene.refGene``) contains a requested gene symbol and writes them to
``multianno_<GENE>.csv`` / ``.xlsx`` (one file per gene), or to a single merged
file for a gene panel.

Run ``python agg.py -h`` for all options.

Seung Won Chae, Moon-woo Seong
Molecular Diagnostics Laboratory, Department of Laboratory Medicine,
Seoul National University Hospital
https://github.com/CleanCham-lab/AGG  (MIT License)
"""

from __future__ import annotations

import argparse
import csv
import glob
import gzip
import importlib
import importlib.util
import io
import logging
import os
import re
import site
import subprocess
import sys
import time
from collections import Counter, OrderedDict
from pathlib import Path

if sys.version_info < (3, 8):
    sys.exit("[ERROR] Python >= 3.8 is required.")

__version__ = "1.0.0"


# --------------------------------------------------------------------------- #
# dependencies
#   Missing packages (pandas, numpy, openpyxl) are installed automatically into
#   the Python that runs this script ("<python> -m pip install ...").
#   Turn this off with -no-install or the environment variable AGG_NO_INSTALL=1.
# --------------------------------------------------------------------------- #
PIP_SPECS = {
    "pandas": "pandas>=1.4",
    "numpy": "numpy",
    "openpyxl": "openpyxl>=3.0",
    "dateutil": "python-dateutil",
    "pytz": "pytz",
    "tzdata": "tzdata",
}
AUTO_INSTALL = (
    not any(arg in ("-no-install", "--no-install") for arg in sys.argv[1:])
    and os.environ.get("AGG_NO_INSTALL", "").lower() not in ("1", "true", "yes")
)


def _setup_msg(msg: str) -> None:
    print(f"[setup] {msg}", file=sys.stderr, flush=True)


def _in_virtual_env() -> bool:
    return sys.prefix != getattr(sys, "base_prefix", sys.prefix) or bool(os.environ.get("CONDA_PREFIX"))


def _manual_command(spec: str) -> str:
    python = f'"{sys.executable}"' if " " in sys.executable else sys.executable
    return f'{python} -m pip install "{spec}"'


def _pip_install(spec: str) -> bool:
    """pip-install `spec` into the running interpreter; try --user outside venv/conda."""
    base = [sys.executable, "-m", "pip", "install", "--disable-pip-version-check", spec]
    attempts = [base] if _in_virtual_env() else [base, base + ["--user"]]
    output = ""
    for cmd in attempts:
        _setup_msg("running: " + " ".join(f'"{c}"' if " " in c else c for c in cmd))
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
        except OSError as exc:
            output = str(exc)
            break
        if result.returncode == 0:
            if "--user" in cmd:
                site.addsitedir(site.getusersitepackages())
            return True
        output = (result.stdout or "") + (result.stderr or "")

    tail = "\n".join(output.strip().splitlines()[-10:])
    if tail:
        print(tail, file=sys.stderr)
    low = output.lower()
    if "no module named pip" in low:
        _setup_msg("pip is not available for this Python. Use a conda environment, or install pip "
                   "(Ubuntu/WSL: sudo apt install python3-pip).")
    elif "certificate_verify_failed" in low or "sslerror" in low or "ssl:" in low:
        _setup_msg("SSL certificate check failed - common behind hospital/company SSL-inspection networks. "
                   "Ask IT for the root CA and run 'python -m pip config set global.cert <ca.pem>', "
                   "or, if your institution allows it, install once with "
                   "'python -m pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org "
                   f"\"{spec}\"'.")
    elif "externally-managed-environment" in low:
        _setup_msg("This Python is managed by the operating system (PEP 668). Run the tool from a "
                   "conda environment or a venv:  python3 -m venv .venv && . .venv/bin/activate")
    elif any(k in low for k in ("temporary failure", "name resolution", "connection", "proxyerror",
                                "no matching distribution")):
        _setup_msg("Could not reach the package index. Check the internet connection / proxy settings.")
    return False


def ensure_package(module: str):
    """Import `module`; if it (or one of its dependencies) is missing, pip-install it and retry."""
    for _ in range(4):
        try:
            return importlib.import_module(module)
        except ModuleNotFoundError as exc:
            missing = (exc.name or module).split(".")[0]
        except ImportError as exc:  # installed but broken (e.g. numpy/pandas binary mismatch)
            sys.exit(f"[ERROR] '{module}' is installed but failed to import: {exc}\n"
                     f"  Repair it with:  {_manual_command(PIP_SPECS.get(module, module))} --upgrade")
        spec = PIP_SPECS.get(missing, missing)
        if not AUTO_INSTALL:
            sys.exit(f"[ERROR] Python package '{missing}' is missing. Install it with:\n  {_manual_command(spec)}")
        _setup_msg(f"'{missing}' is not installed for {sys.executable} - installing {spec} "
                   "(one-time, may take a minute) ...")
        if not _pip_install(spec):
            sys.exit(f"[ERROR] Could not install '{missing}' automatically. Install it manually and re-run:\n"
                     f"  {_manual_command(spec)}")
        importlib.invalidate_caches()
        _setup_msg(f"installed {spec}")
    sys.exit(f"[ERROR] Could not import '{module}' after installing its dependencies.")


pd = ensure_package("pandas")
np = ensure_package("numpy")

LOG = logging.getLogger("agg")

# Multiple genes in one ANNOVAR cell are joined by ';' (current versions) or ','
# (older versions); values that went through a VCF may carry the escaped '\x3b'.
GENE_SPLIT_RE = r"\\x3b|[;,]"
FUNC_SPLIT_RE = r"\\x3b|;"

DEFAULT_PATTERNS = ["*multianno.txt", "*multianno.txt.gz", "*multianno.csv", "*multianno.csv.gz"]
LIST_HEADER_WORDS = {
    "gene", "genes", "symbol", "gene symbol", "gene_symbol", "genesymbol",
    "gene name", "gene_name", "genename", "hgnc", "hgnc symbol", "hgnc_symbol",
}
XLSX_MAX_ROWS = 1_048_576
XLSX_MAX_CELL_CHARS = 32_767
PREFILTER_MAX_GENES = 300  # above this, skip the regex pre-filter and tokenise every row

# Allow very long INFO / Otherinfo fields (sys.maxsize overflows a C long on Windows).
csv.field_size_limit(2**31 - 1)

EXAMPLES = r"""
examples:
  # one gene, one sample
  python agg.py -i S01.hg38_multianno.txt -g RNU4-2 -o results

  # several genes typed on the command line, Excel output
  python agg.py -i S01.hg38_multianno.txt -g DMD BRCA1,BRCA2 -o results -f xlsx

  # gene panel from an Excel list (A = No., B = Gene), every multianno file in a folder
  python agg.py -i annovar_out/ -l panel.xlsx -o results -f xlsx

  # the same panel written as one merged file: results/multianno_panel.xlsx
  python agg.py -i annovar_out/ -l panel.xlsx -o results -f xlsx -merge

  # keep only exonic / splicing / ncRNA_exonic variants
  python agg.py -i S01.hg38_multianno.txt -l panel.xlsx -o results -func exonic splicing ncRNA_exonic

  # Windows paths work as well
  python agg.py -i "Y:\project\Annotated_VCF" -g RNU4-2 -o "C:\Users\me\Desktop\out"
"""


# --------------------------------------------------------------------------- #
# logging
# --------------------------------------------------------------------------- #
class _Formatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        msg = record.getMessage()
        return msg if record.levelno == logging.INFO else f"[{record.levelname}] {msg}"


def _setup_logging(quiet: bool) -> None:
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(_Formatter())
    LOG.handlers[:] = [handler]
    LOG.setLevel(logging.WARNING if quiet else logging.INFO)
    LOG.propagate = False


def _die(msg: str) -> None:
    LOG.error(msg)
    sys.exit(1)


def _check_versions() -> None:
    # An old pandas is reported, not upgraded: replacing it could break other tools in the same environment.
    nums = re.findall(r"\d+", pd.__version__)[:2]
    if len(nums) == 2 and (int(nums[0]), int(nums[1])) < (1, 4):
        _die(f"pandas >= 1.4 is required (found {pd.__version__}). Upgrade with:  "
             f"{_manual_command('pandas>=1.4')} --upgrade")


def _require_openpyxl(why: str) -> None:
    """openpyxl is only needed for .xlsx input/output, so it is installed on first use."""
    if importlib.util.find_spec("openpyxl") is None:
        LOG.info("openpyxl is needed for %s.", why)
    ensure_package("openpyxl")


# --------------------------------------------------------------------------- #
# gene list
# --------------------------------------------------------------------------- #
def column_index(spec: str) -> int:
    """'B' -> 1, '2' -> 1 (0-based)."""
    spec = str(spec).strip()
    if spec.isdigit() and int(spec) > 0:
        return int(spec) - 1
    if spec.isalpha():
        idx = 0
        for ch in spec.upper():
            idx = idx * 26 + (ord(ch) - ord("A") + 1)
        return idx - 1
    raise ValueError(f"invalid column '{spec}' (use a letter such as B or a 1-based number such as 2)")


def _read_text_table(path: Path) -> pd.DataFrame:
    text = None
    for enc in ("utf-8-sig", "cp949", "latin-1"):
        try:
            text = path.read_text(encoding=enc)
            break
        except UnicodeDecodeError:
            continue
    rows = []
    for line in (text or "").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if "\t" in line:
            cells = line.split("\t")
        elif "," in line:
            cells = next(csv.reader([line]))
        else:
            cells = line.split()
        rows.append([c.strip() for c in cells])
    return pd.DataFrame(rows)


def read_gene_list(path: Path, col_spec: str) -> list[str]:
    """Read gene symbols from an .xlsx/.csv/.tsv/.txt list (default: column B)."""
    if not path.is_file():
        _die(f"Gene list not found: {path}")
    if path.suffix.lower() in (".xlsx", ".xlsm", ".xls"):
        _require_openpyxl("reading .xlsx gene lists")
        df = pd.read_excel(path, header=None, dtype=str, sheet_name=0)
    else:
        df = _read_text_table(path)
    df = df.replace("", np.nan).dropna(how="all")
    if df.empty:
        _die(f"Gene list is empty: {path}")

    try:
        idx = column_index(col_spec)
    except ValueError as exc:
        _die(f"-list-col: {exc}")
    if idx >= df.shape[1]:
        if df.shape[1] == 1:
            LOG.info("Gene list has a single column; reading gene symbols from column A.")
            idx = 0
        else:
            _die(f"Gene list has {df.shape[1]} column(s); column '{col_spec}' does not exist.")

    genes, numeric, date_like = [], [], []
    for value in df.iloc[:, idx].tolist():
        if value is None or (isinstance(value, float) and np.isnan(value)):
            continue
        gene = str(value).strip().strip('"').strip("'")
        if not gene or gene.lower() in LIST_HEADER_WORDS or gene.lower() in ("nan", "none"):
            continue
        if re.fullmatch(r"\d+(\.0+)?", gene):
            numeric.append(gene)
            continue
        if re.match(r"\d{4}-\d{2}-\d{2}", gene) or re.fullmatch(r"\d{1,2}-[A-Za-z]{3}", gene):
            date_like.append(gene)
            continue
        genes.append(gene)

    if numeric and not genes:
        _die(f"Column '{col_spec}' of {path.name} holds only numbers. "
             "Point -list-col at the gene-symbol column.")
    if numeric:
        LOG.warning("Ignored %d numeric value(s) in the gene column of %s.", len(numeric), path.name)
    if date_like:
        LOG.warning("Skipped gene-list entries that look like dates (Excel may have converted symbols "
                    "such as MARCH1/SEPT9 - fix them in the list): %s", ", ".join(date_like))
    return genes


def build_targets(genes: list[str], case_sensitive: bool) -> "OrderedDict[str, str]":
    """Map normalised key -> symbol as the user wrote it (first spelling wins)."""
    targets: "OrderedDict[str, str]" = OrderedDict()
    for gene in genes:
        targets.setdefault(gene if case_sensitive else gene.upper(), gene)
    return targets


# --------------------------------------------------------------------------- #
# input files
# --------------------------------------------------------------------------- #
def collect_inputs(items: list[str], patterns: list[str]) -> list[Path]:
    files: list[Path] = []
    for item in items:
        path = Path(item)
        if path.is_dir():
            found = sorted({p for pat in patterns for p in path.glob(pat) if p.is_file()})
            if not found:
                LOG.warning("No file matching %s in directory: %s", " / ".join(patterns), path)
            files.extend(found)
        elif any(ch in item for ch in "*?["):
            found = sorted(Path(p) for p in glob.glob(item) if Path(p).is_file())
            if not found:
                LOG.warning("No file matches pattern: %s", item)
            files.extend(found)
        elif path.is_file():
            files.append(path)
        else:
            _die(f"Input not found: {item}")
    unique, seen = [], set()
    for f in files:
        key = f.resolve()
        if key not in seen:
            seen.add(key)
            unique.append(f)
    return unique


def sample_name(path: Path) -> str:
    """S01.hg38_multianno.txt(.gz) -> S01"""
    name = re.sub(r"\.gz$", "", path.name, flags=re.I)
    name = re.sub(r"\.(txt|tsv|csv)$", "", name, flags=re.I)
    name = re.sub(r"[._]?(?:hg\d+|GRCh\d+|mm\d+)?_multianno$", "", name, flags=re.I)
    return name or path.name


def _delimiter(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return "," if name.endswith(".csv") else "\t"


def _open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace", newline="")
    return open(path, "r", encoding="utf-8", errors="replace", newline="")


def _split_line(line: str, sep: str) -> list[str]:
    line = line.rstrip("\r\n")
    return line.split("\t") if sep == "\t" else next(csv.reader([line]))


def read_column_names(path: Path, sep: str, max_lines: int | None = 5000) -> list[str] | None:
    """Header of a multianno table, padded for unnamed trailing (Otherinfo) columns.

    Older ANNOVAR versions write a single 'Otherinfo' header but several
    tab-separated values behind it when run with -vcfinput.
    """
    with _open_text(path) as fh:
        first = fh.readline()
        if not first.strip():
            return None
        header = [h.strip() for h in _split_line(first, sep)]
        n_fields = len(header)
        for i, line in enumerate(fh):
            if max_lines is not None and i >= max_lines:
                break
            if line.strip():
                n_fields = max(n_fields, len(_split_line(line, sep)))
    return _pad_header(header, n_fields)


def _pad_header(header: list[str], n_fields: int) -> list[str]:
    names = list(header)
    extra = n_fields - len(names)
    if extra > 0:
        match = re.fullmatch(r"Otherinfo(\d*)", names[-1])
        if match:
            if match.group(1):
                start = int(match.group(1)) + 1
            else:
                names[-1], start = "Otherinfo1", 2
            names += [f"Otherinfo{k}" for k in range(start, start + extra)]
        else:
            names += [f"Extra{k}" for k in range(1, extra + 1)]
    seen: Counter = Counter()
    unique = []
    for name in names:
        name = name or "Unnamed"
        unique.append(f"{name}.{seen[name]}" if seen[name] else name)
        seen[name] += 1
    return unique


def func_column_for(gene_col: str) -> str | None:
    """Gene.refGene -> Func.refGene"""
    return "Func." + gene_col[len("Gene."):] if gene_col.startswith("Gene.") else None


# --------------------------------------------------------------------------- #
# scanning
# --------------------------------------------------------------------------- #
class Scanner:
    def __init__(self, targets, gene_col, func_filter, keep_intergenic, case_sensitive, chunksize):
        self.targets = targets
        self.keys = set(targets)
        self.gene_col = gene_col
        self.func_set = {f.lower() for f in func_filter} if func_filter else None
        self.keep_intergenic = keep_intergenic
        self.case_sensitive = case_sensitive
        self.chunksize = chunksize
        self.prefilter = None
        if len(targets) <= PREFILTER_MAX_GENES:
            alts = sorted({re.escape(g) for g in targets.values()}, key=len, reverse=True)
            self.prefilter = "|".join(alts)
        self._warned_func = False

    # -- public ------------------------------------------------------------ #
    def scan(self, path: Path):
        """Return (matched rows with a '_keys' column | None, rows scanned, removed Counter, removed rows)."""
        sep = _delimiter(path)
        names = read_column_names(path, sep)
        if names is None:
            LOG.warning("%s: empty file, skipped.", path.name)
            return None, 0, Counter(), 0
        if self.gene_col not in names:
            gene_cols = [n for n in names if n.startswith("Gene.")]
            _die(f"{path.name}: column '{self.gene_col}' not found. "
                 f"Gene columns in this file: {', '.join(gene_cols) or 'none'} (choose one with -genecol).")
        try:
            return self._scan(path, sep, names)
        except pd.errors.ParserError:
            LOG.info("  %s: rows wider than the header were found further down; re-reading the full file.",
                     path.name)
            return self._scan(path, sep, read_column_names(path, sep, max_lines=None))

    # -- internals --------------------------------------------------------- #
    def _func_col(self, names):
        col = func_column_for(self.gene_col)
        if col in names:
            return col
        if (self.func_set or not self.keep_intergenic) and not self._warned_func:
            LOG.warning("Column '%s' not found; Func filtering (-func / intergenic exclusion) is skipped.", col)
            self._warned_func = True
        return None

    def _func_mask(self, sub: pd.DataFrame, func_col):
        """Boolean Series of rows to keep, or None when nothing is filtered."""
        if func_col is None:
            return None
        func = sub[func_col].fillna("").str.strip().str.lower()
        if self.func_set:
            tokens = func.str.split(FUNC_SPLIT_RE, regex=True).explode().str.strip()
            keep = tokens.isin(self.func_set).groupby(level=0).any()
            return keep.reindex(sub.index, fill_value=False)
        if not self.keep_intergenic:
            return func != "intergenic"
        return None

    def _scan(self, path, sep, names):
        func_col = self._func_col(names)
        # No index_col=False on purpose: a row wider than `names` must raise a
        # ParserError (-> full re-scan of the header width) instead of being truncated.
        kwargs = dict(sep=sep, header=None, skiprows=1, names=names,
                      dtype=str, na_filter=False, chunksize=self.chunksize,
                      encoding="utf-8", encoding_errors="replace", compression="infer")
        if sep == "\t":
            kwargs["quoting"] = csv.QUOTE_NONE

        frames, removed, removed_rows, n_rows = [], Counter(), 0, 0
        with pd.read_csv(path, **kwargs) as reader:
            for chunk in reader:
                n_rows += len(chunk)
                genes = chunk[self.gene_col].fillna("")
                if self.prefilter:
                    candidate = genes.str.contains(self.prefilter, case=self.case_sensitive,
                                                   regex=True, na=False)
                    if not candidate.any():
                        continue
                    chunk, genes = chunk.loc[candidate], genes.loc[candidate]

                tokens = genes.str.split(GENE_SPLIT_RE, regex=True).explode().str.strip()
                if not self.case_sensitive:
                    tokens = tokens.str.upper()
                hits = tokens[tokens.isin(self.keys)]
                if hits.empty:
                    continue

                row_keys: "OrderedDict[int, list]" = OrderedDict()
                for idx, key in zip(hits.index, hits.to_numpy()):
                    keys = row_keys.setdefault(idx, [])
                    if key not in keys:
                        keys.append(key)
                index = list(row_keys)
                sub = chunk.loc[index].copy()
                keys_s = pd.Series([tuple(row_keys[i]) for i in index], index=index, dtype=object)

                keep = self._func_mask(sub, func_col)
                if keep is not None and not keep.all():
                    removed_rows += int((~keep).sum())
                    for keys in keys_s[~keep]:
                        removed.update(keys)
                    sub, keys_s = sub.loc[keep], keys_s.loc[keep]
                if sub.empty:
                    continue
                sub.insert(0, "_keys", keys_s)
                frames.append(sub)

        matched = pd.concat(frames) if frames else None
        return matched, n_rows, removed, removed_rows


# --------------------------------------------------------------------------- #
# output
# --------------------------------------------------------------------------- #
def safe_filename(text: str) -> str:
    name = re.sub(r'[\\/:*?"<>|\s]+', "_", str(text)).strip("._")
    return name or "genes"


def output_path(outdir: Path, label: str, fmt: str, prefix: str | None = None, suffix: str = "") -> Path:
    stem = f"multianno_{safe_filename(label)}{suffix}"
    if prefix:
        stem = f"{safe_filename(prefix)}_{stem}"
    return outdir / f"{stem}.{fmt}"


def _excel_ready(df: pd.DataFrame, numeric: bool) -> pd.DataFrame:
    """Strip characters Excel rejects, cap cell length, optionally turn numeric text into numbers."""
    from openpyxl.cell.cell import ILLEGAL_CHARACTERS_RE

    out = df.copy()
    for col in out.columns:
        series = out[col]
        # text columns are 'object' (pandas < 3) or 'str' (pandas >= 3)
        if pd.api.types.is_numeric_dtype(series) or pd.api.types.is_bool_dtype(series):
            continue
        present = series.notna()
        if not present.any():
            continue
        text = series[present].astype(str).str.replace(ILLEGAL_CHARACTERS_RE, "", regex=True)
        too_long = text.str.len() > XLSX_MAX_CELL_CHARS
        if too_long.any():
            LOG.warning("Column '%s': %d cell(s) over %d characters were truncated in the xlsx output.",
                        col, int(too_long.sum()), XLSX_MAX_CELL_CHARS)
            text = text.where(~too_long, text.str.slice(0, XLSX_MAX_CELL_CHARS))
        values = text.astype(object)
        if numeric:
            blank = text.isin([".", ""])
            filled = text[~blank]
            if len(filled):
                num = pd.to_numeric(filled, errors="coerce")
                is_number = (
                    num.notna().all()
                    and np.isfinite(num.astype(float)).all()
                    and not filled.str.match(r"^-?0\d").any()                               # keep leading zeros
                    and filled.str.replace(r"[^0-9]", "", regex=True).str.len().max() <= 15  # float precision
                )
                if is_number:
                    values.loc[~blank] = num.to_numpy(dtype=object)
        new = series.astype(object)
        new.loc[present] = values
        out[col] = new
    return out


def write_table(df: pd.DataFrame, path: Path, fmt: str, sheet: str = "variants",
                xlsx_numeric: bool = True) -> Path:
    if fmt == "csv":
        df.to_csv(path, index=False)
        return path
    if len(df) + 1 > XLSX_MAX_ROWS:
        alt = path.with_suffix(".csv")
        LOG.warning("%s: %s rows exceed the Excel row limit; wrote CSV instead -> %s",
                    path.name, f"{len(df):,}", alt.name)
        df.to_csv(alt, index=False)
        return alt
    from openpyxl.utils import get_column_letter

    ready = _excel_ready(df, numeric=xlsx_numeric)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        ready.to_excel(writer, index=False, sheet_name=sheet)
        ws = writer.sheets[sheet]
        ws.freeze_panes = "A2"
        if len(ready):
            ws.auto_filter.ref = ws.dimensions
        for i, col in enumerate(ready.columns, start=1):
            ws.column_dimensions[get_column_letter(i)].width = min(max(len(str(col)) + 2, 10), 40)
    return path


def run_label(args, genes: list[str]) -> str:
    if args.name:
        return args.name
    if args.gene_list:
        return Path(args.gene_list).stem
    return "_".join(genes) if len(genes) <= 3 else f"{len(genes)}genes"


def build_summary(targets, groups, combined, removed, func_col, multi, written) -> pd.DataFrame:
    rows = []
    for key, symbol in targets.items():
        idx = groups.get(key, [])
        row = OrderedDict(Gene=symbol, Variants=len(idx))
        if multi:
            row["Samples"] = int(combined.loc[idx, "Sample"].nunique()) if idx else 0
        if idx and func_col and func_col in combined.columns:
            counts = combined.loc[idx, func_col].value_counts()
            row["Func_breakdown"] = "; ".join(f"{k}:{v}" for k, v in counts.items())
        else:
            row["Func_breakdown"] = ""
        row["Removed_by_Func_filter"] = removed.get(key, 0)
        row["Output_file"] = written.get(key, "")
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# command line
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="agg.py",
        description=f"AGG v{__version__} - ANNOVAR Gene Gleaner\n"
                    "Gene-Centric Variant Extraction from ANNOVAR\n\n"
                    "Extract the variants of selected genes (one gene or a gene panel) from ANNOVAR\n"
                    "*_multianno tables and save them as multianno_<GENE>.csv / .xlsx.",
        epilog=EXAMPLES,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
    )
    io_grp = parser.add_argument_group("input / output")
    io_grp.add_argument("-i", "-input", "--input", nargs="+", required=True, metavar="PATH",
                        help="multianno file(s) (.txt/.csv, .gz allowed), a directory, or a glob pattern")
    io_grp.add_argument("-o", "-out", "--outdir", default=".", metavar="DIR",
                        help="output directory, created if missing (default: current directory)")
    io_grp.add_argument("-f", "-format", "--format", choices=["csv", "xlsx"], default="csv",
                        type=str.lower, help="output format (default: csv)")

    gene_grp = parser.add_argument_group("genes (give -g, -l, or both)")
    gene_grp.add_argument("-g", "-gene", "--gene", nargs="+", metavar="GENE",
                          help="gene symbol(s), space- or comma-separated: -g RNU4-2  /  -g DMD BRCA1,BRCA2")
    gene_grp.add_argument("-l", "-list", "--gene-list", metavar="FILE",
                          help="gene list (.xlsx/.csv/.tsv/.txt); layout: column A = No., column B = Gene")
    gene_grp.add_argument("-list-col", "--list-col", default="B", metavar="COL",
                          help="list column holding the symbols, letter or 1-based number "
                               "(default: B; single-column lists fall back to A)")

    out_grp = parser.add_argument_group("output naming")
    out_grp.add_argument("-merge", "--merge", action="store_true",
                         help="write all genes to one file multianno_<NAME>.<fmt> with a Matched_Gene column")
    out_grp.add_argument("-name", "--name", metavar="NAME",
                         help="label for the merged/summary file (default: gene-list file name)")
    out_grp.add_argument("-prefix", "--prefix", metavar="TEXT",
                         help="prefix for every output file, e.g. a sample ID -> TEXT_multianno_<GENE>.csv")

    filt_grp = parser.add_argument_group("matching / filtering")
    filt_grp.add_argument("-genecol", "--gene-col", default="Gene.refGene", metavar="COL",
                          help="gene column to search (default: Gene.refGene; e.g. Gene.ensGene)")
    filt_grp.add_argument("-func", "--func", nargs="+", metavar="FUNC",
                          help="keep only these Func.* categories, e.g. -func exonic splicing ncRNA_exonic")
    filt_grp.add_argument("-keep-intergenic", "--keep-intergenic", action="store_true",
                          help="also keep intergenic rows where the gene is only a flanking gene")
    filt_grp.add_argument("-case-sensitive", "--case-sensitive", action="store_true",
                          help="match gene symbols case-sensitively (default: case-insensitive)")

    misc_grp = parser.add_argument_group("other")
    misc_grp.add_argument("-pattern", "--pattern", nargs="+", default=DEFAULT_PATTERNS, metavar="GLOB",
                          help="file patterns used when -i is a directory "
                               "(default: *multianno.txt *multianno.txt.gz *multianno.csv *multianno.csv.gz)")
    misc_grp.add_argument("-chunksize", "--chunksize", type=int, default=200_000, metavar="N",
                          help="rows read per chunk (default: 200000)")
    misc_grp.add_argument("-xlsx-text", "--xlsx-text", action="store_true",
                          help="in xlsx output keep every cell as text (default: numeric columns become numbers)")
    misc_grp.add_argument("-no-install", "--no-install", action="store_true",
                          help="do not pip-install missing packages automatically (also: AGG_NO_INSTALL=1)")
    misc_grp.add_argument("-q", "--quiet", action="store_true", help="only print warnings and errors")
    misc_grp.add_argument("-v", "--version", action="version", version=f"AGG {__version__} - ANNOVAR Gene Gleaner")

    args = parser.parse_args(argv)
    if not args.gene and not args.gene_list:
        parser.error("give at least one gene with -g or a gene list with -l")
    if args.chunksize <= 0:
        parser.error("-chunksize must be a positive integer")
    return args


def main(argv=None) -> int:
    args = parse_args(argv)
    _setup_logging(args.quiet)
    LOG.info("AGG v%s - ANNOVAR Gene Gleaner", __version__)
    _check_versions()
    if args.format == "xlsx":
        _require_openpyxl("xlsx output")

    # genes ------------------------------------------------------------------
    genes: list[str] = []
    for item in args.gene or []:
        genes.extend(g.strip() for g in item.split(",") if g.strip())
    if args.gene_list:
        listed = read_gene_list(Path(args.gene_list), args.list_col)
        LOG.info("Gene list %s: %d symbol(s).", Path(args.gene_list).name, len(listed))
        genes.extend(listed)
    targets = build_targets(genes, args.case_sensitive)
    if not targets:
        _die("No gene symbols to search for.")
    preview = ", ".join(list(targets.values())[:10]) + (" ..." if len(targets) > 10 else "")
    LOG.info("Searching %d gene(s) in column %s: %s", len(targets), args.gene_col, preview)

    # inputs -----------------------------------------------------------------
    inputs = collect_inputs(args.input, args.pattern)
    if not inputs:
        _die("No multianno input files found.")
    multi = len(inputs) > 1
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    scanner = Scanner(targets, args.gene_col, args.func, args.keep_intergenic,
                      args.case_sensitive, args.chunksize)
    frames, removed, removed_rows = [], Counter(), 0
    for n, path in enumerate(inputs, start=1):
        t0 = time.time()
        matched, n_rows, rem, rem_rows = scanner.scan(path)
        n_hit = 0 if matched is None else len(matched)
        LOG.info("[%d/%d] %s: %s rows scanned, %s matched (%.1fs)",
                 n, len(inputs), path.name, f"{n_rows:,}", f"{n_hit:,}", time.time() - t0)
        removed.update(rem)
        removed_rows += rem_rows
        if n_hit:
            if multi:
                matched.insert(0, "Sample", sample_name(path))
            frames.append(matched)

    combined = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()
    groups: dict = {}
    if not combined.empty:
        exploded = combined["_keys"].explode()
        for key, idx in exploded.groupby(exploded, sort=False).groups.items():
            groups[key] = sorted(set(idx))

    # write ------------------------------------------------------------------
    label = run_label(args, list(targets.values()))
    written: dict = {}
    xlsx_numeric = not args.xlsx_text
    if args.merge:
        if not combined.empty:
            merged = combined.drop(columns="_keys")
            symbols = combined["_keys"].map(lambda keys: ";".join(targets[k] for k in keys))
            merged.insert(1 if multi else 0, "Matched_Gene", symbols)
            path = write_table(merged, output_path(outdir, label, args.format, args.prefix),
                               args.format, xlsx_numeric=xlsx_numeric)
            for key in groups:
                written[key] = path.name
            LOG.info("Wrote %s (%s rows)", path, f"{len(merged):,}")
    else:
        for key, symbol in targets.items():
            idx = groups.get(key)
            if not idx:
                continue
            table = combined.loc[idx].drop(columns="_keys")
            path = write_table(table, output_path(outdir, symbol, args.format, args.prefix),
                               args.format, xlsx_numeric=xlsx_numeric)
            written[key] = path.name
            LOG.info("Wrote %s (%s rows)", path, f"{len(table):,}")

    func_col = func_column_for(args.gene_col)
    summary = build_summary(targets, groups, combined, removed, func_col, multi, written)
    if len(targets) > 1:
        path = write_table(summary, output_path(outdir, label, args.format, args.prefix, "_summary"),
                           args.format, sheet="summary", xlsx_numeric=xlsx_numeric)
        LOG.info("Wrote %s", path)

    # report -----------------------------------------------------------------
    found = int((summary["Variants"] > 0).sum())
    LOG.info("")
    if len(targets) <= 30:
        LOG.info(summary.drop(columns=["Output_file"]).to_string(index=False))
        LOG.info("")
    LOG.info("Done: %d/%d gene(s) with variants, %d input file(s), output in %s",
             found, len(targets), len(inputs), outdir.resolve())
    if removed_rows and not args.func and not args.keep_intergenic:
        LOG.info("Excluded %s intergenic row(s) that list a target only as a flanking gene "
                 "(add -keep-intergenic to include them).", f"{removed_rows:,}")
    missing = [targets[k] for k in targets if k not in groups]
    if missing:
        LOG.warning("No variants for: %s (check the symbol/alias, -genecol, or the filters).",
                    ", ".join(missing))
    return 0


if __name__ == "__main__":
    sys.exit(main())

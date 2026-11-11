# gene_search

Pull the variants of the genes you care about out of ANNOVAR `*_multianno` tables.

After a **GATK → ANNOVAR** (`table_annovar.pl`) run, a single exome or genome yields tens of thousands to millions of annotated rows. `gene_search.py` keeps only the rows whose gene column (default `Gene.refGene`) contains the requested gene symbol(s) and saves them as

```
multianno_<GENE>.csv   or   multianno_<GENE>.xlsx
```

one file per gene, or one merged file for a whole gene panel.

---

## Install

Python ≥ 3.8. Just download and run:

```bash
git clone https://github.com/CleanCham-lab/gene_search.git
cd gene_search
python gene_search.py -h
```

**Missing packages are installed automatically.** If pandas / numpy (or openpyxl, needed only for `.xlsx` input/output) are not found, the script runs `<this python> -m pip install ...` once and continues — inside a conda env or venv they go there, otherwise to the user site (`--user`). An already-installed but too old pandas (< 1.4) is reported, not upgraded, so other tools in the same environment are not affected. Tested with pandas 1.5, 2.2 and 3.0.

To manage packages yourself instead:

```bash
python -m pip install -r requirements.txt
python gene_search.py ... -no-install          # or: export GENE_SEARCH_NO_INSTALL=1
```

If the automatic install fails (no internet, SSL-inspection proxy, OS-managed Python), the script stops and prints the exact `pip install` command and the likely cause.

## Quick start

```bash
# one gene
python gene_search.py -i example/demo.hg38_multianno.txt -gene RNU4-2 -o results
# -> results/multianno_RNU4-2.csv

# a gene panel from Excel, Excel output
python gene_search.py -i example/demo.hg38_multianno.txt -l example/gene_list_example.xlsx -o results -f xlsx
# -> results/multianno_RNU4-2.xlsx, multianno_DMD.xlsx, multianno_BRCA1.xlsx
#    results/multianno_gene_list_example_summary.xlsx
```

More examples

```bash
# several genes on the command line (space- or comma-separated)
python gene_search.py -i S01.hg38_multianno.txt -g DMD BRCA1,BRCA2 -o results

# every multianno file in a folder (adds a "Sample" column)
python gene_search.py -i annovar_out/ -l panel.xlsx -o results -f xlsx

# the whole panel in ONE file -> results/multianno_panel.xlsx (adds a "Matched_Gene" column)
python gene_search.py -i annovar_out/ -l panel.xlsx -o results -f xlsx -merge

# only exonic / splicing / ncRNA_exonic variants
python gene_search.py -i S01.hg38_multianno.txt -l panel.xlsx -o results -func exonic splicing ncRNA_exonic

# prefix output names with a sample ID -> results/S01_multianno_DMD.csv
python gene_search.py -i S01.hg38_multianno.txt -g DMD -o results -prefix S01

# Windows paths
python gene_search.py -i "Y:\project\Annotated_VCF" -g RNU4-2 -o "C:\Users\me\Desktop\out"
```

Every option has a single-dash form (ANNOVAR style, `-gene`) and a double-dash form (`--gene`).

## Gene list format

Excel (`.xlsx`) with **column A = numbering, column B = gene symbol**. A header row is skipped automatically.

| | A | B |
|---|---|---|
| **1** | No | Gene |
| **2** | 1 | RNU4-2 |
| **3** | 2 | DMD |
| **4** | 3 | BRCA1 |

- A list with only one column (genes in column A) also works.
- Genes in another column: `-list-col C` (or `-list-col 3`).
- `.csv`, `.tsv` and `.txt` lists with the same layout, or one gene per line, are accepted too.
- Duplicates and blank cells are ignored. Entries that look like dates (Excel turning `MARCH1` or `SEPT9` into dates) are skipped with a warning.

See `example/gene_list_example.xlsx`.

## Options

| Option | Description |
|---|---|
| `-i`, `-input` | multianno file(s) (`.txt`/`.csv`, `.gz` allowed), a directory, or a glob pattern such as `"out/*_multianno.txt"` |
| `-g`, `-gene` | gene symbol(s) |
| `-l`, `-list` | gene list file (`.xlsx`/`.csv`/`.tsv`/`.txt`) |
| `-list-col` | list column holding the symbols (default `B`) |
| `-o`, `-out` | output directory (default: current directory, created if missing) |
| `-f`, `-format` | `csv` (default) or `xlsx` |
| `-merge` | write all genes to one file `multianno_<NAME>.<fmt>` |
| `-name` | label for the merged / summary file (default: gene-list file name) |
| `-prefix` | prefix for every output file name, e.g. a sample ID |
| `-genecol` | gene column to search (default `Gene.refGene`; e.g. `Gene.ensGene`, `Gene.knownGene`) |
| `-func` | keep only these `Func.*` categories, e.g. `exonic splicing ncRNA_exonic` |
| `-keep-intergenic` | also keep intergenic rows (see *Matching rules*) |
| `-case-sensitive` | case-sensitive symbol matching (default: case-insensitive) |
| `-pattern` | file patterns used when `-i` is a directory (default `*multianno.txt`, `*multianno.txt.gz`, `*multianno.csv`, `*multianno.csv.gz`) |
| `-chunksize` | rows read at a time (default 200000) |
| `-xlsx-text` | keep every xlsx cell as text (default: numeric columns become numbers so they sort/filter correctly in Excel) |
| `-no-install` | do not pip-install missing packages (same as `GENE_SEARCH_NO_INSTALL=1`) |
| `-q` | quiet: warnings and errors only |

## Output

| File | When | Content |
|---|---|---|
| `multianno_<GENE>.csv/.xlsx` | default | all columns of the original table for that gene |
| `multianno_<NAME>.csv/.xlsx` | `-merge` | every target gene in one table + `Matched_Gene` column |
| `multianno_<NAME>_summary.csv/.xlsx` | more than one gene | variants per gene, samples, `Func` breakdown, rows removed by the Func filter, output file |

- With more than one input file a `Sample` column (taken from the file name, e.g. `S01.hg38_multianno.txt` → `S01`) is added as the first column.
- Genes with no variants get no file; they are listed in the summary and in the final warning.
- xlsx output has a frozen header row and auto-filter. Above Excel's 1,048,576-row limit the file is written as CSV instead.

## Matching rules

1. **Exact symbol match.** `Gene.refGene` is split on `;`, `,` and the escaped `\x3b`, and each symbol is compared with the target. `DMD` does not match `DMD-AS1`; a row annotated `BRCA1;NBR2` matches both `BRCA1` and `NBR2`.
2. **Case-insensitive** by default (`dmd` = `DMD`). Use `-case-sensitive` to turn this off.
3. **Intergenic rows are excluded by default.** For `Func.refGene = intergenic`, ANNOVAR lists the two *flanking* genes, which can be tens of kb away, so such a row is not a variant in that gene. `upstream`/`downstream` rows (within 1 kb) are kept. Add `-keep-intergenic` to include intergenic rows; the run log reports how many were excluded.
4. **`-func`** keeps only the listed categories (a row such as `exonic;splicing` is kept if any part matches). When `-func` is given, intergenic rows are kept only if `intergenic` is listed.

## Notes

- **Open CSV results in Excel with care.** Excel converts some gene symbols (`MARCH1`, `SEPT9`, `DEC1`) into dates when a CSV is opened directly. Use `-f xlsx` when the result is meant for Excel.
- Large files are read in chunks, so memory use stays low; a 1-million-row multianno table is processed in a few seconds.
- Old ANNOVAR versions that write a single `Otherinfo` header with several values behind it are handled: the extra columns are named `Otherinfo1`, `Otherinfo2`, …
- Gene symbols must match the annotation database (`refGene` uses current HGNC symbols). If a gene is not found, check for a renamed symbol or alias, or try `-genecol Gene.ensGene`.

## Example data

`example/demo.hg38_multianno.txt` is a small **synthetic** table (made-up variants) in ANNOVAR format for trying the tool; `example/gene_list_example.xlsx` is a 4-gene list in the A = No. / B = Gene layout.

## Author

Seung Won Chae — Molecular Diagnostics Laboratory, Department of Laboratory Medicine, Seoul National University Hospital

## License

MIT — see [LICENSE](LICENSE).

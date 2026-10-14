# Manual curation of pipeline output ("Review" cases)

This document explains what the `Review` label in `table_results.txt` means, how to
decide whether a flagged gene family is a false positive or a false negative, and how
to feed your curated models back into the pipeline.

---

## 1. What `Review` actually means

`Review` is **not** a quality assessment of the gene models. It is purely a
**count-versus-expectation heuristic**. It is assigned by
[`Scripts/table_results.py`](Scripts/table_results.py) by comparing:

- `Step_1_total` — the number of genes recovered in Step 1 (BLAST + domain searches)
- the last Step 2 count (`Total number of identified genes clustering highly identical
  sequences and proteins shorter than 30 aa`) — the number recovered by BITACORA
- the **`Average number of genes`** column of your input `gene_families.xlsx`

If neither step lands near the expected number, the row is labelled `Review`.

So the observation that `Review` rows "coincide with comparatively few or many annotated
genes" is expected: that is the definition of the label, not a correlation.

### The exact rules

**If `Average number of genes` is a range (e.g. `300-400`):**

| Condition | Label |
|---|---|
| Both Step 1 and Step 2 counts inside the range | `Step 2` (or `Step 1` if Step 1 is larger, or Bitacora = `No`) |
| Only the Step 2 count inside the range | `Step 2` |
| Only the Step 1 count inside the range | `Step 1` |
| **Neither count inside the range** | **`Review`** |

**If `Average number of genes` is a single number (e.g. `9`):**

| Condition | Label |
|---|---|
| Expected is `1`, but neither step returned exactly 1 | **`Review`** |
| Closer count differs from the expected number by **more than 10** | **`Review`** |
| Step 2 is closer but the `Bitacora` column is `No` | **`Review`** |
| Otherwise | `Step 1` / `Step 2` |

### Two consequences worth knowing

- The tolerance is an **absolute ±10**, not a percentage. It is strict for large
  families and loose for small ones. If you expect real variation in a family, give a
  **range** rather than a single number in the Excel sheet.
- Any family declared as **single-copy (`1`)** is flagged whenever the count is not
  exactly 1, and any family with **`Bitacora` = `No`** is flagged easily, because its
  Step 2 count is written as `0`. Check these two columns in your spreadsheet first —
  they often account for a large share of `Review` rows.

The `Average number of genes` values distributed with GAGA are **averages across ant
genomes**. A genuine biological difference in your species will produce `Review`, and
the correct outcome of curation may well be "the pipeline is right, accept it".

---

## 2. Which column do I edit?

The decision column is **column 12**, written between `Average number in ant genomes`
and `TRIMMED`. It contains `Step 1`, `Step 2`, `Review`, or `Manual`.

> **Note:** the GAGA documentation refers to this as "column 13". That is off by one.
> [`Scripts/gen_final_output.py`](Scripts/gen_final_output.py) reads index `11`,
> i.e. the **12th** column. Edit column 12.

Also note that `table_results.py` writes a header row with only 12 names for 13 data
columns, so the decision column has no header label. Regenerating the table afterwards
with [`Scripts/table_results_after_running.py`](Scripts/table_results_after_running.py)
produces a correct header in which this column is named **`FINAL RESULTS`**.

---

## 3. Which files to look at

For each flagged family, under
`OUTPUT_DIR/gene_families_pipeline/GENE-FAMILY-NAME/`:

| File | What it tells you |
|---|---|
| `Step1/GENEFAMILY_merged.txt` | Per-protein table: hit start, hit end, protein length, evidence source, and the `complete` / `parcial` / `chimera` label |
| `Step1/GENEFAMILY_annot_genes.gff3` (and `_trimmed.gff3`) | Step 1 models |
| `Step1/temporary_files/*_blast_output.txt`, `*_parsed_blast.txt`, `*_parsed_domain.txt` | Raw BLAST and InterPro/Pfam evidence behind every call |
| `Step2_bitacora/.../GENEFAMILY_genomic_and_annotated_proteins*.gff3` | BITACORA models, including newly found genes absent from the original GFF |
| `Step2_bitacora/GENOME_genecounts_summary.txt` | The counts used to build `table_results.txt` |

### How the `complete` / `parcial` / `chimera` labels are set

From [`Scripts/merge_blast_domain.py`](Scripts/merge_blast_domain.py), using the
`Average protein length` and `Peptide start length` columns of your Excel sheet:

- **`chimera`** — twice the aligned region is still shorter than the protein. The model
  is far longer than the aligned domain, suggesting two genes merged into one.
- **`parcial`** — the aligned region is under **60%** of the family's average protein
  length. A fragment: possible truncation, frameshift, or broken splice site.
- **`complete`** — anything else.

These three counts are columns 2–4 of `table_results.txt`, so a `Review` row with many
`parcial` or `chimera` entries tells you immediately where to look.

---

## 4. Reading the count columns

The most common mistake when triaging `Review` rows is comparing `Step_1_total`
directly with the last Step 2 count. **These two columns do not measure the same
thing**, and comparing them produces discrepancies that look alarming but are
artefacts.

| Col | Name | What it actually counts |
|---|---|---|
| 5 | `Step_1_total` | **Rows in `GENEFAMILY_merged.txt`** — protein entries, not genes. Chimeric proteins split into `_split1`, `_split2`… contribute one row *each* |
| 6 | Number of annotated genes identified | Step 2 count derived from the annotated proteome |
| 7 | Number of putative not annotated genes | Genes BITACORA found in the genome that were absent from the input GFF |
| 8 | Total (Annotated + Genomic) | Step 2 raw total, before any clustering |
| 9 | …clustering identical sequences | After collapsing identical sequences |
| 10 | …clustering highly identical + <30 aa removed | **The gene-level figure. This is the number the `Review` logic uses** |

Two consequences:

- `Step_1_total` is **inflated**: it counts split fragments separately, and when
  `TRIMMED=No` the final GFF3 keeps only `_split1` and discards the rest, so the
  column overstates what Step 1 would actually deliver.
- Column 10 is **deflated** relative to the raw search, because clustering and the
  30 aa filter have been applied.

So the honest comparison is **column 5 against column 10 only after you understand
the collapse in between**. The drop from column 8 to column 10 measures how much
redundancy was in the input. A large drop is normal for proteomes carrying many
isoforms or fragmentary models; a small drop means the input was already clean.

This matters most when the same family is compared across genomes whose annotations
differ in maturity or redundancy. Two genomes can differ severalfold on `Step_1_total`
and agree closely on column 10, because the difference was redundancy in the input
rather than biology. A Step 2 count that looks lower than Step 1 usually means the
input was more redundant to begin with, not that BITACORA missed genes — check the
column 8 → column 10 drop before concluding anything.

**For comparing gene counts across species, use column 10 consistently for every
genome.** It is the only figure that is comparable across assemblies with annotations
of differing maturity. Use `Step_1_total` to understand a single genome's Step 1
output, never to compare genomes.

### Quick checks

Collapse split fragments, to see Step 1's effective gene count:

```
awk '{sub(/_split[0-9]+$/,"",$1); print $1}' GENEFAMILY_merged.txt | sort -u | wc -l
```

Tabulate the evidence source and the per-protein labels (columns 5 and 6 of
`merged.txt`):

```
cut -f5 GENEFAMILY_merged.txt | sort | uniq -c    # evidence source
cut -f6 GENEFAMILY_merged.txt | sort | uniq -c    # complete / parcial / chimera
```

Entries supported only by a domain match, with no blastp support against the query
database, are the usual false positives — especially for families defined by a shared
domain rather than by orthology.

### Caveat: the headers are wrong in `genome` and `protein` mode

When the `Bitacora` column of your spreadsheet is set to `genome` or `protein`,
`run_analysis.py` writes `GENOME_genecounts_summary.txt` itself rather than taking it
from BITACORA, and it fills some fields with constants. The signature to look for,
across the rows of one family:

- **column 7 is exactly `0`** in every genome, and
- **column 8 is identical to column 6** in every genome.

If you see that, those two columns are not measurements. In `genome` mode the figure
under *"Number of annotated genes identified"* is really the count from the genome
search, and the `0` under *"Number of putative not annotated genes"* is a literal
constant. Step 1 and Step 2 are then two **independent** estimates — annotation-based
and genome-based — rather than Step 2 being Step 1 plus newly found genes. Check the
`Bitacora` column before interpreting these two fields.

---

## 5. Deciding false positive vs. false negative

### Fewer genes than expected (possible false negatives)

The gene may be genuinely absent, or it may have been missed. Check:

1. Look for **`parcial`** models in `_merged.txt` — a truncated model often means a
   frameshift, a premature stop, or an incorrect splice site, not a missing gene.
2. Check the **Step 2 / BITACORA** output: BITACORA searches the genome directly and
   recovers genes absent from the input GFF. If Step 2 found more than Step 1, the
   genes were missing from your original structural annotation, not from the genome.
3. **tBLASTn the family database against the genome** to see whether there is sequence
   similarity in a region with no gene model at all. Such regions are the candidates
   for manual annotation.
4. Check whether the hits fall near **scaffold or contig edges**, or in assembly gaps —
   a frequent cause of genuinely unrecoverable genes.

If similarity exists but no acceptable model was built, annotate the gene manually
(GAGA used an Apollo browser) and save a corrected GFF3.

### More genes than expected (possible false positives)

1. Confirm the extra copies are real duplications of the same gene, rather than
   **different genes sharing only a protein domain**. Align the recovered proteins
   (e.g. MAFFT) and build a quick tree — true tandem duplicates cluster tightly;
   domain-only relatives sit well outside.
2. Check the **evidence source** column of `_merged.txt`. Models supported only by a
   domain hit, with no BLAST support against the family database, are the usual
   false positives.
3. Check for **overlapping or redundant models** at the same locus, and for very short
   proteins that survived clustering.
4. Inspect the **`chimera`** entries: a chimera inflates nothing by itself, but it
   signals that neighbouring genes were merged and may need splitting.

Delete incorrect models from the GFF3 and keep the rest.

---

## 6. Feeding curated models back in

For each family you curated:

1. Create a folder named **`Manual`** next to `Step1` and `Step2_bitacora`:
   `OUTPUT_DIR/gene_families_pipeline/GENE-FAMILY-NAME/Manual/`
2. Put your curated GFF3 inside it, named exactly:
   `GENEFAMILY_manualmodels.gff3`
3. In `table_results.txt`, change **column 12** for that family from `Review` to
   **`Manual`**.

For families where you inspected the output and concluded the pipeline was right,
simply set column 12 to `Step 1` or `Step 2` — whichever output you want to keep.

Every row must end up as `Step 1`, `Step 2`, or `Manual`. A row left as `Review`
produces **no output** for that family in the final step.

Then run the final step:

```
python3 Scripts/gen_final_output.py <pipeline_dir> <genome.fasta> <gene_families.xlsx> <pipeline_output> <out_dir>
```

> **Important:** `gen_final_output.py` matches rows of `table_results.txt` to rows of
> your Excel sheet **by position**. Do not add, remove, or reorder rows in either file,
> and check that the table has one row per family in the spreadsheet before running it.
> If a family failed and produced no row, the remaining families will be silently
> shifted and assigned the wrong GFF3.

---

## 7. Chemosensory receptors

ORs and GRs were handled with a separate, dedicated pipeline in GAGA because of their
tandem arrangement and fast evolution. If you are annotating these, see
[GAGA/04_Gene_re-annotation/Chemosensory_gene_families](https://github.com/schraderL/GAGA/tree/main/04_Gene_re-annotation/Chemosensory_gene_families).

## 8. Reference

The GAGA re-annotation protocol this pipeline was built for:
<https://github.com/schraderL/GAGA/tree/main/04_Gene_re-annotation>

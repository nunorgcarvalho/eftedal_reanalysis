# Data

## `eftedal_osf/`

Correlations between pairs of relatives on Norwegian national test scores, as
publicly deposited by Eftedal et al. (2025, *PNAS*, "Beyond additive genetic
effects: Explaining family resemblance in school performance across millions of
pairs of Norwegian relatives"; OSF, https://doi.org/10.17605/OSF.IO/JFW2K).

These are aggregate summary statistics: one correlation (and standard error) per
category of relative, 82 categories per phenotype. `README_osf.txt` is the
authors' own description of the files and columns. The only change made here is
renaming the `remove missing/` folder to `remove_missing/`.

The main phenotype analysed is `MRE8` (average of the 8th-grade math, Norwegian
reading and English tests), as in the main text of Eftedal et al.

### Known issue: duplicated row

In every non-adoptive file, the row for **2nd cousins / `other` / `ssf`** is an exact
copy of the row for **2nd cousins / `paternal` / `ssf`** (identical `corr`, `SE` and
`npairs`). This is almost certainly an export error. The pair counts in Table 1 of
the paper imply the `other`/`ssf` category has roughly 138,000 pairs rather than the
59,735 listed. The row is flagged as `suspect_duplicate` when the data are loaded.
It is kept by default (to reproduce published numbers) and can be dropped as a
sensitivity check.

## `relationship_types.csv`

Our annotation of the 16 relationship types (genetic structure, family degree,
shared-environment class memberships), joined to the correlation files by `rel`.
See the column descriptions at the top of `src/eftedal/data.py`.

# ASHG 2026 results snapshot

The analysis outputs that the ASHG 2026 materials (poster, supplemental slides, Q&A notes) were
built from. The code and the built PDFs are at tag `ASHG2026` (commit 3814faf). `results/` is
gitignored in the main history, so it is archived here, on an orphan branch, with no shared
history.

| path | what |
|---|---|
| `results/` | every CSV and figure under `results/` when the materials were built (2026-10-20). Includes `results/figures/poster/`, the poster and slide figures, and `results/qna/`, the Q&A tables and numbers. |
| `env/python_requirements.txt` | `pip freeze` of the Python environment |
| `env/R_sessionInfo.txt` | R `sessionInfo()` with the plotting packages loaded |

To rebuild the PDFs exactly, check out tag `ASHG2026` in its own worktree, copy this `results/`
into it, and run `bash ASHG/build.sh`. Or regenerate `results/` at that tag with the pipeline in
`ASHG/README.md`. That pipeline is seeded, but package versions may differ from `env/`.

# Data index

## `sequences/`

| File | Origin | Contents |
|---|---|---|
`output_X_0_initial.txt` | **Pinglay lab** | order-5 Markov initialization, 19,968 bp |
`output_X_0_final.txt` | **Pinglay lab** | Ledidi + Borzoi output, 19,968 bp, 84 min runtime |
`gfp_scarlet_halo.fasta` | **Pinglay lab** | eGFP 720 bp, mScarlet 699 bp, HaloTag7 894 bp |
`initial_with_cds.fasta` | derived (this work) | CDS-only control, FASTA |
`output_X_0_initial_with_cds.txt` | derived (this work) | same payload, lab's text format |
`three_arms.fasta` | derived (this work) | initial / control / final in one FASTA |
`designed_arms/*.fasta` | **this work** | the twelve Ledidi ablation arms (A, B, C, E, F, G, H, I, J, K, L, M) |

Note: the third CDS in `output_X_0_final.txt` is **not** HaloTag. It is a 702 bp / 233 aa
RFP-family protein. See [../docs/02-sequence-analysis.md](../docs/02-sequence-analysis.md).

## `predictions/`

| File | Context | Contents |
|---|---|---|
`scores_aavs1.json` | AAVS1 | 15 designs × 19 tracks × 3 genes. **The headline result.** |
`scores_loci.json` | 11 real loci | 3 designs × 11 genomic contexts; tests random-genomic padding |
`arms_aavs1.json` | AAVS1 | 9 perturbation arms (transplant/scramble), 10 tracks |
`arms_N.json` | all-N | same 9 arms |
`arms_orig.json` | random | same 9 arms |
`flank_test.json` | 6 contexts | 3 designs × 6 flanking contexts |
`ism_aavs1_rep0.json` | AAVS1 | 5,400 single-base mutants, 3 readouts each |
`randomflank_means.json.gz` | random | **superseded.** Full 5-track × 6,144-bin profiles, 3 arms |
`randomflank_means_all.json.gz` | random | **superseded.** Same for all 9 perturbation arms |

The two `randomflank_*` files are retained for the record only. They were computed with
composition-matched random padding, which [../docs/08-corrections.md](../docs/08-corrections.md)
shows suppresses the designs up to 437×. Do not draw conclusions from them.

## `context/`

| File | Contents |
|---|---|
`aavs1_flank.txt` | 524,288 bp of hg38 `chr19:54,862,840-55,387,128` (AAVS1 safe harbour). Public domain. |
`context_meta.json` | padding geometry for the original random-flank runs |

The eleven 524 kb random genomic loci used in
[docs/10](../docs/10-multitrack-chromatin.md#6-random-genomic-loci-are-safe-synthetic-random-padding-is-not)
are **not** committed (5.8 MB of redundant hg38). They are reproducible exactly:
`modal run code/borzoi_modal_app.py::main --step fetchloci` samples them from a fixed
seed (20261027) and fetches them from the UCSC REST API.

## `tables/`

Eleven CSV tables, documented in the [README](../README.md#result-tables).
`table8_track_coupling.csv` is regenerated from `table2` by `code/analyze_track_coupling.py`;
`table9_procapnet_motifs.csv` by `code/analyze_procapnet_motifs.py`;
`table11_random_loci.csv` by `code/make_table11_loci.py`.

## Coordinate conventions

All sequence coordinates are 0-based, half-open, relative to the 19,968 bp design.
Output bin $b$ of the design corresponds to design base $32b$, and to Borzoi output bin $2760 + b$.

### A note on `table10`'s columns

`cage_pct_of_target_AREA` and `cage_pct_of_target_PEAK` are the same prediction read
two ways, and they disagree by up to 2× — see
[correction 11](../docs/08-corrections.md). The area reading sums 5 bins from the
designed TSS, matching slide 50's description of the objective. The peak reading
takes the maximum over a 10-bin window, which is what this repository used
originally. **Prefer the area column** unless you are reproducing an older figure.

`initiation_breadth` is area ÷ peak. A value near 1.0 means the entire CAGE signal
sits in one 32 bp bin, which is what the short arms here produce. The lab's `final`
is at 2.14, i.e. a spread-out initiation zone.

`n_edits` counts substitutions against the **CDS-only control**, which is what every
designed arm was initialised from. `initial` shows 1,616 because the Markov backbone
carries different CDS sequence than the control does; it is not an edited design.

`genes_lost` counts genes whose CAGE peak fell below 10% of that gene's own target.

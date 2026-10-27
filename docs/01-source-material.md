# 01. Source material: the draft presentation

54 slides, 720×405 pt, a working draft of a paper from the Pinglay lab. Credits on slide 1 to
Brian Hie, Sud Pinglay and Zach Amador. Rendered for reading with Swift + PDFKit (no pdftoppm,
PyMuPDF or PIL available on the machine used).

## Overall structure

| Figure | Slides | Topic |
|---|---|---|
| 1 | ~8–21 | Establishing Ledidi as a sequence editor |
| 2 | ~22–35 | Steady-state chromatin accessibility |
| 3 | ~36–43 | Dynamic / cell-type-specific accessibility |
| 4 | **44–54** | **Steady-state gene expression (the subject of this audit)** |

Slide 7 frames the programme as a 2×2 "2026 Pinglay lab moonshot": steady-state versus dynamic,
crossed with accessibility versus expression. The dynamic-expression quadrant is untouched.

## Model stack

| Model | Bin size | Role |
|---|---|---|
| Enformer | 128 bp | design oracle inside Ledidi |
| Borzoi | 32 bp | design oracle inside Ledidi |
| AlphaGenome | 1 bp | held out for validation |
| Evo2-7B | — | baseline / initializer |
| NTv3 | — | independent annotator, never in the design loop |

## Figure 4 in detail (slides 44–54)

Slide 46 sets up the logic: CAGE reports the TSS and 5′UTR, RNA-seq reports the transcript body.
Two orthogonal design handles.

### Slide 47, RNA-seq alone
- Target RNA blocks 7 / 5 / 3. After editing, RNA matches almost exactly (~7.0/5.5/3.7, slightly overshooting).
- CAGE (K562, `CNhs11250`) stays near zero, max 0.5–0.8, against a desired 7/5/3.
- DNase (`ENCFF413AHU`) is flat except one spurious spike (~1.25) at ~12,000 bp, **not** at any designed position.

Transcript signal with no initiation and no accessibility.

### Slide 48, CAGE alone
- Edited CAGE hits 28.7 / 21.7 / 12.1 against target 28 / 20 / 12.
- RNA is near zero (max ~1.3, located at bins 0–50, not at the genes); narrow spikes at the TSSs, no transcript body.
- DNase goes **runaway**: broad peaks up to 6.2 around 13–14 kb, mostly off-target.

Initiation without elongation, plus chaotic accessibility.

### Slide 49, CAGE + RNA-seq jointly
- CAGE hits 28.5 / 20 / 12 **and** RNA hits ~7 / 5 / 3.5 at the right bins.
- The never-specified DNase track now aligns with the three designed regions: ~1.5 at 4–5 kb, sharp ~10.4 at ~9.5 kb, ~4.2 at 14–15 kb.

Note the nuance: DNase peaks are correctly *located* but not correctly *scaled*. Against a desired
7 / 5 / 3 the achieved values are roughly 1.5 / 10.4 / 4.2, so the highest-specified position
received the least signal.

### Slides 50–52, the self-consistency machinery

Four GEXP parameters: CAGE area (CAGE summed over the 5′UTR plus 32 bp upstream), 5′UTR length,
3′UTR length, RNA-seq height.

Empirical UTR lengths over 3,000 genes:

| | mean | median | Borzoi bins used |
|---|---|---|---|
| 5′UTR | 163.34 bp | 103 bp | **4 bins (128 bp)** |
| 3′UTR | 2,055.22 bp | 1,346 bp | **33 bins (1,056 bp)** |

CAGE-area-versus-RNA-height is nonlinear and saturating (CAGE plateaus around 57–60 once RNA height
exceeds ~10). A 4-knot natural spline beats a linear fit (10-fold CV RMSE ~13.8 versus ~18, best
AIC/BIC) and is used to derive a consistent CAGE target from a desired RNA height. This prevents
requesting physically impossible readout combinations.

**This target model was independently confirmed in this audit.** See [02](02-sequence-analysis.md):
measured CAGE peaks land at design bins 297–298, 458–460 and 137–140 against predicted TSS bins of
297, 458 and 135 derived from "TSS = ATG − 128".

### Slides 53–54, validation
Borzoi RNA+CAGE design validated on held-out oracles (Borzoi DNase, Enformer DNase, AlphaGenome
ATAC/DNase all peak at designed promoters). NTv3 segmentation, never in the design loop, shows crisp
5′UTR / 3′UTR / enhancer / promoter blocks after editing versus almost nothing before.

## Other findings worth retaining

- **Figure 1:** Ledidi takes ~30–45 min versus ~18 h for Evo2-7B beam search (~40× faster), and
  roughly halves L2 on all six target patterns with much tighter variance.
- **Initializer ranking:** Markov order-5 (3,072 parameters) is best at ~0.150, beating Evo2-7B
  (7×10⁹ parameters, ~0.170). Uniform random is worst.
- **"The genome is a graveyard":** Markov O5 adds ~+41 motifs. Genome-like initializations carry
  dead motifs and the optimizer must install fresh ones.
- **Attribution maps show the oracles are essentially local** (dominant diagonal), in tension with
  the multi-kb ambition. The necessity/sufficiency result in [04](04-context-robustness.md) pushes
  back on this.
- **Experimentally, positions transfer but heights do not.** The large staircase essentially fails
  and the dominant error mode becomes "level".
- **Counterintuitive negative correlation** between specified peak height and attribution/seqlet
  density ($\rho$ from −0.66 to −0.94).
- **Figure 3 headline positive:** dynamic, cell-type-specific orthogonal on/off accessibility works
  experimentally (mESC versus neuron).

## Unresolved metric conflicts in the draft

Three places where two metrics disagree and the deck does not adjudicate:

1. Which initializer is best: RMSE favours Markov O5, Pearson favours Evo2.
2. Which oracle is best: correlation favours one, z-space error the other.
3. Whether tall or short peaks fail worse: absolute versus height-normalized loss disagree.

## Draft artifacts noted

- Slide 4: content clipped off-slide.
- Slides 25 and 26 share a title.
- Slide 51: "borozi" misspelled.
- Open placeholders: slide 21 (global motif analysis), slide 24 (DAF-seq versus ATAC correlation),
  slide 43 (five analyses).

The slide 21 placeholder is partially filled by [05](05-ism-mechanism.md).

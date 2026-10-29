# What CAGE is actually doing: it installs CpG islands

**The question this document answers**, as posed by the authors: the design goal is
accessible promoters, that is DNase peaks at each TSS. CAGE is in the objective
because the sequence elements that drive CAGE signal are thought to be the same
ones that create those peaks. *Why* that works was unknown, and the stated ideal is
to specify only output readouts in the design process — RNA abundance — and have
promoters and accessibility emerge.

**The answer: CAGE's effect on accessibility is almost entirely mediated by CpG
dinucleotide density.** Optimising CAGE installs CpGs, and CpG density is what
Borzoi reads as an accessible promoter. CAGE contributes almost nothing beyond
that.

Everything below comes from 96 promoter windows across the 32 designs already in
this repository. No new model runs. Reproduce with
[`code/analyze_cpg_mediation.py`](../code/analyze_cpg_mediation.py), which needs no
GPU.

## 1. CpG density is the sequence feature that predicts accessibility

The 600 bp upstream of each ATG, 96 windows:

| Sequence feature | Pearson vs DNase | Spearman |
|---|---|---|
| **CpG dinucleotide count** | **+0.592** | **+0.608** |
| CpG observed/expected | +0.442 | +0.533 |
| SP1/KLF GC-box | +0.440 | +0.375 |
| GC % | +0.410 | +0.382 |
| CREB/ATF | +0.251 | +0.268 |
| ETS | +0.298 | +0.229 |
| NFY/CCAAT | +0.065 | +0.190 |
| NRF1 | +0.000 | +0.166 |
| YY1 | −0.010 | +0.118 |
| **TATA box** | **−0.233** | **−0.176** |

CpG count beats every individual motif. The TATA box is **negatively** associated,
which is the expected contrast between CpG-island promoters, constitutively
nucleosome-depleted, and TATA promoters, which are regulated and closed by default.

## 2. CAGE in the objective is what installs the CpGs

Change in CpG count in the promoter window, against the unedited control:

| Objective class | n | mean CpG gain | range |
|---|---|---|---|
| no CAGE in the objective | 7 | **+4.0** | 0 to +9 |
| CAGE in the objective | 23 | **+8.4** | +4 to +17 |
| lab's `final` | — | **+17.0** | — |

## 3. The mediation is near-total

| | |
|---|---|
| r(CAGE, DNase) | +0.460 |
| r(CpG, DNase) | +0.592 |
| r(CAGE, CpG) | +0.674 |
| **r(CAGE, DNase) with CpG held fixed** | **+0.102** |
| **r(CpG, DNase) with CAGE held fixed** | **+0.429** |

**Hold CpG constant and CAGE explains almost nothing about accessibility (+0.10).
Hold CAGE constant and CpG still explains a great deal (+0.43).** CAGE acts on
accessibility *through* CpG content, not independently of it.

## 4. The dose-response, and why asking for DNase directly fails

| Objective | CpG in window | DNase range |
|---|---|---|
| no design at all | 8 | 0.07 – 0.07 |
| RNA only | 8 | 0.11 – 0.15 |
| **RNA + DNase in the objective** | **10** | **0.08 – 0.11** |
| RNA + 6 marks including DNase | 11 | 0.06 – 0.72 |
| CAGE + RNA | 14 | 0.15 – 0.80 |
| lab's `final` | **25** | **0.37 – 5.65** |

CpG 8 → 25 tracks DNase 0.07 → 5.65, monotonically.

**Asking for DNase directly does not work.** Arm C put DNase in the objective and
reached 0.08–0.11, no better than RNA alone and barely above the untouched control,
because it added only 2 CpGs. DNase is a broad, low-amplitude target with weak
gradients; CAGE is a sharp, high-amplitude one whose cheapest solutions happen to
be CpG-rich. **You cannot get accessibility by requesting accessibility. You get it
by requesting something whose cheapest solution is a CpG island.**

This also reconciles the single-nucleotide result in
[docs/05](05-ism-mechanism.md), the 5,400-mutant scan: per-base CAGE and DNase sensitivity correlate at
only +0.351 / −0.041 / +0.169, and the DNase-sensitive bases sit 60–150 bp further
upstream than the CAGE-sensitive ones. CAGE and DNase are not driven by the same
bases. They are driven by the same *regional composition*.

## 5. What this implies for the design loop, and the test that refuted it

CpG density is a **property of the sequence**, not a predicted readout. It can be
constrained directly, for free, with no model in the loop. The obvious implication
was that the objective could be `RNA abundance` + a direct CpG constraint, with
CAGE dropped, leaving exactly one predicted readout.

**That was tested and it does not work.** Full numbers in
[correction 17](08-corrections.md). An RNA-only objective plus a one-sided CpG
hinge targeting 25 CpGs per promoter window:

- **The constraint was met and accessibility still did not follow at two of three
  promoters.** mScarlet reached 29 CpGs, more than the 20 in the lab's `final`,
  and its DNase stayed at 0.13 against `final`'s 0.37. CpG count is not
  sufficient.
- **But eGFP went from 0.07 to 3.16 with no CAGE in the loss at all**, against 3.70
  for the lab's design at that promoter. No CAGE-free objective had previously
  exceeded 0.86. So an accessible promoter does not *require* CAGE in the
  objective.
- **Accessibility arrived without initiation.** CAGE stayed at 10.7–14.8% of
  target. The result is an open CpG-rich region with no TSS.

So the mediation result above is correct about what it measures and was read too
strongly. CpG density is the channel through which CAGE acts *on accessibility*.
CAGE is separately specifying where initiation happens and in which direction, and
that has no CpG substitute. **Accessibility and initiation are two requirements,
not one, and only the first has a free sequence-level proxy.**

## Caveats

- The 32 designs are not independent; many share most of their objective.
- CpG count, CpG o/e, GC% and SP1 GC-boxes are collinear by construction — a
  GC-box contains CpGs. The claim is that *regional CpG density* is what matters,
  not that one specific motif was identified.
- **This describes Borzoi's internal logic, not necessarily the cell.** Borzoi may
  have learned "CpG island ⇒ accessible" as a statistical regularity. That
  regularity is real in genomes, since CpG islands are genuinely
  nucleosome-depleted, so the model is probably right for a defensible reason, but
  only a wet-lab test of a CpG-constrained design would settle it.
- H3K4me3 predicts DNase better than any sequence feature (+0.902 Spearman), but it
  is another Borzoi output, so that is internal consistency rather than independent
  evidence.

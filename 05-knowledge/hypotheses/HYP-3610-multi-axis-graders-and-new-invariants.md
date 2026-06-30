---
id: HYP-3610
title: MORE GRADERS / NEW INVARIANTS -- the apex problemscape is NOT one axis: a family of spread/concentration graders on the Z_p cores splits into >=3 INDEPENDENT axes -- (Axis 1) the GAP g (worst-mode, the THM-590 certificate; concentrated pole = the DOUBLET = the odd cycle C_p); (Axis 2) spectral FLATNESS F = GM/AM ~ -DIFFERENCE-SET-DEFECT Var(a(d)) (global difference-set-ness; concentrated pole = the ARC/interval) -- only 0.50 correlated with g, picking a DIFFERENT concentrated core; (Axis 3) IPR/entropy/CAYLEY-ODD-GIRTH ~ -|O| (effective support size). New invented metrics: difference-set defect, Cayley odd-girth (on-theme: doublet->girth p, diff-set->3, transitive->inf), spectral flatness, binding-mode k*. The difference set (Paley/Fano) is the common RANDOM/spread pole of all axes. Universal-grader candidates (any spectrum/associated graph): flatness, the gap, odd-girth/non-bipartiteness
status: VERIFIED (Spearman rank-correlation over all 127 Z_7 cores; extremal cores per axis). The multi-axis structure and the new metrics are computational facts; the existence/measure interpretation of the two axes is a CONJECTURE.
source: klein-2026-06-29-S22
depends_on:
  - THM-590    # the gap (Axis 1, the certificate)
  - HYP-3611   # the apex core atlas (gap as the first concentration index)
related:
  - HYP-3604   # the doublet = signless Laplacian of the odd cycle (Axis-1 pole)
  - HYP-3602   # mac-mini: intransitivity / odd cycle (the odd-girth grader's theme)
  - HYP-3599   # existence vs measure (the conjectured meaning of the two axes)
results:
  - 04-computation/grader_family_cores_klein.py
  - 05-knowledge/results/grader_family_cores_klein.out
  - 05-knowledge/results/grader_axes_extremal_klein.out
  - 05-knowledge/results/grader_universality_tournaments_klein.out
---

# HYP-3610 — more graders, new invariants: the apex is a multi-axis space

## The family (some standard, some invented here)
For a core `O subset Z_p`, autocorrelation `a(d)=#{x:x,x+d in O}`, spectrum `lam_k=|Ohat(k)|^2`:
- `g` = **gap** = `min_{k!=0} lam_k` (THM-590; the worst-mode concentration / the certificate).
- `F` = **spectral flatness** = `GM(lam_{1..p-1})/AM(...)` in `[0,1]` (1 = difference set).
- `D` = **difference-set DEFECT** = `Var_{d!=0}(a(d))` (**NEW**: distance from a perfect difference set; 0 = diff-set).
- `IPR` = `(sum lam)^2/(p sum lam^2)`, `Hs` = spectral entropy (effective support).
- `godd` = **Cayley ODD-GIRTH** of `Cay(Z_p, supp a)` (**NEW, on-theme**: shortest odd cycle).
- `k*` = **binding mode** = `argmin_{k!=0} lam_k` (which frequency is worst).

## The finding: >= 3 INDEPENDENT axes (Spearman over all 127 cores)
- **Axis 1 = `g`** (the gap). Correlations with the others are all weak-to-moderate (`g`-`F` = 0.50,
  `g`-`IPR` = 0.15, `g`-`D` = -0.44). Concentrated pole = the **DOUBLET** (`g=0.198`, = the odd cycle `C_p`,
  HYP-3604); spread pole = the **difference set** (`g=2`).
- **Axis 2 = `F` ~ -`D`** (`F`-`D` = -0.86: flatness and defect are the same axis). Only **0.50** correlated
  with `g` -- a GENUINELY DIFFERENT grader. Concentrated pole = the **ARC/interval** (`F=0.5`, `D=0.667`),
  NOT the doublet. Spread pole = the difference set (`F=1`, `D=0`).
- **Axis 3 = `IPR` ~ `Hs` ~ `godd` ~ -`|O|`** (`IPR`-`|O|` = -0.96): the effective-support / size axis.

## The key disagreement (why it is multi-axis, not one)
`g` and `F` pick DIFFERENT "most concentrated" cores:
```
doublet {0,1}   : g=0.198 (MIN, Axis-1 pole)   F=0.600   D=0.222
arc-4 {0,1,2,3} : g=0.308                       F=0.500 (MIN, Axis-2 pole)   D=0.667 (MAX)
```
The doublet is the gap-concentrated core (worst single mode = the certificate); the arc is the
flatness-concentrated core (most uneven spectrum). The **difference set is the common spread/random pole**
of every axis. So the apex is a multi-dimensional grading, not a single line.

## Interpretation (CONJECTURE)
The two genuine axes may be the two halves of the sigma-split (HYP-3599):
- **`g` = the EXISTENCE / certificate axis** (worst-mode, discrete, the odd-cycle/non-bipartiteness gap;
  binding = the doublet = the odd cycle `C_p`).
- **`F` = the MEASURE / equidistribution axis** (global flatness = how equidistributed the safe set is;
  binding = the arc/interval).
The difference set (flat on both) is the maximally-equidistributed, certificate-safe, "random" core.

## Universal graders (act on EVERY object)
By construction these need only a spectrum or an associated graph, so they grade cores, coverings (via the
lonely-set or descent spectra), even graphs `E_n`, and the tournament metagraph `G_n` alike:
- **spectral flatness** `F` (any spectrum) -- the structured<->random axis.
- **the gap / least eigenvalue** `g` (any spectrum) -- the certificate / worst-mode axis.
- **odd-girth / non-bipartiteness** (any associated graph) -- the intransitivity / cycle-length-scale axis:
  transitive/orderable = `inf` (the cusp), generic = `3` (the Condorcet triangle), the apex binding core =
  `p` (the full odd cycle). This is mac-mini's intransitivity (HYP-3602) as a graded invariant.
  **VERIFIED cross-object** (grader_universality_tournaments): on TOURNAMENTS the same odd-girth grader
  gives the n=4 Klein-four classes `(odd-girth) = {inf, 3, 3, 3}` -- the TRANSITIVE class is the orderable
  `inf` cusp (the gap=0 analog), the rest are intransitive `3` (Condorcet); the finer **cyclicity**
  (#3-cycles) grades them `{0,1,1,2}` (n=4), `{0,1,1,1,2,2,3,3,4,4,4,5}` (n=5). So the cusp/concentrated/
  spread grading structure of the apex `Z_7` cores recurs verbatim on the tournament side -- one invariant,
  all objects; the apex cusp (`g=0`, full `Z_7`) and the transitive tournament (orderable, odd-girth `inf`)
  are the SAME degenerate pole.

## Net
The apex problemscape is not one concentration line but a multi-axis space: the gap (certificate), the
flatness/defect (equidistribution), and the size/odd-girth (support). The difference set is the common
"random" pole; the doublet (odd cycle) and the arc (interval) are the two distinct "concentrated" poles.
New invariants: difference-set defect, Cayley odd-girth, flatness, binding-mode. Universal graders:
flatness, gap, odd-girth -- each defined on any object's spectrum/graph.

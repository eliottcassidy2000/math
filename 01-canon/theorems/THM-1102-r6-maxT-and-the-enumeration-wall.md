---
id: THM-1102
title: r=6 — MAX T COMPUTED FIRST, AND THE ENUMERATION WALL LOCATED. (I) Following the rule I got wrong at r=5, max T was computed BEFORE any finite-horn run: over all 792 seven-speed cores, scanning every quintuple of removed killers in a width-16 window, **max R = 1.85794** at core [1,2,4,7,9,11,12] with killers (158,160,162,164,166), and **max T = 308.4**, giving KB = 333. The worst case sits at offset 9 inside a window of width 16, so the window is wide enough for the answer to be trustworthy — a check I now run explicitly. (II) THE R-LADDER, extended: **0.51852 (r=2) → 0.73375 (r=3) → 0.98453 (r=4) → 1.28495 (r=5) → 1.85794 (r=6)**. (III) THE r=6 FINITE HORN IS INFEASIBLE, and this is the session's real finding: at KB = 333 roughly **3.64 × 10¹²** sextuples pass the covering-necessary condition — about 140 days of compute, and **13,783×** the r=5 count. The prune that made r=4 and r=5 possible has stopped working, for a structural reason: with six killers the condition Σ frac ≥ 1 is easy to satisfy (mean kill-fraction ≈ 0.13, so 6 × 0.13 ≈ 0.78 with a 6–9% tail), so it discards only ~92% of sextuples while the raw count exploded. (IV) A CORRECTION to my own earlier statements: I had said the method "dies at r ≥ 7 because the union bound needs 7 − r > 0". That describes the SUPERSEDED crude formulation of THM-1051. The current horn removes r−1 killers *exactly* and bounds only the last, so **no r appears anywhere** and there is no structural r-cap — only R < 1 matters. The wall at r=6 is therefore **computational, not mathematical**. (V) STATUS: r = 2,3,4,5 closed; **r = 6 open**, with the obstruction located precisely and quantified
status: (I) COMPUTED and window-checked — max T = 308.4 over all 792 cores, in two chunks (0–300 and 300–792), each run to a printed summary rather than read off a progress line. (II) measured. (III) MEASURED, not merely predicted: the 3.64 × 10¹² figure is an extrapolation from exact per-core counts on an 8-core sample, so it is order-of-magnitude, and the conclusion (infeasible by this method) is robust to a factor of 10. (IV) is a correction to my own prior claims in THM-1051/1093. (V) honest — **r=6 is NOT closed and I make no claim that it is**
source: kind-pasteur-2026-07-18-S128 (cont.64; owner: run the r=6 finite horn, computing max T first)
depends_on:
  - THM-1093    # the r=5 closure and the max-T rule this follows
  - THM-1081    # the R-ladder this extends
script: 04-computation/r6_maxT_kps_S128c64.py, r6_maxT_chunk_kps_S128c64.py, r6_feasibility_kps_S128c64.py (+ .out)
---

# THM-1102 — r=6: max T first, and where the enumeration wall is

## (I) Max T, computed first and window-checked

At r=5 I set the finite-horn bound from max k_removed and it was wrong; the bound is set by
**max T**, since the last killer is certified by the measure horn only once it exceeds T.
So this time T came first.

Over all **792 seven-speed cores**, every quintuple of removed killers in a width-16 window:

> **max R = 1.85794**, at core [1,2,4,7,9,11,12], killers (158,160,162,164,166), T = 308.4
> **max T = 308.4** over the R ≥ 1 region ⟹ **KB = 333**

24,598 + 16,298 quintuples have R ≥ 1, and the largest killer among them is 172. The worst
case sits at **offset 9 inside a window of width 16** — comfortably interior, so the window
is not truncating the answer. I now run that check explicitly rather than assuming it.

Both chunks (cores 0–300 and 300–792) were run to a printed summary. The first attempt at a
single 792-core run was killed mid-way at core 300 with max T still *rising* (129 → 225 →
308), and I did not use those numbers — a partial scan cannot set a bound.

## (II) The R-ladder

> **0.51852 (r=2) → 0.73375 (r=3) → 0.98453 (r=4) → 1.28495 (r=5) → 1.85794 (r=6)**

## (III) The enumeration wall

Last session I predicted r=6 was "where the enumeration finally becomes the binding
constraint rather than the mathematics". It is, and now it is measured:

| r | KB | sextuples/quintuples/… passing the prune | runtime |
|---|---|---|---|
| 4 | 400 | 1.43 × 10⁸ | ~25 min |
| 5 | 235 | 2.64 × 10⁸ | ~9 min |
| **6** | **333** | **≈ 3.64 × 10¹²** | **≈ 140 days** |

**13,783× the r=5 count.** The prune has stopped working, and the reason is structural
rather than incidental: a sextuple can only be uncertified if its six kill-sets cover the
core's safe (q,a) set, requiring Σ frac ≥ 1 — but with *six* killers and a mean kill-fraction
of ≈ 0.13, the sum sits at ≈ 0.78 with a 6–9% tail, so the condition discards only ~92% of
sextuples. At r=4 and r=5 the same condition discarded well over 99%.

The prune's power comes from needing a *large* deviation above the mean; as r grows the
required deviation shrinks toward the mean and the prune dissolves.

## (IV) A correction to my own earlier statements

I wrote in THM-1051 and repeated in THM-1093 that the method "dies at r ≥ 7 because the
union bound needs 7 − r > 0". That describes the **superseded** crude formulation, in which
all r killers were union-bounded together. The current horn removes r−1 killers *exactly*
and bounds only the last, so **r appears nowhere in the estimate** — the threshold is
1/(3L) regardless of r, exactly as noted in THM-1061 and then forgotten by me two sessions
later. There is no structural r-cap. The wall at r=6 is **computational**.

## (V) Status

- r = 2, 3, 4, 5: **closed**.
- r = 6: **open**. max T and KB are known, the measure horn's failure region is mapped
  (max R = 1.858), and the finite horn is infeasible by this enumeration.

## Named next
- r=6 needs a better certificate than enumeration. Three candidates, in order of promise:
  (a) **strengthen the prune** — Σ frac ≥ 1 is weak because it ignores overlap; a bound
  using pairwise |kill(kᵢ) ∩ kill(kⱼ)| would cut the tail far harder, and the positive
  correlation measured in THM-1071(III) says those overlaps are large;
  (b) **shrink KB** by improving the measure horn on the failing region — max T = 308.4 is
  driven by a handful of near-consecutive killer quintuples, and a special argument for
  clustered killers would collapse it;
  (c) **quotient by symmetry** — killers enter only through their residues mod q ≤ 40, so
  many sextuples are certificate-identical; deduplicating on the residue vector could cut
  the count by orders of magnitude.
- (a) and (c) are compatible and both look worth more than raw compute.

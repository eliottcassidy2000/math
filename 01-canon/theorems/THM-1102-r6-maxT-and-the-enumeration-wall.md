---
id: THM-1102
title: r=6 bounded-window R telemetry and the historical candidate-box enumeration wall; no uniform max-T is proved, while THM-1135 later gives a different finite reduction
status: PARTIAL / SUPERSEDED TAIL COORDINATE — the width-16 census over all 792 cores and the candidate-box feasibility estimate are retained, but interiority inside a fixed window does not make its max T uniform. THM-1145 gives finite and infinite covering refutations of the sequential max-T extrapolation. THM-1135 independently supplies a harmonic-discrepancy tail and a genuine finite mixed-scale box, but does not close that box. Uniform r=5 and r=6 remain open; only r<=4 is uniformly closed in this clustered hierarchy
source: kind-pasteur-2026-07-18-S128 (cont.64; owner: run the r=6 finite horn, computing max T first)
depends_on:
  - THM-1101    # bounded r=5 telemetry; its former uniform closure is withdrawn
  - THM-1081    # the R-ladder this extends
related: [THM-1097, THM-1145, THM-1135, MISTAKE-164]
script: 04-computation/r6_maxT_kps_S128c64.py, r6_maxT_chunk_kps_S128c64.py, r6_feasibility_kps_S128c64.py (+ .out)
---

# THM-1102 — r=6: max T first, and where the enumeration wall is

> **Audit correction (codex-S73; MISTAKE-164).**  Every quintuple below was
> restricted to a width-16 bottom window.  Having the maximizer at offset 9
> rules out truncation by *that particular window edge*; it does not reduce
> arbitrary larger or independently shifted killers to the window.  Thus
> `308.4` is the scanned maximum and `KB=333` is a candidate cutoff only.
> The exact all-high r=5 gap in THM-1101 shows why this distinction is real.

> **Tail resolution (codex-S67; THM-1145/1135).**  The sequential coordinate
> is genuinely nonuniform: a covering finite row has `T=1043/3>338`, and an
> infinite covering progression has `T/k_5 -> 28/27>1`.  THM-1135 repairs the
> all-scale reduction by keeping the core safe set intact and charging all six
> killers harmonically.  It proves an explicit finite mixed-scale box, not the
> old `KB=333` box and not an `r=6` closure.

## (I) Max T, computed first and window-checked

At r=5 I set the finite-horn bound from max k_removed and it was wrong; the bound is set by
**max T**, since the last killer is certified by the measure horn only once it exceeds T.
So this time T came first.

Over all **792 seven-speed cores**, every quintuple of removed killers in a width-16 window:

> **max R = 1.85794**, at core [1,2,4,7,9,11,12], killers (158,160,162,164,166), T = 308.4
> **scanned max T = 308.4** over the scanned R ≥ 1 region ⟹ **candidate KB = 333**

24,598 + 16,298 quintuples have R ≥ 1, and the largest killer among them is 172. The worst
case sits at **offset 9 inside a window of width 16**.  This verifies that the
reported bank maximum is not an edge artifact inside that bank.  It does not
prove that the infinite failure region is contained in the bank.

Both chunks (cores 0–300 and 300–792) were run to a printed summary. The first attempt at a
single 792-core run was killed mid-way at core 300 with max T still *rising* (129 → 225 →
308), and I did not use those numbers — a partial scan cannot set a bound.

## (II) The R-ladder

> **0.51852 (r=2) → 0.73375 (r=3) → 0.98453 (r=4) → 1.28495 (r=5) → 1.85794 (r=6)**

## (III) The candidate-box enumeration wall

Last session I predicted r=6 was "where the enumeration finally becomes the binding
constraint rather than the mathematics". It is, and now it is measured:

| r | KB | sextuples/quintuples/… passing the prune | runtime |
|---|---|---|---|
| 4 | 400 | 1.43 × 10⁸ | ~25 min |
| 5 | 235 | 2.64 × 10⁸ | ~9 min |
| **6** | **333** | **≈ 3.64 × 10¹²** | **≈ 140 days** |

**13,783× the r=5 count.** Within this candidate box the prune has stopped working, and the reason is structural
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

## (V) Status after the all-scale audit

- r = 2, 3, 4: **uniformly closed**.
- r = 5: **open**.  The below-235 horn is finite-exact, but THM-1101 has an
  exact covering row above 235 missed by both sides of its former split.
- r = 6: **open but uniformly finite-reduced by THM-1135**.  The width-16
  region remains valid telemetry, while THM-1145 proves that neither its
  numerical `max T` nor the sequential ratio horn is global.  THM-1135 gives
  the valid replacement box
  `k1<=513, k2<=950, k3<=19000, k4<=313500, k5<=4514400, k6<=58687200`.

## Named next
- First close the four-removal `r=5` all-scale bridge; the same endpoint-owner
  overlap/self-similarity information is prerequisite for any honest `r=6`
  finite reduction.
- Inside THM-1135's finite box, r=6 needs a better certificate than raw enumeration. Three reusable candidates are:
  (a) **strengthen the prune** — Σ frac ≥ 1 is weak because it ignores overlap; a bound
  using pairwise |kill(kᵢ) ∩ kill(kⱼ)| would cut the tail far harder, and the positive
  correlation measured in THM-1071(III) says those overlaps are large;
  (b) **shrink KB** by improving the measure horn on the failing region — max T = 308.4 is
  driven within the scanned bank by a handful of near-consecutive killer quintuples, and a special argument for
  clustered killers would collapse it;
  (c) **quotient by symmetry** — killers enter only through their residues mod q ≤ 40, so
  many sextuples are certificate-identical; deduplicating on the residue vector could cut
  the count by orders of magnitude.
- (a) and (c) are compatible and both look worth more than raw compute.

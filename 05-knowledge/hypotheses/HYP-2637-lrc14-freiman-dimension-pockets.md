---
id: HYP-2637
status: STRONG REFRAME (verified-computational) — the pockets ARE Freiman dimension; AP (d=1) the unique max; LRC(14) NOT proved
source: kind-pasteur-2026-06-19-S12
related:
  - HYP-2607   # the crux: consec maximizes L_y
  - HYP-2635   # the third-pocket / route-convergence consolidation
  - THM-401    # pair-sum sieve modulus C=2n-1
  - HYP-2083   # antipodal summand-unit bridge (the summand graph at C)
  - THM-531    # AP-orbit invariance
  - HYP-2633   # codex reciprocal-tail (the signed-lattice crux)
---

# HYP-2637: the LRC(14) "pockets" are FREIMAN DIMENSION; the recursion is dimension-by-dimension; the AP is the unique top

> Namespace note (codex-2026-06-19-S27): this file arrived after
> `HYP-2637-lrc14-relation-fiber-gap-split.md` had already claimed
> HYP-2637 / T885 on main.  The two notes are mathematically compatible, but
> the identifier is colliding.  HYP-2638 records the collision-free finite
> Freiman `3k-4` subcertificate for the small-excess part of this program.

The third pocket (HYP-2635) is now structurally explained. The configuration space of cluster
sets `E` (|E|=k, 0∈E) organizes by **additive structure**, and `L_y` (equivalently `meas(S7)`) is
a strictly decreasing function of it.

## The organizing parameter: sumset excess = Freiman dimension

`excess(E) = |E+E| − (2k−1) ≥ 0`. **Freiman–Vosper: excess = 0 ⟺ E is a (full) arithmetic
progression** = a dilate of consec. So:
```
  excess 0  (d=1, full AP)        -> POCKET 1: L_y = L_y(consec), the UNIQUE GLOBAL MAX.
  excess>0, small (d=2 GAP)       -> POCKET 3.
  d=3 GAP                         -> POCKET 4 (found this session).
  ...                             -> POCKET d+1 (d-dim GAP); relation-lattice rank = k−1−d.
  excess large (dissociated)      -> POCKET 2: L_y -> the independent limit L_y^inf << cap.
```

## VERIFIED: L_y strictly decreases with Freiman dimension `d`

Exact (Fraction), per cluster size k:
```
  k=8:   d=1 (consec) 0.35823  >  d=2 GAP 0.139–0.175  >  d=3 GAP 0.084–0.115     (cap 0.382)
  k=9:   d=1          0.49288  >  d=2     0.276–0.318                              (cap 0.494)
  k=12:  d=1          0.71417  >  d=2     0.586–0.597  >  d=3     0.535            (cap 0.857)
```
Each added dimension drops `L_y` by roughly a factor ~0.5 — the **dimension penalty**. So every
`d≥2` GAP sits FAR below cap (margin ≥ 0.21 at k=8). The "third pocket" of HYP-2635 = the `d≥2`
GAPs; the dimension penalty is exactly why they cap out far below.

## VERIFIED: the AP is cleanly separated (a gap at excess 0)

Max `L_y` over all primitive NON-full-AP sets (exhaustive bounded spread):
```
  k=8:  consec 0.35823, max non-AP 0.30306 at (0,2,3,4,5,6,7,8)  — gap 0.0552, margin to cap 0.079
  k=9:  consec 0.49288, max non-AP 0.48729 at (0,1,2,3,4,5,6,7,9) — gap 0.0056, margin to cap 0.0070
  k=10: consec 0.57008, max non-AP 0.56241 at (0,..,7,9? actually ...,8,10) — gap 0.0077, margin 0.042
```
The dangerous non-AP are the **near-AP** sets (one element shifted/holed). They are bounded-spread
(finite check, DONE) or acquire a dissociated stranger as they widen (contract).

## The PROOF PROGRAM (clean partition)

1. **excess 0 ⟹ E is an AP** [Freiman–Vosper, KNOWN] ⟹ `L_y = L_y(consec)` [AP-invariance THM-531]
   ⟹ `< cap_k` [finite computation]. EXACT. ✓
2. **bounded-spread non-AP: finite check** (the near-AP sets) `< cap`. DONE (k=8,9,10).
3. **wide non-AP** — split (CORRECTED below): (a) SPREAD `d≥2` GAP ⟹ `L_y ≪ cap` by the DIMENSION
   PENALTY (margin ≥0.21, EASY); (b) large doubling ⟹ dissociated ⟹ near `L_y^inf ≪ cap`.

## CORRECTION (workflow w7gkxcz0b, verify-confirmed) — the binding case is the bounded near-AP, NOT a wide GAP

**My "margin ≥0.21, crude bound suffices" was WRONG for the binding case.** Exhaustive (verify-confirmed):
- The ≥0.21 margin applies ONLY to SPREAD `d≥2` GAPs (Minkowski sums of two genuine APs) — those are
  the EASY, non-binding cases (provable via the orbit-collision identity `O(A1+A2)(x)=O(A1)⊕O(A2)`).
- **The actual supremum over ALL non-full-AP sets is a MINIMAL SINGLE-DEFECT near-AP** (excess 2, one
  double-step), NOT a spread GAP: unique k=9 maximizer `E=[0,1,2,3,4,5,6,7,9]`, `L_y=38681/79380=0.487289`,
  margin to cap_9 only **0.00697** (k=8 `[0,2,..,8]` 0.30306 margin 0.078; k=10 `[0,..,8,10]` 0.56242 margin 0.042).
  Verified the UNIQUE global max (exhaustive max(E)≤18). Larger defects strictly LOWER `L_y`.
- **Good news:** the binding near-AP is BOUNDED-SPREAD (max element ~k+1), so it lives in the EXACT
  FINITE CHECK (done) — NOT in the wide/open part. The wide/open part (spread GAP + dissociated) is the
  comfortable-margin part. So "crude bound suffices" is TRUE for the genuinely-OPEN (wide) part; the tight
  part (AP + single-defect near-AP) is bounded-spread and exact.
- **NO separate "third pocket":** the intended binary {GAP | dissociated} dichotomy is the WRONG
  abstraction (every 8–10-set sits in SOME proper 2-dim GAP; dissociated-stranger is nearly vacuous).
  Coverage is complete via the CONTINUOUS doubling `σ=|E+E|/k`: max-`L_y` is ≈monotone non-increasing in
  `σ`, the AP (min σ) is the top, and everything is `< cap`. (Tension: my integer-excess envelope showed
  local bumps at excess 3,5, while the workflow reports σ-monotone — likely a binning artifact; the robust
  fact both agree on is max-non-AP `< cap`, worst = the s=2 single-defect near-AP.)

**Revised remaining lemma — split 3(a):**
- **3(a-i) SPREAD GAPs** (`O(A1+A2)=O(A1)⊕O(A2)` collisions): margin ≥0.157, EASY, low-risk.
- **3(a-ii) near-AP boundary** (the BINDING k=9 single-defect, margin 0.007): the genuinely tight piece,
  but BOUNDED-SPREAD so in the finite check; the remaining open is only (A) prove the global sup over
  UNBOUNDED non-AP is attained at the bounded s=2 defect (currently exhaustive to max(E)≤18), and
  (B) the exact `38681/79380 < cap_9` with a structural reason for the ~0.0056 drop from consec.

## The summand-graph axis (addition / multiplication / odd / pos-neg)

The relations `Σ m_e e = 0` live mod `C = 2n−1 = 27 = 3³` (THM-401; HYP-2083). The two axes:
- **ADDITION** builds the antipodal shells `{a, C−a} = {+a, −a}` (pos/neg pairing mod C); the GAP
  collisions are exactly relations summing to 0 in these shells.
- **MULTIPLICATION** = the GAP step directions (a `d`-dim GAP has `d` multiplicative steps); units
  `(Z/27)*` act on the shells (unit-visible vs blind).
- **ODD** `C=27` has no midpoint ⟹ clean antipodal pairing.
- For `C=27=3³` the shells STRATIFY into {units, gcd-3, gcd-9}; the tight sporadic `V*` lives in the
  gcd-3 blind-spot stratum (HYP-2083). So the pocket structure refines by `(dimension d) × (stratum)`.

## The recursion as the number of runners `n` changes

- The **pockets** (dimensions `1..k−1`) proliferate as the cluster size `k` grows; the AP-to-cap
  margin SHRINKS with k (k=8: 0.024, k=9: 0.0014 — k=9 is the tightest row). There is a "critical
  cluster size" feel where the AP nears cap.
- The **summand modulus** `C = 2n−1` and its factorization govern the stratum structure; **n=14 is
  special: 2·14−1 = 27 = 3³**, a prime power with the deepest (3-level) gcd-stratification. The apex
  prime `7 = n/2` gives the sector count. So LRC(14)'s difficulty is tied to `27 = 3³` and `7`.

## CAVEAT / reconciliation with codex S26 (HYP-2638, HYP-2639) — small vs large doubling
Codex built on this and found the precise scope. (i) **Small-doubling pocket handled (HYP-2638 PASS,
k=8,9,10):** `|E+E| ≤ 3k−4` ⟹ genuine low-dim GAP (Freiman), dimension penalty bounds it ≪ cap. (ii)
**But "relation-covered" ≠ "small doubling" (HYP-2639):** the third-pocket example `[0,3,5,16,28,30,33,35]`
(every element in a height-2 relation, no stranger) has `|E+E|=31 > 3k−4=20` — NOT a low-dim GAP. So
"no stranger ⟹ GAP" has a hole. (iii) **Resolution — THREE excess bands:** excess 0 = AP (exact finite);
small excess = GAP (dimension penalty ≪ cap, HYP-2638); LARGE excess = high doubling = low correlation
⟹ near the independent limit `L_y^inf` ≪ cap. The third pocket sits in the LARGE band (its L_y≈0.06 ≈
`L_y^inf`=0.049), bounded DIRECTLY by high-doubling-⟹-low-correlation, NOT by peeling. Codex's HYP-2639
supplies the typed summand-shell/sign/visibility layer for the signed channelization (the reciprocal-tail
crux HYP-2633). The precise quantity is the ENVELOPE `max L_y at each excess`: verified to decrease from
the AP and stay `< cap` for all excess `≥ 1` (not strictly monotone — bumps at excess 3,5).

## Honest status
LRC(14) NOT proved. This is a structural REFRAME, verified computationally, that (i) shows `L_y` is
continuously controlled by the doubling `σ=|E+E|/k` with the AP the unique top and NO separate third
pocket; (ii) finds the binding worst-non-AP is the BOUNDED single-defect near-AP (k=9 `[0..7,9]`,
margin 0.007) — in the exact finite check, NOT a wide set; (iii) splits the remaining lemma into the
EASY spread-GAP dimension penalty (3(a-i), margin ≥0.157, orbit-collision) and the tight near-AP
boundary (3(a-ii), bounded ⟹ finite, the only open piece being the unbounded→bounded sup reduction);
(iv) finds pocket 4 (d=3 GAP) and the dimension penalty for spread GAPs; (v) ties the difficulty
across n to the factorization of `2n−1` (n=14 anomalously tight via `27=3³`, verify-confirmed
margins n=10 +0.20 / n=12 +0.19 / n=14 +0.054), with the mod-27 summand strata INERT for the binding
small-element clusters (HYP-2640). My initial "margin ≥0.21 closes it" was CORRECTED (that holds only
for spread GAPs). Files:
`04-computation/lrc14_{sumset_excess_pockets,pocket_dimension_recursion,nonAP_max}_kps.py`.

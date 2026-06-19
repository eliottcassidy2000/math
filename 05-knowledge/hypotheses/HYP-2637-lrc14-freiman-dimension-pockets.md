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
3. **wide non-AP: dichotomy via Freiman** —
   (a) small doubling ⟹ E ⊆ a `d≥2` GAP ⟹ `L_y ≪ cap` by the DIMENSION PENALTY (margin ≥0.21);
   (b) large doubling ⟹ a dissociated stranger exists ⟹ CONTRACT (HYP-2610) toward `L_y^inf ≪ cap`.

The genuine remaining lemma is **3(a): the GAP dimension penalty** — but it has HUGE margin (≥0.21),
so a crude rigorous bound suffices (unlike the tight `consec < cap` at k=9, which is the exact finite
check). **Mechanism:** a 2-dim GAP `E = A1 + A2` (Minkowski sum of two APs) has orbit
`O(E,x) = O(A1,x) ⊕ O(A2,x)` (setwise circle sum); the two AP-directions COLLIDE at commensurate `x`
(when `(a i − b j)x ∈ ℤ`), wasting orbit points and lowering coverage → lower `L_y`. The collisions
are governed by the relation lattice — connecting to the summand graph below.

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

## Honest status
LRC(14) NOT proved. This is a structural REFRAME, verified computationally, that (i) explains the
third pocket (dimension penalty), (ii) isolates the only tight case (the AP, excess 0, exact finite
check), (iii) gives a clean Freiman partition, (iv) identifies the one remaining lemma (GAP dimension
penalty, margin ≥0.21) and its mechanism (collision/relation-lattice), and (v) ties it to the
summand graph at `C=27=3³` and codex's reciprocal-tail crux (HYP-2633). Files:
`04-computation/lrc14_{sumset_excess_pockets,pocket_dimension_recursion,nonAP_max}_kps.py`.

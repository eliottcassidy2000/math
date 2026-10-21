---
id: HYP-9350
title: "Polygonal-order growth and the rotation-constrained colouring bridge; the recalled fractional-chromatic equality is refuted"
status: OPEN growth and quantitative questions; PROVED unbounded finite family in corrected THM-4737; REFUTED recalled fractional-chromatic equality by D={1,3,4,7}
source: mac-mini-2026-10-20-kepler
depends_on:
  - 01-canon/theorems/THM-4737-kepler-polygon-reformulation-of-the-lonely-runner-conjecture-observer-angular-momentum-vertex-times-and-the-relation-torus.md
related:
  - 01-canon/theorems/THM-639-hamiltonian-path-frame-for-runner-families.md
  - 01-canon/theorems/THM-4566-hilbert-class-field-planes-are-the-idoneal-planes.md
  - 05-knowledge/reference/CORE-PAPERS.md
scripts:
  - 04-computation/experiments/kepler_lrc_20261020/polygonal_order_bound.py
  - 04-computation/experiments/kepler_lrc_20261020/pair_one_m_growth.py
---

# HYP-9350 -- polygonal order, Kepler measure, and the colouring bridge

**Correction, 2026-10-20:** finite polygonal order is conditional on a
lonely witness; strict area-rate positivity loses boundary witnesses;
and the fractional-chromatic equality recalled in(c) is false.
[Exact counterexample and repaired bridge](../results/runner_area_clock_20261020j/boundary_and_colouring.md).

## (a) The polygonal order is unbounded but slow: its growth law

`N_min(D)` = least `N` such that some vertex time `j/N` of the regular `N`-gon is lonely, or infinity if none exists
(THM-4737 (K3)). A provisional guess that `N_min <= 12` for two runners (true for speeds
`<= 60`) is **refuted by the session's own larger census**: over coprime pairs with speeds
`<= 200` the maximum is `15` (pairs `(1,72)`, `(144,163)`, `(144,187)`), and on the pairs
`(1, m)` the record values are

```text
m      : 2  6  12  18  72  3960  8640  16560
N_min  : 3  5   8  11  15    16    17     19        (m <= 30,000; pair_one_m_growth.py)
```

with record `m` highly composite (`72 = 2^3 3^2`, `3960 = 2^3 3^2 5 11`, `8640 = 2^6 3^3 5`,
`16560 = 2^4 3^2 5 23`). For `(1, m)` the `N`-gon works iff some `j` with `||j/N|| >= 1/3` has
`||mj/N|| >= 1/3`; every `N` dividing `m` fails, and `N | m - 1` or `N | m + 1` with `N >= 3`
succeeds. Thus an odd prime factor of `m^2 - 1`, when one exists, supplies
an upper bound; no such prime exists at m=3. Divisibility by many small
denominators is one explicit obstruction, not a proved characterization.
The resulting comparison with Jacobsthal-type quantities is a research
direction. **Claims.**
(a1) **PROVED in corrected THM-4737:** finite N_min is unbounded for every
fixed k>=2, using D={1,m,2m,...,(k-1)m}, m=lcm(2,...,L), and the explicit
lonely time1/2+1/(km). No extra prime-factor condition is needed.
(a2) **OPEN:** `N_min(D) = O(log max D)` for fixed `k`,
with the `(1, m)` family extremal up to constants. (a3) For the tight sets `N_min = n` exactly
(observed on all six tight sets found). Larger `k`: maxima `18` (`k = 3`, speeds `<= 45`, at
`(19,24,43)`), `17` (`k = 4`, speeds `<= 24`, at `(8,9,10,13)`). The growth law is the
quantitative form of "which polygon catches a lonely moment"; a proof of (a2) would bound
the modulus a divisibility sieve must reach, which is exactly the resource the
computer-assisted proofs for `k <= 12` consume.
If asserted for every D, a finite upper bound in(a2) would itself prove
integer LRC. A bound only on sets already known to have finite N_min is
a different, conditional statement.

## (b) The Kepler measure of lonely times

For integer speeds, `m(D)` is the measure in[0,1) of the **closed** lonely
set. It equals the measure of its strict interior, but those sets can
differ at every actual witness in a tight example. It is positive iff
kappa(D)>1/n, and zero when kappa(D)<=1/n; equality is the tight case and
strict inequality would be a counterexample. For Q-independent real
speeds the **long-time density**, not generally the measure in[0,1), is
`(1 - 2/n)^k`. **Questions.** The minimum of `m(D)` over
non-tight integer sets of size `k` (the "barely lonely" sets of Kravitz give candidates), and
whether `m(D)` is bounded below by a function of the relation rank for real speeds (the
zonotope `M B` of THM-4737 (K4) has volume `(1-2/n)^k |det|`-type lower bounds by Minkowski
when `rank L(D)` is small).

## (c) The colouring bridge toward Hadwiger-Nelson

**REFUTED recalled identity.** The original memory-only attribution to
Chang-Liu-Zhu asserted chi_f(G(Z,D))=1/kappa(D). For D={1,3,4,7}, the graph
has a K4 on{0,3,4,7} and a proper periodic four-colouring
c(z)=2floor((z mod8)/4)+(z mod2), so chi_f=chi_c=chi=4. But kappa=1/5,
by the exact five-interval proof in the linked correction. Thus the
asserted equality is false; no claim is made that the named literature
actually states it.

The survivor is chi_f<=chi_c<=1/kappa. The last quantity is exactly the
optimum when colourings are restricted to the additive rotation form
z -> r{zt}. Arbitrary colourings can lose this common-clock constraint.
Indeed chi(G(Z,D))<=|D|+1 follows already from finite-subgraph degeneracy
and compactness, so that unrestricted bound is not equivalent to LRC.
The new question is whether a restricted, source-preserving conversion
to rotation colourings can meet the LRC bound; the example shows that
optimal unrestricted colourings need not convert without loss.

On the regular `N`-gon the vertex-time criterion (V) is a
statement about the circulant `C_N(D)` (chords of the polygon in the lengths `D`): the
observer is lonely at `j/N` iff the dilation by `j` maps every chord length into the middle
window `[N/n, N - N/n]`, i.e. iff `j D` avoids the "short-chord" set; so a lonely vertex time
is a multiplier `j` sending the chord graph into the long-chord circulant, whose independence
number controls `chi_f(C_N(j D)) = N / alpha`. **Question.** Does the idoneal-plane / residue
colouring toolkit of the HN lane (THM-4558, THM-4566: colourings of field planes by residue
characters) give **rotation** colourings with controlled N? Arbitrary
circular colourings are insufficient for this equivalence. Residue and
character partitions from HYP-9349(c) must retain their literal modulus,
quotient, and chosen multiplier; whether they yield such certificates
remains untested.

## (d) What a polygon-based general proof would need

THM-4737 shows the `n`-gon alone settles every speed set with no multiple of `n`, and the
Kepler/Haar reduction settles relation rank `<= 1`. A general proof along these lines must
handle speeds divisible by `n` through larger polygons: for integer speeds this is the
classical view-obstruction problem in the quotient, and the tight families show the lonely
set can be finite, so no measure argument (Kepler averaging) can close it; the needed input is
arithmetic (which multipliers `j` and moduli `N` exist), i.e. (a). No claim beyond `k <= 12`
is made.

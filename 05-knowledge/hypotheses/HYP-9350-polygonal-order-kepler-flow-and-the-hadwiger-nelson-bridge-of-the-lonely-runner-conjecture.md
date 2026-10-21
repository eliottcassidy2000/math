---
id: HYP-9350
title: "After THM-4737: the polygonal order N_min(D) is unbounded but grows only logarithmically (records 3, 5, 8, 11, 15, 16, 17, 19 on the pairs (1, m) up to m = 30000); the Kepler-flow measure of lonely times, the distance-graph colouring bridge to Hadwiger-Nelson (kappa(D) = 1 / fractional chromatic number of G(Z,D), CITED from memory), and what a polygon-based general attack would have to prove"
status: OPEN (pre-audit); finite evidence and first tests stated per item; one citation flagged as recalled, not verified
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

## (a) The polygonal order is unbounded but slow: its growth law

`N_min(D)` = least `N` such that some vertex time `j/N` of the regular `N`-gon is lonely
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
succeeds, so `N_min(1,m)` is at most the least prime factor `>= 3` of `m^2 - 1` and is driven
up only by `m` divisible by many small `N` -- a Jacobsthal-type quantity. **Claims.**
(a1) `N_min` is unbounded for every `k >= 2` (through `(1, m)` with `m` divisible by all
`N <= L` and `m^2 - 1` free of small primes). (a2) `N_min(D) = O(log max D)` for fixed `k`,
with the `(1, m)` family extremal up to constants. (a3) For the tight sets `N_min = n` exactly
(observed on all six tight sets found). Larger `k`: maxima `18` (`k = 3`, speeds `<= 45`, at
`(19,24,43)`), `17` (`k = 4`, speeds `<= 24`, at `(8,9,10,13)`). The growth law is the
quantitative form of "which polygon catches a lonely moment"; a proof of (a2) would bound
the modulus a divisibility sieve must reach, which is exactly the resource the
computer-assisted proofs for `k <= 12` consume.

## (b) The Kepler measure of lonely times

`m(D)` = measure of `{t : all ||v_i t|| >= 1/n}` (the common strict-increase set of the
observer area clocks). It is `(1 - 2/n)^k` for Q-independent real speeds, `0` exactly on
tight integer sets, and positive otherwise. **Questions.** The minimum of `m(D)` over
non-tight integer sets of size `k` (the "barely lonely" sets of Kravitz give candidates), and
whether `m(D)` is bounded below by a function of the relation rank for real speeds (the
zonotope `M B` of THM-4737 (K4) has volume `(1-2/n)^k |det|`-type lower bounds by Minkowski
when `rank L(D)` is small).

## (c) The colouring bridge toward Hadwiger-Nelson

**Recalled (CITED from memory; verify before use):** for a finite set `D` of positive integers
the fractional chromatic number of the distance graph `G(Z, D)` (vertices `Z`, edges at
differences in `D`) is `chi_f(G(Z,D)) = 1/kappa(D)` (Chang-Liu-Zhu, *Distance graphs and
T-coloring*, JCTB 1999; used by Zhu and Liu-Zhu for circular chromatic numbers of distance
graphs). If so, the lonely runner conjecture is exactly

```text
chi_f(G(Z, D)) <= |D| + 1   for every finite D,
```

a fractional-colouring statement about one-dimensional distance graphs, the integer cousins
of the plane's unit-distance graph. On the regular `N`-gon the vertex-time criterion (V) is a
statement about the circulant `C_N(D)` (chords of the polygon in the lengths `D`): the
observer is lonely at `j/N` iff the dilation by `j` maps every chord length into the middle
window `[N/n, N - N/n]`, i.e. iff `j D` avoids the "short-chord" set; so a lonely vertex time
is a multiplier `j` sending the chord graph into the long-chord circulant, whose independence
number controls `chi_f(C_N(j D)) = N / alpha`. **Question.** Does the idoneal-plane / residue
colouring toolkit of the HN lane (THM-4558, THM-4566: colourings of field planes by residue
characters) give circular colourings of `C_N(D)` with `N` an idoneal modulus, i.e. lonely
vertex times with controlled `N`? The eight one-class genera of `-420` partition the units
mod `420` into six-element register classes (HYP-9349 (c)); whether such a partition can
serve as a colouring certificate for `C_420(D)` is untested.

## (d) What a polygon-based general proof would need

THM-4737 shows the `n`-gon alone settles every speed set with no multiple of `n`, and the
Kepler/Haar reduction settles relation rank `<= 1`. A general proof along these lines must
handle speeds divisible by `n` through larger polygons: for integer speeds this is the
classical view-obstruction problem in the quotient, and the tight families show the lonely
set can be finite, so no measure argument (Kepler averaging) can close it; the needed input is
arithmetic (which multipliers `j` and moduli `N` exist), i.e. (a). No claim beyond `k <= 12`
is made.

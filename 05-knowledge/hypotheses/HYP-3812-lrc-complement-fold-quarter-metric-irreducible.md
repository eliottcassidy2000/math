---
id: HYP-3812
title: THE LRC COMPLEMENT-FOLD + QUARTER-FOLD -- the covering-min obstruction is METRIC-IRREDUCIBLE at the composite Phi6 (does NOT dissolve, unlike the tournament SC-cover), which EXPLAINS why LRC is harder + confirms THM-503 (no Euler product) from the covering side. Applying the tournament fold-thinking (S78: cube folds by <sigma=complement, flip>=Klein-four; an obstruction can DISSOLVE in fold-adapted coords) to LRC-14. (1) COMPLEMENT-FOLD: the antipode iota:t->1-t (S55) is the LRC complement; M(S)=max_t min_v||vt||, the danger pattern, and the lonely measure are ALL iota-invariant (VERIFIED G(t)=G(1-t)), so LRC folds to the half-circle [0,1/2] (the 2 binding atoms {t*,1-t*} fold to 1). This HALF-fold is genuinely useful: lonely measure iota-even => Verblunsky REAL (S66); the parity lemma (S55, odd D => #lonely even) is the iota fixed-point count. (2) QUARTER-FOLD GROUP: the phase space is Z/Phi6, Phi6=3*61 (Eisenstein primes); the complement -1 mod Phi6 = 182 CRT-FACTORS: 182=(-1 mod3, -1 mod61)=62*121, the PRODUCT of the two PARTIAL complements; <62,121>={1,62,121,182}=Klein-four in (Z/Phi6)* folds Z/Phi6 into a QUARTER exactly as <sigma,flip> folded the tournament cube. (3) THE TEST (does the obstruction dissolve, S78?): NO. M(construction) restricted to modulus q is 0 at q=2,3,6,7,14 (COVERED), 4/61 at q=61 (shallow), and =M_C=n/Phi6 ONLY at the composite q=Phi6=183. So the construction COVERS the prime moduli 3,61 and binds ONLY at the composite Phi6 => the covering-min is METRIC-IRREDUCIBLE at Phi6, it does NOT factor over the CRT quarter, so the obstruction does NOT dissolve (unlike the tournament SC-cover which was a coordinate artifact, HYP-3811). This is the analytic non-factorization of THM-503 (L not an Euler product) seen from the covering side: the covering-min lives at the DEEP composite modulus Phi6 (the 2nd continued-fraction convergent [0;n-1,n]), invisible to the prime CRT factors. ABSTRACT LESSON: fold-thinking DISTINGUISHES coordinate-artifact obstructions (dissolve => easy in the right basis, the tournament) from metric-irreducible ones (don't dissolve => genuinely hard, the LRC); the LRC covering-min is the latter, so its proof must engage the FULL Phi6 (the S73 Chebyshev 2-point dual on Z/Phi6), with no prime/CRT shortcut
status: MIXED (exact computation + established-structure synthesis). VERIFIED (n=14, exact Fraction): iota-symmetry G(t)=G(1-t); the Klein-four {1,62,121,182}=(-1 CRT-partials) in (Z/183)*; M(construction) restricted to modulus q = 0 (q=2,3,6,7,14), 4/61 (q=61), 14/183=M_C (q=183) -- the construction binds ONLY at the composite Phi6. HONEST: the complement-fold's utility (halving, real Verblunsky, parity lemma) is established (S55/S66); the quarter-fold Klein-four is an exact group fact; the 'does not dissolve' conclusion is the correct DIS-analogy to HYP-3811 (metric-irreducible vs coordinate-artifact) and RE-DERIVES THM-503's no-Euler-product from the covering side -- a structural explanation of the difficulty, NOT a new bound. Only n=14 computed (the identity M_C=n/Phi6 with binding at Phi6 is general via S68).
source: klein-2026-07-01-S79
depends_on:
  - HYP-3811   # tournament quarter-fold + the coordinate-artifact dissolution (this is the LRC analogue + the DIS-analogy)
  - HYP-3800   # phase-residue p(w)=nw mod Phi6, the binding modulus q*=Phi6 (S68)
related:
  - HYP-3811-opus # opus-S20 CORRECTS the quarter-tiling framing: 'the fold is a SINGLE terminal step (NO quarter-tiling)'. CONSISTENT with this: for LRC the quarter-fold (CRT Klein-four) is a GROUP symmetry that does NOT reduce the problem (obstruction metric-irreducible) -- so the complement-fold is the real (single) reduction, the 'quarter' is not a problem-reducing fold. Both agree the quarter is not a genuine terminal reduction.
  - THM-503    # L is NOT an Euler product -- re-derived here from the covering side (binding at composite Phi6, not primes)
  - HYP-3806   # Chebyshev 2-point dual on Z/Phi6 (the covering-min proof must engage the full Phi6)
  - THM-515    # lonely measure = singular series; complement iota-symmetry
  - HYP-3715   # t*=n/Phi6 hexagonal binding; Eisenstein primes of Phi6
external: Klein four-group; CRT; Eisenstein primes (=1 mod 6); continued fraction [0;n-1,n]
results:
  - 04-computation/lrc_complement_fold_quarter_klein.py
  - 05-knowledge/results/lrc_complement_fold_quarter_klein.out
---

# HYP-3812 — LRC complement-fold + quarter-fold; the covering-min is metric-irreducible at Phi6

## (1) The complement-fold (useful)
The antipode `iota: t -> 1-t` (S55) is the LRC complement. `M(S) = max_t min_v ||vt||`, the danger pattern
`D_v`, and the lonely measure are all `iota`-invariant (`||v(1-t)|| = ||vt||`, VERIFIED `G(t)=G(1-t)`), so
the whole problem folds to the half-circle `[0,1/2]`; the two binding atoms `{t*, 1-t*}` fold to one. This
HALF-fold is genuinely useful: the lonely measure is `iota`-even, so its Verblunsky coefficients are REAL
(S66); and the parity lemma (S55: odd `D` => `#lonely` even) is exactly the `iota` fixed-point count. This
is the direct analogue of the tournament complement-fold (`sigma`).

## (2) The quarter-fold group (Eisenstein CRT)
The phase space is `Z/Phi6`, `Phi6 = 3*61` (Eisenstein primes: `3` ramified, `61 = 1 mod 6`). The complement
`-1 mod Phi6 = 182` **CRT-factors**: `182 = (-1 mod 3, -1 mod 61) = 62 * 121`, the product of the two
**partial complements**. So `<62, 121> = {1, 62, 121, 182}` is a **Klein-four** subgroup of `(Z/Phi6)*`, and
it folds `Z/Phi6` into a **QUARTER** exactly as `<sigma, flip>` folded the tournament cube. The antipode
`iota` is the diagonal (`= 62*121`); the two partial complements are the two axes of the fold.

## (3) The test: does the obstruction dissolve? NO -- it is metric-irreducible at Phi6
S78/HYP-3811's lesson: an obstruction may DISSOLVE in the fold-adapted coordinates (a coordinate artifact).
Testing the LRC covering-min: `M(construction)` restricted to each modulus `q` (the loneliness available at
`q`):
| `q` | 2 | 3 | 6 | 7 | 14 | 61 | **183 = Phi6** |
|---|---|---|---|---|---|---|---|
| `M|_q` | 0 | 0 | 0 | 0 | 0 | 4/61 | **14/183 = M_C** |
> The construction **COVERS** the prime moduli `3, 61` (loneliness `0` and shallow `4/61 < M_C`) and binds
> **only at the composite `Phi6 = 3*61`**. So the covering-min is **METRIC-IRREDUCIBLE at `Phi6`**: it does
> NOT factor over the CRT quarter, and the obstruction does **NOT dissolve** under the quarter-fold.
This is the DIS-analogy to the tournament: there the SC-cover excess was a coordinate artifact (dissolved in
half-address coords); here the covering-min excess is genuinely metric, living at the deep composite `Phi6`.

## The abstract lesson (why this is useful for LRC)
Fold-thinking **distinguishes two kinds of obstruction**:
- **Coordinate-artifact** (tournament SC-cover): dissolves in the fold-adapted basis => the problem is easy
  in the right coordinates.
- **Metric-irreducible** (LRC covering-min): does NOT dissolve, because it lives at the deep composite
  modulus `Phi6` (the 2nd continued-fraction convergent `[0; n-1, n]`), invisible to the prime CRT factors.
This is the covering-side re-derivation of **THM-503** (`L` is not an Euler product): the covering-min binds
at `Phi6`, not its prime factors, so there is **no CRT/prime shortcut**. The covering-min proof must engage
the FULL `Phi6` — precisely the S73/HYP-3806 Chebyshev **2-point dual on `Z/Phi6`** (the alternation `{1,
killer}` at `{t*, 1-t*}`). The fold-thinking does not give a new bound, but it **classifies the crux**: the
half-fold (`iota`) is a genuine reduction (halving, real Verblunsky, parity), while the quarter-fold (CRT)
is a group symmetry that cannot reduce the metric-irreducible core — telling us exactly where the difficulty
is irreducible and where the proof must live.

## Net
Applying the tournament fold-thinking to LRC: the complement-fold (`iota`) halves the problem usefully; the
quarter-fold is the Klein-four `<62,121>` of CRT partial complements over `Phi6`'s Eisenstein primes; but
the covering-min obstruction is **metric-irreducible at the composite `Phi6`** (the construction covers the
prime moduli and binds only at `Phi6`), so it does NOT dissolve — the productive DIS-analogy to the
tournament. This explains why LRC is harder (analytic non-factorization, THM-503) and localizes the proof to
the full-`Phi6` Chebyshev 2-point dual, with no prime shortcut.

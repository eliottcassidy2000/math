---
id: THM-1960
title: "TOURNAMENTS COMPOSE FROM REGULAR SEEDS -- the spectral substitution law, the seed census, and the octonion object C3[C3]. Reading the owner's three recursion modes (full A+B+C-D-E-F+G, even-half A+B-C, odd-half A+B-C+D-E-F+G; THM-442/553/801/830) as the inclusion-exclusion of a RECURSIVE COMPOSITION: every tournament is a substitution tree T[S_1..S_k] over PRIME (indecomposable) SEEDS (Gallai/modular decomposition). (1) THE SEED CENSUS: # prime (no nontrivial module) tournaments = 1,1,1,0,3,15 for n=1..6 -- C3 is the smallest nontrivial seed, the transitive tournament is the fully-LINEAR tree (no prime seeds), and STRIKINGLY n=4 has ZERO seeds (every 4-tournament is decomposable). The linear special case (transitive quotient) is exactly the SCC/order-join composition (kps THM-1860 H=prod H(SCC), boxeph THM-1855). (2) THE SPECTRAL SUBSTITUTION LAW: for the uniform blow-up T[S^m] (each of T's k vertices replaced by a copy of S, |S|=m), the skew matrix is the Kronecker sum S_T (x) J_m + I_k (x) S_S; when the SEED S is REGULAR (skew row-sums 0, so the block-mean vector 1 is in ker S_S), the spectrum SPLITS EXACTLY: nonzero-spec(T[S^m]) = [nonzero-spec(S), each with multiplicity k] UNION [m * nonzero-spec(T)] (the seed spectrum k times, plus the block-tournament spectrum scaled by the block size m). VERIFIED for regular seeds C3, C5 (and any block T); FAILS for irregular seeds (transitive T3) where the block-mean mixes into the perp space. So char_S, var(lambda^2), SC4, and every even moment tr(S^{2j}) of a regular-seed substitution object are CLOSED FORMS in the seed moments. (3) THE OCTONION OBJECT C3[C3] (n=9 = the octonion Cayley-Dickson level, the substitution-SQUARE of the smallest seed): a REGULAR tournament with char_S = x(x^2+3)^3(x^2+27), spectrum lambda^2 in {0, 3 (x6), 27 (x2)} (= seed sqrt3 mult 3 UNION block m^2*3=27), var(lambda^2)=104 (< transitive 168, > Paley-like 0), SC4=81=3^4, H=3159=3^5*13. Because its spectrum is closed-form, C3[C3] and the family C3[C_{2j+1}] give EXACT degree-8 (tr(S^8)) test objects for the octonion wall (THM-1935/1940)"
status: PROVED/VERIFIED. Seed census exhaustive n<=6 (1,1,1,0,3,15; n=4=0 is the classical 'no prime tournaments on 4 vertices'). Spectral substitution law: the Kronecker-sum structure is exact; the regular-seed splitting verified for C3,C5 (fails for irregular seed, mechanism = block-mean not in ker). C3[C3] invariants computed exactly. The recursion-mode reading is the composition/inclusion-exclusion interpretation.
author: opus-2026-07-20-S444
depends_on: [THM-1862-order-join (order-join = substitution into the TRANSITIVE 2-quotient; strong components = order-join atoms), THM-1860 (H=prod H(SCC) = linear/transitive-quotient case), THM-1830-unstable-non-transitive + THM-1925 (char_A = prod char_A(SCC), spec = disjoint union = the linear-quotient spectral case; my skew regular-seed law is the CYCLIC-quotient generalization), THM-1920 (char_S subleading = #arcs = triangular), THM-1940 (var = SC4-census; octonion wall), THM-1935 (the n=5 quaternion wall), THM-442/THM-553/THM-801/THM-830 (the three recursion modes = Mobius/Legendre/Eisenstein characters), kps THM-1880 (transitive Chebyshev-Pell)]
coordinates_with: THM-1955-which-iso-classes-come-from-smaller (boxeph-S196, CLAIMED/stub, same project via the circulant-character angle -- this THM fills the MODULAR-prime + spectral-substitution part, complementary to its reducible/circulant census; do not override)
distinguishes: the repo's order-join/SCC "atoms" (strongly-connected, THM-1862) from MODULAR-prime "seeds" (Gallai) -- the latter is STRICTLY smaller
cite_by_filename: true
---

# THM-1960 — Tournaments compose from regular seeds

> **OPEN Q1 (cyclic-H) RESOLVED by THM-1975 (opus-S446).** `H` under a cyclic quotient composes via
> the **path-cover polynomial** `pc(S,c)`: `H(C₃[S₁,S₂,S₃]) = Σ K(c₁,c₂,c₃)∏pc(Sᵢ,cᵢ)` — scalar `H`
> is not compositional, but `pc` (its refinement, `H=pc(·,1)`) is. The `13` in `H(C₃[C₃])=3159` is the
> cyclic interleaving of the blocks' path-systems.

Owner: consider the three signed recursion modes and think of tournaments as **recursively
composed of smaller subtournament seeds**; what iso-class seeds correspond to larger tournaments?
The frame is **modular (substitution) decomposition**, and the seeds are the **prime** tournaments.

## 1. The seed census

A **module** of `T` is a subset every outside vertex is uniform to (beats all / loses to all). A
tournament is **prime** (indecomposable) iff its only modules are trivial. Every tournament is a
**substitution tree** `T[S_1,…,S_k]` over primes (Gallai/modular decomposition) — the primes are
the **seeds**.

> **Seed census (# prime tournaments), `n=1..6`: `1, 1, 1, 0, 3, 15`.**

- `C_3` (the 3-cycle) is the **smallest nontrivial seed**; the **transitive** tournament is the
  fully-**linear** tree (no prime seeds).
- **`n=4` has ZERO seeds** — every 4-tournament is decomposable (the classical "no prime tournaments
  on 4 vertices").
- The **linear** special case (transitive quotient) is exactly the **SCC / order-join** composition:
  `THM-1860` (`H = ∏ H(SCC)`), `THM-1862` (order-join). Substitution generalizes it to **cyclic**
  quotients (e.g. `C_3` as the top).

> **CRUCIAL DISTINCTION (flagged in the S444 prior-work map).** The repo's "atoms" (order-join /
> SCC, `THM-1862`) are the **strongly-connected** tournaments (counts `1,1,6,35,353` for `n=3..7`).
> **Modular-prime seeds are STRICTLY stricter**: the unique strongly-connected 4-tournament
> `(1,1,2,2)` is an order-join atom **but has a module** (hence `0` modular seeds at `n=4`, vs `1`
> SCC atom). So `modular-primes ⊊ strongly-connected`, and the substitution decomposition is
> **finer** than the SCC/order-join one — it carves inside strong components. This is the new
> content; it **coordinates with** the boxeph stub `THM-1955` ("which iso-classes come from
> smaller"), which reads the three modes as circulant character atoms — complementary to the
> modular-prime + spectral-substitution axis here, not a duplicate.

## 2. The spectral substitution law

For the uniform blow-up `T[S^m]` (each of `T`'s `k` vertices replaced by a copy of `S`, `|S|=m`),
the skew matrix is the **Kronecker sum**

```
S_{T[S^m]}  =  S_T ⊗ J_m  +  I_k ⊗ S_S       (J_m = m×m all-ones).
```

`J_m` has eigenvalues `m` (on the block-mean `𝟙`) and `0` (on `𝟙^⊥`). **If the seed `S` is
regular** (skew row-sums `0`, so `S_S 𝟙 = 0`, i.e. `𝟙 ∈ ker S_S`), the two subspaces decouple and:

> **`nonzero-spec(T[S^m]) = [ nonzero-spec(S), mult k ]  ∪  [ m · nonzero-spec(T) ]`.**

The seed spectrum appears `k` times; the block-tournament spectrum is scaled by the block size `m`.
**Verified for regular seeds `C_3, C_5` (any block `T`); it FAILS for an irregular seed** (transitive
`T_3`), where `S_S 𝟙 ≠ 0` mixes the block-mean into `𝟙^⊥`. Consequently **`char_S`, `var(λ²)`,
`SC4`, and every even moment `tr(S^{2j})` of a regular-seed substitution object are closed forms in
the seed moments.**

## 3. The octonion object `C_3[C_3]` (n = 9)

`n=9` is the octonion Cayley–Dickson level; `C_3[C_3]` is the **substitution-square of the smallest
seed**. It is a **regular** tournament (out-degree 4) with

```
char_S = x·(x²+3)³·(x²+27),     λ² ∈ {0, 3 (×6), 27 (×2)},
```

= the seed `√3` (multiplicity `k=3`) `∪` the block `m²·3 = 27`. Its invariants:
`var(λ²) = 104` (between transitive `2·C(9,3)=168` and the polystable `0`), `SC4 = 81 = 3⁴`,
`H = 3159 = 3⁵·13`. Because the spectrum is closed-form, `C_3[C_3]` and the family `C_3[C_{2j+1}]`
give **exact degree-8 (`tr(S⁸)`) test objects for the octonion wall** (THM-1935/1940): the correctly-
posed second-wall question can be probed on these self-similar seeds without enumerating all `n=9`.

## 4. The recursion modes as the composition's inclusion–exclusion

The three modes (`THM-442/553/801/830`) — full `A+B+C−D−E−F+G`, even-half `A+B−C`
(`=2h_{n−1}−h_{n−2}`), odd-half `A+B−C+D−E−F+G` — are the **inclusion–exclusion for counting the
staircase cell-data under this recursive composition** (the "seven letters as owned overlaps",
literal sub-staircase restriction/face maps `A,B,C` size `n−1`, `D,E,F` overlaps size `n−2`, `G`
triple size `n−3`). The prior-work map pins their meaning as **three characters** with sign
patterns: **Möbius `+++−−−+`** (full), **Legendre `χ₇` `++−+−−+`** (odd-half), **Eisenstein `χ₃`
`++−`** (even-half) — equivalently the cyclotomic factors `(x−1)^{depth}·Φ_d`. And `H(T)=∏H(atoms)`
is the **Euler-product analogue** (`φ`-multiplicativity ↔ `H`-multiplicativity). The even-half
`A+B−C` is the Chebyshev-Pell (`kps THM-1880`) recurrence. Reading the modes as the **module-type
restrictions** of the substitution tree (transitive vs cyclic-prime quotients) is the unification
this frame offers; the character `++−` is exactly the `C_3` base seed.

## Open

1. **`H` under substitution.** `THM-1860` gives `H = ∏ H(SCC)` for the *linear* quotient; find the
   law for a *cyclic* prime quotient (e.g. `H(C_3[S_1,S_2,S_3])`) — the general seed-composition law
   for the Hamiltonian-path count. (`H(C_3[C_3]) = 3159 = 3⁵·13`; the `13` is the cyclic-quotient
   correction.)
2. **The octonion wall on seeds.** Compute `tr(S⁸)` closed-form for `C_3[C_{2j+1}]` via the spectral
   law and test whether a degree-8 invariant decouples from `≤4`-vertex data at `n=9`.
3. **The seed census sequence.** `1,0,3,15` (`n≥3`) — identify/verify against OEIS (indecomposable
   tournaments) and extend to `n=7`.

## Verification

`04-computation/seed_tournaments_and_substitution_opus_S444.py` (+ `.out`) — the seed census `n≤6`,
`C_3[C_3]` invariants. `04-computation/spectral_substitution_law_opus_S444.py` (+ `.out`) — the
regular-seed spectral splitting (holds `C_3,C_5`; fails `T_3`).

---
source: monad-explorer-2026-06-07-S6 (built directly on opus-S703 THM-436 §5, the unused Klein/
  icosahedral handle). Pushes the "5-point cause" of Abel–Ruffini into verified combinatorics and
  ties its central number 15 four ways: overlapping cyclic-triangle pairs = involutions of A_5 =
  two-fold icosahedral axes = degree of Klein's edge-form. PROVED/VERIFIED counts; classical geometry.
status: SYNTHESIS with a PROVED core. (1) #overlapping cyclic-triangle pairs on [n] = 15·C(n,5);
  oriented = 60·C(n,5); the 5-point cause is LOCAL — one copy of A_5 per 5-set. (2) Canonical bijection
  overlapping-pairs ⟷ involutions of A_5 (shared vertex=fixed point; off-pairs=transpositions),
  verified as signature set-equality. (3) Commutator UNIFORMITY: all 60 oriented overlapping pairs give
  3-cycles, 3-to-1 onto the 20 three-cycles — the full content of "A_n perfect". (4) The icosahedral
  dictionary: A_5 class sizes {15,20,24}={#2axes·1,#3axes·2,#5axes·4}={15·1,10·2,6·4}; invariant degrees
  {6,10,15}=axis counts; Klein's T²+H³=1728f⁵ is the quintic-solver. (5) C_5 generates the whole
  icosahedron: 5-cycles=rotation, 3-cycles=cyclic triangles, involutions=overlapping pairs.
tags: [icosahedron, A5, klein, quintic, abel-ruffini, solvability, derived-series, involutions,
  two-fold-axes, cyclic-triangle, round-tournament, C5, overlapping-triangles, invariant-theory,
  j-invariant, 2-3-5, lonely-runner, hadwiger-nelson, monodromy, tournament-analysis]
---

# The icosahedral fifteen: where Abel–Ruffini's "5" becomes the icosahedron

THM-436 (opus-S703) proved the quintic-unsolvability threshold IS the appearance of **two cyclic
triangles sharing one vertex** (`3+3−1=5`), realized as the round tournament `C₅`, and §5 flagged a
*concrete classical handle the repo had not used*: Klein's icosahedral invariant theory. This session
picks up that handle and finds that the central number of the whole story — **15** — is the same `15`
in four guises, all on five points. The counts are computed and exact
(`04-computation/icosahedral_fifteen_monad_s6.py`); the geometry is classical (Klein, 1884).

## 1. The 5-point cause is *local*: `15·C(n,5)`

THM-436's witness counts `0, 0, 15, 90, 315, 840` (n=3..8) for overlapping cyclic-triangle pairs are
exactly `15·C(n,5)`; the oriented count is `60·C(n,5)`. So every overlapping pair lives in a **unique
5-subset**, and each 5-subset contributes exactly `15` unoriented (`60 = |A₅|` oriented) of them. The
engine that makes `A_n` perfect is assembled **five points at a time** — one copy of `A₅` per 5-set,
nothing genuinely larger ever needed. Abel–Ruffini's wall is not a global degree-`n` fact; it is `C(n,5)`
identical local copies of the single fact at `n=5`.

## 2. The `15` is the icosahedron — four times

On any 5-set the number `15` appears as:

| guise | object | why 15 |
|---|---|---|
| combinatorics | overlapping pairs `{X,Y}`, `|X∩Y|=1` | choose shared `v` (5) × match other 4 into 2 pairs (3) |
| group theory | involutions `(ab)(cd)` of `A₅` | choose fixed point (5) × match other 4 (3) |
| geometry | **two-fold axes** of the icosahedron | 30 edges / 2 |
| invariant theory | degree of Klein's edge-form `T` | the 15 axes / 30 edge-points |

The first two are connected by a **choice-free bijection** (verified as set-equality of signatures
`(v, {pair, pair})`): the shared vertex is the fixed point, the two "off-pairs" are the two
transpositions. The second and third are the classical fact that the involution about a 2-fold axis is
the `π`-rotation. So the `15` of THM-436(2) literally *is* the `15` edge-axes of the icosahedron.

And the rest of `A₅` fills in the rest of the solid. The computed class sizes are
`{15, 20, 24}` and they factor as **axis-count × rotations-per-axis**:
```
   15 involutions  = 15 two-fold axes  × 1     (30 edges)
   20 three-cycles = 10 three-fold axes × 2     (20 faces)
   24 five-cycles  =  6 five-fold axes × 4     (12 vertices)
```
so the three nontrivial **icosahedral invariant degrees `{6, 10, 15}`** (equivalently Klein's binary
forms of degree `{12, 20, 30} = {vertices, faces, edges}`) are exactly the **axis counts**
`{#5-axes, #3-axes, #2-axes}`. The relation between Klein's forms, `T² + H³ = 1728 f⁵`, is the
icosahedral equation — the `j`-invariant identity that *solves the general quintic* by transcendental
means. The number 1728 = `12³`, the `5` of `f⁵`, the `2,3` of `T², H³`: the icosahedron's `(2,3,5)`.

## 3. Commutator uniformity: the *whole* of "A_n perfect"

The classical proof of `A_n = [A_n, A_n]` exhibits **one** commutator that is a 3-cycle
(`[(123),(345)] = (032)`). The computation shows something cleaner and complete: **every one of the
60 oriented overlapping pairs has a 3-cycle commutator** — never the identity, never a
double-transposition, never a 5-cycle. The map
```
   {60 oriented overlapping pairs}  ──[ , ]──▶  {20 three-cycles}
```
is **onto and exactly 3-to-1**. Because the commutator is supported on the 5-point union `X∪Y`, this is
universal in `n`. So "A_n is perfect" is not one lucky commutator; it is a uniform `3`-to-`1` covering of
*all* the generators (3-cycles) by overlapping triangle pairs. The perfect core is built with perfect
regularity.

## 4. `C₅` *is* the generator — the tournament reading

The round tournament `C₅` realizes only `5` of the `10` 3-subsets as cyclic triangles (the other `5`
are transitive); among those `5`, exactly `5` pairs share one vertex and `5` share two. So `C₅` carries
`5` of the `15` abstract overlapping pairs — its *own* slice of the icosahedron. But all three `A₅`
classes have direct `C₅` meaning:

> **5-cycles = the rotation `(01234)`** (the round tournament's defining cyclic symmetry);
> **3-cycles = the cyclic triangles** (the rock-paper-scissors atoms, OCF's smallest odd cycle);
> **involutions = the overlapping-pair structure** (the 2-fold axes).

The entire icosahedral group — hence the entire obstruction to the quintic — is generated by the
combinatorics of the **single smallest interesting tournament**. The repo's `C₅` (THM-403's cyclotomic
worry-set witness; THM-436's overlap realization) and Klein's icosahedron are the same object on five
points, read once as a tournament and once as a Platonic solid.

## 5. The `(2,3,5)` thread, and the inversion restated

The icosahedron's signature is `(2,3,5)`: axis-orders `2,3,5`, defect `1 − (1/2 + 1/3 + 1/5) = 1/30`,
forms `T², H³, f⁵`. The repo's worry-set / unit-distance program has independently isolated exactly
three carry-prime frontiers: **prime-2** (doubling, THM-404), **prime-3** (THM-407/428, the `3³ = 27`
shell at the hard `n=14`), and **`n=5` / cyclotomic-5** (THM-403/436). The resonance — that these are
the icosahedron's three axis-orders — is recorded as **HYP-2305** (speculative, with one concrete
handle). It also re-frames THM-436 §4's *inversion*: in Galois theory the icosahedral/`A₅` case is the
*hard* (unsolvable) one; in LRC the cyclotomic/abelian (solvable) case is the *hard/tight* one. The
icosahedron sits at the seam — it is simultaneously the *most symmetric* `ℝ³` object (low complexity)
and the bearer of the *first unsolvable* group (high complexity). The two readings of the
symmetry↔complexity axis cross exactly here, on the 5-point set both problems keep returning to.

## Honest status

- **Proved/verified** (`…fifteen_monad_s6.py`): the `15·C(n,5)`/`60·C(n,5)` counts (n=3..8); the
  overlapping-pair ⟷ involution bijection (signature set-equality); commutator uniformity (all 60 → 20
  three-cycles, 3-to-1); the `A₅` class-size ↔ axis-count factorization; the `C₅` shadow.
- **Classical** (cited, not re-derived): icosahedron `V,E,F = 12,30,20`; invariant degrees `{2,6,10,15}`
  / binary forms `{12,20,30}`; the relation `T²+H³=1728 f⁵` and Klein's icosahedral solution of the
  quintic.
- **Speculative** (HYP-2305): the `(2,3,5)` ↔ (prime-2, prime-3, `n=5`) frontier identification and its
  LRC/HN transfer. Resolves no open case.

**Artifacts:** `04-computation/icosahedral_fifteen_monad_s6.py` (+`05-knowledge/results/…out`); canon
addendum on `THM-436`; new **HYP-2305**. Builds on THM-436/THM-403/THM-406/HYP-2303 and the reflection
`the-solvability-tower-galois-lrc-icosahedron-s703.md` (opus-S703).

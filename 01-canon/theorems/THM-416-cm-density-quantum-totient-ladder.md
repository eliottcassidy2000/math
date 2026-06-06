---
id: THM-416
title: The unit-distance density quantum is the number of roots of unity, and its spectrum is the Euler-totient ladder
status: PROVED
source: monad-explorer-2026-06-06-S706
depends_on:
  - THM-412
related:
  - HYP-2262
  - HYP-2263
  - HYP-2267
  - HYP-2274
  - THM-403
  - THM-414
  - T755
  - T757
  - T760
---

# THM-416 — The Density Quantum is `#(roots of unity)`, and its Spectrum is the Totient Ladder

## Summary

THM-412 proved that a **2D** lattice quantizes the achievable unit-distance
density `r_Q(D)/2` in steps of `w/2`, where `w` = number of proper automorphs =
number of roots of unity in the imaginary quadratic order, `w ∈ {2, 4, 6}`,
capped at `6` (triangular/Eisenstein). This theorem shows the mechanism is **not
special to dimension 2**: it holds in **every CM field**, with the same one-line
free-action proof. Consequently the spectrum of achievable density quanta in
dimension `N = 2d` is exactly the **Euler-totient ladder**

```text
   { w/2 : w even, φ(w) | N },     maximal quantum  M(N)/2,  M(N) = max{m : φ(m) | N}.
```

The "`2D group spectrum between the triangular lattice (κ=6) and the CM field`"
(the dispatched probe; S702/HYP-2263) is therefore a **discrete totient ladder
with a hard cap of `3` in 2D**, not a continuum: every quantum above `3`
requires a **dimension (rank) jump**. This rigorously confirms and *explains*
S702's "no 2D bridge group; the real axis is rank/degree."

## Setting

Let `K` be a CM field of degree `N = 2d` (a totally imaginary quadratic
extension of a totally real field), `O ⊆ K` an order, and `μ(K)` its group of
roots of unity, of order `w = #μ(K)` (always even). Embed `O` as a full-rank
lattice in `R^{2d} ≅ C^d` by the **canonical (Minkowski) embedding**

```text
   ι(g) = (σ_1(g), …, σ_d(g)) ∈ C^d,     ‖ι(g)‖² = Σ_k |σ_k(g)|²,
```

one `σ_k` from each conjugate pair of complex embeddings. For `D > 0` write the
**shell count** `r(D) = #{ g ∈ O : ‖ι(g)‖² = D }`. Choosing the unit distance
to be `√D` makes the unit-distance graph on a large patch of the lattice
`r(D)`-regular, so the interior **unit-distance density** is `r(D)/2`.

## Statement

**(I) Quantization.** For every `D > 0`,  `w | r(D)`.  Hence the unit-distance
density `r(D)/2` is an integer multiple of the **quantum** `w/2 = #μ(K)/2`.

**(II) Spectrum.** A value `w` occurs as `#μ(K)` for some CM field of degree
dividing `N` iff `w` is even and `φ(w) | N`. Therefore the set of achievable
density quanta in dimension `N = 2d` is `{ w/2 : w even, φ(w) | N }`, and the
**maximal quantum** is `M(N)/2` with `M(N) = max{ m : φ(m) | N }`.

**(III) The 2D cap and the dimension jump.** `M(2) = 6` (φ(m) | 2 ⟺ m ∈ {1,2,3,4,6}),
so the largest 2D quantum is `3` — the Eisenstein/triangular value, the maximal
number of roots of unity in *any* imaginary quadratic field. Quanta `> 3`
(equivalently `w > 6`) are **unrealizable by any 2D lattice** and first appear
one dimension up: `M(4) = 12` (quantum `6`, via `Q(ζ_12)`). In particular the
quanta `4` and `5` — which the *special* 2D lattices skip (square forces even
quanta, triangular forces multiples of 3) — become **rigid** quanta of the
degree-4 CM fields `Q(ζ_8)` (`w = 8`, quantum `4`) and `Q(ζ_5)` (`w = 10`,
quantum `5`). This resolves the S703 handoff (4) on the role of density 5.

## Proof

**(I).** `μ(K)` is cyclic of order `w`; let `ζ` generate it. Multiplication by
`ζ` is an `O`-module automorphism, and under the canonical embedding it is an
**isometry**: `σ_k(ζg) = σ_k(ζ)·σ_k(g)` and `|σ_k(ζ)| = 1` because every
conjugate of a root of unity is a root of unity (modulus 1), so
`‖ι(ζg)‖ = ‖ι(g)‖`. Thus `μ(K)` acts on each shell `R_D = {g ∈ O : ‖ι(g)‖² = D}`.
For `D > 0` the action is **free**: `ζ^j g = g` with `0 < j < w` gives
`(ζ^j − 1) g = 0` in the field `K`, and `ζ^j ≠ 1` is a unit, so `g = 0`,
excluded. A free action of a group of order `w` partitions `R_D` into orbits of
size `w`, hence `w | r(D)`. ∎

This is exactly the THM-412 argument with "planar rotation through `2π/w`"
replaced by "isometry of the canonical embedding"; dimension 2 played no role
beyond `|σ(ζ)| = 1`, which holds for all archimedean places.

**(II).** `Q(ζ_m) ⊆ K` iff `K` contains a primitive `m`-th root of unity, which
requires `φ(m) = [Q(ζ_m):Q]` to divide `[K:Q] = N`. Conversely, if `φ(m) | N`
then `K = Q(ζ_m)·F` for a suitable totally real `F` is a CM field of degree `N`
containing `μ_m`. Roots of unity in a CM field form `μ_m` for a single even `m`
(taking `m` even WLOG, since `μ` of odd order `k` equals `μ_{2k}`). Hence the
realizable `w = #μ(K)` are exactly the even `m` with `φ(m) | N`, and the maximum
is `M(N) = max{m : φ(m) | N}` (finite because `φ(m) → ∞`). ∎

## Corollary (the LRC worry-set field) — the cyclic rosette IS a density quantum

For the **LRC modulus** `C = 2n − 1` (always odd), the worry-set witnesses are the
`C`-th roots of unity (THM-401/403), living in the cyclotomic field `Q(ζ_C)`. Its
roots of unity number `w = #μ(Q(ζ_C)) = 2C` (for odd `C`, `Q(ζ_C) = Q(ζ_{2C})`),
so by Part (I) the cyclotomic lattice `Z[ζ_C]`, of degree `φ(2n−1)`, has

```text
   unit-distance density quantum  =  w/2  =  C  =  2n − 1.
```

The **Euclidean** density quantum of the LRC's *own* cyclotomic field equals the
order of the **cyclic** worry-set rosette. This is the exact realization of the
"all `(2n−1)`-th roots at once" intuition: inside `Q(ζ_{2n−1})` the full
`(2n−1)`-rosette is present as a free isometry group, so the dimension is no
longer a constraint — `φ(2n−1)` is precisely large enough to hold all of `μ_C`.
Verified `C = 5, 7, 9` directly (the `Z[ζ_5], Z[ζ_7], Z[ζ_9]` rows below give
quanta `5, 7, 9 = C`); the identity `w/2 = C` is immediate from Part (I) and
`#μ(Q(ζ_C)) = 2C`. The campaign's critical case `n = 14`, `C = 27 = 3³`, sits in
degree `φ(27) = 18`, where it is the **maximal-quantum** field (`M(18) = 54`,
quantum `27`); the `3`-richness of `C = 27` that drives the V* split (THM-413,
HYP-2262) is the same `3`-richness that makes `Q(ζ_27)` rosette-maximal in its
dimension.

## Verification

`04-computation/cm_density_quantum_totient_ladder_s706.py`
(`05-knowledge/results/cm_density_quantum_totient_ladder_s706.out`), **0
failures**:

* **(A)** 2D orders disc `-3, -4, -7, -8, -11, -15, -23`: `w | r(D)` re-confirmed
  (matches THM-412).
* **(B)** CM cyclotomic lattices, the genuinely `>2D` cases, Gram matrix =
  Ramanujan-sum matrix `c_m((i−j) mod m)` (integral):
  `Z[ζ_5]` (deg 4, `w=10`, quantum **5**), `Z[ζ_8]` (deg 4, `w=8`, quantum **4**),
  `Z[ζ_12]` (deg 4, `w=12`, quantum 6), `Z[ζ_7]` (deg 6, `w=14`, quantum 7),
  `Z[ζ_9]` (deg 6, `w=18`, quantum 9). Every trusted shell count is a multiple of `w`.
* **Rigorous capture:** only shells with `D ≤ λ_min(G)·B²` are checked. A vector
  with `cᵀGc = D` has `|c| ≤ √(D/λ_min)`, so `D ≤ λ_min·B²` guarantees the whole
  shell lies inside the box `[−B,B]^dim` — no truncation artifact (a fully
  captured ball is closed under the isometry group).
* **(C)** Totient ladder `M(N)/2` for `N = 2..24`: maximal quanta
  `3, 6, 9, 15, 11, 21, 3, 30, 27, 33, 23, 45` (witnessed by `Q(ζ_{M(N)})`).

## Remarks

**Non-monotone ladder / nontotients.** `M(N)` *drops* at `N` that are
**nontotient** in the relevant sense — e.g. `N = 14`: `φ(m) = 14` has no
solution, so `φ(m) | 14 ⟺ φ(m) ∈ {1,2}`, giving `M(14) = 6` (quantum back down
to 3). The density-quantum spectrum inherits the irregular arithmetic of the
totient image; it is *not* "the more dimensions, the bigger the rosette."

**Why the cyclic LRC floor uses the full `(2n−1)`-rosette.** On the Euclidean
side the rosette is **capped by the ambient dimension** via the totient (a root
of unity `ζ_m` is realizable as an isometry only if `φ(m) | N`). The LRC
worry-set floor (THM-403) instead lives in the cyclic group `Z/(2n−1)`, whose
order-`(2n−1)` rosette is realized **intrinsically** — the cyclic metric carries
*all* `(2n−1)`-th roots regardless of any ambient dimension. This is the precise
sense in which the LRC (cyclic) side is the "`degree → ∞`, all-roots-at-once"
limit of the Euclidean ladder: same root-of-unity rosette skeleton (THM-412/S702
unification), but the Euclidean one is totient-capped and the cyclic one is not.
(Complements THM-414: the cyclic *degree-2 additive* face `r_+(s)` is
matching-capped, escaping only at degree 3.)

**Honesty.** Part (I) is elementary and proved (the free action). Part (II) is
elementary algebraic number theory (`Q(ζ_m) ⊆ K ⟺ φ(m) | N`). The contribution
is the **role** in the unit-distance / LRC density program: it pins the
density-quantum spectrum to the totient ladder, gives the exact 2D cap (3) and
its first escape (dimension 4, quanta 4,5,6), and frames the Euclidean-vs-cyclic
rosette asymmetry. The arithmetic resonance "`quantum 21 = forbidden H-value 21`,
at degree 12 via `Q(ζ_42)`" is recorded as an observation, **not** a theorem
(HYP-2274). This is **not** a result about the regular LRC at any specific `n`.

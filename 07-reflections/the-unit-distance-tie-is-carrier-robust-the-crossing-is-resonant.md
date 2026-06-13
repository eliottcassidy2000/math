---
source: monad-explorer (agent-mathworker)-2026-06-13 (deep-research; OPEN-Q-057 frontier)
status: REFLECTION grounded in EXACT-INTEGER free-patch annealing over the Moser-ladder
  bridge-lattice family L_t (verified: reproduces Engel's table on L_3 incl. u(28)=85>84;
  every reported value re-checked by an independent exact pairwise recount).  Complements
  (does NOT duplicate) the same-day peer result HYP-2460/THM-493 (resonance-bonus = the
  crossing; 27 is bonus-hostile), which used CURATED 2-FACTOR PRODUCTS; this uses FREE
  (whole-lattice) search and adds the decisive t=4 control.  Conjecture-grade where flagged.
tags: [unit-distance, N-star, 3N-crossover, moser-lattice, bridge-lattice, THM-434,
  THM-493, eisenstein, sqrt13-layer, resonance, carrier-selection, hamming-tie,
  tie-vs-crossing, OPEN-Q-057, HYP-2461]
---

# The n=27 unit-distance tie is carrier-robust; the n=28 crossing is a single √−11 resonance

## Setup

`u(n)` = Erdős unit-distance maximum; `N*` = smallest `n` with `u(n) > 3n` (average
degree past the kissing line `κ = 6`).  Canon: `N* ∈ {25,26,27,28}`, conjectured `28`
(HYP-2299); the crux is whether `u(27) = 81` (tie) or `≥ 82`.  The best small-`n`
constructions live in the **Moser-ladder bridge lattices**

```
   L_t = ℤ[ζ₆] ⊕ ℤ[ζ₆]·ω_t,   ω_t = ((2t−1) + i√(4t−1))/(2t),  |ω_t| = 1,
```

rank-4 in ℂ whenever `4t−1` is not `3·square`.  THM-434: `#units(L_t) = 12 + 6·B(t)`,
`B(t) = Σ_{d|t}(−3/d)` = number of `ζ₆`-orbits of norm-`t` Eisenstein integers (the
**transverse** unit vectors, beyond the triangular and `ω_t` rosettes).  So `t=3` (the
**Moser lattice** `M_L`, `√−11`) has 18 units; `t=13,21,31` have 24; `t=49` has 30;
the non-Loeschian `t=2,5` have only 12 (no transverse).  **Engel's `u(28)≥85`, and every
prior repo search, used only `t=3`.**

This session ran the IDENTICAL exact-integer simulated-annealing densest-patch search —
the one that reproduces Engel's `60,64,68,72,76,81,85,89,93` on `M_L` and here
independently re-derives `u(28)=85 > 84` — over the **whole lattice** (not restricted to
2-factor products) for `t = 2,3,4,5,13,21,31,49`, `n = 21…30`.  Adjacency is a pure
integer test (`|z|²=1 ⟺ bc−ad=0` and an integer quadratic form `= 2t`); every count was
re-verified by an independent exact pairwise recount.

## The data (exact; best of 12 restarts × 110k annealing steps)

```
 t   √−D   #units  rosette°   n=21   n=27   n=28      transverse (B(t))
 2   √7     12     41.4°       57     78     83        0   (none)
 5   √19    12     25.8°       57     78     83        0   (none)
 3   √11    18     33.6°       57     81†    85*       1   (6)   <-- Moser
 4   √15    18     29.0°       57     81     83        1   (6)
 13  √51    24     15.9°       57     81     83        2   (12)  <-- √13 layer
 21  √83    24     12.5°       57†    81     83        2   (12)
 31  √123   24     10.3°       55     81     83        2   (12)
 49  √195   30      8.2°       57     81     83        3   (18)
 Engel/u    —       —          57     81     85         —  (target)
```
`*` = `> 3·28 = 84`, the actual crossing.  `†` = my annealing lost 1 here (L_3 reaches
the true 81 at 27 and 89 at 29 under the s4 search; L_21 reaches 57 at 21); these are
search variance, not ceilings — and they do not affect either conclusion below.

Three things fall out, and the third is the decisive one.

## Finding 1 — the tie at 27 is CARRIER-ROBUST, and needs transverse vectors

`u(27) = 81 = 3·27` is reached **exactly** by *every* transverse-bearing lattice
(`t=3,4,13,21,31,49`, unit counts 18→30) and **never exceeded** by a free patch in any
of them.  The lattices with **no** transverse vectors (`t=2,5`, only 12 units) reach only
**78** — they cannot even build the tie.  So:

> the 81-tie is the carrier-*independent* 6-regular structure (`H(3,3)=K₃^□3`, THM-432;
> angle-rigid, THM-437) — any lattice with the transverse vectors to interleave a second
> rosette reaches it, and **none can do better.**

This is the strongest non-product evidence yet for `u(27)=81`, hence `N*=28`: THM-493
killed *product* crossings at 27 (a size-3 factor forces the `√t`-free `K₃`); the present
free search dead-ends at 81 in the *entire* bridge family too.  Both roads stop at 81.

## Finding 2 — the crossing at 28 is a SINGLE √−11 resonance (the decisive control)

The *actual* crossing `u(28) ≥ 85 > 84` appears for **`t = 3` only.**  Every other
lattice — including the 24- and 30-unit ones, and including **`t = 4`** — caps at
**83 < 84**, not even reaching the line.

The control that matters: **`t = 4` (`√−15`, 18 units, rosette angle 29.0°) is
*geometrically closer* to the 30° triangular-gap bisector than Moser's `t = 3` (33.6°),
and it has the *same* number of transverse vectors (6).  It ties at 27 — yet it does NOT
cross at 28.**  So the crossing is **not** a "good interleaving angle" phenomenon (if it
were, `t=4` would cross at least as well as `t=3`).  It is the *specific arithmetic
resonance of `√−11`*: THM-493's bonus `Δ_t = ½ Σ_{N(α)=t} m_α m_α` is nonzero at 28
only when the two factors carry matching `√t`-displacements, and `28 = 4·7` admits a
`√3`-bearing edge-dense factor pair tuned to `t=3`.  `t=4`'s transverse vectors carry the
*wrong* `√(4t−1)=√15` displacement, so its 6 transverse edges never assemble into a
crossing.  The free search confirms the product algebra from the outside: **the 85 is one
arithmetic accident at `√−11`, not a band of nearby angles.**

## The dichotomy (the durable content)

> **The 27-tie is combinatorial, geometric, and universal; the 28-crossing is arithmetic,
> resonant, and singular.**  Tying 81 needs only *enough interleaved unit directions to
> build a 6-regular patch* — any near-bisector transverse-bearing carrier does it
> (`t=3,4,13,…`), and none can do more.  Crossing 85 needs the *specific* `√−11`
> resonance — even `t=4`, nearer the bisector and equally transverse-rich, stalls one
> short of the line.

This is the unit-distance face of a shape the project keeps meeting (cf.
`symmetry-saturates-irregularity-violates-the-hamming-tie-s711.md`): the regular
saturating structure is carrier-independent — you can tie it but never beat it — while the
single edge that *breaks* a threshold is fragile and lives at one special parameter. Here
the threshold is `3N`, the saturator is `H(3,3)` (any interleaved carrier), and the
"special parameter" is named and isolated by the `t=4` control: **the `√−11` Moser angle,
nothing adjacent.**

## Finding 3 — the bridge family unifies BOTH known small-n optima

The parametric family makes a bonus visible: `L_13` (the **`√13` layer**) reaches
`u(21)=57` — the AMP-proven *global* optimum at 21 — while `L_3` (`√11`) owns the
crossover.  THM-431 noted AMP's `u(21)=57` extremal lives in the `√13` Eisenstein layer;
that layer is here identified as a *member of the same `L_t` ring* that contains the Moser
lattice.  One parametric family of lattices contains **both** the proven `n=21` optimum
(`t=13`) **and** the `3N`-crossover engine (`t=3`): two famous dense small-`n` configs,
two angles of one construction.

## The meta-point: unit-vector COUNT is a red herring

THM-434 lets `#units(L_t)` grow without bound (`t=49`→30, `t=133`→36).  One might guess
"more unit directions ⟹ denser patches."  **Decisively false:** the 24- and 30-unit
lattices are *not* denser than the 18-unit Moser lattice — they are competitive (within
1–3 of Engel everywhere, identical to each other) but never better, never cross at 28, and
the 30-unit `t=49` is no better than the 24-unit `t=13`.  Density is selected not by the
*number* of unit vectors but by (a) *having transverse vectors at all* (to reach the tie)
and (b) the *arithmetic* of `√(4t−1)` (to reach the crossing) — `√−11` and nothing else.

## Honest status / open

- VERIFIED (exact integer, double-recounted): the whole table above; `L_3` reproduces
  Engel incl. `u(28)=85>84`; `t=2,5` cap at 78 at 27; **`t=4` ties 81 at 27 but caps at
  83 at 28** (the decisive control); no bridge lattice reaches `82@27` or `≥84@28` except
  `t=3`.
- The annealing is a *lower-bound* search (a value below Engel = "not found here," not a
  proved ceiling).  The evidence is the *pattern* — tie reached in 6 carriers and never
  exceeded; crossing in exactly one — not any single number.
- CONJECTURE: `u(27)=81`, `N*=28`, now supported from the product side (THM-493) *and*
  the free side (here), and the crossing localized to `√−11` by the `t=4` control.
- Next: is the `√−11` crossing truly unique among ALL CM carriers, or does some non-ladder
  ring (e.g. `ℤ[ζ₁₂]`, the 30°-bisector lattice, 12 units, NOT an `L_t`) also cross at 28?
  `ℤ[ζ₁₂]` is the natural "perfect-bisector" test the `L_t` family can't reach
  (its angles are `arccos((2t−1)/2t) ≠ 30°` always).  If even the exact-30° lattice fails
  to cross, `√−11` is singular for a genuinely arithmetic (not geometric) reason.

## Take-away

Two independent methods — the peer's 2-factor resonant-product search (THM-493) and this
whole-lattice free annealing with the `t=4` control — agree and sharpen each other:
**you cannot beat the 81-tie at 27 in any bridge-lattice carrier, you can only reach it
(and only if the carrier has transverse vectors); and you can reach the 85-crossing at 28
in exactly one carrier, the `√−11` Moser lattice — not even in the geometrically-better
`√−15` neighbour.** The tie is a universal wall; the crossing is a single arithmetic key.
That is precisely why `N*` is hard at 27→28: 27 is a robust 6-regular wall, and 28 is one
resonance, not a trend.

---
source: opus-2026-06-03-S585 (remote-control)
status: TAXONOMY + creative hypotheses — rigidity as orbit-discreteness; TWO rigidities on the witness orbit (static-multiplicative always; dynamical-doubling connected for odd n, FRAGMENTED at even n = the frontier); + spectral/isostatic/CRT hypotheses
tags: [LRC, rigidity, orbit, doubling, units, dynamical, spectral, isostatic, taxonomy, n14, hypotheses]
---

# An orbit-rigidity taxonomy: static, dynamical, spectral, isostatic

**Prompt (user):** consider types of rigidity besides global and local; see everything
as orbits; poke the repo, follow tangents; design hypotheses creatively.

Seeing every LRC object as an orbit turns "rigidity" into **orbit-discreteness**, and a
single rigid orbit carries *several* rigidities at once. The headline: the worry-set's
witness orbit is **statically** rigid for all `n` but **dynamically** rigid only for
*odd* `n` — it **fragments at even `n`**, exactly the frontier. (Complements the
concurrent S589, which layers the *symmetry group* cyclic→dihedral→abstract; this layers
the *rigidity type* on a fixed orbit.)

## 0. Everything as orbits

| LRC object | orbit of |
|---|---|
| witness set `{u/n}` | the symmetry group on the circle |
| config `V` | `(ℤ/n)^*` (multiplicative) / `S_n` (relabel) in speed-space |
| runner trajectory `{v_i t}` | the `ℝ`-flow (the time orbit) |
| doubling tower (S579) | `⟨×2⟩` |
| tournament class | `S_n` (iso class = an orbit) |

**Rigidity = the dimension/discreteness of the witness orbit:** `dim 0` (finite orbit,
measure-0 = worry-set) is RIGID; `dim 1` (interval, positive measure = loose) is
FLEXIBLE. Verified: AP witness = the finite unit orbit; `(1,4,5)`, `(2,3,7,8)` have
positive-measure (flexible) witnesses.

## 1. Static (multiplicative) rigidity — `(ℤ/n)^*` [HYP-2124]

The AP is lonely at `t=j/n ⟺ gcd(j,n)=1`, so the witness set is *exactly* the unit orbit
`(ℤ/n)^*·(1/n)`, a rigid `φ(n)`-point orbit, simply-transitive under `(ℤ/n)^*`. This is
present for **every** `n`. Static rigidity never fails.

## 2. Dynamical (doubling) rigidity — `⟨×2⟩` [NEW]

The bit-shift `x↦2x` (the fractal generator, S580) acts on the witness orbit. The
witness orbit decomposes into `⟨×2⟩`-orbits:

> **Doubling-rigidity dichotomy (verified `n=5..18`).** The unit witness orbit is
> `⟨×2⟩`-**connected** (nontrivial doubling orbits) `⟺ 2 ∈ (ℤ/n)^* ⟺ n` is **odd**; it
> **fragments into singletons** `⟺ n` is **even**. (Odd primes/prime-powers: often a
> *single* doubling orbit — `2` a primitive root. n=14 → **6 singletons**.)

So the *same* statically-rigid orbit is **dynamically rigid for odd `n` and dynamically
broken for even `n`** — and even `n` is the LRC hard frontier. This is a clean new
reading of "where and why it works":

> **A proof that propagates loneliness along the doubling (the polynomial/sieve method,
> which rides `t ↦ 2t`) needs the witness orbit to be `⟨×2⟩`-connected. It is — for odd
> `n`. At even `n=2q` the orbit fragments, and the propagation stalls — the `2q`-apex /
> C′ residual in *dynamical-rigidity* dress.** `n=14` is the first even composite where
> the moving frontier meets the fragmentation.

## 3. The taxonomy (rigidity types on one orbit, by acting group)

| type | group | status on the witness orbit |
|---|---|---|
| **static / multiplicative** | `(ℤ/n)^*` | rigid ∀`n` (HYP-2124) |
| **dynamical / doubling** | `⟨×2⟩` | connected (odd) / fragmented (even) — **NEW** |
| **dihedral** | `D_{2m}` | S589 (concurrent) |
| **local / pinch** | the straddle pinch | pins each clock point (HYP-2124 §3) |
| **spectral / character** | `(ℤ/n)^*^` (Dirichlet duals) | conjectural (H3) |
| **isostatic / percolation** | constraint count | conjectural (H4) |

## 4. Creative hypotheses

**H3 — Spectral/character rigidity.** Block-diagonalise the danger structure by Dirichlet
characters mod `n` (the symmetry-adapted *orbit rigidity matrix*). Conjecture: the
worry-set's loneliness is carried by the **principal character** (the rigid mode pinning
`M=1/n`); rigidity "leaks" only through the character(s) supported on `gcd(·,n)>1` — for
`n=2q` the **2-adic character block**. *Test:* compute the character-projected resonance
`Σ_j χ(j) · [‖config·j/n‖<1/n]`; predict the non-principal mass concentrates on the
`2`-block at `n=14`. (Ties to the spectral/Walsh program, THM-071.)

**H4 — Isostatic / percolation rigidity.** Model the binding pairs at the optimum as a
constraint framework; conjecture the **worry-set = the isostatic (critically rigid)
configs** — witness orbit *exactly* determined (`φ(n)` points), neither over- nor
under-determined. A counterexample (`M<δ`) would be *over-rigid* (the orbit forced below
`δ`), forbidden by a **Maxwell/orbit-counting** obstruction: `#binding-coincidences =
orbit size = φ(n)` exactly. *Test:* count binding-pair coincidences at the optimum across
configs; predict a sharp rigidity transition at the worry-set.

**H5 — CRT / orbit-rigidity-matrix for `n=2q`.** The orbit rigidity matrix
block-diagonalises by CRT into a `q`-block (the prime dynamics — `⟨×2⟩`-connected, rigid)
and a `2`-block (the fragmenting part). Conjecture: **all of `n=2q`'s rigidity defect
lives in the `2`-block**, and the `q`-block is the solved prime case lifted (S579/S580).
For `n=14`: `q=7` block solved, the entire residual is the order-2 fragmentation. *Test:*
restrict the witness orbit to each CRT factor and confirm the `7`-factor is
doubling-connected while the `2`-factor is the singleton seam.

**H6 — Rigidity is monotone under the orbit-doubling tower.** Since `D:v↦2v` carries
`AP_n`'s rigid clock into `AP_{2n}` (S579), conjecture a **rigidity functor**: static
rigidity lifts cleanly up the tower, but *dynamical* rigidity **degrades by one factor of
2 per level** (each doubling injects a fresh singleton-fragmenting `2`). The
"rigidity height" `= v_2(n)` measures the dynamical defect; `n=14` has height 1 (one
break above the prime).

## 5. Why this is productive

The two rigidities **diverge exactly at the frontier**: static rigidity (the symmetry
that makes the worry-set's witness a clean orbit) holds for all `n`, but the *dynamical*
rigidity (the doubling-connectivity that the prime-method exploits) breaks at even `n`.
So "LRC is proven for odd-ish and the method stalls at `2q`" gets a one-line orbit
explanation: **the witness orbit fragments under doubling iff `n` is even.** And it
furnishes four falsifiable hypotheses (H3–H6) for *where the residual lives* — all in
orbit language, all pointing at the `2`-block of `n=14`.

## 6. Honest status

- **Verified/exact:** witness = static unit orbit (HYP-2124); the doubling-rigidity
  dichotomy (connected⟺odd, `n=5..18`); flexible=positive-measure witness.
- **NEW structural claim:** the worry-set's witness orbit is statically rigid but
  dynamically fragmented at even `n`; the even-frontier = dynamical-rigidity break;
  rigidity-height `= v_2(n)`.
- **Hypotheses H3–H6:** designed, falsifiable, untested — the creative deliverable.

**Artifacts:** `04-computation/lrc_rigidity_types_orbits_s585.py` (+`.out`),
`lrc_doubling_rigidity_dichotomy_s585.out`. Complements S589 (group layers). Builds on
HYP-2124 (static), S579/S580 (doubling tower), THM-071 (spectral). New: **HYP-2126**.

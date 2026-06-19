---
source: kind-pasteur-2026-06-19-S12
status: reflection on a verified structural reframe (HYP-2637)
tags: [lonely-runner, lrc14, freiman, additive-combinatorics, sumset, gap-dimension, summand-graph, 2n-minus-1, recursion, pockets]
---

# The LRC "pockets" are Freiman dimension, and the difficulty is the factorization of 2n−1

**Prompt (user):** think about the third pocket, look for a fourth, see the recursive structure
that implies these pockets as the number of runners changes; "every element in a small relation =
high additive energy = Freiman GAP"; revisit the summand/multiplicand graph and how
addition/multiplication relate to even/odd and positive/negative.

## The native language was additive combinatorics all along

For a year the project circled the Lonely Runner extremal problem with measure-theoretic and
Fourier tools, and kept hitting the same wall: *the AP is extremal, but no monotone/coupling/order
argument sees it.* This session the wall dissolved into a sentence:

> **The configuration pockets are Freiman dimension.** The extremal functional `L_y` (≈ `meas(S7)`)
> is a strictly decreasing function of the additive dimension of the speed set, and the AP is the
> unique dimension-1 object.

Concretely: `excess(E) = |E+E| − (2|E|−1) ≥ 0`, and by **Freiman–Vosper, excess 0 ⟺ E is an AP.**
The data then shows `L_y` falling by ~½ per Freiman dimension — a clean **dimension penalty**:
```
   d=1 (AP)      L_y = 0.358   POCKET 1   the unique maximum
   d=2 (GAP)     L_y ≈ 0.15    POCKET 3
   d=3 (GAP)     L_y ≈ 0.10    POCKET 4   ← the fourth pocket the prompt asked for
   ...           ...           POCKET d+1
   dissociated   L_y → 0.05    POCKET 2   (the independent limit)
```
The "third pocket" that blocked the proof — wide sets where *every element sits in a small relation*
— was exactly **high additive energy = small doubling = a `d≥2` GAP**, and those sit *far* below the
cap (margin ≥0.21). The thing that looked like an uncontrolled residue is the most additively
structured region, and structure is what bounds it. The only genuinely tight case is the AP itself
(dimension 1), which is a finite, exact check.

This is why every *order* argument failed: stochastic dominance, convex order, monotone descent.
They impose a 1-parameter ladder. But the extremality is **dimensional**, not ordinal — `L_y`
descends a lattice of dimensions, not a line. The AP is not the bottom (or top) of an order; it is
the apex of a *dimension* hierarchy, the unique point of maximal additive rigidity (minimal sumset).

## Addition and multiplication are the two graphs, and they are the dimension

The summand/multiplicand dictionary the project built (S560, S571, HYP-2083) turns out to BE the
dimension story:
- **Addition** = the sumset `E+E`; its excess is the additive-energy deficit; the AP minimizes it.
  The relations `Σ m_e e = 0` live mod `C = 2n−1` as antipodal shells `{a, C−a} = {+a, −a}` — the
  **positive/negative** pairing is literally `a ↔ −a mod C`.
- **Multiplication** = the GAP *step directions*. A `d`-dimensional GAP `{Σ x_i a_i}` has `d`
  multiplicative generators `a_i`; dimension = number of independent multiplicative directions =
  number of independent additive relations' complement. The units `(Z/C)^*` act on the shells,
  deciding visibility.
- **Odd** = `C = 2n−1` is odd, so there is no midpoint shell — clean antipodal pairing.
- **Even** reappears only at the floor/apex degeneracy.

So "addition makes the shells, multiplication makes them visible" (the old slogan) is the statement
that the additive relation lattice (addition) and its GAP-generator directions (multiplication) are
the two coordinates of the Freiman dimension. The **GAP collisions** — the mechanism of the dimension
penalty, where two GAP directions `a_i x, a_j x` become commensurate and orbit points pile up — are
exactly the nonzero relations in the shells.

## The recursion across n: difficulty is the factorization of 2n−1

The deepest thing the dimension lens reveals is *why n=14*. The summand modulus is `C = 2n−1`, and
its prime factorization sets the **stratum depth** of the relation graph:
```
   n=10: 2n−1 = 19      prime         1 stratum  (units only)
   n=12: 2n−1 = 23      prime         1 stratum
   n=14: 2n−1 = 27 = 3³ prime POWER   3 strata   (units, gcd-3, gcd-9)  ← deepest
   n=16: 2n−1 = 31      prime         1 stratum
   n=18: 2n−1 = 35 = 5·7              2 strata
```
**n=14 is anomalous: 27 = 3³ is the unique nontrivial prime power among `2n−1` in range**, and the
AP-to-cap margin on the binding row confirms it QUANTITATIVELY (workflow recursion engine, exact):
`n=10: +0.204, n=12: +0.192, n=14: +0.054 (tightest by ~4×), n≥16: larger again`. The difficulty
spikes exactly at the prime-power `C`. (Honest caveat: the tight sporadic `V*` lives in `C=27`'s
gcd-3 blind-spot stratum, but that mod-27 antipodal structure is INERT for the *binding* reduced
clusters k=8,9,10 — their elements are all `< C/2 ≈ 13.5`, so they carry zero antipodal `{a,27−a}`
relations; the binding danger is in genuine *integer* relations among small elements, and the mod-27
shell structure only activates for the large-element regime like `V*`'s 24. The summand-graph
modulus thus diagnoses *which n is hard* (its factorization) more than it carries the binding tail.) The signed cancellation that
makes `|corr|` diverge while `corr` stays small (MISTAKE-078) is the destructive interference of the
antipodal `±a` shells *within* the gcd-3 stratum — which is codex's reciprocal-tail/residue-lift crux
(HYP-2633) seen from the additive side. The two active LRC(14) threads are addition and multiplication
on the same `C = 3³` graph.

This predicts a recursion: LRC(n) should be governed by `2n−1`'s factorization, with prime-power
`2n−1` the hard cases (single prime, deep strata, no second clock to cross-check). It also explains
the apex prime `7 = n/2` (the sector count) as the *other* modulus — the two clocks of HYP-2101, mod
`7` and mod `27`, are `n/2` and `2n−1`. For LRC(n) generally these are the natural pair, and their
*coprimality and factorization* govern how the additive dimension interacts with the sector cover.

## What is left, and why it is now smaller

The crux is one lemma: **the GAP dimension penalty** — a proper `d≥2` GAP has `L_y ≪ cap`. It carries
a margin of ≥0.21, so a *crude* additive-combinatorial bound suffices, unlike the razor-thin
`consec < cap` at k=9 (margin 0.0014) which is the *exact finite check* and needs no estimate. The
problem has inverted: the part that was "uncontrolled" (the third pocket) is the part with the most
structure and the most room; the part that is tight is finite and already verified. That is the
signature of a problem that has found its right category — here, the Freiman dimension of the speed
set.

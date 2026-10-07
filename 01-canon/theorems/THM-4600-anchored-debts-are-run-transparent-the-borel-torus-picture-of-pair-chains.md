---
id: THM-4600
title: "Anchored debts are run-transparent: if odd u, v satisfy u - c_a = 3^k (v - c_a) with c_a = 1/(2^a - 3) the fixed point of the U-letter a (c_1 = -1, c_2 = +1, c_3 = 1/5, ...), then u and v have U-letter a simultaneously (iff v2(v - c_a) >= a+1) and U(u) - c_a = 3^k (U(v) - c_a); so common runs of any length are invisible to anchored pair-chain states. Anchored states are exactly the affine maps commuting with the run map (a torus in the Borel group), and the translation coordinate relative to the anchor is multiplied by 3/2^a per letter. a = 1 is THM-4555's collision at -1; a = 2 (the trivial cycle) is new and carries the residual first-reset-2 branch (THM-4601)"
status: >
  PROVED (elementary 2-adic algebra; random exact checks of the lemma for a = 1..8, k <= 12, run lengths <= 6 in
  04-computation/experiments/twoanchor_20261007/twoanchor_core.py part A). The general periodic-anchor form (2) is PROVED by the same
  computation; the Borel/torus statement (3) is PROVED; the comparison with the sextic torus-packet paper is ANALOGY.
session: mac-mini-2026-10-07-twoanchor
source: 05-knowledge/results/twoanchor_reset2_friezes_20261007.md
scripts:
  - 04-computation/experiments/twoanchor_20261007/twoanchor_core.py (part A)
related:
  - THM-4555 (uniform switches = collisions at -1: the case a = 1)
  - THM-4581 (the pair chain u = 3^k v + e; Haar absorption)
  - THM-4552/4553 (PSL(2,7) > Borel > torus; the Collatz Frobenius)
  - THM-4601 (application to the residual first-reset-2 branch)
---

# THM-4600 — anchored debts are run-transparent

## Setting

* `U(x) = (3x+1)/2^(v2(3x+1))` on odd x. The *letter* of x is `v2(3x+1)`.
* `F_a(x) = (3x+1)/2^a`.
* For a letter `a ≥ 1` the anchor is `c_a = 1/(2^a − 3)`, the fixed point of F_a: `3c_a + 1 = 2^a c_a`. It is a 2-adic unit.

## Statements

1. **Run transparency (PROVED).** Let u, v be odd integers with `u − c_a = 3^k (v − c_a)`, `k ∈ Z`.
   * (a) v has letter exactly a ⟺ `v2(v − c_a) ≥ a + 1`. The same holds for u, and `v2(u − c_a) = v2(v − c_a)`. So u has letter a iff v does.
   * (b) In that case `U(u) − c_a = 3^k (U(v) − c_a)` and `v2(U(v) − c_a) = v2(v − c_a) − a`.
   * (c) Hence if v has the U-word `a^L`, so does u, the relation persists, and the two runs end at the same step.

2. **Periodic anchors (PROVED).**
   * Let w be a word, `F_w(x) = (3^|w| x + B_w)/2^(Σw)`, and let `c_w = B_w/(2^(Σw) − 3^|w|)` be its rational cycle point.
   * If `u − c = 3^k (v − c)` with c a point of that cycle, and `v2(v − c) > Σw`, then u, v and c share the word w.
   * Moreover `F_w(u) − F_w(c) = 3^k (F_w(v) − F_w(c))`: the anchor moves around the cycle.
   * Integral examples: −1 (w = 1), +1 (w = 2), −5 (w = 12), −17 (w = 1112114).

3. **Borel–torus form (PROVED).**
   * Write the pair-chain state as the affine map `A(x) = 3^k x + e`, with u = A(v).
   * A common letter conjugates A by `τ_a(x) = c_a + λ_a (x − c_a)`, where `λ_a = 3/2^a`.
   * The fixed points of this conjugation are exactly the dilations about c_a: `e = c_a (1 − 3^k)`, the centralizer of τ_a (a torus).
   * In general the translation coordinate `ε = e − c_a(1 − 3^k)` satisfies `ε ↦ λ_a ε`. It is:
     * archimedean-expanding for a = 1;
     * archimedean-contracting for a ≥ 2;
     * 2-adically expanding by `2^a` for every a.

     That is why a non-anchored state leaves a run after about `v2(ε)/a` letters.

## Proof

* `3v + 1 = 3(v − c_a) + 2^a c_a` with c_a a unit. So `2^a | 3v+1` iff `v2(v − c_a) ≥ a`. Given that, `(3v+1)/2^a` is odd iff `v2(v − c_a) ≥ a + 1`.
* `v2(u − c_a) = v2(v − c_a)` because `3^k` is a unit.
* Then `U(u) − c_a = (3u + 1 − 2^a c_a)/2^a = 3(u − c_a)/2^a = 3^(k+1)(v − c_a)/2^a = 3^k (U(v) − c_a)`.
* (2) applies the one-letter computation along the cycle: `F_b(u) − F_b(c) = 3(u − c)/2^b`.
* (3) is the conjugation formula `τ A τ^(−1)(x) = 3^k x + c(1 − 3^k) + λ(e − c(1 − 3^k))`. ∎

## Reading

* **a = 1** is the mechanism of THM-4555. Deleting D ones relates n and `(n+1)/2^D − 1` by `x + 1 = 3^D (y + 1)` at the run ends, an anchored state at −1. That is why the tiling compiler of `collatz_reset2_rules_20261007.md` is uniform in the run length K.
* **a = 2** anchors at the trivial cycle: `u − 1 = 3^k (v − 1)` survives arbitrarily long runs of the letter 2. THM-4601 uses this to make certificates uniform in the second run length of the first-reset-2 branch.
* **Transitions between anchors** are finite word pairs (THM-4601 (iii): `(2,2,2) ~ (4,1,1)` takes "−1, debt 3" to "+1, debt 3").
  * The Terras map preserves positivity of rational limit points.
  * So a positive anchor (+1) can never pass to a negative one (−1, −5, −17) for a positive limit child.
* **ANALOGY** (primitive sextic torus packets, §7):
  * a run is a long diagonal segment;
  * ε is the cross-root coordinate;
  * an anchored state is an inactive cut, the only kind of state that survives arbitrarily long segments.

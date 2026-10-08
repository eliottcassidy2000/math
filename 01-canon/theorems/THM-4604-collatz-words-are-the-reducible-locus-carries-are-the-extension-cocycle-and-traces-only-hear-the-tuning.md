---
id: THM-4604
title: "Collatz words are the reducible locus: the steps x -> (3x+1)/2 and x -> x/2 are upper-triangular matrices, so U-words form a reducible representation G_w = [[3^|w|, B_w],[0, 2^(Sum w)]] whose carry B_w is a twisted 1-cocycle (B_uv = 3^|v| B_u + 2^(Sum u) B_v), trivialized on each cyclic submonoid <w> by the cycle point c_w; every trace function depends only on (|w|, Sum w) and equals 2cosh((|w| ln 3 - Sum w ln 2)/2) after SL2 normalization (the tuning error of HYP-9230); every pair of Collatz words lies on the Cayley cubic x^2 + y^2 + z^2 - xyz = 4 (Fricke invariant tr[A,B] = 2, the reducible characters), not on the Markov/Conway-Coxeter cusp tr[A,B] = -2; merges are cocycle identities, invisible to all trace (Fricke, Markov, frieze) coordinates"
status: >
  PROVED (elementary; the reducible-locus and Fricke facts are standard for representations of free groups into the Borel subgroup).
  The reading - trace-based cluster/frieze coordinates cannot detect merges; only the extension cocycle does; positive minors of carries
  (frieze_collatz_reset2_20261007.md) are minors of the cocycle configuration - is DICTIONARY. Numerical sanity check of the Cayley-cubic
  identity on random word pairs in fricke_check.py.
session: mac-mini-2026-10-07-twoanchor (continuation)
source: 05-knowledge/results/runcompress_orphans_cayley_20261007.md
scripts:
  - 04-computation/experiments/runcompress_20261007/fricke_check.py
related:
  - THM-4600 (anchored states = tori = centralizers in the Borel group), THM-4603 (ladder completeness)
  - THM-4602 (positive friezes and Conway-Coxeter on the irreducible side; Legendre chirotopes)
  - HYP-9230 (musical resonance: the same tuning error |m ln 3 - L ln 2|)
  - 05-knowledge/results/frieze_collatz_reset2_20261007.md (positive marked minors; the ear-move obstruction)
---

# THM-4604 — Collatz words are the reducible locus

## Setting

The odd and even Terras steps act on the projective line as

    O = [[3, 1], [0, 2]]   (x -> (3x+1)/2),     E = [[1, 0], [0, 2]]   (x -> x/2).

* A U-letter a is `G_a = E^(a−1) O = [[3, 1], [0, 2^a]]`.
* Words compose chronologically: `G_(uv) = G_v G_u`.

## Statements

1. **Reducible representation and cocycle.**
   * `G_w = [[3^|w|, B_w], [0, 2^(Σw)]]`, where the carry satisfies

         B_(uv) = 3^|v| B_u + 2^(Σu) B_v.

   * So `w ↦ G_w` is a reducible representation with diagonal characters `χ_1(w) = 3^|w|` and `χ_2(w) = 2^(Σw)`.
   * B is a 1-cocycle twisted by (χ_1, χ_2). Its coboundaries are `c(χ_2 − χ_1)`, i.e. conjugation by a translation.
   * B_w > 0 for every nonempty word.

2. **Cycle points trivialise the cocycle.**
   * For `2^(Σw) ≠ 3^|w|`, on the cyclic submonoid ⟨w⟩ the cocycle is the coboundary of `c_w = B_w/(2^(Σw) − 3^|w|)`: `B_(w^j) = c_w (2^(jΣw) − 3^(j|w|))`.
   * c_w is the fixed point of F_w, i.e. the rational cycle point.
   * The centraliser of G_w in the Borel group is the torus fixing c_w and ∞. These are the anchored states of THM-4600.

3. **Trace blindness.**
   * `tr G_w = 3^|w| + 2^(Σw)`.
   * After SL2 normalisation, `t(w) = tr G_w / √det G_w = 2 cosh(δ_w/2)`, with `δ_w = |w| ln 3 − Σw ln 2`.
   * So every polynomial in traces of words depends only on the pairs (|w|, Σw), i.e. on slopes and the tuning errors δ.
   * Two words with the same counts have the same traces but in general different carries. Collisions (merges) are identities between carries, not between traces.

4. **The Cayley cubic.** For any two words u, v, the triple `(x, y, z) = (t(u), t(v), t(uv))` satisfies

       x² + y² + z² − xyz = 4,

   i.e. the Fricke invariant `tr[ρ(u), ρ(v)] = 2`.
   * This is the reducible locus of the SL2 character variety of the free group, where `x = 2cosh(α/2)`, `y = 2cosh(β/2)`, `z = 2cosh((α+β)/2)`.
   * The Markov / Conway–Coxeter world lives elsewhere: on the punctured-torus cusp `tr[A,B] = −2`, i.e. `x² + y² + z² = xyz` (Markov triples), and on SL2(Z) quiddity products equal to −I.
   * No nonempty Collatz word has a product equal to ±I, since its diagonal is `(3^|w|, 2^(Σw))`. So no Collatz word closes a frieze.
   * Example: the odd/even pair `(O, E)` gives `(5/√6, 3/√2, 7/√12)`, and `25/6 + 9/2 + 49/12 − 105/12 = 4`.

## Proof

1. Matrix multiplication; positivity is by induction on the cocycle identity.
2. `G_w (c_w, 1)^T = 2^(Σw) (c_w, 1)^T`.
3. The diagonal of an upper-triangular product is the product of diagonals.
4. Write `x = 2cosh(a)`, `y = 2cosh(b)`, `z = 2cosh(a+b)` and expand, using `cosh(a+b) = cosh a cosh b + sinh a sinh b`. ∎

## Reading (DICTIONARY)

* **Cluster and frieze machinery is trace-based.**
  * Fricke coordinates, Markov mutations and frieze entries are trace or minor functions on the irreducible locus.
  * On the Collatz (reducible) locus every trace collapses to the characters, which only record slopes and tuning errors.
  * The Markov-type mutation `(x, y, z) ↦ (x, z, xz − y)` still acts on the Cayley cubic. On Christoffel bases it walks the Farey tree of slopes m/L with values `2cosh(δ/2)`: the musical-tuning data of HYP-9230.
* **The dynamics lives in the cocycle.**
  * Merges are carry identities.
  * The positive minors of the tiling session's carry configuration are minors of the cocycle data: `det((3^i, B_i), (3^j, B_j)) = 3^i 2^(A_i) B_(w[i+1..j])`.
  * Anchors are coboundary trivialisations; runs are transparent where the cocycle restricted to the run is trivialised at the state's anchor.
* **Why the frieze route gives coordinates but not rewrites.** This explains the ear-move obstruction (frieze_collatz_reset2_20261007.md §3). Any frieze or cluster identity is a statement on the irreducible side or about traces. Collatz merging is an extension-class phenomenon on the reducible side.

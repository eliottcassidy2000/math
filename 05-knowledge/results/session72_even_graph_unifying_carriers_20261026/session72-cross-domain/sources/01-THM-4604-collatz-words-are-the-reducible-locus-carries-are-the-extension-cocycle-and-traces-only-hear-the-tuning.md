---
id: THM-4604
title: "Collatz words lie on the reducible locus: the steps x -> (3x+1)/2 and x -> x/2 are upper-triangular matrices, so U-words give a faithful reducible representation G_w = [[3^|w|, B_w],[0, 2^(Sum w)]] whose carry B_w is a twisted 1-cocycle (B_uv = 3^|v| B_u + 2^(Sum u) B_v), trivialized on each cyclic submonoid <w> by the cycle point c_w; every trace function depends only on (|w|, Sum w) and equals 2cosh((|w| ln 3 - Sum w ln 2)/2) after SL2 normalization; every pair of words lies on the Cayley cubic x^2 + y^2 + z^2 - xyz = 4 (tr[G_u, G_v] = 2), not on the Markov level tr[A,B] = -2 of the cusped one-holed torus; merges are incidences F_u(n) = F_v(h) governed by the carries, which every trace (Fricke, Markov) coordinate forgets and the frieze minors of the carry configuration record (det(v_i, v_j) = 3^i 2^(A_i) B_(w[i+1..j]); carry exchange relation B_xy B_yz = B_y B_xyz + 3^|y| 2^(Sum y) B_x B_z)"
status: >
  PROVED (elementary). KNOWN in substance: the affine form and c_w = B_w/(2^(Sum w) - 3^|w|) (Boehm-Sontacchi 1978; Lagarias 1985, 1990);
  reducible <=> tr[A,B] = 2 (Goldman 2009, Prop. 2.3.1; Culler-Shalen 1983); the level -2 as the cusped one-holed torus with integer points
  3*(Markov triples) (Cohn 1955; Goldman 2003); the real Cayley cubic and the linear GL(2,Z)-action on its R^2 component (Goldman 2003).
  PROVED one-liners from the cocycle identity: the minor identity, the carry exchange relation, faithfulness of w -> G_w. The reading is
  DICTIONARY. Checks: fricke_check.py; independent audit D (about 1.8*10^5 checks; corrections applied: the reading about friezes, the
  Conway-Coxeter level, faithfulness, KNOWN typing; MISTAKE-587).
session: mac-mini-2026-10-07-twoanchor (continuation)
source: 05-knowledge/results/runcompress_orphans_cayley_20261007.md
scripts:
  - 04-computation/experiments/runcompress_20261007/fricke_check.py
  - 04-computation/experiments/runcompress_20261007/audit_D/thm4604_audit.py (+ .out, REPORT.md)
related:
  - THM-4600 (anchored states = elements of the centraliser torus with multiplier 3^k), THM-4603 (ladder completeness)
  - THM-4602 (positive friezes as a slice of Gr+(2,n), Conway-Coxeter friezes as its integer points; Legendre chirotopes)
  - HYP-9230 (musical resonance: the tuning error)
  - 05-knowledge/results/frieze_collatz_reset2_20261007.md (marked carry coordinates; the ear-move obstruction; example (15))
---

# THM-4604 — Collatz words lie on the reducible locus

## Setting

The odd and even Terras steps act on the projective line as

    O = [[3, 1], [0, 2]]   (x -> (3x+1)/2),     E = [[1, 0], [0, 2]]   (x -> x/2).

* A U-letter a is `G_a = E^(a−1) O = [[3, 1], [0, 2^a]]`.
* Words compose chronologically, `G_(uv) = G_v G_u`. This is an anti-homomorphism for concatenation; `w ↦ G_w^T` is the corresponding homomorphism.
* `A_i = a_1 + ⋯ + a_i` is the prefix exponent sum.

## Statements

1. **Reducible, faithful, with a cocycle.**
   * `G_w = [[3^|w|, B_w], [0, 2^(Σw)]]`, where

         B_(uv) = 3^|v| B_u + 2^(Σu) B_v.

   * So w ↦ G_w is reducible, with diagonal characters `χ_1(w) = 3^|w|` and `χ_2(w) = 2^(Σw)`.
   * B is a 1-cocycle twisted by (χ_1, χ_2). Its coboundaries `c(χ_2 − χ_1)` are conjugation by a translation.
   * B_w is positive, odd and prime to 3 for every nonempty w.
   * B is not a coboundary on the whole monoid (c_1 = −1 ≠ c_2 = +1), so the extension is non-split.
   * w ↦ G_w is injective, even projectively: (|w|, Σw, B_w) determines w, since the prefix sums are successive 2-adic valuations of B_w after its known leading terms are removed.

2. **Cycle points trivialise the cocycle.**
   * For nonempty w (then 2^(Σw) ≠ 3^|w| automatically), on ⟨w⟩ the cocycle is the coboundary of the cycle point `c_w = B_w/(2^(Σw) − 3^|w|)`: `B_(w^j) = c_w (2^(jΣw) − 3^(j|w|))`.
   * The centraliser of G_w in GL2 is the torus fixing c_w and ∞, because the eigenvalues `3^|w| ≠ 2^(Σw)` are distinct.
   * The anchored states of THM-4600 are its elements `x ↦ c_w + 3^k(x − c_w)`.

3. **Trace blindness.**
   * `tr G_w = 3^|w| + 2^(Σw)`, and after SL2 normalisation `t(w) = 2 cosh(δ_w/2)` with `δ_w = |w| ln 3 − Σw ln 2`.
   * So every polynomial in traces of words, including words with inverses, depends only on the pairs (|w|, Σw).
   * δ_w/ln 2 = |w| log2 3 − Σw. When Σw = round(|w| log2 3) this is HYP-9230's θ_f with f = |w|.
   * Distinct words with equal counts have equal traces but always different carries.
   * A merge F_u(n) = F_v(h) is an incidence of two affine maps at integers, governed by the carries. For equal counts it reads `3^|u|(n − h) = B_v − B_u`.
   * Example: 483 and 469 have U-words (1,7) and (7,1), with equal counts and equal traces, and they merge at 17.

4. **The Cayley cubic.** For any two words u, v, the triple `(x, y, z) = (t(u), t(v), t(uv))` satisfies

       x² + y² + z² − xyz = 4,     i.e.   tr[G_u, G_v] = 2   (the commutator is unipotent).

   * This is the reducible locus of the SL2 character variety of the free group.
   * Examples: (O, E) gives `(5/√6, 3/√2, 7/√12)`, and the U-letters 1 and 2 give `(5/√6, 7/√12, 17/√72)`.
   * Markov triples live elsewhere: at tr[A,B] = −2, the cusped one-holed torus, `x² + y² + z² = xyz`, whose positive integer points are 3·(Markov triples).
   * Conway–Coxeter friezes are a different object again: positive quiddities with `M(a_1)⋯M(a_n) = −I`, where `M(a) = [[a, −1], [1, 0]]`. Since `tr[M(a), M(b)] = 2 + (a − b)²`, their transfer pairs are irreducible for a ≠ b and never at −2.
   * No nonempty word has G_w = ±I, or even scalar, since its diagonal is `(3^|w|, 2^(Σw))`. Collatz transfers never close up like a frieze monodromy.
   * Read as quiddities, however, valuation words can close friezes. (1,1,1) and (1,2,2,2,2,2,2,1,7) have M-product −I. They are the actual U-words of 7183 and 2583211, which merge at 24245 (frieze_collatz_reset2_20261007.md, example (15)).

5. **Minors and the carry exchange relation.**
   * For the marked configuration `v_i = (3^i, B_i)` (prefix carriers), `det(v_i, v_j) = 3^i 2^(A_i) B_(w[i+1..j])`.
   * For any three words x, y, z the three-term Plücker relation becomes

         B_(xy) B_(yz) = B_y B_(xyz) + 3^|y| 2^(Σy) B_x B_z.

## Proof

1. Matrix multiplication. Positivity, oddness and coprimality to 3 follow by induction on the cocycle identity. Faithfulness: decode the letters from the 2-adic valuations of B_w − 3^(|w|−1) and its successors.
2. `G_w (c_w, 1)^T = 2^(Σw) (c_w, 1)^T`. An element commuting with a matrix of distinct eigenvalues preserves its eigenlines.
3. The diagonal of an upper-triangular product is the product of the diagonals.
4. Write x = 2cosh a, y = 2cosh b, z = 2cosh(a+b) and expand. The commutator of two upper-triangular matrices is unipotent.
5. Expand `det(v_i, v_j)` with the cocycle identity. The exchange relation is the Plücker relation of v_0, v_|x|, v_|xy|, v_|xyz| for the word xyz, divided by common factors. ∎

## Reading (DICTIONARY; the displayed identities are PROVED)

* **One alphabet, two representations.**
  * The frieze recurrence a ↦ M(a) ∈ SL2(Z) is irreducible as soon as two distinct letters occur. It has relations: the ear identity `M(a)M(b) = M(a+1)M(1)M(b+1)`, and polygon closure.
  * The Collatz recurrence a ↦ G_a is reducible and faithful, so it has no relations.
* **Traces versus minors.**
  * Trace coordinates, such as Fricke traces and Markov-type mutations, are functions on the character variety. That variety records only the semisimplification χ_1 ⊕ χ_2, so on Collatz words they collapse to (|w|, Σw), i.e. to the tuning error.
  * The mutations `(x, y, z) ↦ (x, z, xz − y)`, i.e. `(u, v) ↦ (u, uv)`, and `(x, y, z) ↦ (z, y, yz − x)`, i.e. `(u, v) ↦ (uv, v)`, act on the Cayley cubic. On the Christoffel tree of O/E words grown from (O, E) they walk the Farey tree of slopes with values `2cosh(δ/2)`. This is the linear GL(2,Z)-action on (δ_u, δ_v) of Goldman 2003.
  * Minor coordinates (frieze entries, Plücker and type-A cluster variables) of the carry configuration record the cocycle itself, via statement 5.
* **Why the frieze route gives coordinates but not rewrites.**
  * Coordinates: the minors of the carry configuration are cocycle values. Together with Q_r they are lossless (frieze_collatz_reset2 §2), and the exchange relations are exact carry identities.
  * No rewrites: ear moves are relations of a ↦ M(a), while a ↦ G_a has none. The ear move changes (|w|, Σw) by (+1, +3) and flips the endpoint label mod 3 (frieze_collatz_reset2 §3).
  * A common future of two words is an incidence F_u(n) = F_v(h) at integers. It is decided by the carries and the ternary labels `B_w 2^(−Σw) mod 3^|w|`, which traces do not see.
  * The carry exchange relation holds for every word triple, so on its own it does not detect merges.
* **Anchors.** B restricted to ⟨w⟩ is the coboundary of c_w. So runs of w are transparent for states anchored at c_w, the elements of the centraliser torus with multiplier 3^k (THM-4600).

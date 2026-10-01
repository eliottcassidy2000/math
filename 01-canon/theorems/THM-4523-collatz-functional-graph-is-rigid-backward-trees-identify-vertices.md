---
id: THM-4523
title: "The Collatz functional graph is rigid: the unfolded backward tree of a vertex not divisible by 3 determines it (a 3-adic separation lemma), so the graph of n -> n/2, 3n+1 has no automorphism but the identity, on N, Z and the rationals with denominator prime to 6; the only point where a fibre can be symmetric is the 3-adic number -1/5, which is 2-adically odd"
status: >
  PROVED (elementary) + INDEPENDENTLY AUDITED (2026-10-01, blind re-derivation with
  explicit truncations and an exact child-matching isomorphism search: SOUND; every
  numerical claim reproduced; corrections applied: the fixed point 0 on Z and R_6,
  pruned fibres in the c-family corollary on N); FINITE-EXACT companions (separation
  exponents s(D) for D <= 26; integer collision checks n <= 6000; the
  near-balance law for w < 2*10^5; the shortcut graph to depth 22).
  Setting: T(n) = n/2 (n even), 3n+1 (n odd); B(n) = unfolded backward tree
  (children 2n, and (n-1)/3 when n = 4 mod 6); X ~_D Y = isomorphic
  depth-D truncations (rooted, unordered children).
  (Lemma) For x in Z_3 let e(x) have the children e(2x) and, when
  x = 1 mod 3, o((x-1)/3); let o(x) have the single child e(2x). Then
  B(n) = e(n) (n even), o(n) (n odd). For u in U = 1 + 3Z_3 put A(u) = e(u)
  and let s(D) be the largest k with A(a) ~_D A(b) => a = b mod 3^k on U.
  Then s(D) >= 2 (D >= 4) and s(D) >= 1 + s(D-5) (D >= 9), so
  s(D) >= 2 + floor((D-4)/5) and A is injective on U.
  (Theorem) (i) B(m) and B(n) are non-isomorphic for distinct positive
  integers m, n prime to 3. (ii) Aut(N, T) = 1. (iii) The same on Z and on
  R_6 = Z_(2) cap Z_(3). (iv) A(u) has a non-trivial automorphism iff an even
  vertex of it carries the label -1/5, the unique 3-adic point where the two
  preimage maps 2x and (x-1)/3 coincide (fixed point of x -> 6x+1); -1/5 is
  2-adically odd, so no even rational vertex sits there. (v) For
  w = 4 mod 6 and v = v_3(5w+1) >= 2 the two preimage subtrees of w agree
  exactly to some depth d with 2v - 1 <= d <= 5v - 6.
  Consequences: a functional graph isomorphic to (N, T) has exactly one
  isomorphism to it; (Z, T_{3,c}) = (Z, T_{3,c'}) iff c' = +-c and
  (N, T_{3,c}) = (N, T_{3,c'}) iff c = c' (gcd(cc', 6) = 1; on N the finitely
  many fibres over even n < c with n = c mod 3 are pruned, and the proof uses
  vertices off their forward orbits).
  Extension (note, Theorem R_a, audited): for n/2, an+1 with a an odd prime,
  backward trees separate the vertices prime to a iff 2 is a primitive root mod
  a^2; whole-graph rigidity for every non-Wieferich a is an audit SKETCH.
  Non-consequences: nothing about connectivity (the Collatz conjecture); the
  rigidity is not first-order and gives no grip on a divergent component.
source: opus-2026-10-01-S15 (collatz-functional-uniqueness-20261001), seventeenth note
dependencies: none (elementary 3-adic arithmetic and rooted-tree isomorphism)
---

# THM-4523 -- the Collatz functional graph is rigid

**Status: PROVED (elementary) + INDEPENDENTLY AUDITED (2026-10-01); the audit record is
in the note, section 13.**
Full statement, proof, tables and context:
[`05-knowledge/results/collatz_functional_uniqueness_20261001.md`](../../05-knowledge/results/collatz_functional_uniqueness_20261001.md),
sections 2–5. Script:
`04-computation/experiments/collatz_functional_uniqueness_20261001.py`
(output beside it, `ALL CHECKS PASSED`).

## Statement

Let `T(n) = n/2` for even `n` and `3n + 1` for odd `n`, and let `B(n)` be the
unfolded backward tree of `n` (the children of `n` are `2n` and, when
`n ≡ 4 (mod 6)`, the odd number `(n-1)/3`).

1. **Separation lemma.** With `e`, `o`, `A`, `U` and `s(D)` as in the status
   block, `s(D) ≥ 2 + ⌊(D-4)/5⌋`. Hence `A(a) ≅ A(b)` implies `a = b` for
   `a, b ∈ U`.
2. **Rigidity.**
   * Two distinct positive integers prime to 3 have non-isomorphic backward
     trees.
   * The functional graph `(N, T)` has no automorphism but the identity.
   * The same holds on `Z` and on the rationals with denominator prime to 6
    (there the vertex 0 is the unique fixed point, hence fixed).
3. **The Eckmann–Hilton point.**
   * `A(u)` has a non-trivial automorphism iff an even vertex of it carries
     the label `-1/5`.
   * No even vertex of a rational graph does, since `-1/5` has odd numerator.
   * For `w ≡ 4 (mod 6)` with `v = v_3(5w+1) ≥ 2`, the subtrees `B(2w)` and
     `B((w-1)/3)` agree exactly to a depth `d` with `2v - 1 ≤ d ≤ 5v - 6`.

## Proof (summary)

**Lemma.** The preimage rule is the `e/o` rule. The depth-`d` truncation of
`e(x)` depends on `x mod 3^⌈d/2⌉`.

**Shapes.** `A(u)` has the children `node[A(4u)]` and `Y(u)`, where `Y` is the
bare path, `node[node[A(4z)]]` or `node[A(2z)]` according as `u ≡ 1, 4, 7`
(mod 9), with `z = (u-1)/3`. The three root types are visible at depth 4.

**Separation.** Let `A(a) ~_D A(b)`, so the root types agree.

* Type L: the match is forced straight and gives `a ≡ b` modulo
  `3^(s(D-3)+1)`.
* Type E: a straight match gives `3^(s(D-2)+1)`. A crossed match gives
  `12a ≡ 2b - 2` and `2a - 2 ≡ 12b` modulo `3^(s(D-2)+1)`, and adding gives
  `14(a - b) ≡ 0`.
* Type P: the match is straight into type L, giving `3^(s(D-5)+1)`.

**Rigidity.**

* The shape (first branching depth `0, 1, 2`) and the label `u` are
  invariants, and they return `n = u, u/2, u/4` in `Z_3`, hence `n`.
* An automorphism therefore fixes every `n` prime to 3.
* For `n = 2^j m` with `3 | m` odd, the vertex `w = 3m + 1` is fixed.
  Its preimages are `2w` (fixed) and `m`, so `m` is fixed, and so is the
  chain of unique preimages above it.

**The Eckmann–Hilton point.** `2x = (x-1)/3` iff `x = -1/5`. A moved sibling
pair must be isomorphic, which by separation forces its branch label to be
`-1/5`.

**Near balance.** The labels `4w` and `2(w-1)/3` differ by `2(5w+1)/3`.
Truncations agree while the labels agree to the needed precision (lower
bound), and separation gives the upper bound.

## Consequences and limits

* One labelling: a functional graph isomorphic to `(N, T)` has exactly one
  isomorphism to it.
* The `3n + c` family on `N` is separated by isomorphism type; on `Z` the
  family is separated up to the sign of `c`, by negation.
* The theorem does **not** touch connectivity, i.e. the Collatz conjecture.
  The note's Theorem C makes the conjecture equivalent to `(N, T)` being
  isomorphic to its trivial component. Its Proposition L shows that the
  divergence half is not first-order, so no first-order property of the
  graph can settle it.

---
id: THM-4743
title: "Gluck's trivial set-stabilizers with the converse symmetry added: call a tournament T exceptional if no vertex subset has trivial stabilizer in Aut+-(T), the group of automorphisms and anti-automorphisms. (1) T is exceptional iff it is self-converse and one (equivalently every) anti-automorphism fixes every regular Aut(T)-orbit on subsets. (2) If Aut+-(T) is primitive, T is exceptional iff T is C_3, the circulant Z_5{1,2} or the Paley tournament P_7 (Gluck's Theorem 1 plus |Aut+-| = 2 mod 4). (3) Hence the Paley tournament P_q is exceptional iff q = 3 or 7. (4) Substitution: for A prime and B strong and exceptional, A[B] is exceptional iff Aut+-(A) cannot be distinguished with r(B) colours, r(B) the number of regular Aut(B)-orbits on subsets. (5) r(C_3[...[C_3]]) = 2, so every X[C_3^[k]] (X = C_3, Z_5{1,2}, P_7; k >= 0) is exceptional: an infinite family. (6) FINITE-EXACT: these are the only exceptions among all tournaments with n <= 9 and all circulant tournaments with n <= 15"
status: >
  PROVED ADDENDUM (Session24, independently audited): the family in statement 5 is complete, relative to the cited
  Gluck and Feit–Thompson theorems; see the new structural proof linked below. This is a new necessity argument,
  not retroactive validation of the original synthesis's unsupported completeness sentence.
  PROVED: statements 1-5 (statement 2 rests on Gluck, Can. J. Math. 35 (1983) 59-67, Theorem 1, and on Feit-Thompson
  for the solvability of Aut+-; statement 4 on the modular decomposition of tournaments). FINITE-EXACT: statement 6
  (converse_gluck_tournaments.py (G3), (G5)), the n = 15 substitution checks (G4), the Paley counts (G1), P_27 (G2).
  KNOWN: Gluck's Corollary 1 for tournaments (D(T) <= 2: Albertson-Collins conjecture, Meslem-Sopena arXiv:1608.04866,
  Lozano arXiv:1707.05549). CONCURRENT (session23b, 2026-10-21): every tournament on >= 2 vertices has a regular
  3-colouring under Aut+- (so exceptional = needs exactly 3 colours), a Gluck subset leaves at most one anti-involution,
  P_7's 84 regular subsets, the order-42 failure, and distinguishing numbers 2, 3, 4 for the order-21, -42, -168 actions
  on seven points (05-knowledge/results/signed_tournament_palette_rigidity_20261021/, fano_tournament_move_quotients_20261021/).
  New here: which tournaments need the third colour (statements 1-6). Session opus-2026-10-21-S49; originally pre-audit. Session24 independently verifies the primitive seed and adds full classification.
session: opus-2026-10-21-S49
source: 05-knowledge/results/seven_point_pencil_fano_and_converse_gluck_20261021.md
scripts:
  - 04-computation/experiments/seven_point_spheres_20261021/converse_gluck_tournaments.py (+ .out)
related:
  - HYP-9353 (PROVED by the new independently audited Session24 structural argument)
  - THM-2450 (historical finite-scope grammar); THM-2477 supersedes its frontier with all-order path-rigid atoms {1,C3} and ordered stacks; rigid means H=|Aut|, not Aut=1
  - THM-211 (P_7's fourteen cyclic triples are two Fano planes)
  - 05-knowledge/results/fano_tournament_move_quotients_20261021/session23b-fano/README.md
---

# THM-4743 — Gluck's theorem with the converse added

**Session24 addendum (independently audited).** The [new structural proof](../../05-knowledge/results/converse_exception_classification_20261021/session24-structure/STRUCTURE.md) closes HYP-9353: for every finite tournament on at least two vertices, exceptionality holds iff `T ≅ X[C_3^[k]]`, with `X ∈ {C_3, Z_5{1,2}, P_7}` and `k ≥ 0`. Orbit toggling forces ordinary vertex-transitivity; majority quotients and regular palettes force minimal blocks to be triangles; universal marker coverage forces their full ternary kernel, hence genuine modules. Induction supplies necessity. The original statements and proof history below are retained; the bounded censuses are not all-order premises.

## Setting

* For a tournament `T`, an **anti-automorphism** is an isomorphism `T → T^conv`.
* `Aut±(T)` is the group of automorphisms and anti-automorphisms. A subset `S ⊆ V(T)` is *regular* if its setwise stabilizer in the group is trivial.
* `T` is **exceptional** if `Aut±(T)` has no regular subset.
* `r(T)` is the number of regular `Aut(T)`-orbits on subsets.
* For a group acting on `V(A)`, `D` is the least number of colours admitting a colouring with trivial stabilizer.
* `A[B]` is the lexicographic product: replace every vertex of `A` by a copy of `B`.
* `C_3^[k]` is `C_3[C_3[…]]` with `k` copies; `X[C_3^[0]] = X`.
* `T_5` denotes the circulant `Z_5{1,2}`, the unique vertex-transitive 5-tournament.

## Statements

1. **Orbit pairing (PROVED).**
   * If `T` is not self-converse, `Aut± = Aut(T)` has odd order and `T` is not exceptional (Gluck, Corollary 1).
   * If `T` is self-converse, `Aut(T)` has index 2 in `Aut±(T)`, and any anti-automorphism `σ` acts on the set of regular `Aut(T)`-orbits as an involution.
   * The regular `Aut±`-subsets are exactly the members of the orbits that `σ` moves. Their number is `|Aut(T)|` times the number of moved orbits.
   * So `T` is exceptional iff `σ` fixes every regular `Aut(T)`-orbit.
2. **Primitive classification (PROVED).** If `Aut±(T)` is primitive on `V(T)`, then `T` is exceptional iff `T ≅ C_3`, `T_5` or `P_7`.
3. **Paley tournaments (PROVED).** For every prime power `q ≡ 3 (mod 4)`, `P_q` is exceptional iff `q ∈ {3, 7}`.
4. **Substitution (PROVED).**
   * Let `A` be a prime (indecomposable) tournament with at least 3 vertices, and `B` a strong exceptional tournament. Then `A[B]` is exceptional iff `D(Aut±(A)) > r(B)`.
   * If `A` is not self-converse, `D(Aut±(A)) ≤ 2`, by Gluck.
   * Also `r(A[B]) = N_A(r(B))/|Aut(A)|`, where `N_A(c)` counts the `c`-colourings of `V(A)` with trivial `Aut(A)`-stabilizer.
5. **An infinite family (PROVED).**
   * `r(C_3^[k]) = 2` for all `k ≥ 1`.
   * `D(Aut±(X)) = 3` for `X = C_3, T_5, P_7`.
   * Hence every `X[C_3^[k]]` (`k ≥ 0`) is exceptional, for example `C_3[C_3]` (`n = 9`), `T_5[C_3]` (`n = 15`), `P_7[C_3]` (`n = 21`), `C_3^[3]` (`n = 27`).
   * By contrast `C_3[T_5]` is not exceptional: `r(T_5) = 6 ≥ 3`.
6. **Censuses (FINITE-EXACT).**
   * **All tournaments, `n ≤ 9`.** Among all 191,536 classes with `n = 9`, and all classes with `n < 9`, the exceptional ones are exactly `C_3`, `T_5`, `P_7` and `C_3[C_3]`.
   * **Self-converse classes.** These number 1, 2, 2, 8, 12, 88, 176, 2752 for `n = 2…9`. Every one has `|anti| = |Aut|`, and the orbit-pairing count of 1 holds on each.
   * **Circulants, odd `n ≤ 15`.** Among all `2^((n−1)/2)` circulant tournaments on `Z_n`, the exceptional ones are, up to isomorphism, exactly `C_3`, `T_5`, `P_7`, `C_3[C_3] = Z_9{1,3,4,7}` and `T_5[C_3]`.
   * **`n = 15` directly.** `T_5[C_3]` has no regular `Aut±`-subset, while `C_3[T_5]` has 15,000.

## Proofs

**1.**
* An anti-automorphism composed with an anti-automorphism is an automorphism, so `Aut` has index 2 and is normal. Hence `σ` permutes the regular `Aut`-orbits, and `σ²` acts trivially on them.
* Suppose `S` is `Aut`-regular and `σ(O_S) ≠ O_S`. If an anti-automorphism `τ` fixed `S`, then `τ = gσ` for some `g ∈ Aut`, so `σ(S) = g^(−1)S ∈ O_S`, a contradiction. So `S` is `Aut±`-regular.
* Conversely, an `Aut±`-regular `S` is `Aut`-regular. If `σ(O_S) = O_S`, say `σS = gS`, then `g^(−1)σ` fixes `S`, a contradiction.
* Each moved pair of `Aut`-orbits forms one `Aut±`-orbit of size `2|Aut|`. ∎

**2.**
* `Aut(T)` has odd order: an involution would reverse an arc. So `Aut(T)` is solvable (Feit–Thompson), and so is `Aut±(T)`, an extension of `Aut(T)` by `Z_2`.
* If `T` is exceptional, then by 1 it is self-converse, so `|Aut±| = 2|Aut| ≡ 2 (mod 4)`.
* Gluck's Theorem 1 lists the primitive solvable groups with no regular subset, by degree and point-stabilizer order, which gives the group orders:

  | degree | point-stabilizer orders | group orders |
  |---|---|---|
  | 3 | 2 | 6 |
  | 4 | 3, 6 | 12, 24 |
  | 5 | 2, 4 | 10, 20 |
  | 7 | 6 | 42 |
  | 8 | 21 | 168 |
  | 9 | 16, 24, 48, 8 (`D_8`) | 144, 216, 432, 72 |

  Only 6, 10 and 42 are `≡ 2 (mod 4)`. So `Aut±(T) ∈ {S_3, D_10, F_42}` on 3, 5 or 7 points, and `|Aut(T)| ∈ {3, 5, 21}`.
* The normal subgroup `Aut(T)` of a primitive group of prime degree is transitive, since its orbits form blocks. So `T` is vertex-transitive with `Aut(T) ⊇ Z_p`, i.e. a circulant:
  * on `Z_3`, `T = C_3`;
  * on `Z_5`, every circulant tournament is isomorphic to `Z_5{1,2}` by a multiplier;
  * on `Z_7`, `Aut(T) = F_21` contains multiplication by 2, so the connection set is `{1,2,4}` or `{3,5,6}`, giving `P_7` or its converse, which is isomorphic to `P_7`.
* Conversely these three are exceptional: Gluck's Lemma 6, the orders 6, 10, 42, and (G1)/(G3). ∎

**3.**
* `Aut±(P_q)` contains `{x ↦ ax + b : a ∈ F_q^*}`. Squares give automorphisms and non-squares give anti-automorphisms, since `−1` is a non-square. This group is 2-transitive, so `Aut±(P_q)` is primitive.
* Apply 2: `P_q` is one of `C_3, T_5, P_7` iff `q ∈ {3, 7}`.
* (G1) checks the regular-orbit counts by a closed form for primes `q ≤ 83`, and by brute force for `q ≤ 19`. (G2) gives an explicit regular subset for `P_27`, following Gluck's Lemma 5. ∎

**4.**
* **The block structure.** `A` is prime with `|A| ≥ 3` and `B` is strong, so the copies of `B` are the maximal strong modules of `A[B]` (modular decomposition). Modules are the same for `T` and `T^conv`, so every element of `Aut±(A[B])` permutes the copies. It induces some `π ∈ Aut±(A)`, and on the copies it acts by `β_a ∈ Aut±(B)` of the same type as `π`. All such sign-matched tuples occur.
* **Reduction to colourings.** Let `S = ∪_a {a} × S_a`.
  * If some `S_a` is not `Aut(B)`-regular, a nonidentity automorphism `β` of `B` fixes it, and `(id; β at a)` fixes `S`.
  * Otherwise colour `a` by the `Aut(B)`-orbit of `S_a`, using at most `r(B)` colours. Each of these orbits is fixed by every anti-automorphism of `B`, since `B` is exceptional.
* **If `D(Aut±(A)) > r(B)`.** Some nonidentity `π` preserves the colouring. Choose `β_a` of `π`'s type with `β_a(S_a) = S_(π(a))`. For an automorphism type, use an element of `Aut(B)`. For an anti type, compose an anti-automorphism fixing `S_a` with such an element. Then `(π; β)` fixes `S`, so `A[B]` is exceptional.
* **If `D(Aut±(A)) ≤ r(B)`.** Choose `S_a` in distinct orbits according to a distinguishing colouring, and suppose `(π; β)` fixes `S`.
  * If `π ≠ id`, then `π` preserves the colouring, a contradiction.
  * If `π = id`, then `π` is an automorphism, so every `β_a` is an automorphism fixing an `Aut(B)`-regular set, hence the identity.
  * So `S` is regular.
* **The formula for `r(A[B])`.** The same argument for `Aut` alone shows that the regular subsets are exactly those with every `S_a` regular and the colouring `Aut(A)`-regular. Each colouring is realised by `|Aut(B)|^|A|` subsets, and `|Aut(A[B])| = |Aut(B)|^|A| |Aut(A)|`. ∎

**5.**
* **`r(C_3^[k]) = 2`.** `r(C_3) = 2`: the singletons and the pairs. `C_3` is prime and each tower is strong. `Aut(C_3) = Z_3` fixes only the constant 2-colourings, so `N_(C_3)(2) = 2³ − 2 = 6` and `r(C_3[B]) = 6/3 = 2` whenever `r(B) = 2`. Induct on `k`.
* **`D = 3`.** Two colours fail because `C_3`, `T_5` and `P_7` are exceptional. Three colours suffice:
  * all distinct on `C_3`;
  * `{0}, {1}, rest` on `Z_5`, where the only element fixing 0 is `x ↦ −x` and it moves 1;
  * `{0}, {1}, rest` on `Z_7`, where `F_42` is sharply 2-transitive.
* **The family.** Apply 4 with `A = X` and `B = C_3^[k]`. The towers themselves are exceptional by induction, applying 4 with `A = C_3`. ∎

**6.** These are exhaustive computations ((G3)–(G5)). The classes come from the S36 generator, whose counts match A000568. ∎

## What this does and does not say

* At the original checkpoint, statements 1–6 did not classify all exceptional tournaments; HYP-9353 conjectured completeness. The new Session24 proof linked above closes both formerly open cases: nonmodule permutation blocks and unequal canonical substitutions.
* The orbit-pairing lemma reduces every case to the action of a single anti-automorphism on the regular `Aut`-orbits.
* **Relation to THM-2450.** Its own grammar statement is historical and FINITE-EXACT through order ten. Current [THM-2477](THM-2477-rigid-strong-classification-through-sixty-two.md) and the [all-order bypass proof](../../05-knowledge/results/merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md#3-rigid-strong-tournaments-are-completely-classified) prove that path-rigid tournaments are exactly ordered `{1,C3}` stacks, with the self-converse ones their palindromes. Here rigid means `H(T)=|Aut(T)|`, not `Aut(T)=1`. [THM-2454](THM-2454-pure-blue-alphabet-and-center-law.md) separately gives strong pure-blue atoms `{1,C3,T5}`. These are different from exceptional triangle-substitution towers: `C3[C3]` is exceptional but has `H/|Aut|=39`, not one.

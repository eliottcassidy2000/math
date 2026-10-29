---
id: THM-4849
title: "Circulant tournaments are rational Collatz cycles. (A) Dictionary: for a prime q = 2k+1, the connection set S of a circulant tournament on Z_q, read as a q-periodic parity word, is a rational cycle of the Terras map with k odd steps. The cycle's vertices are the cycle points, its denominator divides D_q = 2^q - 3^k, and the Collatz shift is a translation. Multipliers (the isomorphisms, Turner/Babai) are decimations; the converse is word reversal. (B) Two-base-point lemma: a cycle carries at most two tournament base points (for every odd q, by the audit's (B'')), and a tournament cycle carries exactly two iff it is a rotational-tournament (interval) cycle, iff it is reversal-symmetric. So there are 2^k - k tournament cycles against Cat_k cycles of the shape, and every cycle of the shape is a tournament cycle only for q = 3, 5, 7. (C) The five 3x+101 cycles are the two circulant tournaments on Z_7. R7 gives three achiral cycles {19,79,169}, {31,49,97}, {37,53,65}, and Paley gives the chiral pair {23,85,89}, {29,47,121}; the owner's 35/175 pair is two labellings of R7. (D) Attachment law in fixed-point form (session101's centre law): appending the letter v cancels p^e iff the marked point is p^e-close to 1/(2^v-3), which is -1 for v = 1 and the trivial cycle 1 for v = 2. (E) Locking theorem: for p^e | D(wv) with p not dividing 6(2^v - 3), every marking cancels p^e iff every letter is congruent to v modulo ord_(p^e)(2). The interval cycle of R_q is locked under (1) exactly at the primes of gcd(2^m-1, 3^m-1), m = (q+1)/2, which is 5 at q = 7; this is the 35."
status: >
  Independently audited 2026-10-29: SOUND WITH CORRECTIONS, all applied
  (04-computation/experiments/cycle101_tournament_collatz_20261029/audit/REPORT.md, MISTAKE-696). The audit:
  - corrected the scope of the achirality equivalence in (B): it holds for tournament cycles only;
  - closed a step in the proof of (B);
  - proved (B''), which settles HYP-9446 (4);
  - reclassified the 5x+1 identification in (G) as inherited by construction;
  - added three credits.
  PROVED (elementary; mac-mini-2026-10-29-cycle101): (A) given the cited realisation theorem; (B); (B') in both
  directions; (D); (E) in both directions, with the interval corollary; (F), the trace relation and the Gauss-period
  remark; (G). (B'') was PROVED by the independent audit.
  FINITE-EXACT (04-computation/experiments/cycle101_tournament_collatz_20261029):
  - the dictionary tables for q = 3, 5, 7, 11, 13 (every connection set, necklace, decimation orbit, reduced
    denominator and cycle trace);
  - the (3,7) futures for v = 1..6;
  - the fixed-point law on 8055 random prime-power tests (multipliers a = 3, 5, 7);
  - the locking theorem on 26665 (word, letter, prime-power) tests, all agreeing, 153 of them locked;
  - the interval-locking table for every odd prime q < 100;
  - the Gauss-period identity for every prime q = 3 (mod 4) up to 83 (gauss_periods.py; the audit checked q < 300);
  - (B') and (B'') for every pair (S, S + t) at q = 9, 15, 21, 25, 27 (the audit extended this to 33, 35, 39).
  CITED:
  - the realisation of periodic parity vectors by rational cycles (Boehm-Sontacchi 1978; Lagarias, Acta Arith. 56
    (1990) 33-53);
  - isomorphism of circulant digraphs of prime order iff their connection sets are multiplier-related. This is Turner,
    J. Combin. Theory 3 (1967) 136-145, for graphs; the directed case follows by Sylow, since a regular Z_p is a
    Sylow subgroup of the automorphism group, together with Babai's lemma (Acta Math. Acad. Sci. Hungar. 29 (1977)
    329-336);
  - the cycle lemma (Dvoretzky-Motzkin 1947), giving Cat_k necklaces of k ones in length 2k+1.
  PRIORITY ELSEWHERE IN THE REPO, credited where used:
  - the owner's 35/175 test is C3 of collatz_fermat_reframes_20261029. Its words were taken from the packet's own C1
    table of positive 5x+1 cycles;
  - that the five (3,7) necklaces all have denominator 101, so that 3x+b has exactly five (3,7) cycles iff
    101 | b, is C2 of collatz_fermat_reframes_20261029, and follows from the shift criterion of THM-4484. (C) adds
    their identification with the two circulant tournaments on Z_7;
  - the local centre law v_p(d_v) = e - min(e, v_p(C - xi)), the v = 1 fibre {7:1, 25:1, 35:5, 175:8}, the
    rotation witness 175/175/7 and the R7/Paley7 tournament controls are session101_denominator_connections_20261029
    (pushed 11:08-11:11 the same day, read before this was written). (D) is their law rewritten for the cycle point.
  - reversal as an involution on rational cycles is Proposition 5 and the reversal survey of
    collatz_five_mirrors_20260929, with fixed points the cycles whose word is a rotation of its reverse, and with
    genuine reversal pairs. Reversal invariance of the cycle trace is its Proposition 6. (B) and (F) add that on
    tournament cycles reversal is the converse, and that the achiral tournament cycles are exactly the R_q cycles;
  - the Gauss period zeta + zeta^2 + zeta^4 of the trivial cycle's code (the q = 7 case of (F)) is Proposition 1 of
    ramanujan_heegner7_trivial_cycle_20261006;
  - the binary-fraction reading of the trivial cycle as QR_7 is collatz_paley_bridge_20261001;
  - the transport of the positive 5x+1 cycles {1,3}, {13,33,83}, {17,43,27} to the words (1,4), (1,1,5), (1,3,3)
    is the C1 table of collatz_fermat_reframes_20261029. (G) adds their tournament identification and the
    all-letters-odd criterion.
  NEW HERE:
  - the dictionary itself (the bridge note reads the other way round, see (F));
  - the two-base-point lemma and the count 2^k - k;
  - the identification of the 101 family with the circulant tournaments on Z_7;
  - the fixed-point identification of the centres;
  - the locking theorem and its interval corollary;
  - in (G), the tournament identification and the all-letters-odd criterion.
  NEW BY THE AUDIT: (B'') and its count.
session: mac-mini-2026-10-29-cycle101
source: 05-knowledge/results/cycle101_tournament_collatz_dictionary_20261029.md
depends_on:
  - 05-knowledge/results/collatz_fermat_reframes_20261029/README.md
  - 05-knowledge/results/session101_denominator_connections_20261029/continuations/CONTINUATION_STATE_AND_SYMMETRY.md
related:
  - 01-canon/theorems/THM-4848-multiplier-descent-for-run-observers-of-circulant-tournaments-paley-quadratic-norms-and-the-multiplicity-law.md
  - 05-knowledge/results/collatz_paley_bridge_20261001.md
  - 05-knowledge/results/collatz_five_mirrors_20260929.md
  - 05-knowledge/results/session101_denominator_connections_20261029/tournaments/TOURNAMENT_AND_FACE_CONTROLS.md
  - 05-knowledge/hypotheses/HYP-9446-circulant-collatz-dictionary-and-multiplier-descent-open-questions.md
scripts:
  - 04-computation/experiments/cycle101_tournament_collatz_20261029/circulant_collatz.py (dictionary, q = 3..13)
  - 04-computation/experiments/cycle101_tournament_collatz_20261029/family37.py (the (3,7) family and its futures)
  - 04-computation/experiments/cycle101_tournament_collatz_20261029/attachment_lemma.py (fixed-point law)
  - 04-computation/experiments/cycle101_tournament_collatz_20261029/locking_iff.py (locking theorem, interval table)
outputs:
  - 04-computation/experiments/cycle101_tournament_collatz_20261029/out/ (circulant_collatz.out/.json, family37.out/.json, attachment_lemma.out, locking_iff.out)
---

# THM-4849 — circulant tournaments are rational Collatz cycles

## 0. Conventions

- **The map.** `T(x) = x/2` for even `x` and `(3x+1)/2` for odd `x`, on rationals with odd denominator (2-adic
  parity). The accelerated map is `U(x) = (3x+1)/2^v`, with `v = v_2(3x+1)`.
- **Cycles and their formula.** A `U`-cycle with valuation word `w = (v_1, …, v_k)` and `H = Σ v_i` has the marked
  point `x = C(w)/D(w)`, where `C(w) = Σ_j 3^(k−1−j) 2^(s_j)`, `s_j = v_1 + … + v_j` (`s_0 = 0`), and
  `D(w) = 2^H − 3^k`. This is collatz_fermat_reframes' convention.
- **The `T`-word.** As a `T`-cycle the same orbit has length `H`, and its parity word has ones exactly at
  `s_0, …, s_(k−1)`.
- **Necklaces.** Cyclic valuation words of shape `(k, H)` are the same as `k`-subsets of `Z_H` up to translation:
  the set of odd positions.
- **Realisation (CITED).** Every primitive periodic parity word is the parity vector of exactly one rational point
  (Böhm–Sontacchi; Lagarias 1990).
- **The tournament side.** `q = 2k + 1` is prime. A *connection set* is `S ⊂ Z_q` with `S ⊔ (−S) = Z_q^*`. The
  circulant tournament `T(S)` has the arcs `i → j` iff `j − i ∈ S`.

## 1. Statements

**(A) The dictionary.** Let `w_S` be the `q`-periodic word `(w_S)_j = [j ∈ S]`, and `z_S` its rational point.
- *The cycle.* `z_S` lies on a `T`-cycle of exact length `q` with `k` odd points. Its `U`-word is the cyclic gap
  sequence of `S`, of shape `(k, q)`, and its points are `C/D_q` with `D_q = 2^q − 3^k`. Examples:
  `D_3 = 5`, `D_5 = 23`, `D_7 = 101`, `D_11 = 1805 = 5·19²`, `D_13 = 7463 = 17·439`.
- *Vertices.* Label the `q` points of the `T`-cycle by the vertices, so that the point `z_v` has parity vector
  `1_(v+S)`, the out-neighbourhood of `v`. Then:
  - `T(z_v) = z_(v−1)`, so the Collatz shift is the translation automorphism;
  - `z_v` is odd iff `v → 0`, so the `U`-cycle is the in-neighbourhood of a vertex.
- *Multipliers.* For `m ∈ Z_q^*`, `x ↦ mx` is an isomorphism `T(S) → T(mS)`, and `w_(mS)` is the decimation
  `j ↦ (w_S)_(m⁻¹ j)`. By Turner's theorem (graphs) and Babai's (digraphs), circulants of prime order are isomorphic
  only when their connection sets are multiplier-related, and every isomorphism is affine, `x ↦ m′x + b`. So the
  isomorphism classes of circulant tournaments are the decimation orbits of tournament cycles.
- *Converse.* `T(S)^op = T(−S)`, and `w_(−S)` is the reversal of `w_S`, so the converse is word reversal. The
  decimation by `−1` also induces it: every circulant is self-converse as a tournament, but not as a cycle.
- *Tournament base points.* A point `z` of a shape-`(k, q)` `T`-cycle is a *tournament base point* if its parity
  word `u` satisfies `u_0 = 0` and `u_(−j) = 1 − u_j` for `0 < j < q`: read backwards from `z`, the cycle is the
  bit-complement of itself read forwards. A cycle comes from a circulant tournament iff it has a tournament base
  point, and then the connection set is `{j : u_j = 1}`.

**(B) The two-base-point lemma and the count.**
- *The lemma.* A shape-`(k, q)` cycle has at most two tournament base points.
- *Two base points.* A tournament cycle, i.e. one with a base point, has two iff its connection sets are
  `m·{1, …, k}` for some `m`, i.e. it is a cycle of the rotational tournament `R_q`. Equivalently, iff it is
  reversal-symmetric (achiral).
- *Scope.* Cycles without a base point can be achiral from `q = 11` on: `(1,1,1,4,4)` is one. There are 5 such
  necklaces at `q = 11` and 14 at `q = 13` (audit, `a14_achiral.out`).
- *Counts.*
  - There are `2^k − k` tournament cycles: `2^k − 2k` with one base point and `k` (from `R_q`) with two.
  - All shape-`(k, 2k+1)` cycles together number `C(2k+1, k)/(2k+1) = Cat_k`, by the cycle lemma.
  - `Cat_k = 2^k − k` exactly for `k = 1, 2, 3`.
  - So every rational cycle of the shape is a circulant tournament exactly for `q = 3, 5, 7`. **The 101 family is
    the last one.**
- *Cycles at the next primes.*

| `q` | connection sets | tournament cycles | all cycles of the shape | tournament classes |
|---:|---:|---:|---:|---:|
| 11 | 32 | 27 | 42 | 4: decimation orbits of sizes 5 (`R_11`), 10, 10, 2 (Paley) |
| 13 | 64 | 58 | 132 | 6: orbits of sizes 6 (`R_13`), 12, 12, 12, 12, 4 (the `<3>`-cyclotomic class) |

A class whose multiplier stabiliser is `H` contributes `(q−1)/|H|` cycles, except `R_q`, which contributes
`(q−1)/2`.

**(B') Odd composite `q`.** Let `q` be any odd number, `S` a connection set, `t ≠ 0` and `g = gcd(t, q)`. Then
`S + t` is also a connection set iff both of the following hold:
- every coset of `⟨t⟩ = gZ_q` other than `⟨t⟩` lies inside `S` or misses it;
- `S ∩ ⟨t⟩ = −2t·{1, …, (q/g − 1)/2}`.

That is, `T(S)` is the lexicographic product `C[R_(q/g)]` of a circulant tournament `C` on `Z_q/⟨t⟩ ≅ Z_g` with
the rotational tournament on each coset. For prime `q`, `g = 1` and this is (B).
- *Data* (`composite_basepoints.out`): at `q = 9, 15, 21, 25, 27` every pair `(S, S + t)` has this structure.
- *The bound.* By (B'') no cycle has more than two base points, for any odd `q`.

**(B'') At most two base points, for every odd `q`** (PROVED by the independent audit; this settles HYP-9446 (4)).
- *The statement.* For odd `q ≥ 3`, a connection set `S ⊂ Z_q` has at most one `t ≠ 0` with `S + t` a connection
  set.
- *Corollary (count).* The ordered pairs `(S, t)` number `Σ_(n | q, n > 1) φ(n)·2^((q/n − 1)/2)`. Hence

      #tournament cycles = 2^((q−1)/2) − ½ Σ_(n | q, n > 1) φ(n)·2^((q/n − 1)/2).

  Of the doubly-based cycles, `φ(q)/2` are rotational. For prime `q` the count is `2^k − k`.
- *Data.* `q = 9, 15, 21, 25, 27` give 10, 24, 40, 36, 62 pairs and 11, 116, 1004, 4078, 8161 cycles, matching
  `composite_basepoints.out`.

**(C) The 101 family** (FINITE-EXACT, `family37.out`, `circulant_collatz.out`). The two classes on `Z_7` are `R7`
and Paley `P7`.

| necklace | connection sets | tournament | word | 3x+101 cycle | trace | reversal | locked at 5 under (1) | appending (1) gives, per marking |
|---|---|---|---|---|---|---|---|---|
| `(0,1,2)` | `{1,2,3}`, `{4,5,6}` | R7 (interval) | `(1,1,5)` | 19, 79, 169 | 267/101 | self | **yes** | 35, 35, 35 |
| `(0,1,4)` | `{1,4,5}`, `{2,3,6}` | R7, ×3 | `(1,3,3)` | 31, 49, 97 | 177/101 | self | no | 175, 7, 175 |
| `(0,2,4)` | `{1,3,5}`, `{2,4,6}` | R7, ×2 | `(2,2,3)` | 37, 53, 65 | 155/101 | self | no | 175, 25, 175 |
| `(0,1,3)` | `{1,2,4}` = QR | P7 | `(1,2,4)` | 23, 85, 89 | 197/101 | ↔ `(0,1,5)` | no | 175, 175, 35 |
| `(0,1,5)` | `{3,5,6}` = NQR | P7 (converse) | `(1,4,2)` | 29, 47, 121 | 197/101 | ↔ `(0,1,3)` | no | 35, 175, 175 |

So the owner's two cycles, `(1,1,5)` (C = 19) and `(1,3,3)` (C = 31), are the same tournament `R7` in two
labellings, the interval labelling and its ×3 decimation. Equal denominator, equal clocks, equal automorphism group
and equal achirality cannot separate them. What separates them is (E).

**(D) The attachment law in fixed-point form.** Append the letter `v` to the marked word `w`. Then
`C(wv) = 3C(w) + 2^H`, `D(wv) = 2^(H+v) − 3^(k+1)`, and, with `x = C(w)/D(w)`,

    C(wv) = 2^H (1 − (2^v − 3) x) + x · D(wv)        (exact identity).

Let `p^e ∥ D(wv)` with `p ∤ 2·D(w)`.
- If `p ∤ 2^v − 3`, then `min(e, v_p(C(wv))) = min(e, v_p(x − ξ_v))`, where `ξ_v = 1/(2^v − 3)` is the fixed point
  of the appended branch `y ↦ (3y+1)/2^v`. So the reduced denominator has `v_p(d(wv)) = e − min(e, v_p(x − ξ_v))`.
- *The case `p | 2^v − 3` does not occur here.* Then `D(wv) ≡ 3·D(w) (mod p)`, so `p | D(wv)` would force
  `p | D(w)`. Conversely, under the hypotheses of (E) the cycle points are automatically `p`-integral.
- The centres are the one-letter cycles: `ξ_1 = −1`, the cycle `{−1}`; `ξ_2 = 1`, the trivial cycle;
  `ξ_3 = 1/5`, `ξ_4 = 1/13`, `ξ_5 = 1/29`, …

This is session101's local centre law (their centre `ξ = −2^H/3` for the carry `C`), rewritten for the cycle point.
Their v = 1 fibre `{7:1, 25:1, 35:5, 175:8}` now reads as follows:
- the six markings congruent to −1 mod 5 are 19, 29, 49, 79, 89 and 169;
- of these, 49 is also congruent to −1 mod 25, which gives 7;
- the only marking congruent to −1 mod 7 is 53, which gives 25.

**(E) The locking theorem.** Let `w = (v_1, …, v_k)`, let `v ≥ 1`, and let `p^e | D(wv)` with `p ∤ 6(2^v − 3)`. Then
the following are equivalent:
1. every point of the cycle of `w` is `p`-integral and `≡ ξ_v (mod p^e)`;
2. `v_i ≡ v (mod ord_(p^e)(2))` for every letter.

When they hold, every marking cancels `p^e` on appending `v`, so the response at `p^e` does not depend on the
marking.

*Interval corollary.* The `R_q` interval word `(1^(k−1), k+2)` (connection set `{1..k}`) is locked under the
appended letter `1` exactly at the primes dividing `gcd(2^m − 1, 3^m − 1)`, `m = (q+1)/2`.
- *Small cases.*
  - `q = 7`: gcd(15, 80) = 5. This is the owner's 35: `D(w·1) = 2^8 − 3^4 = 175 = 5²·7`, every marking of
    {19, 79, 169} cancels exactly one 5, and none is `≡ −1 (mod 25)`.
  - `q = 11`: 7; `q = 19`: 11; `q = 23`: 5, 7, 13; `q = 31`: 5, 17. The table to `q < 100` is in `locking_iff.out`.
- *The other four (3,7) necklaces are unlocked at every prime.* Their letters are not congruent modulo any
  `ord_p(2) ≥ 3`. The `(1,3,3)` marking 31/101 is ≡ 1 (mod 5), far from −1, so nothing cancels: 175.
- *Agreement with session101.* Their `g_115(v) = 1 + 4[v ≡ 1 mod 4] + 180[v ≡ 25 mod 36]` splits along this line:
  - the `[v ≡ 1 mod 4]` term is the locked prime 5: `ord_5(2) = 4`, and 1 ≡ 1 ≡ 5 (mod 4);
  - the 37-term is marking-dependent: `ord_37(2) = 36`, and the letters are not congruent mod 36.

**(F) Reversal, traces and the other reading.**
- *Chirality.* Reversal swaps the two Paley cycles. A tournament cycle is chiral iff it is not of `R_q` type, by
  (B). The traces agree, `197/101` for both, by five-mirrors Proposition 6. The `T`-trace is
  `2·(U`-trace`) + k` (837/101 for {19, 79, 169}).
- *The bridge note's reading.* collatz_paley_bridge_20261001 reads a parity word as a binary fraction. There the
  ×2 multiplier is the Collatz shift and translations die. Here translations are the shift and multipliers carry one
  cycle to another.
- *The two readings meet at the Paley base points.* With `c_S = Σ_(j∈S) 2^j` (so `Φ(z_S) = −c_S/(2^q − 1)` by the
  bridge note's Proposition 1), the two Paley base points satisfy, for every prime `q ≡ 3 (mod 4)`,

      c_QR + c_NQR ≡ −1,   c_QR · c_NQR ≡ (q+1)/4   (mod 2^q − 1),

  so `c_QR`, `c_NQR` are the two Gauss periods `η = (−1 ± √−q)/2` evaluated at the `q`-th root of unity `2` modulo
  `2^q − 1`. For `q = 7`: `22 + 104 = 126` and `22·104 = 2288 ≡ 2 (mod 127)`.
  The `q = 7` Gauss period of the trivial cycle's code, `ζ + ζ² + ζ⁴ = τ_7 − 1`, is Proposition 1 of
  ramanujan_heegner7_trivial_cycle_20261006.

**(G) The interval marking, and the positive 5x+1 cycles.**
- *The interval marking.* The `R_q` interval word marked at its first letter has
  `C(1^(k−1), k+2) = Σ_j 3^(k−1−j) 2^j = 3^k − 2^k`. So the owner's numbers are
  - `19 = 3^3 − 2^3`;
  - `79 = C(1,5,1) = 3^2 + 3·2 + 2^6`;
  - `169 = C(5,1,1) = 3^2 + 3·2^5 + 2^6`;
  - `101 = 2^7 − 3^3`.

  That `169` is a square is not structural.
- *The positive 5x+1 cycles.* collatz_fermat_reframes C1 transports the three known positive cycles of
  `n ↦ (5n+1)/2^v` by their words: `{1, 3}` is `(1,4)`, `{13, 33, 83}` is `(1,1,5)`, and `{17, 43, 27}` is
  `(1,3,3)`. The shapes are `(2,5)` and `(3,7)`, so by (B) all three are circulant-tournament cycles: `T5`'s
  interval labelling, `R7`'s interval labelling, and its ×3 decimation.
- *What is inherited and what is new.* The owner's C3 pair was taken by collatz_fermat_reframes from its own C1
  table, so it is the pair of positive 5x+1 three-cycles **by construction** (inherited). (G) adds two things:
  - *Both are `R7`.* At `q = 7` three distinct odd letters would sum to at least 9 > 7, so an all-odd word repeats
    a letter. It is then a rotation of its reverse, hence achiral, hence `R7` by (B).
  - *The criterion below,* which says why exactly these two of the five (3,7) words are integral for 5x+1.
- *Which (3,7) words are integral for 5x+1.* On the (3,7) shape, `D_5 = 2^7 − 5^3 = 3`. The 5x+1 cycle of `w` is
  integral iff every letter is odd, i.e. every gap of `S` is odd. This holds for `R7`'s interval and ×3 labellings,
  and fails for `R7`'s ×2 labelling `(2,2,3)` and for both Paley words.

## 2. Proofs

**(A).**
- *The cycle.* The word `w_S` is primitive (`q` prime, `0 < k < q`), so by realisation `z_S` has exact period `q`.
  Its odd positions are `S`, so its `U`-word is the gap sequence of `S`, and `C/D_q` is the cycle formula with
  `H = q`. `D_q > 0` because `4^k·2 > 3^k`.
- *Vertices.* The parity vector of `T(z)` is the shift of that of `z`. Bit `j` of `T(z_v)` is
  `[j + 1 ∈ v + S] = [j ∈ v − 1 + S]`, so `T(z_v) = z_(v−1)`. And `z_v` is odd iff `0 ∈ v + S` iff `−v ∈ S` iff
  `v → 0`.
- *Multipliers.* `j − i ∈ S` iff `m(j − i) ∈ mS`. The bits of `w_(mS)` are `[j ∈ mS] = [m⁻¹ j ∈ S]`.
- *Converse.* `j − i ∈ −S` iff `i − j ∈ S`, and `[j ∈ −S] = [−j ∈ S]` is the reversed word.
- *Base points.* A base point `z` with word `u` gives `S = {j : u_j = 1}`, and the two conditions say exactly
  `0 ∉ S` and `S ⊔ −S = Z_q^*`. Conversely `z_S` is a base point of its cycle. ∎

**(B).**
- *Setting up.* Let `f = 1_S`, so `f(0) = 0` and `f(x) + f(−x) = 1` for `x ≠ 0`. Suppose `S + t` is also a
  connection set, `t ≠ 0`. Using both properties, `f(x − t) = f(x + t)` for `x ∉ {0, −t}`, i.e.
  `f(y) = f(y + 2t)` for `y ∉ {−t, −2t}`.
- *Walking the cycle.* `2t` generates `Z_q`. Walk `y_j = 2tj`. The equalities break only at `y_(q−1) = −2t` and at
  `y_((q−1)/2) = −t`. So `f` is constant on `{y_0, …, y_((q−1)/2)}`, where it is 0 because `f(0) = 0`, and constant
  on the remaining `k` points, where it is 1. Hence `S = −2t·{1, …, k}`.
- *Conversely.* `−mI = mI + mk` with `I = {1, …, k}`, so every `mI` has a second base point.
- *At most two.* Three base points would give `−2tI = −2uI` with `t ≠ u`, so a multiplier `≠ 1` would fix `I`.
  None does: for `2 ≤ m ≤ k`, the element `x = ⌊k/m⌋ + 1 ∈ I` has `k < mx ≤ 2k`; for `k < m ≤ 2k`, already
  `m·1 = m ∉ I`.
- *Chirality.* A tournament cycle is reversal-symmetric iff `−S` is a translate `S + t` (`t ≠ 0`, as `−S ≠ S`).
  That means two base points.
- *Counting.* The `R_q` class has `q − 1` connection sets (trivial stabiliser), paired into `k` cycles. Every other
  connection set is alone on its cycle, which gives `2^k − 2k + k`.
- *The Catalan comparison.* `gcd(k, 2k+1) = 1`, so there are `C(2k+1, k)/(2k+1) = Cat_k` necklaces. `Cat_k = 2^k − k`
  holds for `k ≤ 3`. From `k = 4` on, `Cat_k` grows by a factor ≥ 3 per step, so it stays above. ∎

**(B').**
- *The cosets.* The same identity `f(y) = f(y + 2t)` holds for `y ∉ {−t, −2t}`. Since `q` is odd, `⟨2t⟩ = ⟨t⟩`, and
  both breaks lie in `⟨t⟩`. On every other coset the walk `y ↦ y + 2t` meets no break, so `f` is constant there.
- *Inside `⟨t⟩`.* The walk `y_j = 2tj` meets `−t` at `j = (q/g − 1)/2`, because `t(2j+1) ≡ 0` iff
  `q/g | 2j + 1`. As in (B), this gives `S ∩ ⟨t⟩ = −2t·{1, …, (q/g − 1)/2}`.
- *Conversely.* Under the structure, `f(y + 2t) = f(y)` off the two breaks, `f(−t) = 0` and `f(−2t) = 1`. These are
  exactly the conditions for `S + t` to be a connection set:
  - the antisymmetry at `x ∉ {0, −t}`;
  - the antisymmetry at `x = −t`, which says `f(−2t) + f(0) = 1`;
  - `0 ∉ S + t`. ∎

**(B'')** (independent audit, 2026-10-29).
- *Setting up.* Suppose `S + t` and `S + u` are both connection sets, with `0, t, u` distinct. Put `f = 1_S`,
  `A = ⟨t⟩` and `B = ⟨u⟩`, of odd orders `n_A, n_B ≥ 3`. By (B'), `f` is constant on every coset of `A` other than
  `A`. In the coordinate `j ↦ j·(−2t)`, `f|_A` is the indicator of `I_A = {1, …, (n_A − 1)/2}`. The same holds
  for `B`.
- *Step 1: `A ≠ B`.* Otherwise `u = mt` with `m` a unit mod `n_A` and `m ≢ 1`, and `m·I_A = I_A`. As in (B), that is
  impossible: for `2 ≤ m ≤ (n_A − 1)/2` take `x = ⌊(n_A − 1)/(2m)⌋ + 1`, and for larger `m` take `x = 1`.
- *Step 2: `A ∩ B = 0`.* Let `E = A ∩ B`, of order `d`, and suppose `0 ≠ E ≠ A`.
  - Each `x ∈ A ∖ E` lies outside `B`, so `f` is constant on `x + B ⊇ x + E`. So `f|_A` is constant on every coset
    of `E` in `A` other than `E`.
  - In the coordinate of `A`, `E = (n_A/d)Z` with `3 ≤ d < n_A`. The coset `1 + E` contains `1 ∈ I_A` and also
    `n_A + 1 − n_A/d`. The latter lies strictly between `(n_A − 1)/2` and `n_A`, so it is not in `I_A`, which is a
    contradiction.
  - Hence `E = 0` or `E = A`. Symmetrically `E = 0` or `E = B`, so by Step 1, `E = 0`.
- *Step 3: contradiction.* Take `a ∈ A ∖ 0` and `b ∈ B ∖ 0`. Then `f(a + b) = f(a)`, by constancy on `a + B`
  (`a ∉ B`), and `f(a + b) = f(b)`, by constancy on `b + A`. So `f` is constant on `A ∖ 0`. This contradicts
  `f(a) + f(−a) = 1`.
- *The count.* For each `t` of order `n` there are `φ(n)` such generators. By (B') the set `S` is then determined on
  `⟨t⟩` and free on the `(q/n − 1)/2` pairs `±` of other cosets. Each doubly-based cycle carries exactly the two
  ordered pairs `(S, t)` and `(S + t, −t)`. ∎

**(D).**
- *The identity.* `C(wv) = 3C(w) + 2^H` is the cycle formula with one more letter. Substituting `C = xD(w)` and
  `3^(k+1) = 2^(H+v) − D(wv)` gives the identity.
- *Cancellation.* Modulo `p^e`, where `x` is `p`-integral, `C(wv) ≡ 2^H(1 − (2^v − 3)x)`, with `2^H` a unit. If
  `p ∤ 2^v − 3`, then `1 − (2^v − 3)x = −(2^v − 3)(x − ξ_v)`. ∎

**(E).**
- *(2) ⇒ (1): integrality.* Put `λ = 3·2^(−v)`. Modulo `p^e`, every branch `y ↦ (3y+1)/2^(v_i)` equals
  `β(y) = (3y+1)/2^v`, and `2^H ≡ 2^(kv)`. Then `p^e | D(wv)` gives `λ^(k+1) ≡ 1`. If `p | D(w)` we would have
  `λ^k ≡ 1 (mod p)`, hence `λ ≡ 1`, i.e. `p | 2^v − 3`, which is excluded. So the cycle is `p`-integral.
- *(2) ⇒ (1): the congruence.* `β^k(x_0) = x_0` and `β(y) − ξ_v = λ(y − ξ_v)` give `(λ^k − 1)(x_0 − ξ_v) ≡ 0`. Here
  `λ^k − 1 ≡ (1 − λ)/λ` is a unit mod `p`, so `x_0 ≡ ξ_v (mod p^e)`, and every `x_i = β^i(x_0)` follows.
- *(1) ⇒ (2).* `x_(i+1) = (3x_i + 1)/2^(v_i)` with `x_i ≡ x_(i+1) ≡ ξ_v` gives `(3ξ_v + 1)(2^(−v_i) − 2^(−v)) ≡ 0`.
  Here `3ξ_v + 1 = 2^v ξ_v` is a unit, so `2^(v_i) ≡ 2^v (mod p^e)`.
- *Interval corollary.* With `v = 1`, condition (2) reads `ord_p(2) | k + 1 = m`. Then `p | D' = 2^(2m) − 3^m` iff
  `p | 3^m − 1`. ∎

**(F).**
- *Traces.* The even points between `x_i` and `x_(i+1)` are `2^r x_(i+1)` for `1 ≤ r < v_i`. They sum to
  `(2^(v_i) − 2)x_(i+1) = 3x_i + 1 − 2x_(i+1)`. Summing over `i` gives `Σ_even = Σx + k`.
- *Gauss periods.* `S = QR_q` is a `(q, k, (q−3)/4)` difference set and `N = −S`. In `Z[x]/(x^q − 1)`,
  `c_S(x)c_N(x) = (q−1)/2 + ((q−3)/4)(J − 1)`, with `J = Σ x^t`. Also `c_S + c_N = J − 1`. Evaluate at `x = 2`
  modulo `2^q − 1`, where `J ≡ 0`. ∎

**(G).**
- *The interval marking.* The geometric sum is `3^k − 2^k`.
- *The 5x+1 criterion: the composite map.* For multiplier 5 the cycle map of `w` is
  `g(x) = (5^3 x + C_5(w))/2^7`. Modulo 3 this is `x − C_5(w)`, so the cycle is 3-integral, and hence integral
  since `D_5 = 3`, iff `3 | C_5(w)`.
- *Branch by branch.* Modulo 3 the branch `y ↦ (5y+1)/2^v` is the translation `y ↦ y − 1` for odd `v` and the
  reflection `y ↦ 1 − y` for even `v`.
- *Checking the shape.* Three odd letters give the translation by `−3 ≡ 0`. The five necklaces were checked
  directly: `C_5 = 39, 51` for the two all-odd necklaces, and `61, 43, 67` for the others, all `≡ 1 (mod 3)`. ∎

## 3. Typing of the owner's numbers (no promotion of coincidences)

The four numbers 19, 79, 169, 101 sort as follows.
- **STRUCTURAL.**
  - 19, 79, 169 are the three markings of the interval cycle of `R7`, the unique locked cycle of the family.
  - 101 is `D_7`, and 7 is the last prime `q` at which every cycle of the shape is a tournament.
- **NUMEROLOGY, recorded without promotion.**
  - `169 = 13²` is a Camion pseudoprime (THM-4746, mod-169 law of period 156).
  - `H(R7) = 175 = 2^8 − 3^4` equals the appended denominator. session101 already flagged this.
  - The 79 in Paley's three-run Smith form is a prime splitting in `Q(√−7)` (THM-4848), not the cycle's 79.
  - `19² | D_11`: the P11 Paley cycle has reduced denominator 361. The cycle {5, 7, 17, 31, 35} of `3x+19` is a Z_11
    tournament cycle. These follow from the residues of `C` modulo 5 and modulo `19²` (i.e. modulo `D_11`)
    only. At `q = 11` the reduced denominators are 1805 (22 cycles), 361 (3), 95 (1) and 19 (1).
- The results note has the full table.

# Ogg's fifteen primes along the levels 1 | 2 | 6 (M, B, Fi22), the square-sum graph between 12 and 25, and what {7, 21} can and cannot share with the Monster

**Session:** opus-2026-10-22-ssj, wave 2 (owner prompt in this session), 2026-10-22. It follows the
[supersingular note](supersingular_divisor_shapes_qseries_20261022.md) and reuses its machinery: supersingular points with level structure and genera of Atkin–Lehner quotients.

**Status:** pre-audit at the checkpoint; the audit record is in §8. Labels are per item:
* PROVED (elementary);
* CITED: Aricheta, SIGMA 15 (2019) 007, arXiv:1809.07421 (Theorems 1.1, 1.3 and Appendix A); Nakaya (as cited there); Ogg 1975; the ATLAS for group orders and centralizers; Gerbicz for square-sum Hamiltonicity at n ≥ 25;
* FINITE-EXACT: scripts in [`04-computation/experiments/supersingular_fus_qseries_20261022/`](../../04-computation/experiments/supersingular_fus_qseries_20261022/README.md);
* ANALOGY or NUMEROLOGY, each with a hostile test.

No HYP was opened.

**The owner's prompt, in short.**
* Analogies between {7, 21} and the baby monster and the monster.
* The number 2^5·3^7·5^3·7^4·11·13^2.
* The fifteen supersingular primes split as {17,19,23,31,47}, {29,41,59,71}, {2,3,5,7,11,13}.
* The cycle 4 → 2 → 5 → 3 → 7 → 13 = 4.
* The square-sum graph: its components up to 13; the 3-4-5 triangle at 12; the new-edge counts per vertex; the picture at 25 with 16, 17, 18 attached; the "double-connection" vertices {15, 23, 25} against {13, 14, 19, 20, 21, 22, 24}.

## 0. The answers

1. **The three sets are a level filtration of Ogg's primes, not a coincidence.**
   * Together with the primes dividing the level, the primes with Ogg's rationality property are:
     * at level 1 (X_0(1)): all fifteen, the primes of M (Ogg);
     * at level 2 (X_0^+(2)): the primes of the baby monster B;
     * at level 6 (X_0^+(6)): exactly {2, 3, 5, 7, 11, 13}, the primes of Fi22.
   * The owner's sets are the successive differences: M ∖ B = {29, 41, 59, 71}, B ∖ Fi22 = {17, 19, 23, 31, 47}, and Fi22 = {2, …, 13}.
   * This is Aricheta's theorem (CITED). It is re-verified here at all thirteen levels N ∈ {1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 21, 22} by an independent count (FINITE-EXACT, §2).
   * Group theory reads the same chain as M ⊃ C_M(2A) = 2.B ⊃ C_M(6A). The last contains Fi22 (ATLAS, CITED).
2. **The number is the Fi22-part of [M : B].**
   * Exactly, |M|/|B| = 2^5·3^7·5^3·7^4·11·13^2 · 29·41·59·71, and the number in the prompt is 39,046,339,332,000, the first factor.
   * The B ∖ Fi22 primes cancel in the index; the M ∖ B primes survive to the first power.
   * The cycle 4 → 2 → 5 → 3 → 7 → 13 = 4 is the map p ↦ v_p([M : B]) on {2, 3, 5, 7, 13}, with 4 read as 13 (the exponent of 7 is 4).
   * Typed NUMEROLOGY (§3). It needs a post-hoc identification and has no mechanism. The same reading fails at once for every other centralizer index checked.
3. **Square-sum graph.** Every stated fact checks (§4). Three of them have exact explanations:
   * **(a) The new-edge counts.** s(n) = ⌊√(2n−1)⌋ − ⌊√n⌋. The runs change exactly where a square enters, at n = ⌈(j^2+1)/2⌉, or leaves, at n = k^2. The two events coincide exactly at n = k^2 with k a Pell number: n = 9, 25, 289, 841, … At n = 25 the square 25 leaves as 49 enters, since 7^2 − 2·5^2 = −1.
   * **(b) Why {15, 23, 25} are "interesting".** If Q_n has a Hamiltonian path and Q_{n−1} has none, then n is interior to every Hamiltonian path of Q_n. So s(n) ≥ 2 (PROVED).
     * Hence every onset is a double-connection vertex. The onsets are exactly 15, 23 and 25.
     * The others, {13, 14, 19, 20, 21, 22, 24}, are double-connection vertices at which no path exists.
   * **(c) 16, 17, 18.** They are the 25-partners of 9, 8, 7 (the reflection x ↦ 25 − x).
     * 16 and 17 append to the two ends 9 and 8 of the unique 15-chain, which is why paths survive at 16 and 17.
     * 18's only partner below 31 is 7, an interior vertex. This is the loss at 18.

   The 3-4-5 component sizes at 12 are exact. Their only structural content is that the three components are welded by the square 25 = 3^2 + 4^2. Typed NUMEROLOGY for the triangle itself; a Pythagorean size triple occurs in no other sum-graph family tested.
4. **{7, 21}.**
   * No map was found from the tournament H-gaps 7 and 21 (THM-1370) to M or B.
   * Three honest bridges exist, each typed (§5):
     * **ANALOGY:** both are finite exceptional sets cut out by a growth law, the strong floor f(m) against the genus of X_0^+(Np).
     * **Moonshine levels 7 and 21:** inside moonshine these behave like levels 1 and 2. The lists are nested: He's primes {2, 3, 5, 7, 17} at level 7, and {2, 3, 5, 7} at level 21. The one prime lost from 7 to 21 is 17, as {29, 41, 59, 71} are lost from 1 to 2.
     * **Inside tournaments:** 7 and 21 are |V| and |Aut| of the Paley tournament QR_7, and F_21 acts freely on its 189 = 9·21 Hamiltonian paths.
   * The shared window (f(6), f(7)) = (15, 25) ∋ 21 coincides with the square-sum period [15, 25]. Typed NUMEROLOGY: the next windows (25, 45) and (45, 75) carry no square-sum event.

## 1. The three sets and the number (exact; group orders CITED)

The orders come from the ATLAS. `level_ogg_lists.gp` factors the decimal |M| and |B|, and computes the quotient, the valuation classes and the exponent table of §3:
* |M| = 2^46·3^20·5^9·7^6·11^2·13^3·17·19·23·29·31·41·47·59·71;
* |B| = 2^41·3^13·5^6·7^2·11·13·17·19·23·31·47.

So:
```text
|M| / |B| = 2^5 · 3^7 · 5^3 · 7^4 · 11 · 13^2 · 29 · 41 · 59 · 71,
2^5 · 3^7 · 5^3 · 7^4 · 11 · 13^2 = 39,046,339,332,000.
```

| set | v_p(\|M\|) | v_p(\|B\|) | in [M : B] | level description (§2) |
|---|---|---|---|---|
| {2, 3, 5, 7, 11, 13} | ≥ 2 | ≥ 1 | exponents 5, 7, 3, 4, 1, 2 | primes of Fi22 (level 6) |
| {17, 19, 23, 31, 47} | 1 | 1 | cancel | B but not Fi22 |
| {29, 41, 59, 71} | 1 | 0 | exponent 1 | M but not B |

Equivalently, {2, …, 13} are the primes whose Sylow subgroup in M is not cyclic (p^2 | |M|).

**Element orders.** M has an element of order 2p (p odd) iff p | |B|. Proof: an element of order 2p is a commuting product t·x of an involution and an element of order p, so p divides |C_M(t)|.
* The involution centralizers are 2.B (2A) and 2^{1+24}.Co1 (2B) (CITED).
* The primes of Co1, {2, 3, 5, 7, 11, 13, 23}, lie inside those of B.
* Hence {29, 41, 59, 71} are exactly the Monster primes p with no element of order 2p.

## 2. The level filtration (CITED: Aricheta 2019; FINITE-EXACT re-verification)

**Definitions.**
* X_0^+(N) is the quotient of X_0(N) by all its Atkin–Lehner involutions.
* A prime p ∤ N has the *rationality property* for X if every supersingular point of X mod p is defined over F_p.

**Aricheta's Theorem 1.1.** For X = X_0(N)^{+e,f,…} and p ∤ N,
* ½ Q_p(X) = genus(X^p) − genus(X), where Q_p counts the non-F_p supersingular points;
* so for genus-zero X_0^+(N), rationality holds iff X_0^+(Np) has genus 0.

For N = 1 this is Ogg's theorem. **Theorem 1.3:** for monstrous X, p has the property iff the centralizer of an element of the class C(X) contains a Fricke element (class pA) of order p. Nakaya found the B and Fi24' cases independently; Aricheta's appendix confirms Nakaya's HN and He conjectures.

**This session's independent count** (`level_ogg_lists.gp`).
* For squarefree level M, the W-invariant subspace of S_2(Γ_0(M)) gets exactly one vector from the oldspace of each newform whose Atkin–Lehner signs are all +1. So genus(X_0^+(M)) = Σ_{d | M} #{newforms of level d with all signs +1}.
* This is checked against PARI's full-space Atkin–Lehner traces at all 91 squarefree levels M ≤ 150, with 0 mismatches.

Results (primes p ∤ N):

| N | rationality list | ∪ primes of N | equals the primes of |
|---|---|---|---|
| 1 | all 15 | — | M (Ogg) |
| 2 | 3, 5, 7, 11, 13, 17, 19, 23, 31, 47 | + 2 | B |
| 3 | 2, 5, 7, 11, 13, 17, 23, 29 | + 3 | Fi24' |
| 5 | 2, 3, 7, 11, 19 | + 5 | HN |
| 6 | 5, 7, 11, 13 | + 2, 3 | **Fi22** |
| 7 | 2, 3, 5, 17 | + 7 | He |
| 10 | 3, 7, 11 | + 2, 5 | HS |
| 11 | 2, 3, 5 | + 11 | M12 |
| 13 | 2, 3 | + 13 | L3(3) |
| 14, 15, 21, 22 | {3, 5}, {2, 7}, {2, 5}, {3, 5} | | (not matched) |

* All twelve lists agree with Aricheta's Appendix A.
* The sporadic group of each matched row is the simple group in the ATLAS centralizer of NA (CITED; the 6A and 10A identifications are recollection). The prime-set equalities are exact against the cited orders.

**Where primes drop out.**
* At level 2, genus(X_0^+(2p)) = 1 for p = 29, 41, 59, 71: exactly one newform of level 2p has both Atkin–Lehner signs +1.
* At level 6, genus(X_0^+(6p)) = 1, 1, 1, 2, 3 for p = 17, 19, 23, 31, 47, while genus(X_0^+(2p)) = 0 for all five.

**Link to the first wave.** The first wave found the non-Fricke level-2 list (X_0(2) itself: all supersingular points F_p-rational) for p ≥ 5 to be {5, 7, 11, 23}. It matches Aricheta's row "2−" = {3, 5, 7, 11, 23}. Passing to the Fricke quotient X_0^+(2) enlarges the list to B's primes. The Fricke involution is what turns the level-2 list into the baby monster.

## 3. The exponent "cycle" (NUMEROLOGY)

The exponents of [M : B] at 2, 3, 5, 7, 11, 13 are (5, 7, 3, 4, 1, 2).
* Following prime ↦ exponent gives 2 → 5 → 3 → 7 → 4. With 4 read as 13 (e.g. 13 ≡ 4 mod 9), it continues 13 → 2.
* 11 ↦ 1 sits outside the cycle.

Each exponent is a difference v_p(|M|) − v_p(|B|) of Sylow exponents. No mechanism is known that would make it a permutation of the primes.

Hostile probe: the same reading for the other centralizer indices fails at the first step, because the exponent at 2 is large.

| index | exponents at 2, 3, 5, 7, 11, 13 |
|---|---|
| [M : Fi24'] | (25, 4, 7, 3, 1, 2) |
| [M : Fi22] | (29, 11, 7, 5, 1, 2) |
| [M : Th] | (31, 10, 6, 4, 2, 2) |
| [M : HN] | (32, 14, 3, 5, 1, 3) |
| [M : He] | (36, 17, 7, 3, 2, 3) |
| [M : Co1] | (25, 11, 5, 4, 1, 2) |
| [B : Fi22] | (24, 4, 4, 1, 0, 0) |

* [M : B] is the only index with all six exponents small, because B is the only subgroup this large.
* The recurring tail (…, 1, 2) at 11, 13 is forced by v_11(|M|) = 2 and v_13(|M|) = 3.

## 4. The square-sum graph Q_n (PROVED / FINITE-EXACT; `square_sum_structure.py`)

**Prior work this inherits.**
* [Wave six](collatz_mod6_20260922_w6_square_sum_hamiltonicity.md):
  * components and leaves to n = 40;
  * Hamiltonian paths exactly at n = 1, 15, 16, 17, 23, 25–40;
  * the cut and leaf certificates for 18–22;
  * Q_n connected for all n ≥ 14;
  * Gerbicz's n ≥ 25 theorem, CITED.
* The [third reflection](collatz_mod6_20260922_lead_reflection3.md): 25 is where the three chains merge, where the scar edge 18 + 7 lies, and where the root 16 + 9 = 3^2 + 4^2 of the Pythagorean tree appears.
* [THM-4510](procgen_sumgraph_20260926_reflection_orbits.md): sum graphs as reflection orbits.

**4.1 Components (checks the owner's list exactly).**
* The number of components for n = 1..14 is 1, 2, 2, 3, …, 3, 2, 1: three components for 4 ≤ n ≤ 12, two at 13, connected from 14.
* Q_12 has the components {2, 7, 9}, {4, 5, 11, 12}, {1, 3, 6, 8, 10}, of sizes 3, 4, 5. All three are paths, and they use only the squares 4, 9, 16.
* Vertex 13 welds two of them by 13 + 3 = 16 and 13 + 12 = 25. Vertex 14 welds the third by 14 + 2 = 16 and 14 + 11 = 25.
* So the components of sizes 3, 4, 5 are joined by the square 25 = 3^2 + 4^2 = 5^2. This is the only structural link found to the "3, 4, 5 right triangle".
* Hostile probe, other sum graphs (a + b ∈ S, n < 120):
  * a Pythagorean triple of component sizes occurs only for squares (once, at n = 12);
  * for twice-squares at n = 24, 25 the sizes contain {3, 4, 5} among five components;
  * never for squares ± 1, triangular, pentagonal, cubes or primes.
* Typed NUMEROLOGY for the triangle as such.

**4.2 New edges per vertex (PROVED).** Vertex n gains one edge per square in [n+1, 2n−1] (the diagonal 2n = k^2 is excluded), so
```text
s(n) = floor(sqrt(2n-1)) - floor(sqrt(n)),        s(n) - s(n-1) = [n = ceil((j^2+1)/2)] - [n = k^2].
```
* Both are checked: the count against brute force for n ≤ 5000, the formula for n ≤ 10^6, and the boundary law for n < 3000.
* The runs are (first vertex, last vertex, s):
  * (5–12, 1), (13–15, 2), (16–18, 1), (19–32, 2), (33–35, 3), (36–40, 2), (41–48, 3), …
  * This is the owner's description: up to 12 one edge each (except 1, 2, 4, which gain none); 13–15 two; 16–18 one; 19–32 two; 33–35 three; 36 two.
* The gains (square j^2 enters) are at 3, 5, 9, 13, 19, 25, 33, 41, …; the losses (square k^2 leaves) at 4, 9, 16, 25, 36, …
* A gain and a loss coincide, leaving an invisible run boundary, exactly at n = k^2 with k^2 − 2t^2 = 1 or 2k^2 − m^2 = 1. Checked by brute force to 10^7, these are n = 1, 9, 25, 289, 841, 9801, 28561, 332929, 970225, for k = 1, 3, 5, 17, 29, 99, 169, 577, 985, the √2-convergent denominators and numerators.
* At n = 25 the window [26, 49] loses 25 and gains 49 (7^2 − 2·5^2 = −1). That is why the run 19–32 passes through 25 unchanged.
* Hostile probe for "Pell coincidences mark Hamiltonicity": n = 9 is also a coincidence, and Q_9 is disconnected. Typed NUMEROLOGY for that link; the coincidence set itself is PROVED.
* The pair (41, 29) that produces the next coincidence 841 = 29^2 also sits in the owner's set {29, 41, 59, 71}. Typed NUMEROLOGY: 59 and 71 are not Pell numbers.

**4.3 The onset lemma (PROVED) and the extension mechanism.**
* *Lemma.* If Q_n has a Hamiltonian path P and Q_{n−1} has none, then n is an interior vertex of P, so s(n) ≥ 2.
* *Proof.* If n were an endpoint, deleting it would leave a Hamiltonian path of Q_{n−1}. Every neighbour of n in Q_n is below n. ∎
* Onsets for n ≤ 40 are 15, 23 and 25, each with s = 2.
* Among the double-connection vertices up to 25 (13, 14, 15, 19–25), the onsets are {15, 23, 25}. The rest, {13, 14, 19, 20, 21, 22, 24}, are vertices where no path exists. That is the owner's split.
* *Extension and loss.*
  * The unique 15-chain is 8,1,15,10,6,3,13,12,4,5,11,14,2,7,9, with ends 8 and 9.
  * 16 = 25 − 9 and 17 = 25 − 8 attach to the two ends, so the chain extends: 17,8,…,9,16 is the unique 17-chain.
  * 18 = 25 − 7 has 7 as its only neighbour below 31, and 7 is interior. So no path can absorb 18: the loss at 18.
  * Independent DFS witnesses were found at 15, 16, 17, 23 and 25; DFS confirms there is no path at 18, 19 or 24.

**4.4 The owner's picture at 25.** Induced edges in Q_25:

| set | edges | shape |
|---|---|---|
| {10, 6, 3, 1, 8, 17} | 5 | a path |
| {12, 4, 5, 11} | 3 | a path |
| {2, 7, 18, 9, 16} | 4 | a tree, 7 of degree 3 |

These are exactly the owner's counts. The three sets are the Q_12 components plus the 25-partners 17 (of 8), 16 (of 9) and 18 (of 7).

## 5. {7, 21} and the Monster

**Inheritance.**
* [THM-1370](../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md): no tournament has exactly 7 or 21 Hamiltonian paths. The proof uses multiplicativity over strong components, the insertion lemma, and the strong floor f(m) = 3, 5, 9, 15, 25, 45, 75 for m = 3..9.
* [THM-200](../../01-canon/theorems/THM-200-h7-impossibility.md) and [THM-079](../../01-canon/theorems/THM-079-h21-component-reduction.md).
* The [Camion/Busch note](camion_busch_gaps_polyhedra_collatz_20261002.md): 7 ∈ (f(4), f(5)) = (5, 9) and 21 ∈ (f(6), f(7)) = (15, 25).

| # | analogy | map / mechanism | test | verdict |
|---|---|---|---|---|
| A | {7, 21} and Ogg's or B's lists are both finite exceptional sets forced by a growth law: the strong floor f(m) increasing past 21; the genus of X_0^+(Np) growing like Np/24 | none between the sides; the shared shape is the repo's "two growth laws" principle (reflection3 §4) | — | ANALOGY (methodological) |
| B | 7 → 21 (multiply by 3, the H of the 3-cycle) as level 1 → 2 (multiply by 2, the order of 2A) | tournaments: H(T) is multiplicative over strong components, 21 = 3·7; Monster: C_M(2A) = 2.B | — | ANALOGY (no map) |
| C | 7 and 21 as moonshine levels | level-7 list ∪ {7} = He's primes {2, 3, 5, 7, 17}; level-21 list ∪ {3, 7} = {2, 3, 5, 7}. The nesting is the Fricke-centralizer nesting (Aricheta Thm 1.3, with C_M(21A) inside C_M(7A) and C_M(3A); ATLAS power maps, recollection). The prime lost from 7 to 21 is 17, as M ∖ B is lost from 1 to 2 | the lists computed (§2) | STRUCTURAL inside moonshine; the tournament meaning of 7 and 21 is not used |
| D | 7 = \|V(QR_7)\|, 21 = \|Aut(QR_7)\| = \|F_21\| | the affine maps x ↦ ax + b (a a square mod 7) act freely on the 189 Hamiltonian paths of the Paley tournament, in 9 orbits, so 21 \| H(QR_7) | computed (`square_sum_structure.py` part 6) | PROVED (tournament side only) |
| E | the gap 21 lies in (15, 25), the square-sum period [15, 25] | none | next windows (25, 45) and (45, 75): no square-sum event (cycle onset is 32, permanence already from 25) | NUMEROLOGY |
| F | prime-set sizes 6, 5, 4 (Fi22, B ∖ Fi22, M ∖ B) against component sizes 3, 4, 5 | none | other level chains give other patterns: 1 \| 3 \| 6 gives 6, 3, 6; 1 \| 7 \| 21 gives 10, 1, 4; 1 \| 2 \| 10 gives 4, 6, 5 | NUMEROLOGY |

## 6. Connection ledger (summary)

| claim | verdict |
|---|---|
| the three sets = M ∖ B, B ∖ Fi22, Fi22 | STRUCTURAL: Aricheta's theorem (CITED), re-verified |
| the number = Fi22-part of [M : B] | exact |
| no element of order 2p in M for p ∈ {29, 41, 59, 71} | PROVED from CITED centralizers |
| the cycle 2 → 5 → 3 → 7 → 4 ≡ 13 → 2 | NUMEROLOGY |
| s(n) formula, run boundaries, Pell coincidences | PROVED |
| onsets need a double connection | PROVED |
| 16, 17 extend at the ends, 18 hits the interior | FINITE-EXACT + mechanism |
| 3-4-5 component sizes at 12 | exact; NUMEROLOGY beyond the welding square 25 = 3^2 + 4^2 |
| {7, 21} versus M and B | analogies A, B; structural only as moonshine levels (C) or inside tournaments (D); NUMEROLOGY for E, F |

## 7. Not claimed, and next questions

**Not claimed.**
* No novelty for §2, which is Aricheta's and Nakaya's.
* No mechanism for the exponent cycle.
* No map from tournament H-gaps to sporadic groups.
* No Collatz or LRC consequence.
* The 6A/10A–Fi22/HS identification and the 21A power maps are ATLAS recollections. The prime-set equalities are what is checked.

**Next questions.**
1. The owner's sets are the differences of the chain 1 | 2 | 6. Do the chains 1 | 2 | 10 (M, B, HS) and 1 | 3 | 6 (M, Fi24', Fi22) also have an index whose small-prime part is a meaningful number?
2. Is there a level-N reading of the square-sum onsets? None is known; a test would need a map from Q_n to modular data, and none is proposed.

## 8. Reproduction and audit

Run from `04-computation/experiments/supersingular_fus_qseries_20261022/`:
* `gp -q -f level_ogg_lists.gp > level_ogg_lists.out` (about 40 s, 124 MB): the cross-check, the level lists, the comparison with Aricheta's table, the sporadic prime sets, the filtration, the drop-out genera, |M|/|B| and the owner's number, the valuation classes, the exponent table and the chain sizes;
* `python3 -B square_sum_structure.py > square_sum_structure.out` (28 s).

Both end with ALL CHECKS PASSED. Audit: pending.

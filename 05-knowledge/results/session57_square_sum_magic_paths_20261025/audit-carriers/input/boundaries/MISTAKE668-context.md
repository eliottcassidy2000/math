# Ogg's fifteen primes along the levels 1 | 2 | 6 (M, B, Fi22), the square-sum graph between 12 and 25, and what {7, 21} can and cannot share with the Monster

**Session:** opus-2026-10-22-ssj, wave 2 (an owner prompt in this session), 2026-10-22. It follows the [supersingular note](supersingular_divisor_shapes_qseries_20261022.md) and reuses its machinery: supersingular points with level structure and genera of Atkin–Lehner quotients.

**Status:** independently audited: PASS WITH CORRECTIONS ([report](../../04-computation/experiments/supersingular_fus_qseries_20261022/audit_wave2/REPORT.md)). All corrections are applied (MISTAKE-668).

Labels, per item:
* PROVED (elementary).
* CITED:
  * Aricheta, SIGMA 15 (2019) 007, arXiv:1809.07421, Theorems 1.1 and 1.3 and Appendix A;
  * Nakaya, arXiv:1809.02363, levels 2 and 3;
  * Ogg 1975;
  * the ATLAS, for group orders, centralizers and involution classes;
  * THM-4510 and wave six, for the square-sum graph;
  * Gerbicz, for Hamiltonian paths at every n ≥ 25 and cycles at every n ≥ 32.
* FINITE-EXACT: scripts in [`04-computation/experiments/supersingular_fus_qseries_20261022/`](../../04-computation/experiments/supersingular_fus_qseries_20261022/README.md).
* ANALOGY or NUMEROLOGY, with a hostile test wherever one exists. Rows A and B of §5 are analogies with no decisive test.

No HYP was opened.

**The owner's prompt, in short.**
* Analogies between {7, 21} and the baby monster and the monster.
* The number 2^5·3^7·5^3·7^4·11·13^2.
* The supersingular primes split as {17,19,23,31,47}, {29,41,59,71} and {2,3,5,7,11,13}.
* The cycle 4 → 2 → 5 → 3 → 7 → 13 = 4.
* The square-sum graph: components up to 13; the 3-4-5 triangle at 12; new-edge counts per vertex; the picture at 25 with 16, 17, 18 attached; the "double-connection" vertices {15, 23, 25} against {13, 14, 19, 20, 21, 22, 24}.

**Concurrent work.** Session32 worked on the same prompt at the same time ([synthesis](session32_square_monster_20261023/SYNTHESIS.md)). Its results agree with this note wherever the two overlap; §7 consolidates them with credit.

## 0. The answers

1. **The three sets are a level filtration of Ogg's primes, not a coincidence.**
   * A prime p ∤ N has the rationality property for a curve if every supersingular point of it mod p is F_p-rational. Taking the rationality lists together with the primes dividing the level:
     * X_0(1) gives all fifteen primes (Ogg's theorem);
     * X_0^+(2) gives the primes of the baby monster B;
     * X_0^+(6) gives exactly {2, 3, 5, 7, 11, 13}, the primes of Fi22.
   * So the owner's sets are the successive differences: M ∖ B = {29, 41, 59, 71}, B ∖ Fi22 = {17, 19, 23, 31, 47}, and Fi22 = {2, …, 13}.
   * Sources and checks:
     * The lists are Aricheta's (Appendix A; CITED). Nakaya found the levels 2 and 3 independently.
     * The equalities with the prime sets of Fi22, HS, M12 and L3(3) are this note's comparisons against ATLAS orders. The group-theoretic reading goes through Aricheta's Theorem 1.3 plus the ATLAS.
     * All thirteen lists, at N ∈ {1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 21, 22}, are re-verified by three independent routes: here, by a newform count (FINITE-EXACT); and by the audit, both by Riemann–Hurwitz and by direct computation of supersingular points.
   * In centralizer terms, C_M(2A) = 2.B (a double cover of B) contains C_M(6A), and C_M(6A) involves Fi22 (recollection: C_M(6A) ≅ 6.Fi22).
2. **The number is the Fi22-part of |M|/|B|.**
   * Exactly, |M|/|B| = 2^5·3^7·5^3·7^4·11·13^2 · 29·41·59·71. The number in the prompt, 39,046,339,332,000, is the first factor.
   * Since B is not a subgroup of M, the class index is [M : 2.B] = |M|/(2|B|), the number of 2A involutions, and its small-prime core is half the number.
   * The cycle 4 → 2 → 5 → 3 → 7 → 13 = 4 is the map p ↦ v_p(|M|/|B|) on {2, 3, 5, 7, 13}, with 4 read as 13 (the exponent of 7 is 4). Typed NUMEROLOGY (§3).
3. **The square-sum graph.** Every stated fact checks (§4). Three of them have exact explanations.
   * **(a) New-edge counts.**
     * s(n) = ⌊√(2n−1)⌋ − ⌊√n⌋.
     * The runs change exactly where a square enters, at n = ⌈(j^2+1)/2⌉, or leaves, at n = k^2.
     * The two events coincide exactly at n = k^2 with k ∈ OEIS A079496: the odd-index Pell numbers 1, 5, 29, 169, … and the even-index half-companion Pell numbers 3, 17, 99, …. So n = 9, 25, 289, 841, …, with n = 1 degenerate.
     * At n = 25 the square 25 leaves as 49 enters, since 7^2 − 2·5^2 = −1.
   * **(b) Why {15, 23, 25} are "interesting" (PROVED).**
     * If Q_n has a Hamiltonian path and Q_{n−1} has none, then n is interior to every Hamiltonian path of Q_n. Hence s(n) ≥ 2.
     * So every onset is a double-connection vertex, and the onsets are exactly 15, 23 and 25.
     * The others, {13, 14, 19, 20, 21, 22, 24}, are double-connection vertices at which no path exists.
   * **(c) 16, 17, 18.**
     * The repo already proves the structural 3-4-5 content (THM-4510, Corollary A3''). The whole first window {15, 16, 17} is the single tight triple of squares (9, 16, 25) = (3^2, 4^2, 5^2), and it closes at 18. Corollary A3''' extends this to every primitive Pythagorean triple.
     * In endpoint language: 16 = 25 − 9 and 17 = 25 − 8 attach to the two ends of the unique 15-chain. 18 = 25 − 7 attaches to an interior vertex.
   * **The 3-4-5 sizes.**
     * The component sizes 3, 4, 5 at n = 12 are exact, and the components are welded by the squares 16 and 25.
     * Typed NUMEROLOGY: among eight sum-graph families, exactly three components of Pythagorean sizes occur only for the squares at n = 12. Twice-squares contain {3, 4, 5} among five components at n = 24, 25.
4. **{7, 21}.**
   * No map was found from the tournament H-gaps 7 and 21 (THM-1370) to M or B.
   * What does connect them, typed (§5):
     * two analogies with no decisive test;
     * the moonshine levels 7 and 21. Their lists nest (He's primes {2, 3, 5, 7, 17}, then {2, 3, 5, 7}), but the nesting is automatic by a monotonicity lemma. Only the prime-set equalities are structural;
     * the Paley tournament QR_7: |V| = 7 and |Aut| = 21, and F_21 acts freely on its 189 = 9·21 Hamiltonian paths. That freeness holds in every tournament.
   * Session32 found a genuine group action carrying 7 and 21 together: N_M(7A)/He ≅ AGL(1,7) acts as the full orientation and converse symmetry group of QR_7 (§7).
   * The shared window (f(6), f(7)) = (15, 25) ∋ 21, against the square-sum period [15, 25], is NUMEROLOGY.

## 1. The three sets and the number (exact; group orders CITED)

The orders come from the ATLAS. `level_ogg_lists.gp` factors the decimal |M| and |B|, and computes the quotient, the valuation classes and the exponent table of §3.
* |M| = 2^46·3^20·5^9·7^6·11^2·13^3·17·19·23·29·31·41·47·59·71;
* |B| = 2^41·3^13·5^6·7^2·11·13·17·19·23·31·47.

So:
```text
|M| / |B| = 2^5 · 3^7 · 5^3 · 7^4 · 11 · 13^2 · 29 · 41 · 59 · 71,
2^5 · 3^7 · 5^3 · 7^4 · 11 · 13^2 = 39,046,339,332,000,
[M : 2.B] = |M| / (2|B|) = |2A| = 97,239,461,142,009,186,000.
```

| set | v_p(\|M\|) | v_p(\|B\|) | in \|M\|/\|B\| | level description (§2) |
|---|---|---|---|---|
| {2, 3, 5, 7, 11, 13} | ≥ 2 | ≥ 1 | exponents 5, 7, 3, 4, 1, 2 | primes of Fi22 (level 6) |
| {17, 19, 23, 31, 47} | 1 | 1 | cancel | B but not Fi22 |
| {29, 41, 59, 71} | 1 | 0 | exponent 1 | M but not B |

**Sylow subgroups.** For M, the primes with non-cyclic Sylow subgroups are exactly {2, …, 13}. p^2 | |M| alone does not imply this. For 11 and 13 it uses the ATLAS: the Sylow 11-subgroup is 11^2 and the Sylow 13-subgroup is 13^{1+2}.

**Element orders.** M has an element of order 2p (p odd) iff p | |B|.
* M has exactly two classes of involutions, 2A and 2B, with centralizers 2.B and 2^{1+24}.Co1 (ATLAS, CITED).
* An element of order 2p is a commuting product t·x of an involution t and an element x of order p, so p divides |C_M(t)|. Conversely, an element x of order p in C_M(t) gives t·x of order 2p.
* The primes of Co1, {2, 3, 5, 7, 11, 13, 23}, lie inside those of B.
* Hence {29, 41, 59, 71} are exactly the Monster primes p with no element of order 2p.

## 2. The level filtration (CITED: Aricheta 2019; FINITE-EXACT re-verification)

**Definitions.**
* X_0^+(N) is the quotient of X_0(N) by all its Atkin–Lehner involutions.
* A prime p ∤ N has the *rationality property* for X if every supersingular point of X mod p is defined over F_p.

**Aricheta's Theorem 1.1.** For X = X_0(N)^{+e,f,…} and p ∤ N, ½Q_p(X) = genus(X^p) − genus(X), where Q_p counts the non-F_p supersingular points. So for a genus-zero X_0^+(N), rationality holds iff X_0^+(Np) has genus 0. For N = 1 this is Ogg's theorem.

**Theorem 1.3.** For a monstrous X, p has the property iff the centralizer of an element of the class C(X) contains a Fricke element (class pA) of order p. On attribution:
* Aricheta names B and Fi24', and records Nakaya's HN and He conjectures, which his Appendix A confirms.
* Nakaya's abstract covers levels 2 and 3.
* The comparisons with Fi22, HS, M12 and L3(3) below are this note's.

**The genus count used here** (`level_ogg_lists.gp`).
* For squarefree M, the W-invariant subspace of S_2(Γ_0(M)) gets exactly one vector from the oldspace of each newform whose Atkin–Lehner signs are all +1. So genus(X_0^+(M)) = Σ_{d | M} #{newforms of level d with all signs +1}.
* Here: checked against PARI's full-space Atkin–Lehner traces at all 91 squarefree levels M ≤ 150, with 0 mismatches.
* The audit: confirmed by Riemann–Hurwitz with CM fixed points at 237 levels up to 1562.

**Monotonicity lemma** (PROVED, from the count; audit F6). For squarefree N | N' and p ∤ N', genus(X_0^+(N'p)) ≥ genus(X_0^+(Np)). Hence L(N') ⊂ L(N) ∖ π(N').
* This makes the nestings L(6) ⊂ L(2) ⊂ L(1) and L(21) ⊂ L(7) automatic.
* It also justifies restricting the search to Ogg's primes, since genus(X_0^+(p)) > 0 forces genus(X_0^+(Np)) > 0.
* Checked on all 24 divisibility pairs among the thirteen levels.

Results, for primes p ∤ N:

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

* All twelve lists with N > 1 agree with Aricheta's Appendix A.
* The audit re-derived all thirteen lists by Riemann–Hurwitz for every prime p < 2000, and by directly computing the supersingular points on explicit Hauptmoduln of X_0^+(N) for p < 200. In the direct route, the number of non-F_p points is always 2·genus(X_0^+(Np)).
* The groups in the last column are the simple groups involved in the ATLAS centralizers of NA. The prime-set equalities are exact against the cited orders; the 6A and 10A identifications are recollection.

**Where primes drop out.**
* At level 2: genus(X_0^+(2p)) = 1 for p = 29, 41, 59, 71. Exactly one newform of level 2p has both Atkin–Lehner signs +1.
* At level 6: genus(X_0^+(6p)) = 1, 1, 1, 2, 3 for p = 17, 19, 23, 31, 47, while genus(X_0^+(2p)) = 0.
* Session32 located the level-2 defect geometrically. For each of 29, 41, 59, 71 there is one Frobenius pair of non-rational coarse points on X_0^+(2), and the pair is two parallel 2-isogeny edges between rational supersingular j-values (§7).

**Link to the first wave.** The first wave found the non-Fricke level-2 list for p ≥ 5 to be {5, 7, 11, 23}: the primes where every supersingular point of X_0(2) itself is F_p-rational. This matches Aricheta's row "2−" = {3, 5, 7, 11, 23}. Passing to the Fricke quotient X_0^+(2) enlarges the list to B's primes.

## 3. The exponent "cycle" (NUMEROLOGY)

The exponents of |M|/|B| at 2, 3, 5, 7, 11, 13 are (5, 7, 3, 4, 1, 2).
* Following prime ↦ exponent gives 2 → 5 → 3 → 7 → 4. With 4 read as 13, it continues 13 → 2.
* 11 ↦ 1 sits outside the cycle.
* Each exponent is a difference v_p(|M|) − v_p(|B|) of Sylow exponents. No mechanism is known that would make these a permutation of the primes.

The decisive test is Session32's ([valuation quotient](session32_square_monster_20261023/VALUATION_CYCLE_AND_QUOTIENT.md)). The identification "13 = 4" cannot come from reduction mod 9: 11 ≡ 2 (mod 9), but 11 and 2 have different output residues. So the exponent map is not a well-defined map mod 9, and the restricted cycle is not induced by any integer polynomial mod 9.

A weaker probe, computed here: the same reading fails at once for the other order quotients |M|/|X|, because the exponent at 2 is large.

| quotient | exponents at 2, 3, 5, 7, 11, 13 |
|---|---|
| \|M\|/\|Fi24'\| | (25, 4, 7, 3, 1, 2) |
| \|M\|/\|Fi22\| | (29, 11, 7, 5, 1, 2) |
| \|M\|/\|Th\| | (31, 10, 6, 4, 2, 2) |
| \|M\|/\|HN\| | (32, 14, 3, 5, 1, 3) |
| \|M\|/\|He\| | (36, 17, 7, 3, 2, 3) |
| \|M\|/\|Co1\| | (25, 11, 5, 4, 1, 2) |
| \|B\|/\|Fi22\| | (24, 4, 4, 1, 0, 0) |

* |M|/|B| is the only quotient with all six exponents small, because 2.B is the largest maximal subgroup of M. B itself is not a subgroup.
* The recurring tail (…, 1, 2) at 11 and 13 is forced by v_11(|M|) = 2 and v_13(|M|) = 3.

## 4. The square-sum graph Q_n (PROVED / FINITE-EXACT; `square_sum_structure.py`)

**Prior work this inherits.**
* [Wave six](collatz_mod6_20260922_w6_square_sum_hamiltonicity.md):
  * components and leaves up to n = 40;
  * Hamiltonian paths exactly at n = 1, 15, 16, 17, 23, 25–40;
  * cut and leaf certificates for 18–22, and the endpoint-aware proof for 24;
  * Q_n connected for every n ≥ 14;
  * Gerbicz's theorems, CITED.
* The [third reflection](collatz_mod6_20260922_lead_reflection3.md): 25 is where the chains merge, where the scar edge 18 + 7 sits, and where the root 16 + 9 of the Pythagorean tree appears.
* [THM-4510](procgen_sumgraph_20260926_reflection_orbits.md), sum graphs as reflection orbits:
  * Corollary A3'': the first window {15, 16, 17} is the single tight triple (9, 16, 25), closing at 18 when 7 gains a third partner.
  * Corollary A3''': Pythagorean zigzags for every primitive triple.

**4.1 Components.** This checks the owner's list exactly.
* The component counts for n = 1..14 are 1, 2, 2, 3, …, 3, 2, 1: three components for 4 ≤ n ≤ 12, two at 13, and connected from 14 on.
* Q_12 has the components {2, 7, 9}, {4, 5, 11, 12} and {1, 3, 6, 8, 10}, of sizes 3, 4, 5. All three are paths, and they use only the squares 4, 9 and 16.
* The welding uses both 16 and 25:
  * vertex 13 joins {1,3,6,8,10} by 13 + 3 = 16 and {4,5,11,12} by 13 + 12 = 25;
  * vertex 14 joins {2,7,9} by 14 + 2 = 16 and {4,5,11,12} by 14 + 11 = 25.
* The structural 3-4-5 content is THM-4510 A3'' above, not the component sizes.
* Hostile probe, over eight sum-graph families (a + b ∈ S, n < 120):
  * exactly three components of Pythagorean sizes occur only for the squares at n = 12;
  * under the weaker criterion "the sizes contain a Pythagorean triple", twice-squares also hit at n = 24 and 25 ({3, 4, 5} among five components);
  * no hits for squares ± 1, triangular, pentagonal, cubes or primes.
* Typed NUMEROLOGY for the sizes 3, 4, 5.

**4.2 New edges per vertex (PROVED).** Vertex n gains one edge per square in [n+1, 2n−1] (the diagonal 2n = k^2 is excluded), so
```text
s(n) = floor(sqrt(2n-1)) - floor(sqrt(n)),        s(n) - s(n-1) = [n = ceil((j^2+1)/2)] - [n = k^2].
```
* Checks: against brute force for n ≤ 5000; the formula for n ≤ 10^6; the boundary law for n < 3000.
* The runs (first vertex, last vertex, s) are (5–12, 1), (13–15, 2), (16–18, 1), (19–32, 2), (33–35, 3), (36–40, 2), (41–48, 3), …. This is the owner's description, with the exceptions 1, 2 and 4, which gain no edge.
* The gains are at 3, 5, 9, 13, 19, 25, 33, 41, …; the losses at 4, 9, 16, 25, 36, ….
* **Coincidences.** A gain and a loss coincide, leaving an invisible run boundary, exactly at n = k^2 with k^2 − 2t^2 = 1 or 2k^2 − m^2 = 1.
  * These are k ∈ OEIS A079496 = 1, 3, 5, 17, 29, 99, 169, 577, 985, ….
  * In terms of the √2-convergents p/q: k is the denominator q of those with p^2 − 2q^2 = −1 (k = 1, 5, 29, 169, 985), and the numerator p of those with p^2 − 2q^2 = +1 (k = 3, 17, 99, 577).
  * The scan of this condition over k up to √(10^7) gives n = 9, 25, 289, 841, 9801, 28561, 332929, 970225.
  * n = 1 satisfies the condition but is degenerate, since its window [2, 1] is empty.
  * At n = 25 the window [26, 49] loses 25 and gains 49 (7^2 − 2·5^2 = −1). That is why the run 19–32 passes through 25 unchanged.
* Hostile probe for "coincidences mark Hamiltonicity": n = 9 is also a coincidence, and Q_9 is disconnected. Typed NUMEROLOGY for that link; the coincidence set itself is PROVED.
* In the owner's set {29, 41, 59, 71}: 29 is a Pell number, 41 a half-companion Pell number, and 59 and 71 neither. NUMEROLOGY.
* Session32 proves the same coincidence law in the form m^2 − 2k^2 ∈ {−1, −2}, and gives the full degree-two birth list 13–15, 19–32, 36–40, 49–50 ([Pell clocks](session32_square_monster_20261023/PELL_ARRIVAL_CLOCKS.md)).

**4.3 The onset lemma (PROVED) and the endpoint reading of THM-4510 A3''.**
* *Lemma.* If Q_n has a Hamiltonian path P and Q_{n−1} has none, then n is an interior vertex of P, so s(n) ≥ 2.
  * *Proof.* If n were an endpoint, deleting it would leave a Hamiltonian path of Q_{n−1}. Every neighbour of n in Q_n is below n. ∎
  * The onsets for n ≤ 40 are 15, 23 and 25, each with s = 2.
  * Among the double-connection vertices up to 25 (13, 14, 15, 19–25), the onsets are {15, 23, 25}. The rest, {13, 14, 19, 20, 21, 22, 24}, are vertices where no path exists: the owner's split.
* *Extension and loss.* This restates A3'' in endpoint language.
  * The unique 15-chain is 8,1,15,10,6,3,13,12,4,5,11,14,2,7,9, with ends 8 and 9.
  * 16 = 25 − 9 and 17 = 25 − 8 attach to those two ends. So the unique 17-chain is 17,8,…,9,16.
  * 18 = 25 − 7 has 7 as its only neighbour below 31, and 7 is interior. This is the loss at 18.
  * Independent DFS: witnesses at 15, 16, 17, 23 and 25; no Hamiltonian path at 18, 19 or 24.

**4.4 The owner's picture at 25.** Induced edges in Q_25:

| set | edges | shape |
|---|---|---|
| {10, 6, 3, 1, 8, 17} | 5 | a path |
| {12, 4, 5, 11} | 3 | a path |
| {2, 7, 18, 9, 16} | 4 | a tree, 7 of degree 3 |

These are the owner's counts. The three sets are the Q_12 components together with the 25-partners 17 (of 8), 16 (of 9) and 18 (of 7). Session32 shows that contracting these blocks loses Hamiltonian feasibility. It also recovers the ten 25-paths by an exact port-state construction ([port transfer](session32_square_monster_20261023/ports/PORT_TRANSFER.md)).

## 5. {7, 21} and the Monster

**Inheritance.**
* [THM-1370](../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) proves that no tournament has exactly 7 or 21 Hamiltonian paths. Its ingredients:
  * multiplicativity over strong components;
  * the insertion lemma;
  * Moon's (n−1)-subtournament theorem;
  * the strong floor f(m), with decisive values f(5) = 9 and f(7) = 25.
* The floor is f(m) = 3, 5, 9, 15, 25, 45, 75 for m = 3..9. Here f(3..7) are exhaustive (the audit re-enumerated every 7-vertex tournament), f(8) = 45 is exhaustive in THM-1370, and f(9) = 75 is Busch's theorem (cited, not used in the proof).
* The THM-1370 number is shared by two canon files. Cite this one by its path.
* Also [THM-200](../../01-canon/theorems/THM-200-h7-impossibility.md), [THM-079](../../01-canon/theorems/THM-079-h21-component-reduction.md), and the [Camion/Busch note](camion_busch_gaps_polyhedra_collatz_20261002.md): 7 ∈ (f(4), f(5)) = (5, 9) and 21 ∈ (f(6), f(7)) = (15, 25).

| # | analogy | map / mechanism | test | verdict |
|---|---|---|---|---|
| A | {7, 21} and Ogg's or B's lists are finite exceptional sets forced by a growth law: the strong floor f(m) passing 21, against genus(X_0^+(Np)) growing like ψ(N)p/(24·2^{ω(N)}) | no map between the sides. The shared shape is the repo's "two growth laws" principle (reflection3 §4) | no decisive test | ANALOGY (methodological) |
| B | 7 → 21 (times 3, the H of the 3-cycle) against level 1 → 2 (times 2, the order of 2A) | tournaments: H is multiplicative over strong components, 21 = 3·7. Monster: C_M(2A) = 2.B | no decisive test | ANALOGY (no map) |
| C | 7 and 21 as moonshine levels | level-7 list ∪ {7} = He's primes {2, 3, 5, 7, 17}; level-21 list ∪ {3, 7} = {2, 3, 5, 7}. The prime lost from 7 to 21 is 17, as M ∖ B is lost from 1 to 2. The nesting is automatic (monotonicity lemma, §2). Only the equality with He's prime set is structural (Theorem 1.3 plus the ATLAS) | the lists are computed (§2) | STRUCTURAL inside moonshine for the prime-set equality. The tournament meaning of 7 and 21 is not used |
| D | 7 = \|V(QR_7)\|, 21 = \|Aut(QR_7)\| = \|F_21\| | Aut(QR_7) = F_21 acts freely on the 189 Hamiltonian paths, in 9 orbits. Freeness is automatic: an automorphism fixing a directed Hamiltonian path fixes every vertex, so \|Aut(T)\| divides H(T) for every tournament T | computed (`square_sum_structure.py` part 6) | PROVED, but generic. The genuine 7/21 group link is Session32's AGL(1,7) action (§7) |
| E | the gap 21 lies in (15, 25), and the square-sum period is [15, 25] | none | (25, 45) and (45, 75) contain no path onset or loss. (25, 45) does contain the cycle onset at 32, and from n = 31 every vertex has degree at least 2 | NUMEROLOGY (post hoc) |
| F | prime-set sizes against component sizes: (4, 5, 6) for (M ∖ B, B ∖ Fi22, Fi22), against 3, 4, 5 | none | other level chains, in the same order: 1 \| 3 \| 6 gives (6, 3, 6); 1 \| 7 \| 21 gives (10, 1, 4); 1 \| 2 \| 10 gives (4, 6, 5), the same multiset {4, 5, 6} | NUMEROLOGY |

## 6. Connection ledger (summary)

| claim | verdict |
|---|---|
| the three sets = M ∖ B, B ∖ Fi22, Fi22 | STRUCTURAL: Aricheta's lists (CITED), re-verified by three routes. The Fi22 identification is this note's comparison |
| the number = Fi22-part of \|M\|/\|B\|, and half of it is the small-prime core of \|2A\| | exact |
| no element of order 2p in M exactly for p ∈ {29, 41, 59, 71} | PROVED from CITED centralizers and involution classes |
| the cycle 2 → 5 → 3 → 7 → 4 ≡ 13 → 2 | NUMEROLOGY (not a mod-9 map; Session32) |
| s(n) formula, run boundaries, coincidences at k ∈ A079496 | PROVED |
| onsets need a double connection | PROVED |
| the window {15, 16, 17} and its closure at 18 | inherited (THM-4510 A3''); the endpoint reading is FINITE-EXACT |
| 3-4-5 component sizes at 12 | exact; NUMEROLOGY |
| {7, 21} against M and B | analogies A, B; STRUCTURAL only as moonshine levels (C); generic inside tournaments (D); Session32's AGL(1,7) action is the genuine group link; NUMEROLOGY for E, F |

## 7. Consolidation with Session32, and what is not claimed

**Session32** ([synthesis](session32_square_monster_20261023/SYNTHESIS.md)) worked the same prompt independently.

Results that agree with this note, each arrived at independently:
* the valuation classes, and the number as the order quotient of M and B divided by 29·41·59·71;
* the Fi22 / B / M prime-support shells ([Fischer support note](session32_square_monster_20261023/FISCHER_SUPPORT_SHELLS.md)). Its pre-audit caveat about this packet's level filtration is now answered by the audit above;
* the s(n) formula, the Pell coincidences, the component counts and the blocks at 25;
* H(QR_7) = 189, |Aut| = 21, and nine orbits.

Results that Session32 adds:
* **A Sylow fixed-point reading of the three sets** for the action of M on its 2A involutions ([Sylow partition](session32_square_monster_20261023/SYLOW_FIXED_POINT_PARTITION.md)):
  * for {17, 19, 23, 31, 47}, the full Sylow subgroups fix a point;
  * for {2, …, 13}, some order-p subgroup fixes a point, but no full Sylow subgroup does;
  * for {29, 41, 59, 71}, no order-p subgroup fixes a point.
* The mod-9 obstruction for the exponent cycle.
* The level-2 edge defect, two parallel 2-isogeny edges between rational supersingular j-values:

  | p | 29 | 41 | 59 | 71 |
  |---|---|---|---|---|
  | edge endpoints | {2, 25} | {3, 32} | {28, 48} | {41, 66} |
* **The Monster–Paley action** N_M(7A)/He ≅ AGL(1,7), acting as the full orientation and converse group of QR_7 ([Paley bridge](session32_square_monster_20261023/paley/PALEY_MONSTER_BRIDGE.md)).
* The port-state construction at 25, and the cycle thresholds at 31 and 32.

**Not claimed.**
* No novelty for §2, which is Aricheta's and Nakaya's.
* No mechanism for the exponent cycle.
* No map from the tournament H-gaps to sporadic groups.
* No Collatz or LRC consequence.
* The 6A/10A identifications with Fi22/HS and the 21A power maps are ATLAS recollections. The prime-set equalities are what is checked.

**Next questions.**
1. Is there a modular (level) characterization of the Sylow regime {2, …, 13} that matches Session32's fixed-point reading? At present the level-6 list matches the set, and the Sylow reading explains it group-theoretically; no common mechanism is stated.
2. Is there a level-N reading of the square-sum onsets? None is known.

## 8. Reproduction and audit

Run from `04-computation/experiments/supersingular_fus_qseries_20261022/`; both scripts end with ALL CHECKS PASSED.
* `gp -q -f level_ogg_lists.gp > level_ogg_lists.out` (about 40 s, 124 MB; stderr not captured). It runs, in order:
  * the cross-check of the genus count against Atkin–Lehner traces;
  * the level lists and the comparison with Aricheta's table;
  * the sporadic prime sets and the filtration;
  * the drop-out genera;
  * |M|/|B| and the number, the valuation classes, the exponent table and the chain sizes;
  * the monotonicity lemma.
* `python3 -B square_sum_structure.py > square_sum_structure.out` (28 s).

**Audit** ([REPORT](../../04-computation/experiments/supersingular_fus_qseries_20261022/audit_wave2/REPORT.md); its scripts are beside it). The verdict is PASS WITH CORRECTIONS. Corrections applied (MISTAKE-668):
* F1: the hostile probe's claim contradicted its own output; both criteria are now reported.
* F2: THM-4510 A3'' and A3''' are now credited, and the welding by 16 and 25 is stated.
* F3: "Pell number" is replaced by A079496, with n = 1 marked degenerate.
* F4: B is not a subgroup; the note now uses |M|/|B| = 2·|2A| and "involves Fi22".
* F5: row E's test sentence was false (the cycle onset at 32).
* F6: the monotonicity lemma is now stated.
* F7: the growth constant is now ψ(N)p/(24·2^{ω(N)}).
* F8: the {4, 5, 6} repeat in row F is recorded.
* F9: the description of THM-1370 is corrected, and its number collision noted.
* F10: the free action is generic, and this is stated.
* F11: attribution is corrected.
* F12: the Sylow and involution-class steps are added.
* F13: missing hostile tests are acknowledged, and Session32's mod-9 test is used.
* F14: outputs are regenerated without stderr.

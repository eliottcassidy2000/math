# Rédei graphs: which shaved tournaments are forced by parity, and u(9) = 14

**Lane:** shave4 (collatz-procgen-20260922, 2026-10-01). It resumes a lane that died in a reboot, with a new scope.

**Builds on:**
- opus S15's note [shaved_tournaments_unavoidable_cores_20261001](shaved_tournaments_unavoidable_cores_20261001.md) and THM-4526 (audit owed);
- the THM-4524 parity machinery.

**Code (independent of S15's code):**
- `04-computation/experiments/procgen_shave4_20261001_engine.c`: modes `hp`, `til`, `chords`, `dag`, `redei`, `rigid`, `witness`, `unav`, `embedall`, `dcheck`;
- `procgen_shave4_20261001_closure.c`: the deletion–reversal Horn closure, with derivation printing;
- `procgen_shave4_20261001_u9.c`: orderly search over a symmetric host, in shaving and parity modes;
- runner `procgen_shave4_20261001_run.py`, which ends ALL CHECKS PASSED (149 checks). Its output is `05-knowledge/results/procgen_shave4_20261001.out`.

**Raw search records:**
- `05-knowledge/results/procgen_shave4_20261001_u9search.out`: the two exhaustive u(9) runs, with the C3[C3] and |Aut| = 27 hosts, and the K = 14 enumeration of the 54 maximum classes;
- `05-knowledge/results/procgen_shave4_20261001_r9search.out`: the r(9) search, giving existence, non-existence and the complete size-11 list.

**Labels:** PROVED, CITED, FINITE-EXACT, VERIFIED, EMPIRICAL, CONDITIONAL, OPEN, REFUTED.

## Status table

| # | Statement | Label |
|---|---|---|
| A1 | THM-4526(A): H_n avoided only by C3 and C3[1,C3,1]; odd count for even n | VERIFIED independently for all classes n ≤ 9; proof re-checked step by step |
| A2 | Theorem A is classical: the cycle type (n−1,1) in Rosenfeld's problem | CITED: Grünbaum 1971 for existence (as reported by El Zein); the exception list (3A;(2,1)), (5C;(4,1)) is due to Havet 2000 / El Zein 2022 |
| B1 | Deletion–reversal (DR): φ_S + φ_{rev_f S} = φ_{S−f} (mod 2) | PROVED |
| B2 | Completeness: DR plus isomorphism generate *all* F_2-linear relations among the φ_S | PROVED |
| B3 | Free-involution, union (Lucas), Forcade formula, halving, twisted-pair, ideal lemmas | PROVED |
| B4 | φ_S rigid with ≥ 1 arc ⟹ Aut(U(S)) contains an involution; Rédei ⟹ also e(S) odd and \|Aut S\| odd | PROVED |
| B5 | Interval heredity of Rédei "path + chords" graphs; deleting a unique source or sink of a Rédei graph leaves a Rédei graph | PROVED (138 cases checked, n ≤ 8) |
| G1 | Goodness propagation: if every U − f is good, all orientations of U are rigid or none is | PROVED |
| G2 | Every orientation of C_n is parity-rigid iff n is even | PROVED (checked n ≤ 8) |
| T | Every orientation of a hereditarily symmetric tree is parity-rigid; no orientation of an asymmetric tree is | PROVED |
| T' | Trees on ≤ 8 vertices: a tree has a non-rigid orientation iff it contains the spider S(1,2,3) | FINITE-EXACT |
| T'' | "Good ⟺ hereditarily symmetric" for all trees | OPEN (true n ≤ 8) |
| D1 | H_n (n even) and D_n = P_n + {(0,n−2),(1,n−1)} (n odd) are Rédei; for even n, φ(D_n) = 1 + Σ_w hc(T−w) | PROVED; identities checked exactly for every class, n ≤ 9 |
| D2 | The Rédei "path + chords" graphs R_n for n ≤ 11 | FINITE-EXACT (n ≤ 9); n = 10, 11 PROVED (B5 + D1 + verified witnesses) |
| D2' | For n ≥ 8: R_n = {∅, {(0,n−1)}} (n even), {∅, {(0,n−2),(1,n−1)}} (n odd) | CONDITIONAL / OPEN beyond 11 (reduces to finitely many explicit refutations per n) |
| D3 | D69: path + span-3 is Rédei for n ≤ 6 (hand proof); it fails for every n ≥ 7, and the reason is the smallest asymmetric tree | PROVED |
| D4 | Doubling: S Rédei on 2^k vertices ⟹ S ⊔ S + ((x,0)→(x,1)) Rédei on 2^{k+1} | PROVED (20/20 instances checked at n = 8) |
| C1 | Census: Rédei classes 1,1,2,5,21,43,156,220 (n = 1..8); rigid classes for n ≤ 7 | FINITE-EXACT, two independent methods |
| C2 | r(n) = max arcs of a Rédei graph = 0,1,2,4,6,8,9,10 (n = 1..8): r = u for n ≤ 7, r(8) = u(8) − 1 | FINITE-EXACT |
| C3 | **r(9) = 11** (so u − r = 0,…,0,1,3 for n ≤ 9) | PROVED by exhaustive computation (orderly parity search, host C3[C3]; maxima verified independently) |
| E1 | **D70: u(9) = 14, κ(9) = 22**; the excess law holds at n = 9; there are 54 classes of maximum 9-shavings | PROVED by exhaustive computation (two hosts; method reproduces S15 at n = 7, 8); class count FINITE-EXACT |
| E2 | **D71: u(n) = n log₂ n − O(n log log n)**, so c = 1 | CITED (Linial–Saks–Sós 1983) |

## Notation

- T is a tournament on [n].
- S is a spanning oriented graph on [n] with m arcs.
- emb(S,T) = #{bijections π : π(S) ⊆ T}.
- **φ_S(T) := emb(S,T) mod 2.**
- U(S) is the underlying graph, and e(S) = emb(S, TT_n) is the number of linear extensions.
- P_n is the directed path 0→…→n−1. "P_n + C" adds forward chords (i,j) with j ≥ i+2.
- |Aut T| is odd for every tournament. So φ_S(T) ≡ m_S(T), the number of completions of S isomorphic to T (S15's multiplicity).

**Definitions.**
- S is **parity-rigid** if φ_S is constant on n-tournaments.
- S is **Rédei** if φ_S ≡ 1.
- A graph U is **good** if every orientation of U is parity-rigid.
- A Rédei graph is contained in every tournament, i.e. it is an S15 "shaving".

## 1. Theorem A (THM-4526): independent verification, proof audit, literature (A1, A2)

Engine mode `hp` uses its own method:
- count all Hamiltonian paths by subset DP;
- count non-closing paths directly as Σ_{s→t} P(s,t);
- check per class that the closing paths number exactly n·hc(T).

| n | classes | no copy of H_n | odd / even counts |
|---|---|---|---|
| 3 | 2 | 1 (C3 = `101`) | 1 / 1 |
| 4 | 4 | 0 | 4 / 0 |
| 5 | 12 | 1 (`1110101111` = C3[1,C3,1], H = 15, hc = 3) | 8 / 4 |
| 6 | 56 | 0 | 56 / 0 |
| 7 | 456 | 0 | 254 / 202 |
| 8 | 6880 | 0 | 6880 / 0 |
| 9 | 191536 | 0 | 100178 / 91358 |

S15 handled n = 9 through one-vertex extensions of the 8-classes; here `gentourng 9` is used directly.

**Proof audit (S15 §1, odd case).** I re-derived every step:
- Lemma 1 holds: words with no factor `io` are o^j i^{r−j}.
- Step 1 needs n ≥ 4.
- Steps 3–4: the non-closing property of Q_1 and Q_2 is exactly the "backward" condition at positions i and i+1.
- Step 5: an HP of a non-strong tournament ends in its terminal component, which lies in R_k; the second bullet uses Moon's lemma and Camion.

No gap.

**Literature (new for the audit; S15 found no reference).**
- H_n is the oriented Hamiltonian cycle of block type (n−1, 1): one directed block of length n−1 plus one reversed arc.
- In the complete classification of tournaments missing some non-directed Hamiltonian cycle (A. El Zein, arXiv:2204.11211, v2 2023, "Oriented Hamiltonian cycles in tournaments: a proof of Rosenfeld's conjecture"), the exceptional pairs (T; C) of this type are exactly A_1 = (3A; (2,1)) and A_9 = (5C; (4,1)).
- El Zein credits exceptions A_1..A_12 to Havet (JCTB 80 (2000) 1–31), and existence for "cycles with a block of length n−1" to Grünbaum (JCTB 11 (1971)).
- So Theorem A is a known special case. S15's proof (constructive, via an insertion lemma) is an independent proof, and the parity statement for even n is a one-line consequence of Rédei.
- I read El Zein's arXiv HTML; Grünbaum's and Havet's papers were not read.

### 1b. Audit summary for THM-4526 (for the orchestrator)

| THM-4526 item | finding of this lane |
|---|---|
| (A) H_n avoided only by C3, C3[1,C3,1]; odd for even n | **Correct.** Independent exhaustive check n ≤ 9; proof sound. **Not new**: the cycle type (n−1,1) case of Rosenfeld's problem (El Zein's exceptions A_1, A_9; Havet 2000; Grünbaum 1971). The phrase "We found no reference" should be replaced by these citations. |
| (B) u(n) = 1,2,4,6,8,9,11 (n = 2..8); 51 classes at n = 7, 1617 at n = 8 | **Reproduced independently:** n ≤ 7 by full completion enumeration (engine `unav`); n = 7, 8 by an orderly host search, which gives 54 → 51 and 2290 → 1617 classes. |
| (B) n ≤ 6: unique maximiser path + span-3 (odd count for every T: "parity curiosity") | Explained and proved (D3): Rédei for n = 4, 5, 6 by short derivations; fails for all n ≥ 7. |
| (C) u(n) = Θ(n log n), c ∈ [1/2, 1] | Superseded: u(n) = n log₂ n − O(n log log n) (Linial–Saks–Sós 1983), so c = 1. |
| D70 u(9) | **u(9) = 14, κ(9) = 22** (§8); this equals the excess-law prediction. |

## 2. The parity calculus (B1–B3)

**B1 (deletion–reversal).** emb(S,T) + emb(rev_f S,T) = emb(S−f,T), because [a→b] + [b→a] = 1. Iterating over the set B of arcs that are backward in a vertex order σ gives the **forward expansion**

φ_S = Σ_{R ⊆ B} φ_{F ∪ rev R}   (F = the σ-forward arcs of S).

**Theorem B2 (completeness).** On isomorphism classes of oriented graphs on n vertices, S ↦ φ_S induces

F_2[classes] / ⟨S + rev_f S + (S−f)⟩ ≅ F_2^{tournament classes}.

*Proof.*
1. On labelled graphs, S ↦ 1_{Q_S}, where Q_S is the subcube of labelled tournaments containing S.
2. The forward graphs map to the ANF monomials, a basis. Every S is DR-equivalent to a sum of the 2^{C(n,2)} forward graphs.
3. So F_2[labelled]/DR ≅ F_2^{labelled tournaments}, S_n-equivariantly.
4. Taking S_n-coinvariants gives F_2[classes]/DR ≅ F_2[orbits], via orbit sums; and the orbit sum of 1_{Q_S} is m_S(T) ≡ φ_S(T). ∎

Hence every parity fact about embedding counts is a finite DR computation. The rest of the note is about *short* certificates.

**B3 (certificate lemmas).**
- **(a) Free involution.** |Aut S| even ⟹ φ_S ≡ 0, because Aut S acts freely on embeddings.
- **(b) Union.** φ_{S_1⊔S_2}(T) = Σ_{|X|=n_1} φ_{S_1}(T[X]) φ_{S_2}(T[X^c]). If both parts are rigid with values c_1 and c_2, this equals C(n,n_1)c_1c_2. By Lucas it is odd iff n_1 and n_2 have disjoint binary digits.
- **(c) Forcade with a closed formula.** For an oriented Hamiltonian path whose backward arcs sit at the positions B ⊆ {1..n−1}:

  φ = #{K ⊆ B : cutting the path at K gives a carry-free composition of n} mod 2.

  *Proof.* Forward expansion along the path; every term is a directed linear forest, whose φ is a multinomial mod 2. ∎ Two consequences:
  - all types are odd when n = 2^k, since no composition of 2^k into ≥ 2 parts is carry-free;
  - *palindromic types are odd* (El Sahili–Abi Aad 2020 [CITED]). The mirror image K ↦ n − K pairs the terms. A self-mirror K ≠ ∅ gives a palindromic composition with a repeated part, which is not carry-free.

  Forcade 1973 [CITED; repo HYP-3545]. Checked against the engine for all 3–8-vertex oriented paths.
- **(d) Halving.** Let ι ∈ Aut(K) be an involution that swaps u and v (non-adjacent in K). Then emb(K+(u→v),T) = emb(K,T)/2 *exactly*. So K+(u→v) is Rédei iff emb(K,·) ≡ 2 (mod 4). (Checked on 180 random pairs.)
- **(e) Twisted pair.** If rev_g S has an automorphism of order 2, then φ_S = φ_{S−g}. If S − g has one, then φ_S = φ_{rev_g S}.
- **(f) Ideal formula.** φ_S(A ⇒ B) = Σ_X φ_{S[X]}(A) φ_{S[X^c]}(B), summed over the down-sets X of S of size |A|.

## 3. Necessary conditions and heredity (B4, B5)

**Theorem B4.** If S has m ≥ 1 arcs and φ_S is constant, then Aut(U(S)) contains an involution. If S is Rédei, then also e(S) is odd and |Aut S| is odd.

*Proof.*
1. Each term of φ_S = Σ_π Π_{(u,v)∈S} [π(u)→π(v)] is a product of literals y or 1+y over the m distinct pairs of π(U(S)).
2. So the degree-m part of the ANF is |Aut U(S)|·Σ_{M ≅ U(S)} y^M.
3. A constant has no degree-m part, so |Aut U(S)| is even (Cauchy). ∎

Consequence: no orientation of an asymmetric graph is ever parity-rigid. Both census engines confirm the condition on every rigid class.

**B5 (interval heredity).** Take T = A ⇒ B. A traceable S has the single down-set {0,…,k−1} of each size, so φ_{P+C}(A⇒B) = φ_prefix(A)·φ_suffix(B). Hence every interval of a Rédei path + chords graph is Rédei.

**B5' (source/sink deletion).**
- Take T = (source z) ⇒ T′. Then the ideal formula gives Σ_{s source of S} φ_{S−s}(T′) = φ_S(T) for every T′.
- So if S is Rédei with a *unique* source s, then S − s is Rédei, and dually for a unique sink.
- Checked on all 138 such cases in the n = 5..8 Rédei lists.

## 4. Goodness propagation, cycles and trees (G1, G2, T)

**Lemma G1.** Let U be a graph with at least one edge, and suppose U − f is good for every edge f. Then either every orientation of U is rigid or none is. They are all rigid, for instance, if U has an involution that reverses no edge.

*Proof.*
- For every orientation S and every edge f, φ_S + φ_{rev_f S} = φ_{S−f} is constant.
- Single flips connect all orientations.
- An involution that reverses no edge admits an invariant orientation, whose φ is 0 by B3a. ∎

**Corollary G2.** Every orientation of C_n is parity-rigid iff n is even.
- *Even n:* C_n − f is a path, which is good by Forcade. A reflection through two opposite vertices reverses no edge.
- *Odd n:* the directed cycle has φ ≡ n·hc ≡ hc. This is not constant: TT_n has hc = 0, and TT_n with its first–last arc reversed has hc = 1.
- So all orientations of an even cycle with an odd number of linear extensions are Rédei (H_n is one). Checked for n = 4..8: C_4, C_6, C_8 have 4, 9, 22 orientation classes, all rigid; C_5, C_7 have 4 and 10, none rigid.

**Theorem T (trees).** Call a tree U *hereditarily symmetric* if every subtree with ≥ 2 vertices has a non-trivial automorphism. Then every orientation of U is parity-rigid. Conversely (B4), no orientation of an asymmetric tree is.

*Proof*, by induction.
1. By induction, the components of U − f are good, so by the union lemma φ_{S−f} is constant. G1 applies.
2. Take g ∈ Aut(U), g ≠ 1, and let v be a moved vertex closest to the centre.
3. **Case 1: the parent w of v is fixed.** This holds whenever the centre is a vertex, or g fixes both ends of the central edge. Then the branches at w rooted at v and at g(v) are isomorphic, and swapping them is an involution that reverses no edge. Orient U invariantly (B3a), and G1 finishes.
4. **Case 2: otherwise.** Every non-trivial automorphism swaps the ends a, b of the central edge, and the rooted halves are rigid. So g² = 1, and ι := g swaps the two halves.
   - Orient the halves compatibly, so that K = S − ab is ι-invariant. By B3d, emb(S) = emb(K)/2 = Σ over unordered splits {X, X^c} of emb(H,X)·emb(H,X^c).
   - This is ≡ c·C(n−1, n/2−1), where c is the constant value of φ on the (good) half H. ∎

| n | oriented trees | rigid | non-rigid |
|---|---|---|---|
| 3–6 | 3, 8, 27, 91 | all | 0 |
| 7 | 350 | 286 | 64 = all orientations of S(1,2,3) |
| 8 | 1376 | 1056 | 320 = the 7 trees containing S(1,2,3) |

- At n = 8, all 128 orientations of the asymmetric tree fail, and 32 orientation classes of each of the six symmetric trees fail. All involution-fixed orientations are rigid.
- So "rigid iff not asymmetric" (my first guess) is **REFUTED at n = 8**. The correct statement is the hereditary one.
- Every oriented tree on ≤ 6 vertices is rigid. This is now PROVED, because the smallest asymmetric tree has 7 vertices.

## 5. Infinite families of Rédei graphs (answer to D72)

| family | arcs | why | label |
|---|---|---|---|
| oriented Hamiltonian paths with odd carry-free count (all palindromes; all types for n = 2^k) | n−1 | B3c | PROVED (Forcade; El Sahili–Abi Aad) |
| H_n = P_n + (0,n−1), n even | n | Rédei + rotations | PROVED |
| D_n = P_n + (0,n−2) + (1,n−1), n odd | n+1 | D1 | PROVED |
| orientations of C_{2k} with e(S) odd | 2k | G2 | PROVED |
| orientations of hereditarily symmetric trees with e(S) odd | n−1 | T | PROVED |
| S_1 ⊔ S_2 with Rédei parts and carry-free sizes | sum | B3b | PROVED |
| doubling S ⊔ S + ((x,0)→(x,1)), S Rédei on 2^k vertices | 2\|S\|+1 | B3d | PROVED |
| (H_4 ⊔ H_4) + twisted pair (the two 10-arc maxima at n = 8) | 10 | B3e + census | FINITE-EXACT |
| path + span-3 for n = 4, 5, 6 | 2n−4 | D3 | PROVED (finite family) |

**Theorem D1.**
- (b) For odd n ≥ 5, D_n is Rédei. Its embeddings number H − A − B + AB, where:
  - A = #HPs with v_{n−2}→v_0, which equals Σ_w hc(T−w)·d⁻(w);
  - B = #HPs with v_{n−1}→v_1, which equals Σ_w hc(T−w)·d⁺(w);
  - so A + B = (n−1)Σ_w hc(T−w) is even;
  - AB is even, because (v_0,…,v_{n−1}) ↦ (v_{n−1},v_1,…,v_{n−2},v_0) is a fixed-point-free involution on the double violations.
- (c) For even n the same computation gives φ(D_n) = 1 + Σ_w hc(T−w).
- These integer identities hold for every class with 5 ≤ n ≤ 9 (engine mode `dcheck`: 0 failures). For even n they give 32/56 odd classes at n = 6 and 3588/6880 at n = 8, so D_n is not rigid there.

**Doubling (D4).** Let S be Rédei on m = 2^k vertices. Then D(S,x) := S ⊔ S + ((x,0)→(x,1)) has, by B3d, emb = emb(S⊔S)/2 ≡ C(2m−1, m−1) = 1 (mod 2). All 20 instances with m = 4 are in the n = 8 census.

None of these infinite families has more than about 1.25·n arcs. Whether Rédei graphs can have Θ(n log n) arcs is **OPEN** (§7).

## 6. Path + chords (D2, D3)

| n | R_n (chord i-j means i→j on the path 0→…→n−1) |
|---|---|
| 3 | ∅ |
| 4 | ∅, {03} |
| 5 | ∅, {03}, {14}, {03,14} |
| 6 | ∅, {03}, {25}, {05}, {03,25}, {05,14}, {03,14,25} |
| 7 | ∅, {03,16}, {05,16}, {03,36}, {05,36} |
| 8 | ∅, {07} |
| 9 | ∅, {07,18} |
| 10 | ∅, {09} |
| 11 | ∅, {09,1-10} |

How each row was obtained:
- n ≤ 8: a zeta transform on the tiling cube. Each class's tilings come from the Hamiltonian paths of its representative. The table is a partition into odd parts (Rédei) with orbit sum 2^{C(n,2)}.
- n = 9: all single chords, all chord pairs, and all subsets of the B5 candidates {07,18,08}, over the 191536 classes.
- n = 10: B5 leaves {}, {09}, {07,18,29}, {07,18,29,09}.
- n = 11: B5 leaves the 8 subsets of {09,1-10,0-10}.
- At n = 10, 11 every non-family candidate is refuted by an explicit tournament, re-counted independently in Python by the runner.

**Inductive step for D2'.** If R_{n−1} and R_{n−2} have the stated form, B5 confines the chords at n to spans ≥ n−3:
- even n: the candidates are ∅, {(0,n−1)}, X_3 = {(0,n−3),(1,n−2),(2,n−1)} and X_3 ∪ {(0,n−1)};
- odd n: the candidates are the subsets of {(0,n−2),(1,n−1),(0,n−1)}.

So D2' reduces to refuting two explicit chord sets (n even) or six (n odd) at each n. It is CONDITIONAL beyond n = 11.
- Some refutations are uniform. H_n with odd n is refuted by TT_n with its first–last arc reversed (hc = 1).
- H_{n−1}·P_1 is refuted by the same tournament when n ≡ 3 (mod 4), because A = (n−1)(n−2)/2 is then odd.
- I did not find uniform witnesses for the remaining cases.

**D69 (Theorem D3).** Path + span-3 is Rédei exactly for n = 4, 5, 6.
- **n = 4** is H_4, and **n = 5** is D_5.
- **n = 6** in three steps:
  1. *Twisted pair:* reversing 1→4 gives a graph invariant under (1 3)(2 4), so φ = φ(P_6 + {03,25}).
  2. *Inclusion–exclusion:* AB is even by the involution (v_0,…,v_5) ↦ (v_4,v_5,v_2,v_3,v_0,v_1), and B(T) = A(T^op).
  3. *A = emb(C_4 + tail of length 2) is even:*
     - DR on a cycle arc splits it into Y_2 and a tree T_1.
     - (c_0 c_2) twists Y_2 down to the oriented path c_0←c_1→…→w_2, with φ = 1.
     - Three DR/twist steps reduce T_1 to (antidirected P_4) ⊔ P_2, with value C(6,2) ≡ 1.
     - So A ≡ 1 + 1 = 0.

  Every step is re-checked on all 56 classes by the runner.
- **n ≥ 7:** P_7 contains no copy at all. By B5 no larger n works either.
- **A much simpler witness, from the polynomial.** In the arc variables y_ij = [i→j] (i < j), the ANF of φ − 1 at n = 7 has 25768 monomials, of degrees 2 to 8, including 13 of degree 2.
  - So no single arc reversal of the transitive tournament changes the parity (all 21 counts are odd), but 13 double reversals do.
  - For example, TT_7 with 0→6 and 3→6 reversed contains exactly **2** copies of P_7 + span-3. (FINITE-EXACT; re-checked in the runner.)
  - For comparison, the ANF of φ − 1 for H_7 and for H_4·P_3 already has the degree-1 term y_06. Reversing the single arc 0–6 of TT_7 creates exactly one Hamiltonian cycle, which flips their parity.
- **Why 7:**
  - The same expansion of H_4·P_3 at n = 7 produces T_1 = S(1,2,3), the smallest asymmetric tree. By B4 it can never be rigid.
  - Indeed φ(H_4·P_3) = 1 + φ(T_1) + φ(Y_3), with φ(Y_3) = 0 and φ(T_1) non-constant (verified).
  - The same 7 marks the first non-rigid oriented tree (§4).

The deletion–reversal Horn closure (`closure`) also certifies path + span-3 for n = 4, 5, 6 automatically; the derivations are printed in the .out.

## 7. Census, densest Rédei graphs, r(n) versus u(n) (C1–C3)

| n | oriented graphs | rigid | Rédei | Horn-closure certifies (all / Rédei) | DAG classes | r(n) | u(n) |
|---|---|---|---|---|---|---|---|
| 1 | 1 | 1 | 1 | 1 / 1 | 1 | 0 | 0 |
| 2 | 2 | 2 | 1 | 2 / 1 | 2 | 1 | 1 |
| 3 | 7 | 5 | 2 | 5 / 2 | 6 | 2 | 2 |
| 4 | 42 | 23 | 5 | 22 / 5 | 31 | 4 | 4 |
| 5 | 582 | 178 | 21 | 174 / 20 | 302 | 6 | 6 |
| 6 | 21480 | 2999 | 43 | 2976 / 43 | 5984 | 8 | 8 |
| 7 | 2142288 | 130574 | 156 | 130345 / 151 | 243668 | 9 | 9 |
| 8 | — | — | 220 | — | 20286025 | 10 | 11 |
| 9 | — | — | — | — | — | **11** | 14 |

Methods:
- Rédei lists come from two engines that agree class for class for n ≤ 7: Gray-code completion parity with a full class table, and backtracking parity with the B4 filters. The n = 8 run uses the second engine over all 20286025 DAG classes.
- The host-based orderly search in parity mode reproduces the n = 7, 8 lists by size (labelg classes).

Horn closure: from the free-involution, linear-forest and union bases, deduce the third of S, rev_f S, S−f when two are known.
- The rigid graphs it misses include C_4 + chord, whose count is 2·hc, and the unions (TT4 − e) ⊔ R, which need linear combinations such as φ_{TT4−e} = 1 + φ_{C3⊔K1}.
- It also misses P_7 + {05,16} = D_7.
- All of these are covered by B2 in principle and by D1/B3 in practice.

**Densest Rédei graphs.**
- n ≤ 6: path + span-3.
- n = 7: two 9-arc classes, so r(7) = u(7). Both are among S15's 51 maximum shavings (re-derived here: 51 classes).
  - They are 0>3 0>4 6>0 1>3 5>1 1>6 2>4 5>2 6>2 and 0>3 4>0 0>6 1>3 5>1 6>1 4>2 5>2 2>6.
  - Both lie on the underlying graph "hexagon 0-3-1-5-2-4-0 plus a hub 6 joined to {0,1,2}".
  - In each, an involution of the underlying graph reverses exactly two hub arcs, a twisted pair: (1 2)(3 4) reverses {1>6, 6>2} in the first, and (0 1)(4 5) reverses {0>6, 6>1} in the second.
- n = 8: two 10-arc classes, one fewer than u(8) = 11. Both are (H_4 ⊔ H_4) plus a twisted pair of cross arcs, under an involution that swaps the two H_4's:
  - 0>4 6>0 0>7 1>4 1>6 2>5 6>2 7>2 3>5 3>7, with (0 2)(1 3)(4 5)(6 7) reversing {0>7, 6>2};
  - 0>4 6>0 0>7 1>4 1>6 5>2 2>6 2>7 5>3 7>3, with (0 7)(1 5)(2 6)(3 4) reversing {0>7, 2>6}.
- n = 9: **r(9) = 11**, three arcs below u(9) = 14.
  - *Search:* the orderly host search in parity mode (C3[C3], pool = the 64 killers of the u(9) search) tests the Rédei property at every node of size 12–14 that survives pruning, and the depth-15 nodes as well. The pre-filters are the proven conditions of B4 (e(S) odd, an involution of U(S)); 561306 nodes failed the first, 270028 the second. Result: REDEI_BY_DEPTH 12:0 13:0 14:0 and none at 15.
  - *Maxima:* a complete depth-11 run (Rédei tests at size 11; DFS to size 12) finds 26 orbit representatives. These form **6 isomorphism classes**, each re-verified on all 191536 classes by the engine's backtracking tester.
    - Four classes are S' ⊔ K_1: 0>1 0>3 0>4 1>3 2>4 2>5 4>5 4>6 5>6 7>0 7>1; 0>1 0>3 0>4 1>4 1>5 3>6 3>7 4>5 6>2 6>7 7>2; 0>1 0>3 1>3 1>4 2>4 3>4 5>6 5>7 6>2 6>7 7>2; 0>1 0>3 1>3 1>4 3>4 3>6 5>6 5>8 6>2 8>2 8>6. Each has e = 315.
    - Two classes are spanning (e = 19): 0>1 0>3 1>2 1>4 2>3 2>5 3>6 4>7 5>8 6>7 7>8 and 0>1 0>3 1>2 1>4 2>5 3>6 4>7 5>6 5>8 6>7 7>8. Each is three parallel directed paths 0→3→6, 1→4→7, 2→5→8, joined by the paths 0→1→2 and 6→7→8 plus one extra arc (2→3, resp. 5→6). Their underlying graph is the 3×3 grid minus its middle row, plus a diagonal.
  - *A new phenomenon:* four of the six classes found first are S' ⊔ K_1, where S' is an 11-arc graph on 8 vertices. Since r(8) = 10, S' itself is not Rédei, but its vertex-deleted sum Δφ_{S'}(T) = Σ_v φ_{S'}(T − v) is identically 1 (as for (TT_4 − e) ⊔ K_1 at n = 5).

So parity certificates are as strong as unavoidability up to n = 7, lose one arc at n = 8 and three at n = 9 (r = 0,1,2,4,6,8,9,10,11 against u = 0,1,2,4,6,8,9,11,14). Asymptotically r(n) ≤ u(n) ~ n log₂ n, and every known infinite family is linear. The growing gap suggests r(n) = o(n log n). **OPEN:** is r(n) = Θ(n log n), or is r(n) linear?

## 8. D70: u(9) = 14 and κ(9) = 22 (E1)

**Lower bound.** The 14-arc graph

S_14 = 0>1 0>3 0>6 1>4 1>5 1>7 2>3 2>5 3>4 3>7 4>7 5>6 5>8 6>8

embeds in every one of the 191536 nine-vertex classes.
- It was found by local search from TT_4 ⊔ (path + span-3 on 5). The check was repeated with an independent plain recursive search (`embedall`: 0 avoiders). Adding 0>2 leaves 217 avoiding classes.
- Its properties: connected, not bipartite (it contains the transitive triangle 3>4>7, 3>7); sources {0,2}, sinks {7,8}; longest directed path of 5 vertices.
- |Aut S_14| = 1 and |Aut U(S_14)| = 2. It has e(S_14) = 80 linear extensions, which is even, so this maximum shaving is **not** a Rédei graph.

**Upper bound, by orderly generation (`u9.c`).**
- Every shaving embeds in the host H = C3[C3]: the unique 9-class with |Aut| = 81, regular, with no TT_5 (checked).
- So it suffices to enumerate acyclic arc subsets of H up to Aut(H), keeping a set iff it is lexicographically minimal in its orbit.
  - *Lemma:* deleting the largest arc of a lex-minimal set leaves a lex-minimal set, since inserting an element into a sorted tuple can only decrease it lexicographically. So orbits are generated from their parents.
- A branch is pruned only for a directed cycle, or for non-embeddability into a pool tournament (exact backtracking). Both are inherited by supersets.
- The pool starts as the 14 classes with |Aut| > 9 and grows by every class that kills a depth-15 candidate.
- Every depth-15 survivor is tested against all classes.

Results:

| host | \|Aut\| | nodes | depth-15 candidates | killed | final pool |
|---|---|---|---|---|---|
| C3[C3] | 81 | 4.6 × 10^6 (re-run by the runner with the final code: 4,604,146 nodes) | 50 | all 50 | 64 |
| the \|Aut\| = 27 class | 27 | 1.6 × 10^7 | 32 | all 32 | — |

Both hosts give **no 15-arc shaving**.

**The maximum shavings on 9 vertices (FINITE-EXACT).**
- A third run (C3[C3], K = 14, every surviving 14-set fully checked) finds 64 orbit representatives, which form **54 isomorphism classes**. S_14 is one of them.
- Each class is re-verified in the runner by `embedall` (0 avoiders), and the classes are pairwise non-isomorphic (labelg).
- So the number of maximum-shaving classes is 1, 1, 1, 51, 1617, **54** for n = 4..9.
- 20 of the 54 have an odd number of linear extensions, but none is Rédei, since r(9) = 11.

**Method validation.** The same code at N = 7, 8 reproduces S15 exactly:
- n = 7 (host P_7): no 10-arc shaving; 54 Aut(P_7)-orbits of 9-arc shavings = **51** classes.
- n = 8 (host P_7 + source): no 12-arc shaving; 2290 orbits of 11-arc shavings = **1617** classes (nauty labelg).

**The excess law holds at n = 9.**
- HYP-3817/3819 predicts excess(n) = κ(n) − ⌈log₂ T(n)⌉ = #{self-converse classes with |Aut| > n}.
- At n = 9, nauty finds 14 classes with |Aut| > 9: |Aut| = 15 (8 classes), 21 (4), 27 (1) and 81 (1). Exactly 4 are self-converse (|Aut| = 21, 21, 27, 81).
- So κ(9) = 18 + 4 = 22, as predicted, and the excess sequence is 0,0,0,1,3,4,4 (n = 3..9).

| n | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|
| u(n) | 1 | 2 | 4 | 6 | 8 | 9 | 11 | **14** |
| κ(n) = C(n,2) − u(n) | 0 | 1 | 2 | 4 | 7 | 12 | 17 | **22** |

S15's lower bound gives u(9) ≥ C(4,2) + u(5) = 12 (every 9-tournament contains TT_4). The truth is 14.

## 9. D71 (E2)

Linial, Saks and Sós, "Largest digraphs contained in all n-tournaments", Combinatorica 3 (1983) 101–104, prove

n log₂ n − c_1 n ≥ f(n) ≥ g(n) ≥ n log₂ n − c_2 n log log n.

Here f(n) is S15's u(n) and g(n) its weakly connected version. So **u(n) ~ n log₂ n and S15's constant is c = 1**.
- The statement is CITED from the abstract (via indexers); the paper itself (PDF) was not read.
- The upper bound coincides with S15's counting bound.

## 10. Relation to THM-4524

- **Two dual identities.** DR is the *pattern-side* identity; THM-4524's shaving lemma H(T − e) = H(T) − c(e) is the *host-side* one. Both are [a→b] + [b→a] = 1.
- **Rédei in derivative form.** c_T(e) ≡ c_{T^e}(ē) (mod 2). Reversing an arc preserves the parity of its HP count, because H(T−e) + c_T(e) = H(T) ≡ 1 ≡ H(T^e) = H(T−e) + c_{T^e}(ē).
- **Cycle counts through THM-4524's B1.** The counts in Theorem A and D1 are Hamiltonian-cycle counts, which THM-4524's fixed-point OCF (B1) controls.
  - For odd N, H(T;0) = Σ_U (−1)^{N−|U|} H(T[U]) ≡ 2·hc(T) (mod 4): only partitions into a single odd cycle carry weight 2; all others carry weight ≥ 4.
  - Hence φ(H_n) ≡ 1 + H(T;0)/2 for odd n, and φ(D_n) ≡ 1 + Σ_w H(T−w;0)/2 for even n.
- **Different "shavings".** THM-4524 deletes arcs of the *host* keeping H odd. Here we fix a *pattern* and ask that every host contain it an odd number of times. Both are Rédei parity read through the arc variables.

## 11. Open problems

- (a) Is r(n) = Θ(n log n)? The data r = 9, 10, 11 against u = 9, 11, 14 for n = 7, 8, 9 suggests that parity certificates fall behind.
- (b) Trees: is "good ⟺ hereditarily symmetric" true for all n? (True for n ≤ 8.)
- (c) D2' for all n ≥ 12: show that X_3 = {(0,n−3),(1,n−2),(2,n−1)} (n even), and the six non-D_n subsets (n odd), are never Rédei.
- (d) A characterization of *good* graphs beyond G1. Which graphs with cycles have all orientations rigid?
- (e) The excess law beyond n = 9 (n = 10 has 9733056 classes).

## 12. Plain-language summary for the owner

- **Your 4-vertex shape is forced by parity.** Every 4-tournament contains A→B→C→D plus A→D an *odd* number of times, so it is never missing. Rédei's 1934 theorem says the same for any Hamiltonian path. I call such shapes *Rédei graphs*.
- **One identity generates every parity statement.** "a beats b" equals "1 minus (b beats a)", so flipping one arc of a shape changes its count by the count of the shape with that arc deleted. These flip/delete identities generate *every* true parity statement about such counts.
- **One extra ingredient makes the proofs short:** a shape with a mirror symmetry is always contained an even number of times. Flips, deletions and mirrors prove:
  - your shape and its 5- and 6-vertex extensions (path + every "skip-two" arc);
  - "path + first→last arc" for every even size;
  - "path + two almost-first→last arcs" for every odd size;
  - every orientation of an even cycle;
  - every orientation of a "hereditarily symmetric" tree.
- **Why the extensions stop at 7 vertices.** A shape can only be parity-forced if its undirected skeleton has a mirror symmetry, and 7 is the first size with a tree that has none (a spider with legs 1, 2, 3). The 7-vertex attempt runs straight into that tree.
- **How dense can parity-forced shapes be?** Up to 7 vertices they are as dense as any shape every tournament must contain. At 8 they fall one arc short (10 versus 11), and at 9 three arcs short (11 versus 14). Whether they stay within a constant factor is open.
- **The S15 numbers, one step further.** At 9 vertices the largest shape that every tournament contains has exactly 14 arcs (proved by computer, two different ways). This confirms the repo's "excess law" prediction. For large n the answer grows like n·log₂ n, a 1983 theorem of Linial, Saks and Sós.
- **Theorem A is older than we thought.** Every tournament has a Hamiltonian cycle with exactly one arc reversed, except the 3-cycle and one 5-vertex tournament. This turns out to be a special case of known results of Grünbaum, Havet and El Zein.

## 13. Reproduction

```bash
python3 04-computation/experiments/procgen_shave4_20261001_run.py            # full run (includes the n = 9 search)
python3 04-computation/experiments/procgen_shave4_20261001_run.py --skip-u9  # without the exhaustive n = 9 search
```

Needs cc, nauty (gentourng, geng, directg, countg, labelg, pickg) and python3. Build files go to `scratch/procgen_shave4/build/`. The full run makes **149 checks in about 8.5 minutes** on this machine; the exhaustive n = 9 search takes about 3.5 minutes of that. With `--skip-u9` it takes about 5 minutes. The output's SHA-256 is 9c0c842493ab24980a54a480cb1b0fb9923544c2a5f0d9055d2a40488ed52dac (raw bytes). The options `--with-k14` and `--with-r9` also re-run the two long n = 9 enumerations recorded in the search files.

**Process notes.**
- No git state was changed.
- All processes were single-threaded and under 300 MB.
- The heavy searches (n = 8 DAG census, n = 9 orderly searches) were run one at a time under `nice`.
- A first exploratory n = 9 run at K = 13 was stopped, because it full-checks every true 13-arc shaving; it is not part of any deliverable.
- Two earlier r(9) parity runs (Rédei tests from size 10) were stopped as too slow and replaced by the runs recorded in the r9 file, which use the B4 pre-filters, a stronger pool and tests from size 11/12. The 11-arc graphs they had found were all among the final 6 classes.
- The search code was optimized during the session (several cached embeddings per pool tournament). Its final version is re-validated against S15's n = 7, 8 results and re-runs the n = 9 search inside the runner.

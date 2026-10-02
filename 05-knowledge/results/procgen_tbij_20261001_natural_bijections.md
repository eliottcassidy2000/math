# procgen_tbij_20261001 — natural bijections for the twisted Mallows–Sloane identities (OPEN-Q-060, bijective form)

Lane `tbij`, session `collatz-procgen-20260922`, 2026-10-01. The lane was resumed after a reboot. This note supersedes
nothing; it extends THM-4524 (E1).

- Runner: `04-computation/experiments/procgen_tbij_20261001_run.py --full`. It re-verifies every computational
  claim below and ends with ALL CHECKS PASSED: 39749 checks in 1256 s, peak RSS about 360 MB, run with nice.
- Output: `05-knowledge/results/procgen_tbij_20261001.out` (sha256 58e27d3ae9f28de41fffa59d9997af6d6f55511926979ec08cac8250ff3ff917;
  it contains timing fields).
- Library: `04-computation/experiments/procgen_tbij_20261001_{lib,nat,twist,hu,blocks}.py`. It calls only the
  nauty binaries geng, gentourng, labelg and dreadnaut, plus sympy (Smith normal form).

## 0. Status table

Notation (all objects on the vertex set [n], with S_n acting by relabelling):

| species | objects | iso-type count |
|---|---|---|
| CLASS | switching classes of tournaments | A049313 |
| EVENE | "even" Euler graphs: no automorphism reverses an odd number of edges (the untwisted Euler graphs of THM-4524 E1) | A049313 |
| TOUR | tournaments | A000568 |
| EVENG | even graphs (Andersson / von Brömssen) | A000568 |
| TWOG | two-graphs, i.e. switching classes of graphs | A002854 |
| EULER | Euler graphs | A002854 |

A *natural map* X → Y is an S_n-equivariant map of labelled objects (a species morphism).

| id | statement | status |
|---|---|---|
| G | Twisted Mallows–Sloane principle. For every S_n-submodule U ⊆ F_2^{E(K_n)}: #(tournaments mod U-reversals)/iso = #(even graphs in U^⊥)/iso, with a per-permutation fixed-point identity | PROVED; VERIFIED for 8 submodules, n ≤ 8 (two independent computations) and per permutation for n ≤ 5 |
| G0 | U = 0: tournaments and even graphs are equinumerous | CITED (Royle–Praeger–Glasby–Freedman–Devillers 2023) |
| G1 | U = cuts: A049313 = #EVENE | PROVED (THM-4524 E1). Here it is identified as the switching analogue of G0 |
| G2 | #self-converse tournaments (= A002785) = Σ_{even graphs} (−1)^{#edges}; #self-converse switching classes = Σ_{EVENE} (−1)^{#edges} | PROVED; VERIFIED (n ≤ 8, resp. n ≤ 9) |
| S | odd n: #self-converse switching classes on n = #self-converse tournaments on n − 1 | PROVED; VERIFIED n = 3, 5, 7, 9 |
| N1 | A natural map X → Y bijective on iso types exists iff the stabiliser-containment bipartite graph on orbits has a perfect matching | PROVED |
| O1 | No natural map EVENE → CLASS exists, for any n ≥ 3 | PROVED |
| O2 | CLASS → EVENE, bijective on types: exists for n ≤ 4; does not exist for 5 ≤ n ≤ 100 | PROVED for 5 ≤ n ≤ 100 (block lemma; two non-isomorphic forced classes per n, certificates checked by nauty; n = 5 entirely by hand). FINITE-EXACT exhaustive matching for n ≤ 10 |
| O3 | EVENG → TOUR: no natural map, n ≥ 2. TOUR → EVENG, bijective on types: exists for n ≤ 4; does not exist for 5 ≤ n ≤ 15, nor for n ∈ {22, 23, 26, 31, 32, 33, 38, 46, 47, 61, 62, 63, 71, 74, 83, 86} | PROVED (block-type Hall certificates; n = 5 by hand) and FINITE-EXACT (exhaustive, n ≤ 9) |
| O4 | Mallows–Sloane at even n ≥ 4: TWOG and EULER have equal cycle index but no natural map in either direction is bijective on types. For odd n Seidel's map is an isomorphism | PROVED (all even n ≥ 4); FINITE-EXACT check n ≤ 8 |
| P1 | Odd n: CLASS ≅ {tournaments with d⁺ ≡ d⁻ (mod 4) at every vertex} as species. So A049313(n) counts tournaments with all scores even (n ≡ 1 mod 4) or all scores odd (n ≡ 3 mod 4) | PROVED; VERIFIED n ≤ 9 (labelled uniqueness n ≤ 7) |
| P2 | Every bipartite Euler graph is even | PROVED; VERIFIED n ≤ 9 |
| P3 | Higashitani–Ueyama, arXiv:2409.10904v4, Ex. 4.10: "We expect s_{4,n} = t_{4,n} for any n, though we do not have a proof." In fact s_{ℓ,n} = t_{ℓ,n} for all ℓ ≥ 2, n ≥ 1 | PROVED (Brauer's lemma for finite abelian groups); VERIFIED ℓ ∈ {2,3,4,6,8,9,12}, n ≤ 6 or 7 |
| B | THM-479 branch integrality for n ≡ 2 (mod 4), n ≥ 6. τ: C ↦ C + M_g reverses the arcs of the perfect matching of an involution g ∈ Aut C. It is a fixed-point-free involution on the iso classes of switching classes with even-order Aut. Hence N_lev(n) = ½·#{such classes} ∈ Z | PROVED; VERIFIED n = 6, 10. n ≡ 0 (mod 4): OPEN |
| C1 | No natural bijection CLASS → EVENE for any n ≥ 5, and none TOUR → EVENG for any n ≥ 5 | OPEN (conjecture) |

Evidence for C1. The deficit (symmetric objects left unmatched) grows: CLASS → EVENE gives 1, 1, 2, 3, 4, 10 for
n = 5..10, and TOUR → EVENG gives 1, 2, 5, 12, 44 for n = 5..9. Certificates cover 5 ≤ n ≤ 100 (CLASS) and
5 ≤ n ≤ 15 plus sporadic n ≤ 86 (TOUR).

**Answer to OPEN-Q-060 (bijective form).** Read "bijective" as "natural", i.e. given by a relabelling-invariant
construction. Then the answer is no, in a strong form:

- No such construction sends even (untwisted) Euler graphs to switching classes, for any n ≥ 3.
- No such construction sends switching classes onto even Euler graphs bijectively on iso types, for any
  5 ≤ n ≤ 100 (conjecturally for all n ≥ 5).
- The same holds for the 2023 theorem "tournaments and even graphs are equinumerous", whose natural-bijection
  question was posed as open by its authors: proved for 5 ≤ n ≤ 15 and for 16 sporadic larger n ≤ 86;
  conjecturally for all n ≥ 5.
- Mallows–Sloane wrote that "for n even, we have been unable to find such a correspondence". That remark is
  explained exactly: none exists, for every even n ≥ 4.

On the positive side, for odd n there is the honest tournament analogue of Seidel's theorem. Every switching
class contains exactly one tournament with d⁺ ≡ d⁻ (mod 4) at every vertex, and this is an isomorphism of
species. So A049313(n) counts these tournaments. It counts no natural family of Euler graphs.
The unifying principle (Theorem G) explains why "even/untwisted" objects keep appearing.

## 1. The twisted Mallows–Sloane principle

Setting. E = E(K_n), V = F_2^E (graphs on [n]), S_n permutes E. Tournaments form a V-torsor: v ∈ V reverses the arcs
on the pairs in v. Fix T0: i → j for i < j, and let c_g ⊆ E be the set of pairs on which g(T0) and T0 differ. For a
graph X and g ∈ Aut(X) put

  sgn_X(g) = (−1)^{#edges {u<v} of X with g(u) > g(v)} = (−1)^{|X ∩ c_g|}.

This is the sign of Royle et al. (DFGPR below) and THM-4524's eps_F. It equals the sign of the permutation g induces
on the 2|E(X)| arcs of X. In cycle form it is (−1)^{#even cycles of g, of length 2k, whose antipodal pairs
{x, g^k x} are edges}. X is *even* if sgn_X ≡ 1 on Aut X.

**Theorem G.** Let U ⊆ V be an S_n-submodule and 𝒯_U the set of tournaments modulo reversal along members of U.
For every g ∈ S_n, #Fix_{𝒯_U}(g) = Σ_{X ∈ U^⊥, gX = X} sgn_X(g). Hence the number of iso classes of 𝒯_U equals the
number of iso classes of even graphs X ∈ U^⊥.

*Proof.* 𝒯_U is a torsor under W = V/U, whose character group is U^⊥ (χ_X(v) = (−1)^{|X ∩ v|}). For X ∈ U^⊥ the
function φ_X(T) = (−1)^{|X ∩ δ(T,T0)|} is well defined on 𝒯_U; here δ(T,T0) is the set of pairs where T and T0
differ. The φ_X form a basis of C[𝒯_U]. From δ(g⁻¹T, T0) = g⁻¹(δ(T, T0) + c_g) we get
g·φ_X = (−1)^{|gX ∩ c_g|} φ_{gX}. So the permutation representation on 𝒯_U is monomial in this basis, and Aut X
acts on the line of φ_X by sgn_X. Therefore sgn_X is a homomorphism on Aut X. The trace of g gives the fixed-point
identity. Frobenius reciprocity gives dim C[𝒯_U]^{S_n} = #{[X] : sgn_X ≡ 1 on Aut X}. ∎

Instances, n = 2..8. Both columns are computed independently: affine Burnside by F_2 linear algebra per cycle type,
and a nauty orbit count of even graphs in U^⊥. They agree in every entry.

| U | 𝒯_U | U^⊥ | counts, n = 2..8 |
|---|---|---|---|
| 0 | tournaments | all graphs | 1, 2, 4, 12, 56, 456, 6880 (G0) |
| ⟨J⟩, J = K_n | tournaments up to converse | graphs with an even number of edges | 1, 2, 3, 10, 34, 272, 3528 |
| C (cuts) | switching classes | Euler graphs | 1, 1, 2, 2, 6, 12, 79 (G1 = E1) |
| C + ⟨J⟩ | oriented two-graphs up to sign | Euler graphs with an even number of edges | 1, 1, 2, 2, 6, 12, 69 |
| Z (cycle space) | score-parity vectors | cuts | 1, 2, 3, 3, 3, 4, 5 |
| Z + ⟨J⟩ | score-parity vectors up to complement (n even; = Z for n odd, as J ∈ Z) | cuts with an even number of edges | 1, 2, 2, 3, 2, 4, 3 |
| C ∩ Z | tournaments mod switching at even-size sets (n even; = tournaments for n odd) | graphs whose degrees all have the same parity (n even; all graphs for n odd) | 1, 2, 3, 12, 8, 456, 156 |
| C + Z | (V for n odd) | Euler cuts δW, \|W\| even | 1, 1, 2, 1, 2, 1, 3 |

**Corollary G2.** Let N_sc be the number of self-converse objects. Orbits on 𝒯_{⟨J⟩} number (N + N_sc)/2.
Subtracting G0 (resp. G1) gives:
- #self-converse tournaments = Σ_{even graphs X} (−1)^{|E(X)|};
- #self-converse switching classes = Σ_{X ∈ EVENE} (−1)^{|E(X)|}.

Values for tournaments: 1, 2, 2, 8, 12, 88, 176 (n = 2..8). Values for classes: 1, 1, 2, 2, 6, 12, 59, 176
(n = 2..9). For example, at n = 9 the 792 even Euler graphs split as 484 with an even and 308 with an odd number of
edges, and 484 − 308 = 176 classes are self-converse.

**Proposition S (self-converse classes, odd n).** For odd n, the number of self-converse switching classes on n
vertices equals the number of self-converse tournaments on n − 1 vertices. The latter is OEIS A002785 (self-complementary
oriented graphs = self-converse tournaments, McKay). Values: 1, 2, 12, 176 for n = 3, 5, 7, 9.

*Proof.*
1. Burnside count. Let c be the converse involution, which commutes with S_n. For an S_n-set X with such an
   involution, #(c-stable orbits) = (1/n!) Σ_g #{x : gx = c(x)}.
2. Anti-automorphisms. Call g an anti-automorphism of C if gC = C^op = C + J. Write gT = T + J + δW for a member T and
   sum this over a pair orbit O of g. This gives the condition: |O| + (c_g-holonomy on O) + |O ∩ δW| ≡ 0
   (THM-479 bookkeeping).
3. No odd cycles of length ≥ 3. A non-diameter orbit inside an odd cycle of length l ≥ 3 has |O| = l odd,
   c_g-holonomy 0 and cut contribution 0, which is impossible.
4. Fixed points. Two fixed points force W to split them, and three fixed points cannot all be split pairwise. For n
   odd the remaining cycles are even, so g has exactly one fixed point v.
5. Reduction to tournaments. Write g = h ⊕ (v fixed). The pointed correspondence C ↦ C_v − v (C_v = the member with v
   a source) is an S_{n−1}-equivariant bijection from classes on [n] to tournaments on [n] ∖ v, and
   (C^op)_v − v = (C_v − v)^op. So #{C : gC = C^op} = #{T : hT = T^op} for every h.
6. Summing over the n choices of v turns (1/n!)Σ_{g ∈ S_n} into (1/(n−1)!)Σ_{h ∈ S_{n−1}}. ∎

So by G2, Σ_{F ∈ EVENE(n)} (−1)^{|E(F)|} = Σ_{X ∈ EVENG(n−1)} (−1)^{|E(X)|} for odd n.

*Remark (the idea in goal 2).* An untwisted F does not give an Aut(F)-invariant class. For F = ∅, Aut F = S_n, and
no class is S_n-invariant once n ≥ 3. What untwisting gives is an Aut(F)-invariant *partition* of the labelled
classes into the two halves {φ_F = ±1}, with no canonical choice of half. For bipartite F the choice is canonical:
count disagreements with any Eulerian orientation of F, whose parity does not depend on that orientation (P2). Even
then this yields a half, not a class.

## 2. Natural maps and the matching criterion

**Lemma N1.** Let X, Y be finite G-sets. A G-equivariant map X → Y inducing a bijection X/G → Y/G exists iff the
bipartite graph on X/G ⊔ Y/G has a perfect matching, where O ~ O' iff Stab(x) ≤ Stab(y) for some x ∈ O, y ∈ O'
(equivalently, O' contains a Stab(x)-fixed point).

*Proof.* Equivariance forces Stab(x) ≤ Stab(f(x)). Conversely, given a matching β choose x_O and y_O ∈ β(O) with
Stab(x_O) ≤ Stab(y_O), and put f(g x_O) = g y_O. This is well defined, because g x_O = g' x_O implies
g⁻¹g' ∈ Stab(y_O). ∎

In species language this is a natural transformation that is bijective on unlabelled structures. "Multi-valued
constructions" (an Aut(x)-invariant set of mutually isomorphic outputs) impose no constraint at all, so the
single-valued notion is the meaningful one. Orbits with trivial stabiliser are adjacent to everything. So when
|X/G| = |Y/G|, a perfect matching exists iff the orbits with nontrivial stabiliser can be matched, and that is
what the computation checks.

*The requested Aut-multiset comparison is trivially negative.* Σ_{[C]} 1/|Aut C| = 2^{C(n−1,2)}/n!, while
Σ_{[F] ∈ EVENE} 1/|Aut F| is strictly smaller for n ≥ 3, because twisted Euler graphs exist. So no bijection preserves
|Aut| (n = 3..9 checked). The same holds for TOUR vs EVENG. For TWOG vs EULER the |Aut| multisets coincide for odd n
(Seidel) and differ for n = 4, 6, 8. Lemma N1 replaces equality of groups by the containment a construction needs.

## 3. Obstructions

### 3.1 Even Euler graphs vs switching classes (E1)

**Theorem O1.** For n ≥ 3 there is no S_n-equivariant map EVENE[n] → CLASS[n].
*Proof.* ∅ is even and S_n-invariant, but no class is S_n-invariant. A transposition is not level for n ≥ 3
(Babai–Cameron Lemma 7.1 = THM-479). ∎

A group H acting on a block B is *twist-rigid* if every nonempty H-invariant graph on B has an automorphism with
sgn = −1.

**Block lemma.** Let [n] = B_1 ⊔ … ⊔ B_r with every |B_i| = b_i odd. Let H = H_1 × … × H_r, where H_i acts on B_i
transitively, with odd order, and twist-rigidly.
- (a) Every H-invariant even graph F is the blow-up of a graph R on [r] in which every block is independent.
- (b) If F is also Euler, R is Euler.
- (c) If r ≤ 2, or r = 3 with two blocks of equal size, the only H-invariant even Euler graph is ∅.

*Proof.*
- (a) For i ≠ j, H_i × H_j is transitive on B_i × B_j, so F contains all or none of the B_i–B_j edges. Suppose
  F_i = F[B_i] ≠ ∅. Take an odd automorphism σ of F_i and extend it by the identity. It is an automorphism of F,
  because every outside vertex sees all or none of B_i. Its cycles lie in B_i, so by the cycle form
  sgn_F(σ) = sgn_{F_i}(σ) = −1. Hence every block of an even F is independent.
- (b) A vertex of B_i has degree Σ_{j ~_R i} b_j ≡ deg_R(i) (mod 2).
- (c) There is no nonempty Euler graph on ≤ 2 vertices, and on 3 vertices only the triangle, which blows up to
  K_{b,b,b'}. Pairing the two equal parts by an involution gives b 2-cycles that are all edges, so the sign is
  (−1)^b = −1. ∎

**Twist-rigid families.**
- b = 1.
- b = q ≡ 3 (mod 4) a prime power, H = {x ↦ ax + c : a a nonzero square}. H has odd order q(q−1)/2 and is
  2-homogeneous, so the only invariant graphs are ∅ and K_q, and a transposition is odd on K_q. This includes
  (3, C_3), (7, F_21), 11, 19, 23, 27, 31, 43, …
- b = q ≡ 5 (mod 8) a prime power, H = {x ↦ ax + c : a ∈ □_odd}, where □_odd is the odd-order index-2 subgroup of
  the squares; H has order q(q−1)/4. Since −1 ∈ □ ∖ □_odd, ±□_odd = □, so the invariant graphs are ∅, the Paley graph
  P_q, its complement and K_q. On P_q the map x ↦ −x has (q−1)/2 2-cycles {x, −x}. Exactly (q−1)/4 of them are
  edges (an odd number), and the same holds for the complement. This includes (5, C_5), 13, 29, 37, 53, 61, …
  The invariant tournament is x → y iff y − x ∈ □_odd ∪ r□_odd, with r a primitive root.
- b = 9, H = Z_3² (regular). The invariant graphs are Cayley graphs on the k = 1..4 directions:
  - k = 1: 3K_3, which has adjacent twins;
  - k = 2: K_3 □ K_3, where a row swap has 3 edge 2-cycles;
  - k = 3: K_{3,3,3};
  - k = 4: K_9.
  All are odd.

All sizes ≤ 100 in these families were also checked by computer: every nonempty invariant graph was enumerated and
an odd automorphism found with dreadnaut.

If Aut(C) ⊇ H with H as in (c), Lemma N1 forces every natural partner of C to be ∅. Such classes exist: put an
H_i-invariant tournament on each block and orient all arcs between two blocks the same way. With 3 blocks the
transitive and cyclic block patterns are switching-equivalent.

**Theorem O2.** A natural map CLASS[n] → EVENE[n] that is bijective on types exists for n ≤ 4 and does not exist for
5 ≤ n ≤ 100.

*Proof for n = 5, by hand.* A049313(5) = 2, and the even Euler graphs are ∅ and C_4 ∪ K_1. The class of the regular
R_5 has an automorphism of order 5. The class of "source → 3-cycle → sink" has one of order 3. The Euler graphs on 5
vertices with an automorphism of order 3 are ∅, K_3 ∪ 2K_1, K_5 − K_3 and K_5. Those with one of order 5 are ∅, C_5
and K_5. Only ∅ is even. Both classes need ∅. ∎

*Proof for 5 ≤ n ≤ 100.* The block constructions of type (c) available at size n are: (n) when n is rigid;
(b_1, b_2) with b_1 + b_2 = n; and (b, b, b') with 2b + b' = n. For every n in the range the runner finds at least
two non-isomorphic ones. Non-isomorphism is checked with the complete class invariant
min_w canon(T_w − w), where T_w is the member with w a source. The minimum number is 2, attained at
n ∈ {5, 6, 8, 9, 28, 100}. Each is forced to ∅, so two classes compete for one graph. ∎

The exhaustive matching (n ≤ 10) confirms the theorem. For n ≤ 10 the classes forced to ∅ are *exactly* the
block-lemma classes:

| n | forced classes (blocks; |Aut C|) | exhaustive: #classes with nontrivial Aut / max matched |
|---|---|---|
| 5 | (5) R_5, order 5; (1,1,3) "source → C_3 → sink", order 3 | 2 / 1 |
| 6 | (1,5), order 5; (3,3), order 18 | 6 / 5 |
| 7 | (7) QR_7, order 21; (1,1,5), order 5; (3,3,1), order 9 | 7 / 5 |
| 8 | (1,7) = the PSL(2,7) Paley class, order 168; (3,5), order 15 | 39 / 36 |
| 9 | (1,1,7), order 21; (9) = (3,3,3), order 81 | 73 / 69 |
| 10 | (1,9), order 81; (3,7), order 63; (5,5), order 50 | 1066 / 1056 |

(The Z_3²-Cayley class and the (3,3,3) class coincide at n = 9.) Going beyond 100 by this method needs additive
representations n = b_1 + b_2 (n even) or n = 2b + b' (n odd) over twist-rigid sizes. That is a Goldbach/Lemoine-type
question, so C1 stays a conjecture.

### 3.2 Even graphs vs tournaments (the DFGPR problem)

DFGPR (J. Algebraic Combin. 57 (2023) 515–524) write: "It is an open problem to find a natural bijection between
the sets of unlabelled even graphs and tournaments on n vertices."

**Theorem O3.**
- (a) For n ≥ 2 there is no S_n-equivariant map EVENG[n] → TOUR[n]: ∅ is S_n-invariant, and no tournament is.
- (b) A natural map TOUR[n] → EVENG[n] that is bijective on types exists for n ≤ 4. It does not exist for
  5 ≤ n ≤ 15, nor for n ∈ {22, 23, 26, 31, 32, 33, 38, 46, 47, 61, 62, 63, 71, 74, 83, 86}.

*Proof of (b) for n = 5, by hand.* The tournaments with a nontrivial automorphism are R_5 (order 5) and four with a
C_3 of type (3,1,1). The four have scores 43111, 42220, 32221 and 33310: a 3-cycle plus two fixed points, with uniform
arcs. Among the 16 graphs invariant under (abc)(d)(e), the even ones are, up to isomorphism, ∅, K_{1,3} ∪ K_1, K_{1,4}
and K_{2,3}. Any graph containing the triangle abc has adjacent twins. The edge de alone, K_{1,1,3}, K_5 − e and the
rest have an odd transposition. The only even graph with a 5-cycle automorphism is ∅. Five tournaments compete for four
graphs. ∎

*General mechanism.* Block lemma (a) holds for arbitrary graphs. Call T0 *forced* if its block group admits only ∅.
This holds for a single rigid block (n rigid) and for two equal blocks (b, b), where K_{b,b} is odd through the swap.
Let b' = n − k be rigid with k ≤ 4. The pattern family (b', 1^k) consists of all arc patterns among the k singletons
and between them and B. It has P iso types, while the H_{b'}-invariant even graphs have G types. The block-level
version of Theorem G suggests P ≈ G, and P = G in most computed cases. P < G also occurs, e.g. n = 6, k = 3, where
non-block isomorphisms merge patterns; such cases give no certificate. Whenever P + 1 > G, the set
S = {T0} ∪ patterns has |S| > G ≥ |N(S)|.

*Block-type Hall search (this is what proves (b)).* Take as left vertices all pattern tournaments of
(b', 1^k) for every rigid b' = n − k, 1 ≤ k ≤ 4, together with the forced tournaments. Give each the bound
N_H(T) ⊇ N(T): the types of H-invariant even graphs, intersected over its constructions. A Hall violator in this
restricted bipartite graph is a Hall violator for TOUR[n] → EVENG[n]. Since the true neighbourhoods are smaller,
the restricted deficit is a lower bound for the true one.
- It finds a violator for every 5 ≤ n ≤ 15 and for n ∈ {22, 23, 26, 31, 32, 33, 38, 46, 47, 61, 62, 63, 71, 74, 83,
  86}.
- The restricted deficits for n = 5..15 are 1, 2, 5, 2, 5, 3, 5, 2, 5, 3, 4. At n = 5, 6, 7 they equal the
  exhaustive deficits.
- Example, n = 8 (no forced tournament exists): the 11 patterns of (5, 1³) plus the two tournaments
  "QR_7 + source/sink" (13 tournaments) have all their bounds inside the 12 H_5-invariant even graph types.
  This is exactly the violator found by the exhaustive search.
- Gaps (16 ≤ n ≤ 21, …) occur where fewer than two rigid sizes lie in [n − 4, n − 1]; richer block families would
  be needed there.

Exhaustive matching (n ≤ 9): the tournaments with nontrivial Aut number 5, 15, 57, 324, 3199 for n = 5..9, and at
most 4, 13, 52, 312, 3155 of them can be matched (deficits 1, 2, 5, 12, 44).

### 3.3 Mallows–Sloane at even n (two-graphs vs Euler graphs)

By Brauer's permutation lemma (Brauer 1941; e.g. Isaacs, *Character Theory of Finite Groups*, Thm 6.32), TWOG and
EULER have equal fixed-point counts for every permutation, i.e. *the same cycle index series* (species in the sense of
Bergeron–Labelle–Leroux 1998). Equal permutation characters do not force isomorphic permutation representations:
the tables of marks can differ, as here. This was checked for all cycle types with n ≤ 8. For odd n they are
isomorphic species (Seidel). Mallows and Sloane (1975): "But for n even, we have been unable to find such a
correspondence."

**Theorem O4.** Let n ≥ 4 be even.
- (a) TWOG[n] ≇ EULER[n] as S_n-sets, and no natural map TWOG → EULER is bijective on types. The two two-graphs
  "all triples" and "no triple" are S_n-invariant and distinct, but the only S_n-invariant Euler graph is ∅,
  since K_n is not Euler for even n.
- (b) No natural map EULER → TWOG is bijective on types. Let B_n = K_{n/2,n/2} if 4 | n and B_n = 2K_{n/2} if
  n ≡ 2 (mod 4). The pairwise non-isomorphic Euler graphs ∅, K_{n−1} ∪ K_1 and B_n have automorphism groups S_n,
  S_{n−1} and S_{n/2} ≀ S_2. Each of these groups fixes only the two trivial two-graphs, so three orbits compete for
  two.

*Proof of (b).* A two-graph is a set of triples meeting every 4-set evenly.
- S_{n−1}, the stabiliser of n, has two orbits on triples (containing n or not). The 4-set {n, a, b, c} contains
  three triples of the first kind and one of the second, so both orbits are in or both out.
- For S_{n/2} ≀ S_2 with parts P, Q, the orbits are "3 + 0" and "2 + 1". A 4-set split 3 + 1 contains one triple
  of the first kind and three of the second, so again both are in or both out. ∎

The runner verifies (a) and (b) symbolically for every even n ≤ 30, and the exhaustive matching for n ≤ 8 agrees:
perfect for odd n, not perfect for n = 4, 6, 8.

## 4. Positive results

**P1 (odd n: the tournament Seidel theorem).** For odd n every switching class contains exactly one tournament with
d⁺(v) ≡ d⁻(v) (mod 4) at every vertex, i.e. all scores ≡ (n−1)/2 (mod 2).
- C ↦ (this member) is an S_n-equivariant bijection, so Aut(C) = Aut(member), which has odd order.
- So A049313(n) is the number of tournaments with all scores even (n ≡ 1 mod 4) or all scores odd (n ≡ 3 mod 4).
- The labelled mod-4-Eulerian tournaments form one coset T* + Z (reversal along Euler graphs), of size
  2^{C(n−1,2)}.

*Proof.* By THM-1470 (Lemmas 1–3), the members of a class realise every score-parity vector of weight
≡ C(n,2) (mod 2) exactly once, and the constant vector ((n−1)/2)·1 has weight ≡ C(n,2). ∎

THM-1470 proved the case n ≡ 1 (mod 4) ("even tournaments"). For n ≡ 3 (mod 4) it recorded only "n tournaments with
a single odd score", whereas the canonical member is the all-odd one. The same row-sum condition appears for skew
matrices over Z/ℓ as Higashitani–Ueyama's "modular Eulerian" matrices (§6); P1 is the odd-entry, labelled ℓ = 4
analogue of their Theorem 4.4. Verified: iso classes 1, 2, 12, 792 for n = 3, 5, 7, 9 (= A049313). For n = 3, 5, 7,
every labelled mod-4-Eulerian tournament has exactly one mod-4-Eulerian member among its 2^{n−1} switchings.

**P2.** Every bipartite Euler graph is even.
*Proof.* Orient each component from one colour class to the other. An automorphism maps components to components and
preserves or swaps the classes. A swap reverses all |E_i| edges of a component, and |E_i| = Σ_{v in one class} deg v
is even. ∎

The converse fails. The (bipartite, even) Euler graph counts are (1, 1), (1, 1), (1, 1), (2, 2), (2, 2), (4, 6),
(5, 12), (12, 79), (19, 792) for n = 1..9.

**P3 (Higashitani–Ueyama for composite ℓ).** Let A = Alt_n(Z/ℓ)/⟨X_1, …, X_n⟩ (their switching classes) and
K = {modular Eulerian matrices}. The pairing ⟨M, N⟩ = Σ_{i<j} m_ij n_ij is perfect and S_n-invariant, and
⟨X_v, N⟩ = −(row sum v of N). So K ≅ Â equivariantly. For a finite abelian group A and an automorphism g,
|Â^g| = |A/(g−1)A| = |ker(g−1)| = |A^g|, so Burnside gives s_{ℓ,n} = t_{ℓ,n} for all ℓ, n. Their Theorem 4.7
(ℓ prime) used vector spaces; the finite-group form of Brauer's lemma removes that restriction.

Verified with Smith normal forms; both sides are computed independently, with the per-permutation equality checked
for every cycle type.
- ℓ = 4 gives 1, 1, 3, 8, 62, 1760, 224248 (n = 1..7). Their table stops at t_{4,6} = 1760.
- ℓ = 3 gives 1, 1, 2, 4, 14, 120, 3222 (Cheng–Wells, A240973).
- ℓ ∈ {6, 8, 9, 12} also agree.
- Brute force agrees for ℓ = 4, n ≤ 4 and ℓ = 6, n ≤ 3.

## 5. The branch split (THM-479)

THM-479 splits A049313(n) = N_odd(n) + N_lev(n) by the 2-adic type of level permutations and observed both parts
integral for n ≥ 3. Equivalently N_odd(n) = Σ_{[C]} ρ(Aut C), where ρ(G) is the proportion of odd-order elements of
G; this was checked for n ≤ 9. The same sum over all Euler graphs (or two-graphs) gives N_odd, because odd-order
elements fix equally many classes and Euler graphs.

**Theorem B.** Let n ≡ 2 (mod 4), n ≥ 6. For a class C with |Aut C| even, pick an involution g ∈ Aut C and let M_g be
its perfect matching {x, gx}. Then τ[C] = [C + M_g] is a well-defined fixed-point-free involution on the iso classes
of switching classes with even-order automorphism group. Consequently N_lev(n) = ½·#{such classes} is an integer, and
so is N_odd(n).

*Proof.* Throughout, γ_C(x,y,z) = f(x,y)f(y,z)f(z,x) is the oriented two-graph of C (Babai–Cameron §2).
1. Sylow structure. A 2-subgroup of Aut C acts semiregularly (BC Thm 5.1), so it has order dividing 2. Hence every
   class with even |Aut C| has Sylow 2-subgroups of order 2, its involutions are fixed-point-free (m = n/2 odd
   transpositions) and conjugate in Aut C, and Aut C has a normal 2-complement. So ρ = ½ and
   N_lev = ½·#{classes with even Aut}.
2. Well defined, and an involution. Conjugate choices of g give isomorphic results. g ∈ Aut(C + M_g), and
   (C + M_g) + M_g = C.
3. No fixed point. Suppose h(C) = C + M_g.
   - Making h commute with g: h⁻¹gh ∈ Aut C is an involution, hence equals kgk⁻¹ for some k ∈ Aut C. Replacing h by
     hk, h commutes with g.
   - Then h(C + M_g) = C, so h² ∈ Aut C.
   - h is a 2-element: replace h by its odd-order power h^q, which still maps C to C + M_g. If h^q = 1, then M_g would
     be a cut, which is impossible for n ≥ 4.
   - h² ∈ ⟨g⟩, since ⟨g, h²⟩ is a 2-subgroup of Aut C.
   - h² ≠ g, since a square root of g pairs its m transpositions into 4-cycles and m is odd. So h is an involution
     and h ∉ ⟨g⟩.
4. The orbit picture. K = ⟨g, h⟩ ≅ C_2² has orbits of size 2 and 4. Since n ≡ 2 (mod 4), the number of size-2 orbits
   is odd. On one of them, {a, ga}, h acts trivially or as g; replacing h by gh if needed, h fixes a. Note:
   (i) h is an isomorphism C → C + M_g, so γ_{C+M_g}(hx, hy, hz) = γ_C(x, y, z);
   (ii) γ_{C+M_g}(x, y, z) = (−1)^{#M_g-pairs in {x,y,z}} γ_C(x, y, z).
5. Case 1: some y has hy ∉ {y, gy}. The triple {a, y, hy} contains no M_g-pair. Then (i) and (ii) give
   γ_C(a, y, hy) = γ_{C+M_g}(a, hy, y) = γ_C(a, hy, y) = −γ_C(a, y, hy), a contradiction.
6. Case 2: h(y) ∈ {y, gy} for all y. Then h is the product of the transpositions of g over a set S of g-pairs, with
   ∅ ≠ S ≠ all. If |S| = 1, replace h by gh; since m ≥ 3, this gives |S| ≥ 2. Take b, d in different pairs of S.
   - (i) and (ii) give γ_C(b, gb, gd) = γ_C(b, gb, d).
   - g ∈ Aut C gives γ_C(b, gb, gd) = −γ_C(b, gb, d).
   This is a contradiction. ∎

Verified for n = 6 (4 classes, which τ pairs) and n = 10. At n = 10 there are 560 classes with even-order Aut, all
with Sylow 2-subgroup of order 2; τ is a fixed-point-free involution on them; and N_lev(10) = 280 = 560/2 (THM-479
value). For n = 2, M_g is a cut and τ = id, matching THM-479's exceptional ½ + ½.

The branch split thus has a combinatorial meaning for n ≡ 2 (mod 4). N_lev(n) is the number of τ-pairs
{C, τC}, and N_odd(n) = #{classes with odd-order Aut} + #{τ-pairs}.

*n ≡ 0 (mod 4) remains open.* At n = 8 the 27 classes with even Aut have Sylow 2-subgroups of types C_2 (19 classes,
an odd number), C_4 (2), V_4 (2: Aut ≅ V_4 and Aut ≅ A_4), C_8 (2) and D_8 (2: one Aut of order 8, and
Aut ≅ PSL(2,7)), with ρ = ½, ¼, ¼ resp. ¾, ⅛ and ⅛ resp. ⅝. Their 1 − ρ contributions sum to 15 = N_lev(8). The 19
classes with Sylow C_2 cannot be paired among themselves, so integrality must mix Sylow types. The converse map does
not pair classes either: at n = 6 all four even-Aut classes are self-converse.

## 6. Literature (goal 1)

Only HTML pages, abstracts and the previous instance's local transcript were read; no PDF was downloaded in this
session.

- **Mallows & Sloane**, "Two-graphs, switching classes and Euler graphs are equal in number", SIAM J. Appl. Math. 28
  (1975) 876–880. Their Theorem 2 is a labelled bijection via the star tree at vertex 1 and fundamental circuits;
  it is not equivariant. For even n: "we have been unable to find such a correspondence." Theorem O4 shows that no
  natural one exists.
- **Babai & Cameron**, EJC 7 (2000) #R38.
  - Switching classes ↔ oriented two-graphs ↔ S-digraphs, naturally and with the same Aut.
  - Aut groups are exactly the groups with cyclic or dihedral Sylow 2-subgroups (Thm 5.2); they act semiregularly
    (Thm 5.1).
  - Level permutations (Lemma 7.1, Thm 7.2); odd-order automorphisms fix a member (Lemma 3.1); Cor 7.7 counts the
    classes with even-order Aut for n ≡ 2 (mod 4).
  - Remark 7.4 ("We cannot do this") asks for that count in general.
  - There is no Euler-graph partner and no bijection discussion.
- **Royle, Praeger, Glasby, Freedman, Devillers** ("DFGPR"), "Tournaments and even graphs are equinumerous",
  J. Algebraic Combin. 57 (2023) 515–524, arXiv:2204.01947.
  - They prove von Brömssen's conjecture (OEIS A000568, A334335) by Cauchy–Frobenius, with sgn_X a homomorphism.
  - They pose the natural-bijection problem: "It is an open problem to find a natural bijection…".
  - This is the closest prior art. THM-4524 E1 is the U = C instance of the same mechanism (Theorem G), and
    Theorem O3 answers their question negatively for the n listed, in the equivariant sense.
  - No statement of E1 was found: OEIS A049313 and A002854 carry no such comment as of 2026-10-01, and OEIS full-text
    searches for "odd automorphism" and "even graphs" return only A334335 and A000568.
- **Higashitani & Ueyama**, "Combinatorics of graded module categories over skew polynomial algebras at roots of
  unity", arXiv:2409.10904 (v4, 2025-03-06).
  - Switching of skew matrices over Z/ℓ and "modular Eulerian" matrices (row sums ≡ 0).
  - Theorem 4.4: if gcd(n, ℓ) = 1, each switching class contains a unique modular-Eulerian iso class.
  - Theorem 4.7 (ℓ prime): s = t. Example 4.10 leaves ℓ = 4 open; P3 settles it.
  - Their switching adds X_v to any skew matrix. Tournament switching adds 2X_v to odd-entry matrices.
  - See also arXiv:2107.12927 (±1-skew projective spaces ↔ graph switching).
- **Cheng & Wells**, "Switching classes of directed graphs", JCTB 40 (1986) 169–186 (abstract read). The ℓ = 3 theory:
  two-digraphs, cohomological invariants, Burnside count.
- **Oh, Yoo & Yun**, "Rainbow graphs and switching classes", SIAM J. Discrete Math. 27 (2013) 1106–1111,
  arXiv:1108.6143. A natural bijection between switching classes of graphs and n-rainbow graphs on 2n vertices, i.e.
  double covers. It does not involve Euler graphs.
- **Hage, Harju & Welzl**, "Euler graphs, triangle-free graphs and bipartite graphs in switching classes",
  Fund. Inform. 58 (2003) 23–37 (ICGT 2002). Polynomial algorithms only.
- **Cameron**, "Cohomological aspects of two-graphs", Math. Z. 157 (1977) 101–119. Citation only: the Springer page
  returned a bot challenge, which was not bypassed. A Brauer-lemma proof of Mallows–Sloane appears, according to
  search snippets, in arXiv:1406.7870 (Multicoloured random graphs, PDF only; the ar5iv HTML timed out).
- **Further, not on bijections:**
  - Cameron & Tarzi, "Switching with more than two colours", Europ. J. Combin. 25 (2004) 169–177.
  - Cameron & Spiga, Australas. J. Combin. 62 (2015) 76–90 (switching classes with primitive groups).
  - Moorhouse, "Two-graphs and skew two-graphs in finite geometries", Linear Algebra Appl. 226–228 (1995) 529–551
    (skew two-graphs as invariants).
  - Cameron, "Counting two-graphs related to trees", EJC 2 (1995) #R4.
- **Searches.** WebSearch, the arXiv API and OEIS turned up no bijective proof of Mallows–Sloane for even n and no
  Euler-graph incarnation of A049313.

## 7. Open

- **C1:** no natural bijection CLASS → EVENE (resp. TOUR → EVENG) for every n ≥ 5. The block method reduces CLASS →
  EVENE to additive representations over twist-rigid sizes. A uniform proof needs a new forcing family available for
  every n (or a counting argument at the top of the symmetry lattice). For single cyclic types the twisted Burnside
  counts match exactly, so a counting argument must use several types at once.
- **B for n ≡ 0 (mod 4):** explain N_lev ∈ Z.
- A natural *graph-type* incarnation of A049313 at even n. P1 covers odd n by tournaments; O2 rules out EVENE.
- Statistics (goal 3).
  - The 3-cycle count of the mod-4-Eulerian member is constant mod 4: c3 ≡ C(n,3) − n·C((n−1)/2, 2) (mod 4).
    Writing s_v = (n−1)/2 + 2u_v with Σ u_v = 0 gives Σ C(s_v, 2) = n·C((n−1)/2, 2) + 2Σ u_v², and Σ u_v² is even.
    So c3 mod 4 carries no information.
  - Neither c3 nor the number of "circular" 4-sets of a class (a switching invariant) is equidistributed with |E(F)|
    on EVENE (n = 5, 7, 9).
  - Two natural constructions from the mod-4-Eulerian member to Euler graphs are not injective on types. One is
    {uv : |N⁺(u) ∩ N⁺(v)| odd}, which is Euler because all in- and out-degrees have the same parity; it has 11
    images for 12 classes at n = 7 and 246 for 792 at n = 9. The other is the S² ≡ 1 (mod 4) graph, which is
    constant.
  - None of this is surprising given O2. The q = −1 shadow is G2; the q-graded version of E1 is the Fourier-degree
    grading of S_n-invariant functions on CLASS.
  - Skew-Seidel spectra (switching invariants of CLASS) were not pursued.

## 8. Process

- No git state was changed. Files: the five library modules, the runner, this note and the runner output. Scratch
  is under `scratch/procgen_tbij/` and is not committed.
- Second implementations.
  - Theorem G: affine Burnside vs direct even-graph count.
  - TOUR → EVENG: two matching codes (`analyse` and `tour_eveng_exhaustive`) agree for n ≤ 8.
  - E1 obstruction: exhaustive matching vs block certificates; for n ≤ 10 the forced classes coincide exactly.
  - P3: both sides by Smith normal form, plus brute force.
  - Runner section L: every automorphism group obtained from nauty generators (S-digraph for classes, colour-
    partitioned incidence graph for two-graphs, dreadnaut for graphs and tournaments) equals the group found by
    independent backtracking (classes n ≤ 8, Euler graphs n ≤ 8, tournaments n ≤ 7, graphs n ≤ 6) or by brute force
    over S_n (two-graphs n ≤ 7). Evenness agrees with the cycle form and the arc-sign form of sgn.
- Resources: every process stayed at or below about 360 MB RSS (the n = 10 class enumeration), run with nice and one
  heavy process at a time.
- The previous instance (before the reboot) had downloaded the Babai–Cameron and Mallows–Sloane PDFs to /tmp; they
  were lost. The quotations above come from that instance's local transcript.

## 9. Suggested status line for OPEN-Q-060 (for the orchestrator)

> Bijective form: ANSWERED NEGATIVELY in the natural (equivariant / species) sense.
> - No relabelling-invariant construction maps untwisted Euler graphs to switching classes (n ≥ 3).
> - None maps switching classes onto untwisted Euler graphs bijectively on types (proved for 5 ≤ n ≤ 100;
>   conjectured for all n ≥ 5).
> - The same holds for DFGPR's tournaments/even graphs (proved for 5 ≤ n ≤ 15 and sporadic n ≤ 86) and for
>   Mallows–Sloane at every even n ≥ 4.
> - Positive: for odd n the switching classes are naturally the mod-4-Eulerian tournaments.
> - The odd/even-level branch split is explained for n ≡ 2 (mod 4) by the involution τ (Theorem B).
>
> See procgen_tbij_20261001_natural_bijections.md.

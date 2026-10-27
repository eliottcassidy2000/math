---
id: THM-4819
title: "Regular tournaments under native moves: (A) every geodesic between the transitive class and the regular classes lies in the slab L + R = (n^2-1)/8, and the slab classes on geodesics are exactly those reachable from TT by R-descending moves; (B) slab censuses at orders 5, 7, 9 (every slab class with R <= m-1 is on a geodesic); (C) every regular class through order 11, all 64 circulants and 1,000 random regular 13-tournaments, all 128 order-15 circulants and 200 random regular 15-tournaments, R5[C3], C3[R5], D(P7) and P19 are exactly (n^2-1)/8 moves from the transitive class; (D) the black endpoint lemma (for every Hamilton cycle of T - v, some rotation closed at v is black) and hence every prescribed vertex of a regular tournament on 2m+1 >= 7 vertices is detached from a strong core in exactly m black native moves, in either native role; (E) the even balanced analogue (forced role, m - 1 black moves) at every even order >= 10, finite-exact at 6, 8, 10, false at 4"
status: >
  PROVED (mac-mini-2026-10-27-slab; independent audit pending at filing): (A) the slab inequality and its equality case,
  confinement of geodesics, the reachability characterization and the source/sink obstruction; (D) the black endpoint
  lemma for every regular tournament on 2m+1 >= 7 vertices and, with Session75's audited named-start construction from
  its second move on, the prescribed-singleton theorem for both native roles. This resolves the prescribed-singleton
  strengthening that HYP-9422 (1) left OPEN. The endpoint lemma also replaces the pure-blue classification as the source
  of the first black path. FINITE-EXACT: (B) the slab censuses; (C) the listed instances (search words replayed and
  checked; the lower bound is THM-4810 (C)); (D') every (class, vertex, role) at orders 7, 9, 11 by an independent
  search (42; 270; 26,906 instances). OPEN: HYP-9426 (1)-(2); HYP-9423 (7) (every regular class at every odd order)
  stays OPEN.
session: mac-mini-2026-10-27-slab
source: 05-knowledge/results/regular_slab_prescribed_black_core_eisenstein_catalyst_20261027.md
scripts:
  - 04-computation/experiments/regular_transitive_slab_20261027/slab_census.py, slab_closure.cpp, forward_closure.cpp (+ slab_closure_n9.out, forward_closure_n7_n9.out)
  - 04-computation/experiments/regular_transitive_slab_20261027/regular_geodesic_search.c, paley_geodesic_search.c, rerun_undecided.sh, special15.py (+ regular_geodesic_n13_*.out, regular_geodesic_n15_*.out, paley19_geodesic.out, special15_rows.txt, special15_geodesic.out)
  - 04-computation/experiments/regular_transitive_slab_20261027/prescribed_isolation.c (+ prescribed_isolation_n7_n9.out, prescribed_isolation_n11.out)
  - 04-computation/experiments/regular_transitive_slab_20261027/start_set_conditions.py, gap_parity.py, black_endpoint_check.py (+ .out)
related:
  - 01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md
  - 01-canon/theorems/THM-4815-dipole-principle-native-move-distance-bounds-endpoint-relation-cut-or-colour-exactness-and-the-p7-parity-correction.md
  - 01-canon/theorems/THM-4814-order-ten-prime-cycle-index-types-order-eleven-ledger-skew-hadamard-nuclear-norm-extremality-and-black-core-distance-of-regular-tournaments.md
  - 01-canon/theorems/THM-4818-anti-automorphism-cycles-are-2-mod-4-order-11-prime-types-by-symmetry-order-12-ledger-dipole-census-condensation-criterion-and-p7-class-distance-six.md
  - 01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md
  - 01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md
  - 01-canon/theorems/THM-4557-doubling-a-doubly-regular-tournament-keeps-its-automorphism-group-hyp-9162.md
  - 05-knowledge/hypotheses/HYP-9422-regular-tournaments-exactly-m-black-moves-from-a-singleton-plus-strong-core-nuclear-norm-maximizers-off-the-hadamard-orders-and-the-order-11-types.md
  - 05-knowledge/hypotheses/HYP-9423-dipole-endpoint-conjectures-two-connected-descent-and-the-prime-tower-predictions.md
  - 05-knowledge/hypotheses/HYP-9426-slab-exactness-black-endpoint-counts-eisenstein-exceptional-set-and-catalyst-thresholds.md
  - 05-knowledge/results/session75_overnight_black_root_spectral_20261027/session75-black-core/BLACK_CORE_GEODESICS.md
  - 05-knowledge/results/session74_eleven_mechanisms_20261027/session74-paley/supplement/REPORT.md
---

# THM-4819 — the transitive slab, exact distances, and the prescribed black core of regular tournaments

**Notation.**
- n = 2m + 1 is odd, and c = m(m+1)/2 = (n² − 1)/8.
- A native move τ_Q = K_Q reverses every arc that is not on the directed Hamilton path Q (THM-4810, Session74). It is an involution, and Q is still a path of the image, so the move graph is undirected.
- Sorted scores are s_(0) ≤ … ≤ s_(n−1). Write L = ½ Σ |s_(i) − i| and R = ½ Σ |s_(i) − m|.
- THM-4810 (C): one move changes L by at most 1, so d(T, TT) ≥ L(T).
- THM-4815 B2 and B8: one move changes R by at most 1, and the distance from T to the set of regular classes is exactly R(T).
- A move is **black** unless the positional reversal of Q is an anti-automorphism, in which case it is **blue** (THM-4815, Session75 §1).

## (A) The slab (PROVED)

**(A1) The inequality.** For every tournament of order n,

    L(T) + R(T) >= c,

with equality iff every sorted score s_(i) lies between i and m inclusive. *Proof.* Termwise,
|s_(i) − i| + |s_(i) − m| ≥ |i − m|, with equality iff s_(i) lies between i and m. Then Σ_i |i − m| = 2(1 + ⋯ + m) = m(m+1);
halve. ∎ The equality set is the **slab** G_n.

**(A2) Geodesics lie in the slab.** The distance from TT to the regular set is exactly c:
- L(regular) = c gives the lower bound;
- B8 gives the upper bound R(TT) = c.

On a word of length c from TT to a regular class, L must rise and R must fall by exactly one at every step, since both are 1-Lipschitz. So every state has L + R = c.

**(A3) Which slab classes are on geodesics.** For T in G_n the following are equivalent:
- (i) T lies on a geodesic from TT to the regular set;
- (ii) d(TT, T) = L(T);
- (iii) T is reachable from TT by a chain of R-descending moves.

*Proof.*
- (i) ⇒ (ii): this is (A2).
- (ii) ⇒ (iii): along a chain of length L(T) from TT, the j-th state has L = j. It has R ≥ c − j by (A1), and R ≤ R(T) + (L(T) − j) = c − j by the Lipschitz bound. So R falls by one at each step.
- (iii) ⇒ (i): every R-descending move from a slab state is L-ascending and stays in the slab, since ΔR = −1 and L + R ≥ c force ΔL = +1. So the chain has length L(T). Continue with R(T) further R-descending moves to a regular class (B8), for a total of c. ∎

**(A4) Sources and sinks.** A tournament with a source or sink has R ≥ m: one deviation is ±m, and the deviations sum to zero. It also has no L-descending move (Session74 supplement, pruning lemma, PROVED). Since moves are involutions, such a class with L ≥ 1 is never on a geodesic. A geodesic must therefore realize, at every step, an endpoint-descending pair (THM-4815 B1) by a Hamilton path, without ever passing through a class with a source or sink.

## (B) Slab censuses (FINITE-EXACT)

For each level (L, R) the table gives the slab classes, the geodesic classes, and the failing slab classes split as (source or sink, strong with λ = 1, other).

| n | classes | slab | on geodesics | levels (L,R): slab / geodesic (failing) |
|---|---|---|---|---|
| 5 | 12 | 10 | 6 | (0,3) 1/1; (1,2) 5/1 (3,1,0); (2,1) 3/3; (3,0) 1/1 |
| 7 | 456 | 389 | 252 | (0,6) 1/1; (1,5) 25/1 (17,7,0); (2,4) 101/17 (33,51,0); (3,3) 162/133 (19,9,1); (4,2) 82/82; (5,1) 15/15; (6,0) 3/3 |
| 9 | 191,536 | 165,342 | 152,276 | (0,10) 1/1; (1,9) 113/1 (81,31,0); (2,8) 1,770/82 (747,941,0); (3,7) 10,525/5,424 (2,496,2,441,164); (4,6) 30,123/26,678 (3,438,7,0); (5,5) 49,881/47,666 (2,215,0,0); (6,4) 45,989/45,484 (505,0,0); (7,3) 21,489/21,489; (8,2) 4,979/4,979; (9,1) 457/457; (10,0) 15/15 |

Four observations:
- At orders 5, 7 and 9, every slab class with R ≤ m − 1 is on a geodesic.
- At order 9, every failure with R ∈ {m, m + 1} has a source or sink. At order 7, R = m = 3 already has 9 strong failures with λ = 1 and one other.
- The slab pinches at L = 1, because TT has a single Hamilton path and hence a single move.
- An independent forward closure (`forward_closure.cpp`: BFS from TT by R-descending moves on canonical forms) reproduces the geodesic column level by level:
  - order 7: 1, 1, 17, 133, 82, 15, 3;
  - order 9: 1, 1, 82, 5,424, 26,678, 47,666, 45,484, 21,489, 4,979, 457, 15.

## (C) Regular classes are exactly c moves from the transitive class (FINITE-EXACT instances)

The lower bound d(T, TT) ≥ L(T) = c holds for every regular T (THM-4810 (C)). The upper bound comes from explicit words found by depth-first search:
- **Moves considered.** The search tries only L-descending moves. By (A2)–(A3) this loses nothing.
- **Candidate moves.** A pair (a, b) is descending when THM-4815 B1's endpoint test says so. Hamilton paths for such pairs are sampled from the subset-DP endpoint relation.
- **Ranking and memo.** Candidates are ranked by lookahead (first no source or sink, then the number of descending pairs). Failures are memoised by nauty canonical form.

Every word is replayed from the input labelled tournament. Each path is checked to be a directed Hamilton path at its own stage, and the final state is checked to be transitive. A search that exhausts its node budget is reported as UNDECIDED, never as a counterexample.

| order | instances | result |
|---|---|---|
| 5, 7, 9 | every regular class (1, 3, 15) | exact, by (B) |
| 11 | every regular class (1,223) | exact (THM-4815 B8', random descending chains from TT); P₁₁ also by Session74's checked certificate |
| 13 | all 64 circulant connection sets; 1,000 random regular tournaments (rotational circulant followed by 60n² random 3-cycle reversals) | 64/64; 999 at node budget 50,000, the last at budget 200,000 |
| 15 | all 128 circulant connection sets | 121 at budget 50,000, the other 7 at budget 300,000 with 20 candidates |
| 15 | R5[C3] and C3[R5] (THM-4812's vertex-transitive non-Cayley classes); the doubling D(P₇) (THM-4557, Aut = F₂₁) | all three exact (28 moves) |
| 15 | 200 random regular tournaments | 189 at budget 50,000, the other 11 at budget 300,000 with 20 candidates (`regular_geodesic_n15_random_rerun.out`) |
| 19 | P₁₉ | exact (45 moves); the word is in `paley19_geodesic.out` |

**The order-11 fact merged.** "Every regular tournament of order 11 is exactly 15 moves from the transitive class" is the statement that the score-transport distance 15 = 1 + 2 + 3 + 4 + 5 is realized. By (A), every geodesic is a monotone flow of score units inside the slab:
- at each step one endpoint pair, one vertex above its sorted target and one below, moves one unit toward the middle;
- the sorted scores never leave the band between i and 5;
- (A4) says no state on the way has a source or sink.

**The dipole law decides the bookkeeping.** What remains hard at higher orders is realizing a descending pair by a Hamilton path at every step while avoiding sources and sinks. This is exactly where the order-13 and order-15 searches spend their nodes.

Prior work merged here:
- **Order-11 certificate.** Session74's P₁₁ certificate fixes the class distance at 15. Its labelled distance to TT_id is at least 16, by THM-4815 (C) parity.
- **Lower bounds and B8.** THM-4810 (C) supplies the lower bound, and THM-4815 B8 and B8' supply the distance to the regular set.
- **Forest frames.** HYP-9420 (4) found forest-frame excess 0 at q = 7 and 11. THM-4810 (C) shows that every regular tournament from order 15 on is more than one complement from any forest frame (28 > 26).
- **Special order-15 classes.** THM-4812 (c) identifies the vertex-transitive 15-tournaments, and THM-4557 gives D(P₇).

## (D) The prescribed black core (PROVED for every order 2m+1 ≥ 7)

**The target.** Following Session75, the target is

    C = {1 -> S, S -> 1 : S strong on n − 1 vertices}.

For a vertex v and a role ρ ∈ {source, sink}, C(v, ρ) is the set of members of C in which v is the singleton, in role ρ.

**Theorem D.** Let T be regular on 2m + 1 ≥ 7 vertices, let v be a vertex and let ρ be a role. There is an all-black native word of length exactly m from T into C(v, ρ). No native word of length less than m makes v a vertex of score 0 or 2m.

The lower bound is THM-4814 (D) and the dipole law of THM-4811 (A): s'(v) = n − 1 − s(v) + [v starts Q] − [v ends Q], so |s(v) − m| changes by at most one per move.

**Lemma B (black endpoint lemma).** Let T be regular on 2m + 1 ≥ 7 vertices and let v be a vertex.
- **The cycle.** T − v is strong (Session75 eq. (4): every cut of T − v is at least m − 1). So it has a Hamilton cycle c_0 → c_1 → ⋯ → c_{2m−1} → c_0 (Camion). Fix any such cycle.
- **The rotations.** For each in-neighbour u = c_j of v, put P_j = (c_{j+1}, c_{j+2}, …, c_j, v).

Then at least one P_j is black. Consequently every vertex of T is the last vertex of a black Hamilton path. Applying this to T^op, every vertex is also the first vertex of one.

*Proof.* Suppose every P_j is blue.

**Setup.** Index by starts: put U = {t : c_{t−1} ∈ N⁻(v)}, so |U| = m. For t ∈ U, the positional reversal of the path starting at c_t is

    σ_t :  v <-> c_t,   c_x <-> c_{2t−x}   (x ≠ t, indices mod 2m).

**The key identity.** Each σ_t is an anti-automorphism, and σ_t(v) = c_t. So N⁺(c_t) = σ_t(N⁻(v)): for x ≠ t,

    c_t -> c_x   iff   c_{2t−x} ∈ N⁻(v)   iff   2t − x + 1 ∈ U.          (★)

Two consequences:
- **(a)** Take x = t − 1. The cycle arc c_{t−1} → c_t excludes c_t → c_{t−1}, so t + 2 ∉ U for every t ∈ U.
- **(b)** For t ≠ t′ in U, exactly one of c_t → c_{t′} and c_{t′} → c_t holds. So exactly one of 2t − t′ + 1 and 2t′ − t + 1 lies in U.

**Counting.** By (a), U contains no two elements that differ by 2. In Z_{2m} the map x ↦ x + 2 has two orbits (evens and odds), each an m-cycle when m ≥ 3. Hence |U| ≤ 2⌊m/2⌋.
- **m odd.** Then 2⌊m/2⌋ = m − 1 < m, a contradiction.
- **m ≥ 4 even.** Then |U| = m forces U to meet each orbit in an alternate class, so U = (a + 4Z) ∪ (b + 4Z) with a even and b odd. Take two distinct t, t′ ≡ a (mod 4); there are m/2 ≥ 2 of them. Both 2t − t′ + 1 and 2t′ − t + 1 are ≡ a + 1 (mod 4), so both or neither lie in U. This contradicts (b). ∎

**Remarks on Lemma B.**
- **Order five.** At m = 2 the conditions (a) and (b) have solutions (U = {0,1} or {0,3}), as they must: the regular 5-tournament is pure blue (THM-4730).
- **What the proof avoids.** It uses neither the automorphism group nor the pure-blue classification. Session75 used the classification only to obtain *some* black first path.
- **An alternative route that fails.** The products σ_t σ_{t′} are automorphisms of an explicit cycle type, computed in `gap_parity.py` and checked against direct composition for m ≤ 40. Write d = t − t′, e = gcd(d, m) and ℓ = m/e:
  - if d/e is odd, the type is [1 + 2ℓ, ℓ^(2e−2)];
  - if d/e is even, the type is [((ℓ+1)/2)², ℓ^(2e−1)].

  An even order is impossible, because tournament automorphism groups have odd order. A full (2m+1)-cycle would make T vertex-transitive and hence pure blue, which is impossible for n ≥ 7. These two facts alone settle every m ≤ 40 except 5, 10, 13, 17, 25, 26, 29, 34, 37. Those are the m at which every even gap is bad, so U can be the set of even positions. In this range they are the m such that m or m/2 is odd with every prime factor ≡ 1 (mod 4); this pattern is observed, not proved. The neighbourhood identity (★) is the uniform argument.
- **Data** (`black_endpoint_check.py`).
  - Instances: 14,862 (tournament, vertex, Hamilton cycle of T − v), at orders 7 to 15.
  - Tournaments: every circulant and 30 random regular tournaments per order, every vertex, three random cycles each.
  - Result: no cycle has more than two blue rotations.

*Proof of Theorem D.*

**Reversal gauge.** Session75 §1 converts routes in the reversal gauge into native routes exactly. In the reversal gauge, R_P reverses only the path arcs. The conversion is T_j = U_j or U_j^op and Q_j = P_j or reversed(P_j) by parity of j, and it preserves colours and the class C.

**Session75 §3, run from a chosen first path.** Take P₁ to be a black Hamilton path ending at v; Lemma B supplies one. Every later step uses only that P₁ is black and ends at v:
- **The second move.** The core U₁ − v is strong, with every cut at least m − 2 ≥ 1. One Hamilton cycle of the core, rotated at the m − 1 ≥ 2 in-neighbours of v, gives a second path ending at v whose start differs from P₁'s. That path is black, because its endpoint scores sum to 2m + 1 ≠ n − 1.
- **Later moves.** For j ≥ 2 every path ending at v is black, by the endpoint sum (m + j) + (m − j + 1) > n − 1.
- **States stay strong** for j < m, by the cut identity c(A) = a(2m + 1 − a)/2 and the one-unit change per move.
- **At j = m**, v is a source of U_m, every other score is at least 1, every nonsingleton cut is positive, and the core is strong.

**The two roles.** The native terminal is U_m if m is even and U_m^op if m is odd. So this route makes v a source in the native gauge when m is even and a sink when m is odd. For the other role, run the same construction on T^op, which is regular, and conjugate by the converse:
- K_{rev Q}(S^op) = (K_Q(S))^op;
- colour and the class C are converse-invariant;
- the converse swaps source and sink.

So both native roles occur for every v. ∎

**(D′) Census (FINITE-EXACT).** `prescribed_isolation.c` is independent of the proof. For every regular class, every vertex a and both roles, it searches:
- **Paths.** The dipole law fixes a's endpoint role at each move (alternating, with the first role set by the parity of m and the target role). The search tries the black Hamilton paths with a at that end.
- **Memo.** Failed states are memoised by a nauty canonical form with a marked.
- **Final check.** It checks that a is a source or sink and that T − a is strong.

Every instance succeeds:

| order | regular classes | (vertex, role) instances | run time |
|---|---|---|---|
| 7 | 3 | 42 | — |
| 9 | 15 | 270 | — |
| 11 | 1,223 | 26,906 | 37 s |

**Consequence.** HYP-9422 (1)'s prescribed-singleton strengthening is PROVED: the detached singleton can be chosen in advance, in either native role. This is the regular-start case of the atlas's "black-core connectivity" row. The even balanced case of Session75 §4 with a prescribed singleton is HYP-9426 (2).

## (E) Even balanced orders (PROVED for n ≥ 10; FINITE-EXACT at n = 6, 8, 10; FALSE at n = 4)

Let n = 2m and let every score lie in {m − 1, m}. Session75 §4 (audited) proves that the black distance to C is m − 1 with an unprescribed singleton.

**Which role is forced.** Put d = s − (n − 1)/2 = ±1/2. One move changes |d| by at most one, so after m − 1 moves a singleton has been an endpoint of every path, with alternating roles. Hence:
- a vertex of score m ends as a source if m is odd and as a sink if m is even;
- a vertex of score m − 1 ends in the opposite role.

Only that role is reachable in m − 1 moves. The other role costs at least m moves (Session75 §2: a named sink needs d⁺(v) moves).

**Lemma B′ (even black endpoint lemma).** Let n = 2m ≥ 10, let v have score m, and let c_0 → ⋯ → c_{2m−2} → c_0 be any Hamilton cycle of T − v. This cycle exists, because T − v is strong by Session75 §4's cut bound. Then some rotation (c_t, …, c_{t−1}, v), with c_{t−1} ∈ N⁻(v), is black.

*Proof.*

**Rotations starting at score m.** A blue path has endpoint scores summing to n − 1 = 2m − 1. So a rotation whose start has score m is black.

**Setup.** Suppose all m − 1 rotations are blue. The derivation of Lemma B carries over verbatim with indices mod M = 2m − 1. It gives (★), (a) and (b) for U = {t : c_{t−1} ∈ N⁻(v)}, where |U| = m − 1.

**Halving.** M is odd, so U′ = U/2 is well defined. By (a), U′ has no two cyclically consecutive elements of Z_M. Since |U′| = m − 1 = ⌊M/2⌋, the m − 1 cyclic gaps are all at least 2 and sum to 2m − 1. So all gaps equal 2 except one, which equals 3. After a translation, which preserves (a) and (b), U′ = {0, 2, …, 2m − 4}.

**The pair (0, 2).** Take t = 0 and t′ = 4, which correspond to 0 and 2 in U′. Using 1/2 ≡ m (mod M), the witnesses of (b) are, in halved coordinates, m − 2 and m + 4. They have the same parity, and for m ≥ 6 both lie in [0, 2m − 2].
- **Odd m ≥ 7.** Both witnesses are odd, so neither is in U′.
- **Even m ≥ 8.** Both are even and at most 2m − 4, so both are in U′.

Either way (b) fails.

**m = 5 and 6.** An exhaustive search over the subsets of Z_9 and Z_11 (`start_set_conditions_even.py`) finds no U satisfying (a) and (b). ∎

At m = 3 and 4, (a) and (b) do have solutions: (0,1) and (0,4) mod 5; (0,1,4), (0,3,4) and (0,3,6) mod 7. So orders 6 and 8 rest on the census. `start_set_conditions_even.py` and an inline check confirm that the progression always has a violating pair for m = 5, …, 60.

**Theorem E.** Let T be balanced on n = 2m ≥ 10 vertices, let v be a vertex and let ρ be the role the dipole law forces. There is an all-black native word of length exactly m − 1 into C(v, ρ), and none shorter.

*Proof.*
- **Score m.** If s(v) = m, take P₁ to be a black Hamilton path ending at v (Lemma B′) and run Session75 §4 from its second move. The core stays strong, a second start is available because v has m − 2 ≥ 2 in-neighbours, and the endpoint sums are at least 2m > n − 1. All of this needs only m ≥ 4.
- **Score m − 1.** Pass to T^op. It is balanced, v has score m there, and the converse conjugation applies as in Theorem D. ∎

**Census** (FINITE-EXACT, `prescribed_isolation.c`, even mode). For each vertex the census tests only the forced role.

| n | balanced classes | instances | result |
|---|---|---|---|
| 6 | 5 | 30 | all succeed |
| 8 | 85 | 680 | all succeed |
| 10 | 13,333 | 133,330 | all succeed (24 s) |
| 4 | 1 (S4 = "100111") | 4 | 2 fail |

At n = 4, the only Hamilton path ending at S4's score-2 vertex is 3 → 0 → 1 → 2, and it is blue. This is the even analogue of the T5 exception.

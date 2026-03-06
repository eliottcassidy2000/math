# Tangents, Rabbit Holes & Novel Ideas

**Purpose:** Dense, scannable index of interesting directions that have emerged. Each entry is 2–4 lines. A new Claude reads this whole file quickly to find a jumping-off point. Add entries with the tag format shown. Do not expand entries here — create a separate file in `03-artifacts/drafts/` if you want to develop a tangent.

**Format:** `#tag1 #tag2 | [certainty: low/medium/high] | [source]`

---

## Combinatorics & Structure

**T001** #ballot-sequences #dyck-paths #distribution | certainty: medium | source: MASTER_FINDINGS F4
The number of internal signature patterns giving exactly k Type-II positions within an L-cycle window is C(L-2, 2k-1). This is suspiciously clean and suggests a connection to ballot sequences or Dyck paths. A combinatorial proof of this formula is missing and would be a clean standalone result.

**T002** #asymptotic #cycle-length #mu-weights | certainty: low | source: MASTER_FINDINGS F4
The average Type-II contribution per L-cycle window is (L-4)/4, growing with L. Does this lead to an asymptotic formula for Σ_C μ(C) as a function of cycle-length distribution? Could give insight into why Claim A holds in aggregate even when per-path identity fails.

**T003** #n5-mystery #5-cycles #per-path-identity | certainty: high (RESOLVED by THM-008) | source: MASTER_FINDINGS F5
RESOLVED: mu(C)=1 trivially for ALL cycles at n<=5. At n=5, 5-cycle C\{v} exhausts all T-v vertices so no independent cycles exist. Per-path identity is vacuously true at n<=5. There was no "delicate balance" -- the identity is degenerate. See THM-008.

**T004** #per-path-formula #generalization | certainty: low | source: MASTER_FINDINGS open questions
Is there a correct per-path formula for ALL n (not just n≤5)? The 3-cycle-only formula undercounts at n≥6 because it misses longer-cycle μ contributions. The natural generalization (sum over all cycles) overcounts. The maximal-embedding-only formula also fails. Finding the right per-path formula would likely imply Claim A.

**T005** #conflict-graph #arc-reversal | certainty: low | source: general
How does the conflict graph Ω(T) change under a single arc reversal? Arc reversals are the natural symmetry of the tiling model. If Ω(T) changes predictably, this could give an inductive handle on Claim A.

---

## Connections to Other Mathematics

**T006** #lattice-gas #statistical-mechanics | certainty: medium | source: LaTeX paper abstract
H(T) = I(Ω(T), 2) connects directly to the **hard-core lattice gas model** at fugacity λ=2. The independence polynomial I(G, x) is the grand partition function of the hard-core model on G. This might import tools from statistical mechanics (transfer matrices, Bethe ansatz, correlation inequalities) into the tournament parity problem.

**T007** #2-adic-tower #higher-redei | certainty: medium | source: LaTeX paper abstract
The OCF formula gives H(T) mod 2 from I(Ω(T), 2) mod 2. But I(Ω(T), x) at x=4 gives a mod-4 invariant, at x=8 a mod-8 invariant, etc. This creates a "2-adic tower" of higher-order Rédei theorems. What is the 2-adic valuation of H(T)? Is there a combinatorial characterization?

**T008** #tsscpp #alternating-sign-matrices | certainty: low | source: LaTeX paper abstract
The self-evacuating SYT count |Fix(σ)| = 2^{m²} for n=2m+1 (verified n=5,7). This connects to TSSCPPs (Totally Symmetric Self-Complementary Plane Partitions) and possibly to alternating sign matrices. The full proof is conditional on a classical reference not yet pinned down.

**T009** #hook-lengths #rsk #representation-theory | certainty: high | source: LaTeX paper
All hook lengths of δ_{n-2} (the staircase Young diagram) are odd. This is proved. What are the representation-theoretic consequences? The RSK correspondence applied to tournaments through the pin grid encoding could be worth exploring.

**T010** #forcade-generalization #f2-invariance | certainty: high | source: LaTeX paper
The paper gives the first purely combinatorial proof of F₂-invariance for ALL k-block decompositions (Forcade 1973). This generalizes the Q-Lemma. The proof method might generalize further to other combinatorial structures with parity constraints.

---

## Proof Strategy Ideas

**T011** #claim-a-strategy #deletion-vertex | certainty: low | source: LaTeX paper §claim_strategies
The paper lists 5 organized proof strategies for Claim A. These need to be detailed in a separate document. Key question: which of the 5 routes is most likely to close the gap between n=5 (proved) and n≥6 (open)?

**T012** #involution #fixed-points | certainty: medium | source: Q-Lemma (Route A)
The Q-Lemma uses a fixed-point-free involution on two-block path decompositions. Can a similar involution be constructed that directly witnesses Claim A? Such an involution would pair up Hamiltonian paths in a way that reveals the 2Σμ(C) structure.

**T013** #mu-recursion | certainty: low | source: general
μ(C) = I(Ω(T-v)|_{avoid C\{v}}, 2). This is itself a restricted independence polynomial. Does it satisfy a useful recursion? Can Claim A be proved by simultaneous induction on n and on the structure of Ω(T)?

---

## Dead Ends (documented to prevent re-exploration)

**T016** #dead-end #cycle-bijection #arc-reversal | certainty: high (confirmed dead end) | source: file.txt
The naive cycle bijection under arc reversal FAILS. When T' flips i->j to j->i, cycles through v containing i->j do NOT biject with cycles containing j->i preserving vertex sets, because reversing directed path segments doesn't preserve validity. See MISTAKE-005. Any proof of arc-reversal invariance must work at the SUM level, not cycle-by-cycle.

**T017** #dead-end #contraction | certainty: high (confirmed dead end) | source: file.txt
Tournament contraction T/(a->b) (merging a,b into single vertex) does NOT always give a valid tournament. The "mixed case" (x beats a, b beats x) leaves the arc between x and the merged vertex undefined. adj_T(a,b) != H(any fixed tournament) in general.

**T018** #dead-end #involution #parity | certainty: high (confirmed dead end) | source: file.txt
The naive involution pairing TypeII liftings with orphans fails because the map (Q, TypeII position k) -> (orphan path with swapped block) may produce an invalid path when the arc Q_{k-1} -> Q_{k+1} doesn't exist. A more sophisticated involution or a completely different approach is needed.

---

**T019** #paley-tournament #H-Tp #1729 #taxicab | certainty: high (computation verified) | source: kind-pasteur-2026-03-05-S2
H(T_11) = 95095 = 5*7*11*13*19 = 55*1729, where 1729 is the Hardy-Ramanujan taxicab number (7*13*19 = 12^3+1^3). H/|Aut| gives the sequence 1, 9, 1729 for p=3,7,11. The formula 3^((p-3)/2) was wrong. Whether 1729's taxicab property is meaningful or coincidental is open. The factorization 95095 = 5*7*11*13*19 (product of five consecutive-ish primes) may have a structural explanation via the automorphism group and its action on Ham paths.

**T020** #paley-tournament #c9 #symmetry | certainty: high | source: kind-pasteur-2026-03-05-S2
h_QR = h_NQR = 201 (directed Ham cycles in T_11\{0,1} and T_11\{0,2}). This equality reflects the anti-automorphism x↦2x of T_11 (which maps QR pair {0,1} to NQR pair {0,2} via multiplication by 2 ∈ NQR). The anti-automorphism reverses arcs, mapping Ham cycles to reverse-Ham cycles, but T_11 is self-complementary (T_11 ≅ T_11^op), giving a bijection between Ham cycles of T_11\{0,1} and T_11\{0,2}. This symmetry can be exploited for future Paley tournament computations.

## Computational

**T014** #n7-verification #random-sampling | certainty: high | source: LaTeX paper
Claim A is verified at n=7 by random sampling (3,500 pairs, 0 failures). Full exhaustive verification at n=7 would be computationally expensive (2^21 tournaments × 7 vertices = ~14M pairs) but could be worth attempting on a cluster for increased confidence before investing heavily in a proof.

**T015** #failure-characterization #n6 | certainty: high (RESOLVED) | source: MASTER_FINDINGS F5, FINAL_FINDINGS
RESOLVED: The per-path identity fails iff some Type-II position has mu > 1, iff V\{v,a,b} contains a 3-cycle in T-v. See THM-009. The ~30% failure rate corresponds to the probability that complement vertices form a 3-cycle.

**T021** #arc-reversal-invariance #sum-equality #key-step | certainty: high (key open problem) | source: file.txt
Arc-reversal invariance D(T,v) = D(T',v) is the key unproved step. Requires showing sum_{C: i->j in C} H(T[V\V(C)]) = sum_{C': j->i in C'} H(T'[V\V(C')]). This is a SUM equality, not a bijection. Both sides can be expressed in terms of path-pairs in T_0 (T with arc (i,j) deleted), since T[V\V(C)] = T'[V\V(C)] when V(C) contains both i and j. See OPEN-Q-009.

---

## Paley Tournaments

**T022** #paley-tournament #sub-tournament #h-QR-h-NQR | certainty: high (RESOLVED) | source: PALEY_T11_c9_ANALYSIS.md, kind-pasteur-S2
c_9(T_11) = (55/2)(h_QR + h_NQR) exactly (LEM-001). RESOLVED: h_QR = h_NQR = 201, c_9 = 11055, H(T_11) = 95095 = 55*1729 (direct enumeration, kind-pasteur-S2). CONJ-002 refuted for p=11.

**T023** #paley-tournament #ratio-coincidence #false-pattern | certainty: high (confirmed false) | source: PALEY_T11_c9_ANALYSIS.md
The ratios c_k / C(11,k) for T_11 are 1/3 (k=3), 9/7 (k=5), 4 (k=7) -- no pattern. Ratio=4 at k=7 that motivated c_9=220 estimate was wishful thinking (MISTAKE-006). Don't use ratio arguments for c_k in structured tournaments.

**T024** #paley-tournament #spectral #gauss-sum | certainty: medium | source: PALEY_T11_c9_ANALYSIS.md
Off-diagonal matrix powers in Paley tournaments: (A^k)_{0,d} = [5^k - Re(omega^k)]/11 +/- [Im(omega^k)*sqrt(11)/11]*chi(d), where omega = (i*sqrt(11)-1)/2 and chi is the Legendre symbol. Gives a closed form for counting walks in terms of chi(d) without full matrix powers.

**T025** #paley-tournament #H-Tp #1729 #taxicab | certainty: high | source: kind-pasteur-2026-03-05-S2
H(T_11) = 95095 = 5*7*11*13*19 = 55*1729, where 1729 is the Hardy-Ramanujan taxicab number (7*13*19). H/|Aut| gives the sequence 1, 9, 1729 for p=3,7,11. The formula 3^((p-3)/2) fails at p=11. Whether 1729 structure is meaningful is open. Next: H(T_19) unknown; see OPEN-Q-013.

**T026** #paley-tournament #c9 #symmetry | certainty: high | source: kind-pasteur-2026-03-05-S2
h_QR = h_NQR = 201. This equality reflects the anti-automorphism x->2x of T_11 (maps QR pair {0,1} to NQR pair {0,2}) combined with T_11 being self-complementary (T_11 ~= T_11^op). Exploitable for all Paley tournament sub-tournament computations.

**T027** #n7-ABD #near-cancellation #a-eq-d #per-path-failure | certainty: high (1050 pairs tested) | source: kind-pasteur-2026-03-05-S3
At n=7, A (TypeII count) and D (Claim A RHS/2) do NOT coincide in general. Only 5.9% exact equality. Mean A-D ~ 0.097 (near zero), but range -39 to 26 (large variance). A-B has mean -73.78, B-D has mean +73.88 -- large but nearly opposite. The near-cancellation (A-D ~ 0 on average) is statistical, not algebraic. The per-path formula does not simplify at n=7 as hoped. Key: the 5-cycle contributions (which are trivially mu=1 at n=7) do NOT collapse A-D to 0. See 03-artifacts/code/test_n7_ABD.py for computation.

---

## Arc-Flip Approach (opus-S2 + kind-pasteur-S5, merged)

**T028** #conflict-graph #complete #n5 #ocf-proof | certainty: high (proved) | source: opus-2026-03-05-S2
At n=5, the conflict graph Ω(T) is ALWAYS complete: any two odd cycles share a vertex. This is because 3-cycles use 3/5 vertices and 5-cycles use all 5. Consequence: I(Ω,2) = 1 + 2·|cycles| (only independent sets of size 0 and 1). The OCF formula ΔH = ΔI reduces to ΔH = 2·(#destroyed - #created), verified 732/732. Combined with OCF for transitive tournament (base case), this proves OCF for all n=5 tournaments.

**T029** #arc-flip-formula #ocf-reduction #n6 | certainty: high (verified 2216/2216) | source: opus-2026-03-05-S2, THM-013
At n=6: ΔI = -2·Σ_x s_x·H(B_x) + 2·(D5-C5) where s_x = 1-T[x][i]-T[j][x], B_x = V\{i,j,x}. This is a closed-form expression for how the independence polynomial changes under a single arc flip. Combined with OCF for transitive base case, proving this = ΔH proves OCF (and hence Claim A) for all n=6 tournaments. The identity involves only sub-tournament H-values, arc directions, and 5-cycle counts. See THM-013 and 04-computation/q009_fixed_formula.py.

**T030** #n8-formula-failure #vd-pairs-35 #general-formula | certainty: high (verified) | source: opus-2026-03-05-S2
The simplified formula DH = -2*sum(s_x*H(B_x)) + 2*sum(DL-CL) FAILS at n=8 (8/30 random flips). The failure is due to VD 3-5 cycle pairs that first appear at n=8 (3+5=8 vertices). The correct general formula is DH = sum_{k>=1} 2^k * Delta(alpha_k) where Delta(alpha_k) recursively depends on alpha_{k-1} of cycle complements. Verified n=4,...,9 including n=9 where alpha_3 is nonzero. See THM-013 updated.

**T031** #swap-involution #proof-strategy #blocking-vertices | certainty: high (verified) | source: opus-2026-03-05-S2, THM-014
For proving ΔH = ΔI, the swap involution (swap i,j positions in a Ham path) creates a perfect matching between adj(i,j) and adj'(j,i) paths. Unmatched T-paths are blocked by s_x=-1 vertices (either predecessor of i doesn't beat j, or successor of j beats i). Unmatched T'-paths blocked by s_x=+1 vertices. This gives adj(i,j)-adj'(j,i) = #U_T - #U_T'. The blocking structure directly connects to the s_x values in the cycle formula. See THM-014.

**T032** #arc-flip #E-invariance #proof-strategy | certainty: high (exhaustive n<=5, random n=6) | source: kind-pasteur-2026-03-05-S5
E(T) = H(T) - I(Omega(T), 2) is invariant under arc flips. Verified exhaustively for n<=5 (10240 flip tests) and by random sampling for n=6 (600+ tests). This gives a NEW PROOF STRATEGY for OCF/Claim A: prove E flip-invariant, then since E(transitive)=0, E=0 for all T. See PROP-001.

**T033** #arc-flip #A-clique #delta-I-formula | certainty: high (algebraic derivation + verification) | source: kind-pasteur-2026-03-05-S5
The change in I(Omega(T), 2) under an arc flip has a clean algebraic formula via the A-clique argument: delta_I = 2*[sum_gained I(R_C, 2) - sum_lost I(R_C, 2)] where R_C = Omega(complement). This is the SAME technique as the Claim B proof, adapted to arc flips instead of vertex deletion. Verified at n=5,6.

**T034** #arc-flip #simple-formula #n-dependence | certainty: high | source: kind-pasteur-2026-03-05-S5
The simple formula delta = 2*(#gained_cycles - #lost_cycles) holds at n<=5 but FAILS at n>=6 (31% failure rate). This is because I(R_C, 2) = 1 for all cycles at n<=5 (complements too small for independent cycles), but I(R_C, 2) can be 3 at n=6 (complement has a free 3-cycle). The correct formula uses the full I(R_C, 2) weights.

**T035** #contiguous-block #dead-end | certainty: high (confirmed dead end) | source: kind-pasteur-2026-03-05-S5
The "contiguous block" decomposition (cycles embedded as contiguous substrings in Ham paths) does NOT give the right counting. f_C = 2*H(comp) per cycle fails (50% at n=4). The "non-cycle" path count f_0 is NOT flip-invariant at n>=5. The correct decomposition is algebraic (A-clique on Omega), not combinatorial (contiguous blocks).

**T036** #paley-tournament #integrality #cycle-counts #structural-theorem | certainty: medium | source: kind-pasteur-2026-03-05-S5 (inbox/other.txt)

**T037** #polynomial-identity #proof #ocf #n4-hand-proof | certainty: high (proved) | source: kind-pasteur-2026-03-05-S6
At n=4, the swap involution identity U_{T'}-U_T = 2*sum(s_x) can be proved BY HAND as a polynomial identity: U_T = pt(r+q-rq) + qr(t+p-tp) + pr + qt, and U_{T'}-U_T simplifies to 4-2(p+q+r+t) = 2*(s_a+s_b). All cross-terms from double-blocking cancel perfectly. Extended to n=5 (512 cases) and n=6 (16384 cases) by exhaustive symbolic verification. See THM-015 and symbolic_proof.py.

**T038** #sign-convention #thm-013 #destroyed-created | certainty: high | source: kind-pasteur-2026-03-05-S6
CAUTION: THM-013 uses D=destroyed, C=created convention. In the formula H(T')-H(T) = 2*sum(s_x*H(B_x)) + 2*(gained_5cycles - lost_5cycles), the 5-cycle term has a PLUS sign. Easy to get wrong — the sign error caused initial n=6 failures.

**T039** #contraction #proof-strategy #digraph | certainty: high (verified) | source: opus-2026-03-05-S3
Contract i,j into single vertex w: adj(i,j) = H(G) and adj'(j,i) = H(G') where G,G' are digraphs on n-1 vertices. Key: G has bidirectional arcs w<->v for s_v=-1 vertices and no arcs for s_v=+1; G' swaps these (bidirectional for s=+1, none for s=-1; s=0 vertices identical). So ΔH = H(G)-H(G') is controlled entirely by the s_v != 0 vertices' connection to w. Paths where w's neighbors are all s=0 cancel. This recovers the swap involution from a graph-theoretic perspective. OCF proof reduces to showing H(G)-H(G') = ΔI for these specific digraphs.
**T040** #even-odd-split #proof-strategy #alternating-sum | certainty: high (verified n=5,...,8) | source: opus-2026-03-05-S4
The adj decomposition delta = sum_S Delta(S, others\S) has a perfect EVEN-ODD SPLIT: the sum over even-sized S equals the sum over odd-sized S. Equivalently, the alternating sum sum (-1)^|S| Delta(S,R) = 0. This means delta = 2*(odd-S sum), which connects directly to the cycle formula (only odd cycles). This is a CONSEQUENCE of OCF (not equivalent — see correction in even-odd-split-lemma.md). The cancellation is global (pairwise Delta(S,R) != ±Delta(R,S) in general). See 03-artifacts/drafts/even-odd-split-lemma.md and 04-computation/q009_even_odd_split.py.

For T_11, the average c_k(T_11)/C(11,k) is an exact integer for ALL k >= 6 = (p+1)/2, but fractional for k < 6. Specifically: c_7/C(11,7)=12, c_8/C(11,8)=45, c_9/C(11,9)=201, c_10/C(11,10)=971, c_11/C(11,11)=5505 are all integers; while c_3/C(11,3)=1/3, c_4/C(11,4)=1/2, c_5/C(11,5)=9/7, c_6/C(11,6)=145/42 are not. Proposed structural theorem: C(p,k) | c_k(T_p) for k >= (p+1)/2, for Paley primes p ≡ 3 (mod 4). The regularity of T_p (each vertex has out-degree (p-1)/2) likely forces this integrality via symmetry. Proving this rigorously is an accessible target -- see OPEN-Q-013.

**T044** #ocf #n7-exhaustive #proof | certainty: high (exhaustive) | source: opus-2026-03-05-S3
OCF PROVED at n=7 by exhaustive verification: 1,048,576/1,048,576 arc variable assignments passed (27min). Script: symbolic_proof_n7.py. Fixes one arc 0->1, varies all C(7,2)-1=20 remaining arcs. Combined with relabeling invariance, proves H(T)=I(Omega(T),2) for ALL n<=7 tournaments. Next: n=8 needs 2^27=134M cases -- hours but feasible with optimization.

**T045** #proof-strategy #dead-ends #per-vertex | certainty: high | source: opus-2026-03-05-S3
DEAD END: Per-vertex decomposition of U_T' - U_T does NOT match per-vertex ΔI. Specifically: nadj(x,i,j)/H(B_x) is NOT constant across vertices. The ratio varies from 2 to 7+ depending on tournament structure. The identity is fundamentally GLOBAL, not per-vertex. This rules out approaches that try to match blocking counts to cycle weights vertex by vertex. The proof must work at the level of the full subset convolution Σ_S [f_i(S)·g_j(R) - f_j(S)·g_i(R)].

**T046** #bijection #2-colored-cycles #interpretation | certainty: medium | source: opus-2026-03-05-S3
I(Omega(T), 2) counts "2-colored vertex-disjoint odd cycle collections" in T. At n=3 with the 3-cycle: 3 paths = 1 (empty) + 2 (two colors). At n=4 with two 3-cycles (H=5): cycle-consecutiveness signatures are {.., C., .C, CC, CC} suggesting paths encode WHICH cycles they "use" and HOW (the color). But the base path (signature ..) uses NO cycles and doesn't correspond to a unique "transitive" path. The bijection, if it exists, is non-trivial and might require a global construction (e.g., Lindström-Gessel-Viennot style).

**T047** #bracket-structure #subset-convolution #proof-strategy | certainty: high (algebraically verified) | source: opus-2026-03-05-S3
KEY ALGEBRAIC INSIGHT for OCF proof. The subset convolution integrand B(u,w) = T[u][i]*T[j][w] - T[u][j]*T[i][w] has a 4-way type structure. Classifying V\{i,j} into M- (s=-1), M+ (s=+1), Z1 (s=0, beats both), Z0 (s=0, beaten by both): the bracket table is {M-: (1,0,0,1), M+: (0,-1,0,-1), Z1: (1,-1,0,0), Z0: (0,0,0,0)}. CRITICAL: Z0 rows and Z1 columns are ALL ZERO. This means vertices beaten by both i,j are invisible in u-position, vertices beating both are invisible in w-position. Delta_H = sum_S sum_{u in S, w in R} B(u,w)*h_end(S,u)*h_start(R,w) has only 6 nonzero bracket types. See bracket_structure.py.

**T045** #transfer-matrix-symmetry #bilinear-form #proof-strategy | certainty: high (verified n=4,...,8, 7500+ tests) | source: opus-2026-03-05-S4b
MAJOR DISCOVERY: The transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S)*B_b(M\S) for a,b in {i,j} is ALWAYS SYMMETRIC. This means sum (-1)^|S| E_i(S)*B_j(R) = sum (-1)^|S| E_j(S)*B_i(R) — each cross-term individually agrees, not just their difference. The Even-Odd Split Lemma (T040/T044) is a COROLLARY: the off-diagonal symmetry implies the difference vanishes. This connects to Feng's Dual Burnside framework where Q=AB is symmetric when the underlying structure satisfies detailed balance. Proving M is symmetric would prove OCF. See 04-computation/symmetry_check.py and 03-artifacts/drafts/paper-connections.md.

**T046** #tournament-representations #flip-classes #rank #proof-strategy | certainty: medium | source: opus-2026-03-05-S4b
Rajkumar et al. (arXiv:2110.05188) show every tournament is in the flip class of an R-cone (vertex beating/losing to everyone). R-cones have constrained structure for OCF: vertex v with no in-arcs means every Hamiltonian path starts or ends at v. STRATEGY: prove OCF for R-cones (where the formula simplifies), then extend via cut-flip phi_S. Also: locally transitive tournaments (= rank 2) have 5-cycles and 7-cycles despite constrained structure (67% have 5-cycles at n=5), so the "only 3-cycles" conjecture fails. Upper bound on representation dimension 2(mu(T)+1) relates to arc-flip distance — could give induction on mu(T). See 03-artifacts/drafts/paper-connections.md.

**T048** #real-roots #claw-free #chudnovsky-seymour #n-geq-9 | certainty: high (proved n<=8, conjectured general) | source: kind-pasteur-2026-03-05-S13
All roots of I(Omega(T), x) are real and negative. PROVED for n<=8 via claw-freeness + Chudnovsky-Seymour theorem (2007). Computationally verified at n=9 (50 random, 3-cycle subgraph) and n=10 (30 random). At n>=9 Omega(T) has claws, so C-S doesn't apply. The Jerrum-Patel (2026) results on subdivided-claw-free graphs (JLMS) may help extend. This is a strong structural conjecture: if true, gives alternative proof of Redei and implies log-concavity of independence polynomial coefficients. See THM-020.

**T049** #comparability #omega-fails-n7 | certainty: high | source: kind-pasteur-2026-03-05-S13
Omega(T) is a comparability graph for n<=6 (exhaustive n=5, 200 random n=6). FAILS at n=7 (1/100 random) and n=8 (58/100 random). So comparability does NOT explain real roots. The hierarchy: Omega is chordal for n<=5 only, perfect for n<=7, comparability for n<=6, claw-free for n<=8.

**T050** #hard-core #fugacity-2 #uniqueness | certainty: high | source: kind-pasteur-2026-03-05-S13
Lambda=2 is OUTSIDE the hard-core uniqueness regime for Omega(T) when max_degree(Omega) >= 4 (which happens already at n=5 where max degree=6). The uniqueness threshold lambda_c(Delta) = (Delta-1)^{Delta-1}/(Delta-2)^Delta gives lambda_c(4)=1.6875 < 2. So standard hard-core correlation decay results don't apply to our setting. The real-rootedness of I(Omega(T), x) is a special property of tournament conflict graphs, not explained by generic hard-core theory.

**T051** #mitrovic #hopf-algebra #redei-berge-noncommuting | certainty: medium | source: kind-pasteur-2026-03-05-S13
Mitrovic (arXiv:2407.18608, 2024) develops Hopf algebra structures on the Redei-Berge function: three CHAs on permutations, posets, and digraphs, with deletion-contraction properties. Follow-up arXiv:2504.20968 (2025) gives a noncommuting variable analogue with deletion-contraction (which the commutative version lacks). These papers don't address our conflict graph/independence polynomial formulation but may provide algebraic tools.

**T052** #paley-oeis #new-sequences | certainty: high | source: kind-pasteur-2026-03-05-S13
The sequences H(T_p)/|Aut(T_p)| = 1, 9, 1729 (for p=3,7,11) and H(T_p) = 3, 189, 95095 are NOT in OEIS. These could be submitted as new sequences once H(T_19) is computed. The Tribonacci connection (A000213 for full transitive tournaments) also has no known tournament interpretation in OEIS.

**T053** #paley-maximizer #A038375 #conjecture | certainty: high (verified p=3,7,11; extended p=19,23) | source: kind-pasteur-2026-03-05-S14, opus-2026-03-05-S9
MAJOR CONJECTURE: Paley tournaments T_p (p = 3 mod 4 prime) achieve the MAXIMUM H(T) among all p-vertex tournaments. Confirmed via OEIS A038375: a(3)=3=H(T_3), a(7)=189=H(T_7), a(11)=95095=H(T_11). P(7) confirmed as GLOBAL maximizer by exhaustive enumeration of all 2^21 tournaments on 7 vertices (240 achieve maximum). NEW: H(P(19))=1172695746915, H(P(23))=15760206976379349. These are candidate values for a(19), a(23). IMPORTANT: a(8)=661 is NOT achieved by the Paley extension T_657 (H=657); the n=8 maximizer has |Aut|=1 and does NOT contain P(7). The conjecture applies only at Paley primes p=3 mod 4. The ratio H(P(p))/(p!/2^{p-1}) converges slowly toward e, consistent with Szele-Alon-Friedland theory.

**T054** #line-graph #heilmann-lieb #omega | certainty: medium | source: kind-pasteur-2026-03-05-S14
Is Omega(T) always a line graph? If so, real roots follow from Heilmann-Lieb (1972) rather than Chudnovsky-Seymour. A line graph must be claw-free AND satisfy Beineke's 9 forbidden subgraph conditions. Testable at n=5,6. If Omega(T) = L(H) for some H, then I(Omega(T), x) = matching polynomial of H, connecting OCF to matching theory.

**T055** #interval-graph #omega-fails-n6 | certainty: high | source: kind-pasteur-2026-03-05-S14 (Tribonacci agent)
Omega(T) is an interval graph for n<=5 (exhaustive) but FAILS at n=6 (13.9% fail exhaustively). Failure is due to non-chordality (induced C4 or longer). The interval property holds for T_full_n (consecutive odd intervals), but not for general tournaments. Interval => Chordal => Perfect => Comparability (no!). Failures: interval n=6, chordal n=6, perfect n=8, comparability n=7.

**T056** #chromatic-redei-berge #mitrovic-stojadinovic | certainty: medium | source: kind-pasteur-2026-03-05-S14 (Redei-Berge agent)
Mitrovic & Stojadinovic (arXiv:2506.08841) prove that the chromatic symmetric function and Redei-Berge function are "almost identical" at the poset level. This connects Stanley's chromatic symmetric function theory to tournament path counting. The Stanley-Stembridge conjecture techniques may thus apply to tournament problems. Also: Grujic & Stojadinovic (arXiv:2402.07606) establish a Hopf algebra structure with deletion-contraction for U_D, potentially providing an algebraic foundation for our vertex-deletion recurrence (Claim A).

**T057** #LGV-lemma #even-odd-split | certainty: medium | source: kind-pasteur-2026-03-05-S14 (Novel hypotheses agent)
THM-016/017 (even-odd split) has strong structural parallels to the Lindstrom-Gessel-Viennot lemma. Talaska (2012) extended LGV to arbitrary digraphs where cycles appear in the denominator. Our alternating sum sum_S (-1)^|S| H(S) h_start(W\S, v) may be an instance of extended LGV. Investigating this connection could yield a proof of the even-odd split independent of OCF.

**T058** #subdivided-claw #bezakova-dichotomy | certainty: medium | source: kind-pasteur-2026-03-05-S14
Bezakova et al. (arXiv:2404.07615, 2024) prove a dichotomy: for bounded-degree H-free graphs, Glauber dynamics mixes in O(n log n) for ALL fugacities iff H is a subdivided claw or path. Check whether Omega(T) at n>=9 is subdivided-claw-free (S_{a,b,c}-free). If so, Jerrum-Patel (2026) zero-free regions extend to cover lambda=2.

**T059** #n8-maximizer #omega-structure #5-cycle-dominance | certainty: high | source: opus-2026-03-05-S9
At n=8, the H-maximizer (H=661, a(8) from OEIS A038375) is a self-converse tournament with |Aut|=1. It does NOT contain P(7) as a vertex-deletion and has non-uniform D_v=(67,63,68,72,72,68,63,67). The full Omega has 78 vertices (20 tri, 50 pent, 8 hept) with density 0.98. The 3-cycle subgraph accounts for only 30% of H: I(Omega_3,2)=201 vs H=661, meaning 70% of Hamiltonian paths come from 5-cycle and 7-cycle contributions. Compare T_657 (Paley extension, H=657): 76 Omega vertices (20,48,8), I(Omega_3,2)=165, 75% from longer cycles. The transition to 5/7-cycle dominance is THE structural phenomenon at n=8.

**T060** #paley-ratio #szele-bound #asymptotic | certainty: high | source: opus-2026-03-05-S9
The ratio H(P(p))/(p!/2^{p-1}) for Paley primes: 2.000 (p=3), 2.400 (p=7), 2.440 (p=11), 2.527 (p=19), 2.557 (p=23). This converges slowly toward e=2.718..., consistent with the Szele-Alon upper bound max H <= O(n^{1/2} * n!/2^{n-1}). If Paley tournaments achieve the asymptotic max, then H(P(p)) ~ e * p!/2^{p-1} as p->infinity. The convergence rate H(P(p))/(p!/2^{p-1}) ~ e(1 - c/sqrt(p)) fits the data well. This connects tournament Hamiltonian path optimization to quadratic residue structure.

**T061** #line-graph-refuted #beineke #heilmann-lieb | certainty: high (REFUTED) | source: kind-pasteur-2026-03-05-S14b
T054 (line graph hypothesis) REFUTED at n=6. K5-e (Beineke forbidden subgraph #2) appears in 53% of n=6 tournaments with >=4 cycles. Omega(T) is NOT a line graph in general. Heilmann-Lieb (1972) cannot explain real roots of I(Omega(T), x). Claw-free at n<=8 is a stronger property that works.

**T062** #subdivided-claw-free #S211-free #real-roots | certainty: high (computational) | source: kind-pasteur-2026-03-05-S14b
Omega(T) (3-cycle subgraph) is S_{2,1,1}-free (no subdivided claw with one path of length 2) at n=9 (0/100 random failures) despite 86% having claws. FAILS at n=10 (92% have S_{2,1,1}). Combined with opus-S9 finding: S_{1,1,1}-free through n=11, fails at n=12. Hierarchy: claw-free (n<=8) -> S_{2,1,1}-free (n<=9) -> S_{1,1,1}-free (n<=11) -> ??? (n>=12). Each subdivision level buys ~3 more vertices. Jerrum-Patel (2026) applies but requires fixed H.

**T063** #quasi-line-fails #chudnovsky-structure | certainty: high (computational) | source: kind-pasteur-2026-03-05-S14b
Omega(T) is NOT quasi-line at n=8 (49%), n=9 (10%), n=10 (0%). A quasi-line graph has N(v) = union of two cliques for every v. Since quasi-line is intermediate between line graphs and claw-free, and Omega has claws at n>=9, this was expected. Not a viable approach for real roots beyond n=8.

**T064** #alon-hamiltonian-maximum #regular-tournaments | certainty: high (literature) | source: kind-pasteur-2026-03-05-S14b
Alon (1990): max H(T) <= c * n^{3/2} * n!/2^{n-1}. Adler-Alon-Ross (2001): max H(T) >= (e-o(1)) * n!/2^{n-1} using random regular tournaments. This means max H(T) = Theta(n!/2^{n-1}). Regular tournaments (including Paley) are near-maximizers. Our conjecture that Paley achieves the EXACT maximum (A038375) is consistent but much stronger than asymptotic bounds.

**T065** #claw-free-characterization #barrier | certainty: high (published) | source: kind-pasteur-2026-03-05-S14b (web search)
Engstrom (arXiv:1610.00805, Alco 2019): A stability-like property of the MULTIVARIATE independence polynomial characterizes claw-freeness. This is a BARRIER RESULT: no broader graph class than claw-free can guarantee real-rootedness of I(G,x) via the same mechanism. For Omega(T), which has claws at n>=9, the real-rootedness conjecture (if true) requires a tournament-specific argument — not a general graph property. The explanation MUST be algebraic, arising from the specific structure of tournament conflict graphs.

**T066** #schweser-redei-revisited #dirac-berge-hierarchy | certainty: high (published) | source: kind-pasteur-2026-03-05-S14b (web search)
Schweser et al. (arXiv:2510.10659, Oct 2025, revised Feb 2026): "The Tournament Theorem of Redei revisited". Exhibits and connects the stronger theorems of Redei, Dirac, and Berge about Hamiltonian paths. Clarifies that Dirac's stronger theorem has two corollaries: one equivalent to Redei's stronger theorem, the other related to Berge's. May sharpen our understanding of the parity structure.

**T067** #knuth-taocp-8a #hamiltonian-enumeration | certainty: medium | source: kind-pasteur-2026-03-05-S14b (web search, OEIS A391999)
Donald Knuth's TAOCP Prefascicle 8a (work in progress, 2025) covers Hamiltonian path/cycle enumeration. Referenced by OEIS A391999. If Knuth has new computational results on tournament Hamiltonian paths, this could include data beyond the 11 terms of A038375 or provide new structural insights.

**T068** #undecidability-tournaments #polynomial-inequalities | certainty: high (published IMRN 2025) | source: kind-pasteur-2026-03-05-S14b (web search)
Chen-Lin-Ma-Wei (arXiv:2412.04972, IMRN 2025): Proving polynomial inequalities in digraph homomorphism densities for tournaments is UNDECIDABLE. This concerns flag algebra / homomorphism density context, not our specific setting, but is a cautionary note about the limits of algebraic approaches to tournament problems.

**T069** #root-gap #independence-polynomial | certainty: high (published) | source: kind-pasteur-2026-03-05-S14b (web search)
Prakash & Sharma (arXiv:2510.09197, FSTTCS 2025): "On The Roots of Independence Polynomial: Quantifying The Gap." Quantifies the gap between the smallest real root and the smallest absolute value among all other roots. Could help analyze the root structure of I(Omega(T), x) even when claw-freeness doesn't apply — e.g., by showing the gap remains positive (forcing real-rootedness) for tournament conflict graphs.

**T070** #perpendicular-maximizer #self-converse #H-maximization | certainty: high (computational) | source: kind-pasteur-2026-03-05-S15
Self-converse tournaments (T ~ T^op) have significantly more Hamiltonian paths: mean H ratio SC/NSC = 1.97 (n=5 exhaustive), 1.12 (n=6 exhaustive), 1.59 (n=7 random). SC tournaments CONTAIN the global H-maximizer at every tested n. Paley tournaments (which maximize H at Paley primes) are always self-converse. In the tiling framework, SC tournaments sit on the hyperplane PERPENDICULAR to the transitive<->full diagonal at the hypercube center. The "perpendicular direction" is where H(T) is maximized.

**T071** #ultra-log-concavity #independence-polynomial #newton | certainty: high (computational) | source: kind-pasteur-2026-03-05-S15
The independence polynomial coefficients alpha_k of Omega(T) satisfy ULTRA-LOG-CONCAVITY: alpha_k/C(m,k) is log-concave where m=|V(Omega)|. Tested at n=5 (exhaustive 784), n=6..9 (random) with 0 failures. This is STRONGER than Newton inequalities. Despite Omega(T) NOT being a matroid (fails at n=6: 62% non-matroid), ultra-log-concavity holds. If I(Omega,x) has all real roots, ULC follows from Newton inequalities for polynomials with real roots. ULC is a necessary condition for real-rootedness.

**T072** #claudes-cycles #odd-even-dichotomy #three-axes | certainty: medium | source: kind-pasteur-2026-03-05-S15
Knuth's "Claude's Cycles" paper (Feb 2026): decomposing K_{m^3} into 3 Hamiltonian cycles using coordinates (i,j,k) with modular serpentine (Gray code) rule s=(i+j+k) mod m. Works for ODD m, fails for EVEN m. Parallel to our Paley maximizer working at ODD primes (p=3 mod 4), failing at even n. The three coordinate axes mirror three "directions" in tournament space: (1) Hamming weight (transitive<->full, Tribonacci), (2) self-converse (perpendicular, H-maximizer), (3) automorphism complexity. Claude's construction is a 3D modular Gray code; our tournament structure lives in a hypercube with analogous directional decomposition.

**T073** #irving-omar-schur #not-positive #barrier | certainty: high (published) | source: kind-pasteur-2026-03-05-S15
Irving-Omar (arXiv:2412.10572, Prop 21): The Schur expansion of U_D is NOT always positive for digraphs with cycles. Example: a 4-vertex tree gives negative s_{(2,2)} coefficient. Schur-positivity holds only for ACYCLIC digraphs (Prop 23). This is a BARRIER: cannot prove real-rootedness of I(Omega(T),x) via Schur-positivity of the Redei-Berge function. However, Corollary 20 CONFIRMS our OCF: ham(D-bar) = sum_{sigma, all odd cycles} 2^{psi(sigma)}.

**T065** #irving-omar #cayley-transform #real-roots #matrix-algebra | certainty: medium | source: opus-2026-03-05-S12
Irving-Omar (arXiv:2412.10572) Corollary 20 proves OCF via matrix algebra: ham(D_bar) = sum over odd-cycle permutations of 2^{psi(sigma)}, derived from exp(sum_k tr(A^k) p_k/k). The ODD-cycle extraction uses arctanh-like splitting: sum_{k odd} tr(A^k)/k = (1/2)*tr(log((I+zA)/(I-zA))). This is the CAYLEY TRANSFORM of zA. For tournaments A = (J-I)/2 + S (S skew-symmetric), so the Cayley transform maps to an orthogonal-like structure. CREATIVE LEAD: If I(Omega(T), x) can be expressed as a determinant involving the Cayley transform of the tournament matrix, real roots would follow from the spectral theory of orthogonal/unitary matrices. The eigenvalues of P(p) are {(p-1)/2, (-1 +/- i*sqrt(p))/2} (Gauss sums), and alpha_3(P(11))/|Aut| = 21 = |Aut(P(7))| — suggesting deep number-theoretic structure in the independence polynomial coefficients.

**T074** #discriminant #turan #real-roots #elementary | certainty: high (PROVED for n<=8) | source: opus-2026-03-06-S15
ELEMENTARY PROOF of real-rootedness via Turán's theorem (THM-021). For n<=8, alpha(Omega(T))<=2, so I(Omega,x)=1+a1*x+a2*x^2. Real roots iff a1^2>=4*a2. At n<=7: only 3-3 disjoint pairs exist, and three pairwise disjoint triples need 9>n vertices, so the disjoint graph is TRIANGLE-FREE. Turán gives a2<=c3^2/4, hence 4*a2<=c3^2<=a1^2. QED. At n=8: add 3-5 pairs bounded by B<=c5, use AM-GM (c5-1)^2>=0. Tight at n=6,7: tournaments with exactly 2 complementary 3-cycles achieve disc=0, I=(1+x)^2. For n>=9 (degree>=3), need Newton inequalities or stronger bounds.

**T075** #perspective-conjecture #vertex-orbits #isomorphism-classes | certainty: medium (fails at n=5->6) | source: human inbox (tournament-tiling-explorer.html)
PERSPECTIVE CONJECTURE from human's Tournament Tiling Explorer: #vertex-orbits at n = #isomorphism-classes at n+1. Holds for n=3->4 (2 perspectives -> 4 classes) and n=4->5 (4 perspectives -> 12 classes). FAILS at n=5->6 (predicts some number, actual is 56). Parent class mapping (each n-vertex class maps to an (n-1)-vertex class via top-vertex deletion) is visualized with purple center labels. Grid transpose = tournament converse (T->T^op). Tool archived at 03-artifacts/visualizations/tournament-tiling-explorer.html.

**T079** #perpendicular-bell-curve #hamming-weight #H-maximization | certainty: high (exhaustive n=5) | source: kind-pasteur-2026-03-05-S16
H(T) forms a perfect symmetric bell curve as a function of tiling Hamming weight at n=5: avg H = 1.0, 3.2, 5.3, 6.9, 8.0, 8.4, 8.0, 6.9, 5.3, 3.2, 1.0. Ratio perpendicular/diagonal = 2.79. NSC tournaments have CONSTANT avg H ≈ 4.5 across all HW levels, while SC tournaments vary from 1.0 to 9.7. The bell curve property likely holds for all n and explains why tournaments at "maximum confusion" (HW = m/2, halfway between transitive and full) have the most Hamiltonian paths.

**T080** #NSC-constrained #SC-dominance | certainty: high (exhaustive n=5) | source: kind-pasteur-2026-03-05-S16
Non-self-converse tournaments have severely constrained H values. At n=5: NSC H ∈ {3, 5} while SC H ∈ {1, 3, 9, 11, 13, 15}. At n=4: NSC H = {3} (constant!) while SC H = {1, 5}. NSC tournaments cannot achieve the global maximum H. This suggests a deeper algebraic constraint: the T ~= T^op symmetry creates constructive interference among Hamiltonian paths.

**T081** #hook-palindrome #irving-omar-proposition26 | certainty: high (exhaustive n=3,4,5) | source: kind-pasteur-2026-03-05-S16
The Irving-Omar hook coefficients [s_{(i,1^{n-i})}]U_T are PALINDROMIC: h_i = h_{n+1-i} for ALL tournaments. This follows from ham(T) = ham(T^op). The transitive tournament has hooks = binomial coefficients C(n-1, k). At n=4, exactly 3 hook patterns (one per isomorphism class): H=1 → [1,3,3,1], H=3 → [3,3,3,3], H=5 → [5,3,3,5]. The hook coefficients distinguish more than just H — they form a finer invariant.

**T082** #J(n,3)-NOT-hereditary #counterexample #tournament-specific | certainty: high (computational) | source: kind-pasteur-2026-03-05-S16
The hereditary real-rootedness conjecture for J(n,3) is FALSE: arbitrary induced subgraphs of the Johnson graph J(n,3) do NOT always have real-rooted independence polynomials. Counterexample at n=9: 7 triples giving alpha=[1,7,5,1] with discriminant -44 (complex roots). About 1.1% of random 7-element subsets of J(9,3) have non-real I.P. roots. BUT: tournament-realizable triple sets (sets of cyclic triples for some tournament) ALWAYS give real roots (0 failures across ~1000 tests at n=5-20). The real-rootedness is NOT a Johnson graph property — it is SPECIFICALLY a tournament structure property. This rules out proving real roots via general intersection graph theory.

**T083** #degree6-real-roots #n20 | certainty: high (computational) | source: kind-pasteur-2026-03-05-S16
First verification of real roots for a DEGREE 6 independence polynomial of Omega_3(T) at n=20. Coefficients [1, 298, 26477, 839457, 9050120, 26503054, 11559099]. All 6 roots real and negative: -1.90, -0.27, -0.076, -0.029, -0.013, -0.006. This is the highest-degree polynomial tested for real-rootedness. The polynomial has enormous coefficients (11M at degree 6) yet roots remain firmly real. Extends verification table to n=20, degree 6.

**T084** #interlacing #vertex-deletion #real-roots #FAILS-n7 | certainty: high (computational) | source: kind-pasteur-2026-03-06-S17, CORRECTED S18
CORRECTED: Interlacing HOLDS exhaustively at n=5 (0/5120) and n=6 (0/196608) but FAILS at n=7 (36/1113 for full Omega, 7/1113 for 3-cycle-only). Failure mode: vertex deletion can dramatically change the root spread — T with [1,20,3] (roots -6.6,-0.05) and T-v with [1,6,1] (roots -5.8,-0.17) fails because the larger root moves the wrong way. NOT a viable proof approach for general n. The original 0/1324 result was an artifact of computing only small Omega graphs. Interlacing works at n<=6 because Omega has small independence number (alpha<=1 at n=5, alpha<=2 at n=6).

**T085** #full-omega-T11 #21169-cycles #degree3 | certainty: high (computational) | source: kind-pasteur-2026-03-06-S17
Full Omega(T_11) has 21169 vertices (all odd cycles), I.P. = [1, 21169, 10879, 1155], degree 3. All 3 roots real. For comparison, Omega_3(T_11) has only 55 vertices. The full Omega is vastly larger but the independence number stays at 3 (same as Omega_3). The ratio I(full Omega, 2)/I(Omega_3, 2) measures how much longer cycles contribute to H.

**T086** #polymer-gas #claw-free #arxiv-2505-22766 | certainty: medium | source: kind-pasteur-2026-03-06-S17 (web research agent)
arXiv:2505.22766 (May 2025): Z_G(z) = (1+z)^|V| * sum_F (-1)^|F| (z/(1+z))^|V_F| for claw-free graphs, connecting independence polynomial to abstract polymer gas. Improved zero-free disk radius r*_Delta = 1/(1+2*Delta) for even Delta. While the zero-free disk doesn't reach x=2, the polymer-gas representation could provide structural insight for Omega(T) at n<=8 where claw-freeness holds. Does NOT help at n>=9.

**T087** #tiling-flip #grid-symmetry #flip-invariant #blueself-blackself | certainty: high (PROVED) | source: opus-2026-03-06-S16
THM-022: Grid-symmetry is a FLIP INVARIANT (proved via commutativity of pointwise ops with permutations). Flip and transpose commute: F(G(t))=G(F(t)). Consequence: blueself tilings come in flip-pairs. Additionally, transpose preserves I(Omega(T),x) exactly (cycle reversal argument). See THM-022.

**T088** #blueself-parity #odd-n-obstruction #score-sequence | certainty: high (PROVED ALL ODD n) | source: opus-2026-03-06-S16, S17
NO BLUESELF at ANY odd n (PROVED): For grid-symmetric tilings, k_0+k_{n-1}=n-2 (endpoint constraint). Flip maps endpoint multisets {1+k_0, n-2-k_0} -> {n-1-k_0, k_0}. These match iff k_0=(n-2)/2 (Case A, non-integer at odd n) or 1=0 (Case B, impossible). Pure algebraic proof — no enumeration needed. Verified: n=3,5,7 exhaustive. Even n: blueself exists (1 class at n=4, 2 at n=6). Script: `04-computation/blueself_odd_n_proof.py`.

**T089** #blueself-blackself #mutual-exclusivity #class-level | certainty: high (verified n<=7) | source: opus-2026-03-06-S16
No isomorphism class contains BOTH blueself and blackself tilings, verified exhaustively at n=4,5,6,7. At odd n this is vacuous (no blueself exists). At even n (4,6), the result is non-trivial. Self-flip members within a class always have uniform grid-symmetry. The mechanism: flip preserves grid-symmetry (T087), so self-flip members come in same-type pairs. Open whether the class-level separation holds at n>=8.

**T090** #flip-not-top #tiling-operation #tournament-level | certainty: high (computational) | source: opus-2026-03-06-S16
The tiling flip F(T) is NOT isomorphic to T^op in general. F reverses non-path arcs only; T^op reverses ALL arcs. Measured: flip(T)≅T^op for only 25% (n=4), 12.5% (n=5), 3.9% (n=6) of tilings. Flip does NOT preserve H or I(Omega,x). The flip is a path-dependent operation with no clean tournament-level description.

**T091** #SC-maximizer-within-score-class #blueself #H-maximization | certainty: high (exhaustive n=4,5,6,7) | source: kind-pasteur-2026-03-06-S18
WITHIN each self-complementary score sequence class, the tournament achieving max H is ALWAYS self-converse. Verified exhaustively at n=4,5,6,7 (all 15 self-complementary score classes at n=7 confirmed). Only self-complementary score sequences can contain SC tournaments (s_i + s_{n-1-i} = n-1). At n=7: global max H=189 at regular score (3,3,3,3,3,3,3) = Paley T_7. All 15 classes had FIRST max-H achiever already SC. The SC advantage comes from more vertex-disjoint cycle pairs (alpha_2 in independence polynomial). Key examples at n=6: (3,3,3,2,2,2) SC max=45 > NSC max=43; (4,3,3,2,2,1) SC max=37 > NSC max=31.

**T092** #score-regularity-H-correlation #variance | certainty: high (exhaustive n=6) | source: kind-pasteur-2026-03-06-S18
Score sequence regularity (low variance) strongly correlates with high H at n=6. Regularity 0.250 (most regular, (3,3,3,2,2,2)) → H=41-45. Regularity 2.917 (transitive, (5,4,3,2,1,0)) → H=1. SC classes cluster at low regularity (regular scores).

**T094** #blueself-maximizer #even-n #global-max | certainty: high (exhaustive n=4,6) | source: kind-pasteur-2026-03-06-S18
At even n, a BLUESELF class always achieves or ties for the global maximum H. n=4: unique blueself class has H=5 (global max). n=6: blueself classes have H=41 and H=45; H=45 ties with a non-self-flip class for global max. The blueself property (grid-symmetric + self-flip) is the strongest structural constraint, and it correlates with maximal H. At odd n, blueself doesn't exist (THM-022 Theorem 5), so the SC maximizer operates through blackself instead.

**T093** #SC-mechanism #anti-automorphism #cycle-orbits #disjoint-pairs | certainty: high (computational n=5,6) | source: kind-pasteur-2026-03-06-S18
KEY MECHANISM: The anti-automorphism sigma of an SC tournament is an involution (sigma^2=id at n=6) that partitions vertex set into orbits {i, sigma(i)}. This creates a natural pairing on 3-cycles: C and sigma(C) are vertex-disjoint iff their vertex sets partition into 3 sigma-orbits. At n=6 score (3,3,3,2,2,2): SC tournament with sigma=(1,0,5,4,3,2) has 8 three-cycles forming 4 disjoint complementary pairs (each pair covers all 6 vertices), giving alpha_2=4 and IP=[1,14,4], H=45. The NSC tournament with same score has only 1 disjoint pair, alpha_2=1, H=43. The 3-cycle vertex sets in SC case form a combinatorial design: {0,1}x{2,3}x{4,5} structure from sigma orbits. IMPORTANT REFINEMENT: SC maximizer only applies within SELF-COMPLEMENTARY score sequences (where s_i + s_{n-1-i} = n-1). Non-SC score sequences have NO SC tournaments, so the comparison is vacuous. Testing at n=7 (20 self-complementary score sequences).

**T095** #anti-aut-involution #disjoint-cycle-pairs #even-n-advantage | certainty: high (proved + verified) | source: kind-pasteur-2026-03-06-S18
THEOREM: Every anti-automorphism of an SC tournament at n=5,6 is an involution (sigma^2=id). At even n=6: sigma is fixed-point-free, decomposing vertices into n/2=3 orbits {i,sigma(i)}. All 3-cycles that select one vertex from each orbit are paired by sigma with their vertex-disjoint complement, creating C(n/2,3)*4 disjoint pairs. At n=6 with 3 orbits: 4 disjoint pairs, yielding alpha_2=4 in I.P. The NSC tournament with same score has only 1 disjoint pair. At odd n=5: sigma has exactly 1 fixed point + 2 two-cycles, limiting disjoint pairs. This explains why even-n SC tournaments have disproportionately high alpha_2 and thus high H. The "involution on orbits" mechanism is the algebraic engine behind the SC maximizer phenomenon.

**T096** #grid-overlap #n-copies #sigma-orbit #induction | certainty: high | source: kind-pasteur-2026-03-06-S18e (user insight)
The n-vertex tiling grid (triangular, C(n,2) cells) contains n copies of the (n-1)-grid (delete one vertex) and C(n,2) copies of the (n-2)-grid (delete two vertices). Each cell is in exactly n-2 copies of (n-1) and C(n-2,2) copies of (n-2). For SC tournaments at even n with orbits {a,sigma(a)}: deleting an entire orbit gives 3 copies of n=4 inside n=6 that inherit anti-automorphism symmetry. The 8 transversals (one vertex per orbit) correspond exactly to the 8 three-cycles, pairing into 4 vertex-disjoint complements. This "grid overlap" view connects the sigma-orbit structure to the physical tiling geometry and may enable an inductive proof.

**T097** #paley-deletion #maximizer-heredity | certainty: high (VERIFIED p=3,7,11) | source: kind-pasteur-2026-03-06-S18e, opus-2026-03-06-S18
PALEY DELETION GIVES MAXIMIZER: T_p − v achieves max H at n=p−1 (= OEIS A038375(p−1)). Verified: T_3−v: H=1=a(2) ✓. T_7−v: H=45=a(6) ✓. **T_11−v: H=15745=a(10) ✓ (NEW, opus-S18).** By vertex-transitivity all deletions are equivalent. The Claim A decomposition at T_7: H(T_7)−H(T_7−v) = 144 = 2×72, with sum_mu = 6×3+30×1+24×1 = 72 (all 3-cycle complements in T_7 have a 3-cycle, so mu=3). CONJECTURE: For all Paley primes p≡3 mod 4, H(T_p−v) = a(p−1). Combined with T053 (Paley achieves a(p)), this says the maximizer chain T_p → T_{p-1} is "hereditary" via vertex deletion. Scripts: `04-computation/paley_deletion_test.py`.

**T099** #blueself-not-SC-maximizer #clique-deletion #interlacing | certainty: high (computational) | source: opus-2026-03-06-S17
BLUESELF ≠ SC MAXIMIZER: At n=6, two blueself classes (H=41, H=45) both in score class (3,3,3,2,2,2). SC max in that class is H=45, but blueself H=41 is NOT the maximizer. Blueself classes are always SC, always have regular scores, but not always the max-H SC class. They have distinctive I.P. structure: H=41 has I.P.=[1,16,2] while H=45 has I.P.=[1,14,4] — the blueself maximizer has MORE disjoint pairs (alpha_2=4 vs 2) but FEWER total cycles (alpha_1=14 vs 16). Script: `04-computation/blueself_sc_maximizer_connection.py`.

**T100** #clique-deletion-interlacing #through-v-clique #omega-structure | certainty: high (proved + computed) | source: opus-2026-03-06-S17
THROUGH-V CYCLES ALWAYS FORM A CLIQUE in Omega(T) (proved: any two cycles through v share vertex v). Deleting vertex v from T = deleting this clique from Omega(T). At n=5: 100% of remaining cycles adjacent to some through-v cycle. The clique deletion recurrence I(G,x) = I(G-u,x) + x*I(G\N[u],x) can be applied sequentially to remove through-v cycles. The high density of Omega(T) (~98%) means the "remainder" G\N[u] is very small at each step. This connects to Heilmann-Lieb for matching polynomials and Chudnovsky-Seymour for claw-free graphs. See `03-artifacts/drafts/interlacing-clique-deletion.md`.

**T101** #omega-quasi-regular #spectral-density #real-roots-mechanism | certainty: high (computational, n=5-15) | source: opus-2026-03-06-S17
OMEGA_3(T) IS QUASI-REGULAR: The ratio lambda_max/avg_degree is ~1.005-1.011 for ALL tested n (5-15). This means the conflict graph is spectrally indistinguishable from a regular graph. Density starts at 1.0 (n=5) and decreases slowly: 0.95 (n=6), 0.83 (n=8), 0.71 (n=10), 0.52 (n=15). Still above 0.5 at n=15! I.P. degree stays very small (max 3 at n=10). 0 real-root failures in 100 random samples per n for n=5-10. The quasi-regularity may explain persistent real-rootedness: for near-regular dense graphs, the Chudnovsky-Seymour mechanism's local structure (neighborhoods are quasi-line) may hold even when global claw-freeness fails. Script: `04-computation/omega_spectral_fast.py`.

**T102** #anti-aut-NOT-always-involution #corrects-T095 #PROVED | certainty: high (PROVED THM-024) | source: opus-2026-03-06-S18
CORRECTS T095/T093 and PROVES existence: Not all anti-auts are involutions (counterexamples at n=6 with |Aut|>1), BUT every SC tournament has ≥1 involution anti-aut (THM-024). PROOF: Moon's theorem gives |Aut(T)| odd. H = ⟨Aut(T), σ⟩ has even order. By Cauchy, H has order-2 element; can't be in Aut(T) (odd), so it's in σ·Aut(T) = anti-auts. At even n the involution is fixed-point-free, creating n/2 orbits. This is the algebraic foundation for the SC maximizer mechanism. Script: `04-computation/anti_aut_involution_test.py`, THM-024.

**T103** #quasi-regularity-explained #johnson-graph #pseudorandomness | certainty: high (proved + verified n=5-20) | source: opus-2026-03-06-S18
WHY Omega_3(T) is quasi-regular: Two 3-cycles are adjacent iff they share ≥1 vertex. This depends on VERTEX SETS, not arc orientations. So Omega_3 is an induced subgraph of J(n,3) (Johnson graph complement). Since all 3-element subsets have identical intersection statistics, the degree of each 3-cycle in Omega_3 concentrates around the mean: E[deg] ≈ (C(n,3) - C(n-3,3))/4 - 1. The CV (coefficient of variation) of degree is O(1/√m) → 0 as n grows, forcing λ_max/avg_deg → 1. Verified computationally: CV drops from 0.05 (n=6) to 0.03 (n=20), matching the empirical λ_max/avg_deg ≈ 1 + CV². The quasi-regularity is a CONSEQUENCE of Johnson graph regularity + pseudorandom 3-cycle distribution. Script: `04-computation/omega_quasireg_proof.py`.

**T104** #hereditary-maximizer #even-odd-dichotomy #vertex-deletion | certainty: CORRECTED (S18g) | source: kind-pasteur-2026-03-06-S18f, S18g
HEREDITARY MAXIMIZER: REGULAR-ONLY AT ODD n (CORRECTED from "all maximizers"). At odd n, only REGULAR maximizers (score (k,...,k)) are hereditary. At n=5: 24/64 hereditary (regular) vs 40/64 non-hereditary (score (1,2,2,2,3)). At even n: 0% hereditary. See MISTAKE-010. Deletion spectrum NOT always constant: n=5 non-regular maximizer has del_hs=[3,5,5,5,3]. Two types at n=6: Type A (del=11) and Type B (del=13), both vertex-homogeneous.

**T105** #AAT-spectral #regularity-H | certainty: high (exhaustive n=5,6) | source: kind-pasteur-2026-03-06-S18f
STRONG SPECTRAL CORRELATE OF H: corr(H, lambda_1(AA^T)) = -0.97 at n=5 and -0.96 at n=6. H-maximizers have the SMALLEST leading eigenvalue of AA^T (most spectrally regular). Equivalently, maximizers minimize the spectral gap of the common-successor matrix. This is consistent with regular/Paley tournaments being maximizers: AA^T has equal diagonal entries (= out-degrees) and the most uniform off-diagonal structure.

**T106** #fano-paley #steiner-triple #combinatorial-design | certainty: high (PROVED for T_7, generalized) | source: opus-2026-03-06-S4
The 14 cyclic triples of Paley T_7 decompose into EXACTLY TWO COPIES of the Fano plane PG(2,2). The 7 vertex-disjoint 3-cycle pairs each take one line from each Fano copy. This perfectly explains alpha_2=7: each omitted vertex uniquely determines one disjoint pair. More generally, the c_3 = p(p^2-1)/24 cyclic triples of T_p form a 2-(p, 3, (p+1)/4) BIBD. Values: p=3: lambda=1 (trivial), p=7: lambda=2 (two Fano planes), p=11: lambda=3 (2-(11,3,3)), p=19: lambda=5. |Aut(T_7)|=21 divides |Aut(Fano)|=168. Transitive action of Aut on cyclic triples confirmed at p=7,11.

**T107** #transfer-symmetry-polynomial #tournament-constraint #INV-001 | certainty: high (PROVED n=4,5,6,7 symbolically) | source: opus-2026-03-06-S4
BREAKTHROUGH for INV-001: The transfer matrix symmetry M[a,b] = M[b,a] is a POLYNOMIAL IDENTITY modulo the tournament constraint T[x,y]+T[y,x]=1. Proved by symbolic computation (sympy) at n=4,5,6,7. With INDEPENDENT arc variables, M[a,b]-M[b,a] is nonzero (12 terms at n=4, 48 at n=5). After tournament substitution T[j,i]=1-T[i,j], all terms cancel EXACTLY. The identity does NOT hold for general digraphs. Combined with path reversal M_{T^op}[i,j]=(-1)^{n-2}M_T[j,i], symmetry is equivalent to M_{T^op}=(-1)^{n-2}M_T. A conceptual proof should exploit the tournament constraint algebraically.

**T108** #paley-deletion-p19 #maximizer-conjecture | certainty: high (computed) | source: opus-2026-03-06-S4
H(T_19 - v) = 117,266,659,317 for all v (Aut transitivity). Self-complementary deletion scores (8^9, 9^9). Extends Paley deletion conjecture: T_3-v->a(2)=1, T_7-v->a(6)=45, T_11-v->a(10)=15745, T_19-v->a(18)=117266659317 (conjectured). Cannot verify against OEIS A038375 (only goes to a(11)=95095). Diff H(T_19)-H(T_19-v)=1055429087598, so sum mu-weighted cycles through v = 527714543799.

**T109** #omega3-complement-matching #real-roots #structural | certainty: high (exhaustive n=6, counterexample n=7) | source: opus-2026-03-06-S18
OMEGA_3(T) COMPLEMENT IS A MATCHING for n<=6 (exhaustive: 31088/31088 at n=6). Meaning: each 3-cycle is disjoint from AT MOST one other 3-cycle. This gives I(Omega_3, x) = 1 + mx + px^2 with discriminant m^2-4p >= 0 (trivially, since p <= m/2). FAILS at n=7: 75.3% of tournaments have a 3-cycle disjoint from 2+ others. At n=8: 96.7% fail. The structural transition at n=7 is because vertex sets {0,1,2} and {3,4,5} can be disjoint 3-cycles while {0,1,2} and {4,5,6} can also be a disjoint pair — this requires n >= 7 to have overlapping disjoint pairs. The degree of I(Omega_3,x) stays <= 2 through n=8, but the complement structure is richer than a matching. Scripts: `04-computation/determinantal_ip.py`.

**T110** #turan-newton #real-roots-n9-n11 #proof-strategy | certainty: high (proved) | source: opus-2026-03-06-S18
TURAN PROVES NEWTON'S FIRST INEQUALITY for n<=11: At n=9-11, alpha(Omega_3) = 3, so the disjoint-pair graph is K_4-free (can't have 4 pairwise disjoint 3-cycles since 12>n). Being triangle-free (which is K_3-free, even stronger), Turan gives a2 <= c3^2/4 < c3^2/3, proving a1^2 >= 3*a2 (Newton 1 for degree 3). This extends THM-021 from n<=8 (where degree <= 2 makes it trivial) to n<=11 (degree = 3, Newton is non-trivial). For n>=12 (alpha>=4), Turan alone fails (ratio > 1). Actual min ratios: 2.14 (n=9), 1.92 (n=10), 1.65 (n=12), 1.35 (n=15), 1.19 (n=18) — always above 1 but approaching it. Real-rootedness tested: 0 failures through n=21, 1470 samples, degrees up to 5. Scripts: `04-computation/omega3_real_roots_fast.py`.

**T111** #real-roots-DISPROVED #n9-counterexample #newton-inequality | certainty: high (PROVED) | source: opus-2026-03-06-S18
REAL-ROOTEDNESS OF I(Omega(T), x) DISPROVED AT n=9. Counterexample: score sequence [1,1,3,4,4,4,6,6,7]. Full Omega has 94 directed odd cycles (12 tri, 40 pent, 36 hept, 6 ham). I.P. = [1, 94, 10, 1]. Newton's 2nd inequality fails: a_2^2=100 < a_1*a_3*3/2=141. Roots: one real (-0.011), two complex (-4.995 +/- 8.303i). H(T)=I(Omega,2)=237 (OCF unaffected). The Omega_3 restriction also fails: I(Omega_3,x) = [1,12,6,1], disc=-1323. The failure is caused by extreme density of Omega (avg degree 92.8/93) leaving very few independent sets. Earlier sampling missed this because: (1) prior full-Omega tests at n=9 used too few samples, (2) Omega_3-only tests miss the structure. The Chudnovsky-Seymour bound (claw-free => real roots, n<=8) is SHARP. THM-020 remains valid for n<=8. See THM-025. Scripts: `04-computation/full_omega_n9_correct.py`, `04-computation/real_roots_n9_verify.py`.

**T112** #omega3-also-fails #discriminant-negative #degree3 | certainty: high (PROVED) | source: opus-2026-03-06-S18
Even the 3-CYCLE-ONLY polynomial I(Omega_3(T), x) fails real-rootedness at n=9. Same counterexample: I(Omega_3, x) = 1 + 12x + 6x^2 + x^3 with disc = -1323 < 0. CRITICAL SUBTLETY: Newton's inequalities BOTH hold (a1^2=144 >= 3*a2=18; a2^2=36 >= a1*a3*3/2=18) yet the discriminant is negative. Newton's inequalities are NECESSARY but NOT SUFFICIENT for real-rootedness of degree-3 polynomials. The discriminant = 18a1a2a3 - 4a1^3a3 + a1^2a2^2 - 4a2^3 - 27a3^2 is a more complex quantity. The structural cause: 11 of 12 three-cycles share vertex 3, making Omega_3 star-like. Only cycle (0,4,6) avoids vertex 3 and is disjoint from 5 others. This creates a2=6 disjoint pairs but a3=1 triple, a "lopsided" structure where few independent sets exist despite many vertices.

**T113** #R-minimization #deletion-ratio #variational #OCF-consequence | certainty: REFUTED at n=7 | source: kind-pasteur-2026-03-06-S18g
DELETION-SUM RATIO FORMULA (THM-025): R(T) = n - E_weighted[|U(S)|] is PROVED. But the conjecture that H-maximizer minimizes R is FALSE at n=7 (tournaments with H=123 achieve R=1.585 < R(max)=5/3). The formula itself remains a clean identity connecting vertex deletion to independent sets. The breakdown at n=7 occurs because some non-maximal tournaments have very efficient vertex coverage despite lower total cycle count. The formula is useful structurally but NOT as a variational characterization of the maximizer.

**T114** #hopf-algebra-comultiplication #subset-convolution #OCF-structure | certainty: high (confirmed) | source: opus-2026-03-06-S5
Grujić-Stojadinović (arXiv:2402.07606) Hopf algebra comultiplication Δ([T]) = Σ_S [T|_S] ⊗ [T|_{V\S}] is EXACTLY our subset convolution. The Hopf algebra deletion property for cycles (U_X = Σ (-1)^{|S|-1} U_{X\S} for edge subsets forming a cycle) could relate H(T) to H(T-v) via Claim A. The Hopf antipode encodes Berge's theorem. See web-synthesis-opus-S5.md.

**T115** #feng-reversibility #tournament-constraint #transfer-symmetry | certainty: medium (connection identified) | source: opus-2026-03-06-S5
Feng's Dual Burnside (arXiv:2510.25202) proves Q=AB is reversible/symmetric under detailed balance. Our transfer matrix has AB structure. The tournament constraint T[x,y]+T[y,x]=1 plays the role of detailed balance condition. This could prove INV-001 (transfer matrix symmetry) for all n by showing tournaments are "reversible digraphs" in Feng's framework.

**T116** #DRT-skew-Hadamard #paley-equivalence #maximizer | certainty: high (classical result) | source: opus-2026-03-06-S5
Doubly regular tournaments ↔ skew Hadamard matrices (Reid-Brown 1972). Paley T_p is the canonical DRT. DRTs have the most balanced cycle distribution. Nozaki-Suda (arXiv:1202.5374) characterize skew Hadamard via tournament spectra of size n-2. This connects our spectral regularity finding (corr(H,λ₁)=-0.97) to DRT theory and potentially explains Paley deletion maximality.

**T117** #adler-alon-ross #asymptotic-e #random-regular | certainty: high (published 2001) | source: opus-2026-03-06-S5
Adler-Alon-Ross proved max H(T) ≥ (e-o(1)) · n!/2^{n-1}. Our ratios H(T_p)/(p!/2^{p-1}) = 2.00, 2.40, 2.44, 2.53, 2.56 converge toward e≈2.718. Paley tournaments are quasi-random, consistent with the random regular tournaments that achieve this asymptotic. Open: does H(T_p)/(p!/2^{p-1}) → e as p → ∞?

**T118** #irving-omar #matrix-algebra #det-per-formula | certainty: high (published 2025) | source: opus-2026-03-06-S5
Irving-Omar (arXiv:2412.10572) give ham(D) = Σ_S det(Ā[S])·per(A[S^c]) and reprove OCF for tournaments as Corollary 20. Their walk generating function W_D(z) = det(I+zXĀ)/det(I-zXA) may connect to our transfer matrix. The det/per decomposition over subsets parallels our subset convolution structure.

**T119** #transfer-matrix-trace #proved #parity | certainty: high (PROVED) | source: opus-2026-03-06-S6
THM-027 (opus): The n×n transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S) B_b(R) satisfies tr(M) = H(T) for odd n, 0 for even n. PROOF: Bijection between (σ through S ending at a, τ through R starting from a) pairs and Hamiltonian paths; each path P contributes (-1)^{pos(a,P)} to M[a,a]. Summing over a gives sum_{k=0}^{n-1} (-1)^k per path. ALSO: the old claim M = [[1,0],[0,-1]] always is FALSE (MISTAKE-011). M values vary widely with tournament; only SYMMETRY M[a,b] = M[b,a] holds universally.

**T120** #bibd-cycle-maximization #paley-maximizer #directed-cycles #alpha1-dominates | certainty: high (PROVED at n=7) | source: kind-pasteur-2026-03-06-S18h
BIBD ARRANGEMENT MAXIMIZES TOTAL DIRECTED CYCLES, NOT DISJOINT PAIRS (THM-027-kind-pasteur). Among regular n=7 tournaments, the BIBD (Paley, uniform lambda=2) MINIMIZES alpha_2=7 (disjoint 3-cycle pairs) but MAXIMIZES total directed odd cycles (80 = 14+42+24). H = 1 + 2*80 + 4*7 = 189. Three rigid classes: (alpha_2=7, H=189, 240 tours), (alpha_2=10, H=171, 1680), (alpha_2=14, H=175, 720). The BIBD forces every 5-vertex subtournament to be the regular T_5 (2 directed Ham cycles each), yielding 42 directed 5-cycles (vs 28-36 for non-BIBD). H-maximization is driven by alpha_1 (linear term 2*alpha_1 >> quadratic 4*alpha_2).

**T121** #noncommuting-redei-berge #deletion-contraction #inductive-OCF | certainty: medium | source: kind-pasteur-2026-03-06-S19 (web search)
Mitrovic (arXiv:2504.20968, Apr 2025) introduces Rédei-Berge function in NONCOMMUTING variables satisfying deletion-contraction: W_X = W_{X\e} - W_{X/e}↑. This is the missing inductive structure: the commutative version lacks deletion-contraction. For tournaments: W_X = sum over odd-cycle permutations of 2^{psi} p_{Type}. Properties: W_X = W_{X^op}, product rule. Could enable an INDUCTIVE proof of OCF by reducing n-vertex tournaments to (n-1)-vertex via edge deletion/contraction. See INV-051.

**T122** #chromatic-redei-bridge #converse-redei #positivity | certainty: medium | source: kind-pasteur-2026-03-06-S19 (web search)
Mitrovic-Stojadinovic (arXiv:2506.08841, Jun 2025): chromatic function and Rédei-Berge function are "almost identical" at the poset level. This enables translating chromatic polynomial results to tournament Hamiltonian path problems. They prove a "converse of Rédei's theorem" and generalize the triple deletion property. The VAST literature on chromatic symmetric functions (Stanley, Shareshian-Wachs, etc.) could be imported to our setting via this bridge. See INV-052.

**T123** #savchenko-cycle-formulas #spectral-cycle-counting #phase-transition | certainty: high (published series) | source: kind-pasteur-2026-03-06-S19 (web search)
Savchenko has EXACT polynomial formulas for c_k(T) in regular tournaments: c5,c6 (JGT 2016), c7 (DM 2017), c8 (arXiv:2403.07629, 2024). KEY: c8(DRT_n) is INDEPENDENT of which DRT is chosen! Phase transition at n=39: DRTs have MORE 8-cycles for n≤35 but FEWER for n≥39 vs locally-transitive tournaments. This suggests our cycle-maximization mechanism (BIBD → max directed cycles → max H) may REVERSE at large n. Critical for understanding whether Paley remains H-maximizer asymptotically. See INV-053.

**T124** #komarov-mackey #5-cycle-formula #score-sequence | certainty: high (published JGT 2017) | source: kind-pasteur-2026-03-06-S19 (web search)
Komarov-Mackey (arXiv:1410.6828): exact formula for c5(T) in terms of edge score sequence. Max c5 ≈ (3/4)*C(n,5) achieved by almost all random tournaments. Could connect our alpha_1 analysis to score-sequence invariants. See INV-054.

**T125** #jerrum-patel-2026 #subdivided-claw #zero-free-regions | certainty: high (published JLMS 2026) | source: kind-pasteur-2026-03-06-S19 (web search)
Jerrum-Patel (JLMS 2026): zero-free regions for I(G,x) on H-free graphs. For H NOT a subdivided claw or path, H-free + max degree 3 can have complex roots. This means our Omega_3 real-rootedness (through n=20, 0 failures) requires tournament-specific structure beyond any fixed forbidden subgraph. The "what breaks it" classification could identify exactly which subgraphs Omega_3 avoids. See INV-056.

**T126** #terwilliger-DRT #n27-classification #algebraic-invariant | certainty: medium | source: kind-pasteur-2026-03-06-S19 (web search)
Herman (arXiv:2404.11560, 2024): Terwilliger algebras classify DRTs up to n=23 but FAIL at n=27 (237 non-isomorphic DRTs). If H(T) is a DRT invariant (constant across DRTs of same order), Terwilliger algebras would suffice. If H varies across non-isomorphic DRTs at n=27, this gives test cases. See INV-057.

**T127** #linial-morgenstern #cycle-minimization #spectral-methods | certainty: high (partially proved) | source: kind-pasteur-2026-03-06-S19 (web search)
Linial-Morgenstern conjecture: fixed c3 density → min c4 by random blowup of transitive. Proved d ≥ 1/36 using spectral methods (Ma-Tang arXiv:2011.14142). Extended to c_ℓ for ℓ ≢ 2 mod 4. This is the "dual" extremal problem to our maximization. Their spectral methods (eigenvalue-based cycle density bounds) could provide tools for proving Paley maximizes H. See INV-055.

**T044** #signed-adjacency #polynomial-identity #proof-strategy | certainty: high (verified n=3,4,5 as polynomial identity) | source: opus-2026-03-05-S4
THE EVEN-ODD SPLIT IS A POLYNOMIAL IDENTITY (not just over {0,1}). D(x) = F(x)-G(x) where F counts T-paths using i->j by position, G counts T'-paths using j->i by position. Then D(-1) = 0 holds for REAL-VALUED arc variables with x=-1 as the UNIQUE universal root. Equivalent to B(L_i, R_j) = B(L_j, R_i) where B is the alternating subset convolution. B is sigma-invariant (p_w->1-q_w, q_w->1-p_w) and EVEN in s-variables. Since max s-degree is 2, OCF reduces to proving all s-degree-1 terms vanish: C_w + D_w = 0 for each w, where C_w = dB/dp_w, D_w = dB/dq_w. Clean proofs at n=3,4; for n>=5 the identity is global. See 03-artifacts/drafts/signed-adjacency-identity.md.

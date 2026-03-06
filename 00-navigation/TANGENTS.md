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

**T044** #signed-adjacency #polynomial-identity #proof-strategy | certainty: high (verified n=3,4,5 as polynomial identity) | source: opus-2026-03-05-S4
THE EVEN-ODD SPLIT IS A POLYNOMIAL IDENTITY (not just over {0,1}). D(x) = F(x)-G(x) where F counts T-paths using i->j by position, G counts T'-paths using j->i by position. Then D(-1) = 0 holds for REAL-VALUED arc variables with x=-1 as the UNIQUE universal root. Equivalent to B(L_i, R_j) = B(L_j, R_i) where B is the alternating subset convolution. B is sigma-invariant (p_w->1-q_w, q_w->1-p_w) and EVEN in s-variables. Since max s-degree is 2, OCF reduces to proving all s-degree-1 terms vanish: C_w + D_w = 0 for each w, where C_w = dB/dp_w, D_w = dB/dq_w. Clean proofs at n=3,4; for n>=5 the identity is global. See 03-artifacts/drafts/signed-adjacency-identity.md.

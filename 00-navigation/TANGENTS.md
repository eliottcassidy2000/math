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
For T_11, the average c_k(T_11)/C(11,k) is an exact integer for ALL k >= 6 = (p+1)/2, but fractional for k < 6. Specifically: c_7/C(11,7)=12, c_8/C(11,8)=45, c_9/C(11,9)=201, c_10/C(11,10)=971, c_11/C(11,11)=5505 are all integers; while c_3/C(11,3)=1/3, c_4/C(11,4)=1/2, c_5/C(11,5)=9/7, c_6/C(11,6)=145/42 are not. Proposed structural theorem: C(p,k) | c_k(T_p) for k >= (p+1)/2, for Paley primes p ≡ 3 (mod 4). The regularity of T_p (each vertex has out-degree (p-1)/2) likely forces this integrality via symmetry. Proving this rigorously is an accessible target -- see OPEN-Q-013.

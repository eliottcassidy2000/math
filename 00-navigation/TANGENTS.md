# Tangents, Rabbit Holes & Novel Ideas

**Purpose:** Dense, scannable index of interesting directions that have emerged. Each entry is 2–4 lines. A new Claude reads this whole file quickly to find a jumping-off point. Add entries with the tag format shown. Do not expand entries here — create a separate file in `03-artifacts/drafts/` if you want to develop a tangent.

**Format:** `#tag1 #tag2 | [certainty: low/medium/high] | [source]`

---

## Combinatorics & Structure

**T001** #ballot-sequences #dyck-paths #distribution | certainty: medium | source: MASTER_FINDINGS F4
The number of internal signature patterns giving exactly k Type-II positions within an L-cycle window is C(L-2, 2k-1). This is suspiciously clean and suggests a connection to ballot sequences or Dyck paths. A combinatorial proof of this formula is missing and would be a clean standalone result.

**T002** #asymptotic #cycle-length #mu-weights | certainty: low | source: MASTER_FINDINGS F4
The average Type-II contribution per L-cycle window is (L-4)/4, growing with L. Does this lead to an asymptotic formula for Σ_C μ(C) as a function of cycle-length distribution? Could give insight into why Claim A holds in aggregate even when per-path identity fails.

**T003** #n5-mystery #5-cycles #per-path-identity | certainty: high (the mystery is real) | source: MASTER_FINDINGS F5
At n=5, 5-cycles through v DO exist (v→a→b→c→d→v uses all 5 vertices). Yet the per-path identity holds at n=5 (0 failures). At n=6 it fails (2,758/9,126 failures). WHY? Hypothesis: at n=5 every 5-cycle embedding is always "accompanied" by exactly the right 3-cycle embeddings to compensate its μ contribution. This delicate balance breaks at n=6. This is **the deepest mystery** in the current work.

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

## Computational

**T014** #n7-verification #random-sampling | certainty: high | source: LaTeX paper
Claim A is verified at n=7 by random sampling (3,500 pairs, 0 failures). Full exhaustive verification at n=7 would be computationally expensive (2^21 tournaments × 7 vertices = ~14M pairs) but could be worth attempting on a cluster for increased confidence before investing heavily in a proof.

**T015** #failure-characterization #n6 | certainty: medium | source: MASTER_FINDINGS F5
The per-path identity fails for 2,758/9,126 (≈30%) of (T,v,P') triples at n=6. What characterizes the failing triples? Is there a simple condition on the tournament structure (e.g., number of 5-cycles, structure of Ω(T)) that predicts failure? This characterization might hint at a corrected formula.

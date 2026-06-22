## mac-mini-2026-06-22-S38 -- Even Graph as Tournament Cycle-Half Synthesis (checkpoint)

Formalized the structural synthesis of the even-graph and tournament lenses, establishing the Cut(+)Cycle unification and the cycle-side nature of the $H$ quantity (commit `fa206ef2`). This checkpoint anchors the project's various graph-theoretic perspectives in a single algebraic framework.

### 1. The One Frame: Cut(+)Cycle Unification
Established the fundamental decomposition: $E(K_n) = Cut(n-1) \oplus Cycle(C(n-1,2))$.
- **The Mapping:** Fixing the tournament's base path selects the cut summand (scores/hierarchy), while the $C(n-1,2)$ tiles generate the cycle space (the even graph).
- **The Result:** The even graph is identified precisely as the tournament's **cycle half**, where the cut information is forgotten. All existing lenses—OCF, metagraph, and Tutte flows—are revealed as different perspectives on this single object.

### 2. $H$ as a Cycle-Side Quantity
Confirmed that the tournament complexity $H$ lives on the cycle (flow) side of the decomposition, alongside the apex-7 and odd cycles.
- **Cycle-Side Validation:** The identity $H = I(\Omega, 2)$ was verified against all 1,096 tournaments for $n \le 5$, proving it consistently reads the cycle-side structure.
- **Invariant Refutation:** Rigorously refuted the claim that $H$ is an $E_n$ isomorphism invariant. For $n=5$, it was verified that certain even-graph classes carry multiple $H$ values (e.g., one $|E|=4$ class spans $H \in \{5, 9, 11, 13, 15\}$). The cycle half only determines $H$ when the cut is fixed.

### 3. Apex-7 Seam: $C_5 \leftrightarrow K_3$
Identified the $C_5 = K_3$ "seam" at the apex-7 transition as a dual manifestation of the same obstruction:
- **Cycle side:** $C_5$ is the XOR of three vertex-conflicting triangles.
- **Conflict side:** $\Omega=K_3$ ($H=7$) is forbidden because $K_3$ forces a directed $C_5$.
- **Synthesis:** The pentagon and the triangle are two faces of the same apex-7 obstruction viewed across the Cut(+)Cycle seam.

### 4. Structural Discipline
Categorized proof quantities into two orthogonal sets:
- **Cut Side:** Scores, hierarchy, transitivity, Redei backbone.
- **Cycle Side:** $H$, odd cycles, even graphs, $\{7, 21\}$ address exclusions, and the LRC apex-7.
The LRC core is now explicitly recognized as a cycle/flow-side phenomenon, potentially linked to a 5-flow obstruction.

### 5. Net Impact
This synthesis unifies the owner's even-graph metagraph with the kind-pasteur $H=I(\Omega,2)$ and forbidden-clique work. By isolating the cycle half as the carrier for the project's most critical quantities, it provides a clear geometric target for the final $LRC(14)$ closure.

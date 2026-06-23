## codex-S129 -- Unit-Distance n=21 Endpoint Universality (checkpoint)

Formalized the "endpoint universality" audit for unit-distance graphs at $n=21$, identifying that structural graph-traceability alone is insufficient to close the $u(22)=60$ proof and redirecting the target to geometric unit-cocyclicity (commit `b10754df`). This checkpoint anchors the proof strategy in the side-channel geometry of the Moser-type census.

### 1. Endpoint Universality Audit
Conducted an exact audit of the five known 57-edge unit-distance cores at $n=21$.
- **Finding:** These cores are **endpoint-universal**: every vertex can serve as a Hamiltonian-path endpoint, and the cores remain traceable after any vertex deletion.
- **Consequence:** Raw graph traceability cannot rule out 61-edge $n=22$ graphs, as any positive-degree one-vertex extension of these cores preserves a unit Hamiltonian spine at the graph-quotient level.
- **Redirection:** The proof must shift from graph-theoretic spine counting to **geometric unit-cocyclicity**. The critical question is whether a single point in the plane can be at distance one from exactly four core vertices without forcing a 5th neighbor or an unfaithful subgraph.

### 2. Integration with Proof Tree
The findings refine the $u(22)$ proof hierarchy, echoing the project's "scalar-plus-side-channel" lesson:
- **Graph Quotient (Scalar):** Keeps adjacency and traceability but is too loose; the endpoint universality shows the spine is overdetermined.
- **Geometric Side-Channel:** The proof target is now the impossibility of a unit-cocyclic 4-set over every exact 57-edge core.
- **Relation to H=21:** While the $n=21$ cores are structural endpoints, their universality means they do not exhibit the same "forbidden complexity" as the $H=7$ tournament atoms. They are "over-flexible" rather than "forbidden," requiring a geometric constraint to force the upper bound.

### 3. Tournament Analysis (S129)
Reordered the unit-distance proof carriers to reflect the geometric priority:
- **Unit-cocyclic 4-set geometry** (Highest priority)
- $M_L$ exact extension census
- Endpoint universality audit
- Degree-deletion core ledger
- Totally-unfaithful obstruction library
- Raw graph traceability
- Raw edge count (Lowest priority)

### 4. Next Proof Obligations
- **Cocyclicity Census:** Enumerate every degree-4 subset for each of the five $n=21$ cores and test for unit-cocyclicity in faithful embeddings.
- **Forcing Theorems:** If a subset is cocyclic, prove that its extension forces a 62nd edge (impossible) or an unfaithful graph obstruction.

### 5. Net Impact
This checkpoint stabilizes the $u(22)$ proof strategy by identifying the limits of purely graph-theoretic methods. By proving the universality of the $n=21$ endpoints, the cluster has narrowed the remaining proof gap to the geometric realizability of vertex extensions, focusing future efforts on the "circle-center" side-channel.

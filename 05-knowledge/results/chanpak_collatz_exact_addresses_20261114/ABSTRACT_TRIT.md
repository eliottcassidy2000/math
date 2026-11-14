# Two abstract realizations of 3H+2 that cannot be literal tournament cycle-conflict graphs

**PROVED, structural.** The readout H=I(Ω,2) has a natural arithmetic move of exactly 3H+2 in ordinary graphs. Its obstruction is not a guessed scarcity of special integers: the required abstract graph fails tournament-specific cycle geometry. The statement is not an obstruction to getting the same scalar via a *different* graph.

Write Ω(T) for the intersection graph whose vertices are all directed odd cycles of T, adjacent when they share a tournament vertex. The [odd-cycle packing identity](https://arxiv.org/html/2412.10572v3#S3.SS4.Thmcorollary20) gives H(T)=I(Ω(T),2), where I is the independence polynomial. For any graph G, adjoining one isolated vertex triples I(G,2), then adjoining one universal vertex adds 2. Let G_0 be empty and let G_k be obtained from G_(k−1) by those two moves. Then

    I(G_k,2)=2·3^k−1,       |V(G_k)|=2k,       alpha(G_k)=k   (k≥1).

The independence statement holds inductively: the new isolated vertex can join an old maximum independent set, and the final universal vertex can occur only alone. G_k is connected. A second abstract graph, K_(k,k), has the same readout, since an independent set lies entirely on one side and the empty set must only be counted once. These two different graphs therefore encode the proposed shifted-ternary itinerary perfectly at the scalar level.

Both fail to be **full** tournament odd-cycle conflict graphs for k≥2. For the complete bipartite graph this follows at once from the [inherited universal-vertex theorem](../chanpak_cross_problem_theorems_20261114/CONFLICT_GRAPHS.md): every connected component of Ω(T) has a universal vertex, whereas K_(k,k) does not. The graph G_k does have a universal vertex, so this first test does not catch it.

The stronger test: if Ω(T) is nonempty and connected and alpha(Ω(T))=k, then

    |V(Ω(T))| ≥ floor((3k−1)^2/4).

Indeed all cycle-bearing vertices are in one nonsingleton strong component, by the inherited component theorem. A packing of k vertex-disjoint odd cycles puts at least 3k vertices in that component, say n in total. [Moon, Corollary 2.1](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on-subtournaments-of-a-tournament.pdf) gives at least n−l+1 directed l-cycles for every 3≤l≤n in a strong tournament. Summing over odd l gives floor((n−1)^2/4), and this is at least the claimed bound. For k≥2 it exceeds 2k, ruling out G_k as well. At k=1 both abstract graphs are K2, physically realized by the strong four-vertex tournament with exactly two triangles and H=5.

This illustrates two distinct ways to lose a putative arithmetic compiler: a chosen-subgraph model can omit the additional odd cycles, and an exact integer formula can hide the tournament-level reason extra terms must appear. The [induced-source theorem](TRIANGLE_TARIFF.md) addresses instead what happens when one tries to retain the old tournament and append a physical ternary block. A successful repair can look for a graph Ω' with the same numerical value but a different packing geometry, not just for a new presentation of either forbidden graph.

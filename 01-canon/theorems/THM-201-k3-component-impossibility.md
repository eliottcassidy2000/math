# THM-201: K₃ Component Impossibility in Ω(T)

**Status:** PROVED
**Proved by:** opus-2026-03-14-S71g
**Dependencies:** THM-200 (H=7 impossibility), directed cycle enumeration

## Statement

For any tournament T, the conflict graph Ω(T) does not contain K₃ (complete graph on 3 vertices) as a connected component.

Equivalently: there do not exist 3 directed odd cycles in T such that
(a) all three pairwise share at least one vertex, AND
(b) no other directed odd cycle in T shares a vertex with any of these three.

## Proof

Suppose for contradiction that C₁, C₂, C₃ are three directed odd cycles forming a K₃ component of Ω(T). Since they are isolated from all other cycles:

**Step 1.** Each Cᵢ must be a 3-cycle.

If any Cᵢ were a 5-cycle (or longer), the 5 vertices of that cycle would contain ≥3 directed 3-cycles (verified exhaustively at n=5). These additional cycles would share vertices with Cᵢ and hence connect to the K₃ in Ω, contradicting component isolation.

**Step 2.** Three pairwise-sharing 3-cycles span at most 7 vertices.

Let S = V(C₁) ∪ V(C₂) ∪ V(C₃). Each Cᵢ uses 3 vertices. With pairwise sharing, |S| ≤ 7 (maximum when all share a common vertex: 1 + 2×3 = 7).

**Step 3.** The subtournament T|_S has exactly 3 directed odd cycles.

The three planted cycles C₁, C₂, C₃ exist on S. Any additional cycle on S vertices would also be a cycle of T, and would share vertices with at least one of C₁, C₂, C₃ (since all cycles on S touch the same bounded vertex set). This contradicts the K₃ component assumption (no other cycle connects to these three).

**Step 4.** No tournament on |S| ≤ 7 vertices has exactly 3 directed odd cycles with Ω = K₃.

- |S| ≤ 6: No tournament has exactly 3 directed odd cycles (verified exhaustively at n = 3, 4, 5, 6). Contradiction.
- |S| = 7: 3360 tournaments have exactly 3 directed odd cycles. In ALL 3360 cases, the conflict graph is P₃ (path on 3 vertices), not K₃. Contradiction.

∎

## Corollary 1: SCC Product Constraint

For non-strongly-connected tournaments, H(T) = ∏ H(SCCᵢ). Since H(SCC) = 7 requires Ω(SCC) = K₃ (connected = the whole graph), which is impossible by this theorem, no SCC can have H = 7. Therefore H(T) is never divisible by 7 through SCC decomposition.

More precisely: the set of achievable H values is closed under products (by SCC decomposition), and 7 is not a "prime factor" in this multiplicative system.

## Corollary 2: Broader Forbidden Products

Any H value that can only be expressed as I(G, 2) where every such G has a K₃ component is forbidden. This potentially extends the H=7 prohibition to other values whose only graph realizations involve K₃ components.

## Key Scripts

- `04-computation/h7_theorem.py` — Main verification
- `04-computation/h7_n7_check.py` — Exhaustive n=7 check
- `04-computation/knacci_simplex_cuboid.py` — H-spectrum analysis showing forbidden values {7, 21, 63}

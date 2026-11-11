# Fourth extension: the universal quantitative defect inequality

**Status: PROVED for every finite tournament and marked vertex**, using the inherited odd-cycle formula and standard tournament endpoint facts. The stronger single-cycle inequality is also proved. **Chan–Pak remains OPEN:** the theorem controls excess; it does not cancel it or realize arbitrary prescribed counts.

This supersedes the OPEN labels for the two inequalities in Session260 `bridges/DEFECT_QUANTITATIVE.md` and for (SC) in this packet's `CHAIN_COMPONENTS.md`. Those earlier records remain untouched. No external novelty claim is made.

## 1. Main statements

Let T be a tournament, z a marked vertex, U=T−z, A=N_T⁺(z), B=N_T⁻(z). Let G be the bipartite graph on A⊔B with edge ab exactly when a→b in U. Let P(U;A,B) count all directed simple paths of positive even order beginning in A and ending in B. Write H for Hamilton-path count and d for dominating-set count. Then

    d(G)−1 ≤ 2P(U;A,B) ≤ H(T)−H(U).                     (1)

Consequently,

    H(T)−d(G) ≥ H(U)−1 ≥ 2t₃(U),                        (2)

where t₃(U) counts directed triangles. Thus both global quantitative targets of Session260 follow. The second inequality in (2) is the inherited positive odd-cycle formula, not a new triangle estimate.

The proof works through an exact positive expansion of domination by induced even components. It neither reserves Hamilton paths in the old domination injection nor requires a universal Hall/exchange lemma.

## 2. Inclusion-exclusion becomes an even-component expansion

For **any finite simple undirected graph** G=(V,E), not necessarily bipartite, let C(W) denote the set of connected components of G[W], κ(W)=|C(W)|, and call W admissible if every component has an even number of vertices. The empty set is admissible with κ=0. Then

    d(G)=Σ_{W⊆V admissible} 2^{κ(W)}.                    (3)

### Proof

Inclusion-exclusion on undominated vertices gives

    d(G)=Σ_{X⊆V}(−1)^{|X|} 2^{|V\N[X]|}.

Expand the power of two by choosing Y⊆V\N[X]. Thus d is the signed sum over ordered disjoint pairs (X,Y) with no graph edge between X and Y, the sign being (−1)^{|X|}. Group these pairs by W=X⊔Y. Having no edge between the parts means that X is a union of connected components of G[W]. Hence the signed sum for this particular W is

    ∏_{C∈C(W)} [1+(−1)^{|C|}].

It vanishes if any component has odd order; otherwise it equals 2^{κ(W)}. This proves (3), including the empty graph. This is an exact cancellation followed by a positive expression, not an estimate of alternating terms. In particular it explains domination parity directly.

## 3. Endpoint lemma: no isolated bipartite vertex is enough

**Lemma.** Let R be a nonempty tournament with partition A₀⊔B₀. Suppose every vertex of A₀ has an outgoing neighbour in B₀, and every vertex of B₀ has an incoming neighbour in A₀. Then R has a Hamilton path beginning in A₀ and ending in B₀.

These hypotheses say that the bipartite graph of A₀→B₀ edges has no isolated vertices. They do not require that graph to be complete, a chain graph, or connected.

### Proof

The strong components of a tournament form a total order with all edges pointing from earlier to later components. The first strong component must meet A₀. Otherwise it lies wholly in B₀, but every one of its vertices has an incoming edge from A₀ outside it, contradicting that it is first. Dually, the last strong component must meet B₀.

If these are different components, choose a Hamilton path of the first component beginning in A₀ and one of the last component ending in B₀; choose arbitrary Hamilton paths in the intermediate components and concatenate. A singleton component causes no difficulty. In a nontrivial strong component, the required one-sided endpoint path is obtained by cutting a Hamilton cycle at the chosen endpoint (the standard Camion consequence already used in Session260).

If there is just one strong component, both colours occur. Take a Hamilton cycle. Traversing it, some edge crosses from B₀ to A₀. Cut that edge: the resulting Hamilton path begins in A₀ and ends in B₀. There is no one-vertex exceptional case under the lemma's hypotheses. This proves the lemma.

## 4. Every even-component fibre supplies enough actual paths

Fix a nonempty admissible W in (3), and let its components be C₁,…,C_c. Every C_i has even order at least two, so its induced bipartite graph is connected and has no isolated vertices. The endpoint lemma supplies a Hamilton path of U[C_i] beginning in A∩C_i and ending in B∩C_i. Choose one deterministically, for example the lexicographically least such path using fixed labels.

For distinct components i,j, there are no graph edges between them. Thus every B∩C_i vertex points to every A∩C_j vertex in U. The selected component paths can therefore be concatenated in **all c! linear orders**. Each concatenation has even total order, begins in A and ends in B, and uses exactly W.

They are distinct: the vertex set recovers W, and the sequence of component labels recovers the order. Consequently, if p(W) is the number of Hamilton paths of U[W] beginning in A and ending in B,

    p(W) ≥ κ(W)!   for every nonempty admissible W.       (4)

The elementary inequality 2c!≥2^c holds for c≥1: equality holds at c=1,2, and multiplying by each next integer at least two proves the rest. Therefore

    2P(U;A,B)
      =2Σ_{W≠∅, |W| even}p(W)
      ≥2Σ_{W≠∅ admissible}κ(W)!
      ≥Σ_{W≠∅ admissible}2^{κ(W)}
      =d(G)−1.

This proves the first inequality in (1). Close each even A-to-B path through z and use that odd cycle as an OCF packing, leaving other vertices unused. These distinct objects have weight two, proving the second inequality in (1).

All counts are of literal paths, directed cycles, or dominating sets. The passage from dominating sets to (3) is a proved cardinal identity; no unsupported bijection between the two object classes is claimed.

## 5. Exact nonnegative quantitative-defect fibres

The inherited marked OCF decomposition is more precise than the second inequality in (1):

    H(T)−H(U)
      =2Σ_{W≠∅, |W| even}p(W)H(U−W).                    (5)

Indeed a packing containing z has one unique z-cycle, equivalently one even A-to-B path on W, and an arbitrary packing on the complementary vertices. The cycle has weight two, and the complementary weights sum to H(U−W).

Put ε(W)=1 if W is admissible and 0 otherwise. Combining (3) and (5) yields the **exact identity**

    H(T)−H(U)−d(G)+1
      =Σ_{W≠∅, |W| even}
          [2p(W)H(U−W)−ε(W)2^{κ(W)}].                  (6)

Every bracket is nonnegative individually. If ε=0 it is immediate. Otherwise (4), H≥1, and 2κ!≥2^κ prove it. This is a new positive refinement of the original domination-cut defect formula: the fibres are indexed by the vertex set of the unique newborn odd cycle rather than the cut of a spanning path at z.

In particular,

    H(T)−H(U)−d(G)+1
      ≥ Σ_{W≠∅ admissible}[2κ(W)!−2^{κ(W)}].             (7)

The previous component and chain-component theorems are now specific selections of these universal fibres.

### Equality boundary

Equation (6) gives a complete, exact fibre criterion. Equality in the first bound of (2) holds if and only if:

1. every nonempty admissible W has κ(W)≤2;
2. for each such W, p(W)=κ(W)! and H(U−W)=1;
3. every even nonadmissible W has p(W)=0.

Necessity follows by making each nonnegative bracket zero. For κ=1,2, the minimum p is respectively1,2, so equality forces complementary H=1 and minimum p. For κ≥3 the factorial inequality is strict. Sufficiency follows by substitution in (6).

This is a precise structural/counting boundary, not a claim of a simple forbidden-subgraph classification. One graph-only consequence is useful: if G contains an induced 3K₂, its six endpoints give κ=3, so the quantitative slack is at least 12−8=4. More generally an induced kK₂ gives the single-fibre lower bound 2k!−2^k and all its matching-subset fibres may also be summed. The condition “no induced 3K₂” alone is not asserted sufficient for equality.

## 6. Fourth random-source wave and the actual turn in the proof

Seed **2611710**. Pool: all `.md` paths under `01-canon/theorems/` in `git ls-tree -r --name-only HEAD` order; size **4,130**. The canonical corpus at the sampling point agrees with f92e070445e43c59f1252f2979db564df88d7b60. Three sampled sources, read through their relevant initial material:

- `01-canon/theorems/THM-2702-polynomial-graph-coordinate-projection-keller-classification.md`
- `01-canon/theorems/THM-1179-six-wall-slow-gap-pair-gcd-obstruction.md`
- `01-canon/theorems/THM-4132-fixed-body-exceptional-scale-two-lrc14-completion.md`

**Coordinate change / nilpotence source → grouped inclusion-exclusion.** The transferable question was whether a signed formula becomes triangular or factorized in another coordinate system. Source map: polynomial graph coordinates adapted to a derivation. Target map: (X,Y)↦(W=X⊔Y, selected induced components). Preserved: exact signed mass. Lost: the original individual dominating-set interpretation, restored by the proved inclusion-exclusion equality rather than an invented bijection. Cheap decisive test: each component contributes 1+(−1)^{|C|}. Outcome: complete positive expansion (3). No Jacobian or locally nilpotent derivation theorem is invoked in the proof.

**Adaptive graph debt → component-rank factorial capacity.** Source uses a pointwise graph/rank inequality before integrating. Target uses the pointwise fibre inequality 2c!≥2^c before summing over W. Preserved: local nonnegative budget implies global comparison. Lost: the original period and gcd data, which have no role here. Cheap test: c=1,2 are exact, c=3 produces slack4. Outcome: universal nonnegative fibres (6), not merely a total-count inequality.

**Multiple physical lifts → endpoint choices.** Source's warning is that no single fixed lift need work everywhere. Target initially seemed to require a fixed global A-to-B Hamilton path choice compatible with all subsets. That is unnecessary: choose a canonical endpoint path separately in each induced component and retain W. Cheap test: examine first and last strong components; their colours are forced by the no-isolates hypothesis. Outcome: the endpoint lemma, which removes the chain-order restriction entirely. No LRC completion theorem is transported.

### Routes abandoned before the successful regrouping

- A determinant encoding would need a cancellation rule distinguishing one cycle through z from complementary cycle packings. Directly taking absolute coefficients would lose the signs and had no justified inequality. It was not pursued as a theorem.
- A direct involution on dominating sets would need one fixed point and a recoverable path label for every other pair. The even-component expansion bypasses that unconstructed involution; it does not assert one was found.
- The source-replacement fallback to chain components would already face the inherited count29 obstruction. It is unnecessary for proving (1), and this pass does not claim a count-preserving replacement.

## 7. Focused controls and reproducibility

Run `python3 05-knowledge/results/session261_creative_chanpak_20261111/bridges/extension4_probe.py`.

The checker imports no earlier author counting modules. It directly enumerates dominating sets and verifies (3) for every labelled simple graph of orders0..5, not just bipartite graphs. For the new endpoint/path fibre statement, it uses literal permutation enumeration on all tournament bases through order4 and every A/B colouring, plus35 seeded bases of each order5,6 and every colouring. It checks (4), the factorial lower bound, and the single-cycle inequality. The receipt states the exact scopes and counts.

The all-order argument is the proof in sections2–5. These controls test its delicate new interfaces; they are not an expanded spectrum census, a substitute for proof, or a re-audit of the inherited odd-cycle formula.

## 8. What remains genuinely open

The domination-to-tournament interface now has a universal positive quantitative defect, not just an injection and a conjectured bound. But positive defect is not exact-count preservation. The same fixed-input one-apex equality obstruction still applies: H(T)=d(G) only when U is transitive. To prove Chan–Pak, one must still change the source count, change the interface, or construct a legal correction with the needed state. This extension closes the requested quantitative bridge and no more.

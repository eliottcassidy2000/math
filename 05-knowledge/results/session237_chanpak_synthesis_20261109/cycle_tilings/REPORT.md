# Session237: universal induced cycle carriers and the completion tax

Status: **PROVED here** (elementary propositions below), **FINITE-EXACT** computational controls. Independent audit PASS WITH TWO MINOR WORDING CORRECTIONS, applied; THM-5090. Chan–Pak completeness remains **OPEN**. No new all-order H-value interval is claimed.

## Outcome first

Every finite simple graph G occurs as an **induced subgraph** of the odd-cycle conflict graph Ω(T) of a tournament. An explicit strengthened construction uses

    |V(T)| = |E(G)| + 5|V(G)| + #{v : deg_G(v) is odd}
           <= |E(G)| + 6|V(G)|.

It simultaneously installs the corona G∘K1 as an induced subgraph of Ω(T), and hence gives the coefficientwise inequality

    I(Ω(T),x) >=_coeff (1+x)^m I(G,x/(1+x)) + m x,   m=|V(G)|.

The right side is a polynomial, interpreted as the independent-set sum below. In particular

    H(T) >= 3^m I(G,2/3) + 2m >= 3^m+2m.

This is a carrier embedding, **not** a count-preserving reduction. The theorem proves that no nontrivial hereditary graph property can hold for every Ω(T). In particular, looking for a forbidden-induced-subgraph characterization of the full family of realizable Ω(T) is doomed: its hereditary closure is the class of all finite graphs. The full realizability obstruction resides in forced extra cycles and their packing contributions.

This sharpens the conceptual lesson of the old induced C5 and claw examples. It also diagnoses why the sharp even-graph domination theorem cannot be imported through a generic graph-intersection representation: embedding supplies only a lower bound, whereas exact spectrum coverage needs control of **all** newly created cycles.

## Inheritance and concept board

Anchor: build a rigorous transfer toward tournament H-spectrum completeness.
Niche: the hereditary closure of odd-cycle conflict graphs.
Wildcard: poset/antichain and tiling encodings, distinguishing exact coordinates from count-preserving maps.

Five concepts kept live:
1. Odd-cycle packings and the positive independence polynomial.
2. The full cycle carrier versus a selected induced subcarrier.
3. Tiling coordinates and completion of a partially oriented complete graph.
4. Domination/neighbourhood-complex inclusion–exclusion and graph offsets.
5. Poset antichains versus general independent sets.

Closest proved mechanism: OCF, H(T)=I(Ω(T),2), and general-base packing laws in THM-4940, `01-canon/theorems/THM-4940-odd-cycle-packing-counts-in-every-base.md`.

Canonical hostile: I(K3,2)=7 although Ω(T)=K3 is impossible. Three selected intersecting cycles do not certify that those are all cycles. The older K3 discussion is grounded in THM-029/THM-343 and the all-order 7/21 exclusion THM-1370, not in a historical reflection's wording.

Corrected near miss: Ω(T) need not be perfect (THM-019, `01-canon/theorems/THM-019-omega-perfectness.md`), nor is I(Ω(T),x) always real-rooted (THM-020/THM-025). The old THM-019 file contains legacy empirical/log-concavity prose; only its explicit induced-C5 construction is used here.

Least-used sidecar: the **list of unselected directed odd cycles**, with their vertex supports and conflicts. A graph induced in Ω does not retain this data. It is exactly what an H-preserving synthesis needs.

## Proposition A: universal induced realization with a corona

Let G be any finite simple graph with vertices 1,...,m and edge set E; m>=1. For each edge e introduce a vertex s_e. For each i introduce 5+(deg(i) mod 2) private vertices

    p_i, q_i, r_i, a_i, b_i [, c_i].

List the incident edge vertices s_e in any fixed order. Form a directed cycle C_i by following this list and then the private list, cyclically. Its length is

    ell_i = deg(i)+5+(deg(i) mod 2),

which is odd and at least five. Also impose r_i -> p_i, giving the private directed triangle Δ_i=(p_i,q_i,r_i).

These constraints are compatible. Two distinct C_i share exactly the one vertex s_{ij} if ij is an edge, and share no vertex otherwise. Thus they cannot share an unordered pair, so no prescribed cycle arcs conflict. The added chord r_i -> p_i is not a cycle edge because each private list has at least five vertices. Complete every remaining unordered pair arbitrarily to a tournament T.

The selected cycles C_i induce G inside Ω(T). The triangles Δ_i are pairwise disjoint, meet C_i, and meet no C_j with j≠i. Consequently the 2m selected cycles induce exactly G∘K1. The vertex count is the formula above. This construction works for **every** completion of the unassigned pairs, not merely the lexicographic completion used by the program.

The empty graph on zero vertices has the trivial empty selected carrier; it can be realized inside any transitive tournament. No nonempty-construction bound is needed for this boundary case.

### Why the cycle vertices really are distinct

Ω has directed cycles, up to cyclic rotation, as vertices; it is not a quotient by cycle length or support. The selected C_i have distinct supports, the Δ_i have private supports, and |C_i|>=5 differs from |Δ_i|=3. Thus all 2m selected vertices are distinct even under the cycle-order-forgetting/support-only convention. No cycle multiplicity has been discarded.

## Proposition B: the positive completion tax

An independent set in the selected corona first chooses an independent set S of G, selecting C_i for i∈S. It may then select any of the m-|S| available private triangles. Therefore its independence polynomial is

    J_G(x) = Σ_{S independent in G} x^|S| (1+x)^(m-|S|)
           = (1+x)^m I(G,x/(1+x)).

Every selected packing is a genuine packing in T, giving I(Ω(T),x)>=_coeff J_G(x).

There are at least m additional cycles. The tournament induced on V(C_i) is strongly connected because it contains the spanning directed cycle C_i. Every vertex of a strongly connected tournament lies in a directed triangle: take a shortest directed cycle through that vertex; if its length exceeds three, the chord to the second vertex after it produces a shorter cycle through it in either orientation. Apply this to a_i. It yields a triangle Γ_i contained in C_i and containing a_i. It differs from Δ_i and from every long selected cycle. For i≠j, Γ_i and Γ_j cannot be the same cycle, because V(C_i)∩V(C_j) has size at most one. None equals any Δ_j. Their m singleton packings contribute at least m x beyond J_G.

All remaining coefficients of the difference are nonnegative because every selected packing injects into the full packing set. This proves the stated coefficientwise inequality. For any real x>=0 it gives a numerical inequality; it is not asserted for negative x. The m x term is a conservative bound, not an optimal completion tax. For m=1 the lexicographic completion actually gives H=9, while J(2)+2=7.

### Exact decomposition that isolates the missing task

For any selected induced cycle carrier A⊆Ω(T), write B=V(Ω(T))\A. Partition independent sets by their intersection R with B. Then identically

    I(Ω(T),x) = Σ_{R independent in Ω(T)[B]} x^|R| I(A-N_A(R),x).

The R=∅ term is I(A,x). All other terms are nonnegative for x>=0. This standard independent-set partition identity gives the exact sidecar: the independent sets of extra cycles, together with their neighbour unions in A. Knowing only A is insufficient; knowing only the number of extra cycles is also insufficient.

For arbitrary positive fugacity, equality with I(A,x) occurs iff B is empty. At x=0 equality is automatic. Thus no cancellation inside the ordinary positive packing model can remove the unselected cycles.

## Corollaries and exact limitations

1. **No hereditary proper class.** If a graph property is inherited by induced subgraphs and every Ω(T) has it, then every finite graph has it. This excludes universal perfectness, claw-freeness, chordality, comparability, and incomparability as structural mechanisms on the *whole* conflict graph. It does not exclude order-bounded versions, restricted tournament classes, or a non-hereditary global constraint.
2. **No induced-obstruction characterization.** No nonempty collection of forbidden finite induced subgraphs characterizes the full realizable Ω family. The family itself is not all graphs (K3 is already impossible as the full carrier).
3. **No real-root conclusion from universality.** An induced subgraph with a non-real-rooted independence polynomial need not force the full polynomial to have non-real roots. The independent old THM-025 counterexample, not this embedding, handles that question.
4. **No minimal-order claim.** The quadratic vertex bound is explicit and convenient, not asserted optimal. In particular, the old C5-on-eight-tournament-vertices construction is much smaller than this uniform construction.
5. **No obstruction to a more clever reduction.** The completion tax belongs to this construction and, more generally, any positive selected-carrier inclusion with extras. It does not prove that an unrelated transformation of domination counts into tournament H-counts is impossible.

## Connecting the historical lanes without merging their predicates

### A. Tiling dictionary

THM-830 (`01-canon/theorems/THM-830-b3-deletion-deck-mirror-current-calculus.md`) identifies fixed-path tilings with labelled tournaments on path gaps. THM-5041, `01-canon/theorems/THM-5041-matrix-tiling-pair-dictionary.md`, identifies tournaments with complementary tiling pairs.

Our partial arc prescription is therefore a partial assignment of these tiling bits. Every completion produces the same selected conflict graph; the remaining bits control extra cycles. This is an exact coordinate transfer, not an independent H-preserving operation on posets or graph domination. The selected carrier forgets those completion bits and the full cycle sidecar restores their effect on H.

### B. Poset approach

Antichains of a finite poset are independent sets in its comparability graph. Since comparability graphs are perfect and Ω(T) is induced-universal, there is no representation of every Ω(T) as a comparability graph using its same cycle vertices. Likewise an incomparability-graph identification cannot hold universally. One may still encode the *family of packings* as a face poset, but counting its elements counts I(Ω,1), while H counts 2-coloured packings, I(Ω,2). The fugacity/colour sidecar cannot be dropped.

This does not preclude a different poset whose number of ideals, antichains, or linear extensions equals H. Such an equality would need its own map and is not supplied by the word “poset.”

### C. Latest even-graph theorem

THM-5040, `01-canon/theorems/THM-5040-sharp-even-graph-chan-pak.md`, and the independently credited THM-5055 prove the even-graph domination spectrum is exactly the odd numbers other than 3,7,13. Here “even” means the automorphism edge-reversal sign condition, **not** even vertex degrees. This terminological distinction is essential when reading old cycle-space reflections.

The proof uses positive domination counts but a complement-neighbourhood model:

    d(G)=2^n-D(complement G).

In the ordinary closed-neighbourhood formulation, inclusion–exclusion also gives

    d(G)=Σ_{S⊆V(G)} (-1)^|S| 2^(n-|N_G[S]|).

Proof: S indexes vertices required to be undominated; a chosen dominating-set candidate must then avoid the union N_G[S]. This formula is signed. Turning its terms into positive independent cycles would discard the signs and therefore change the count. The neighborhood complex and the positive cycle-packing complex are different carriers.

The operation check is decisive: a graph cone satisfies d(K1∨G)=2^n+d(G), whereas adjoining a universal source or sink to a tournament leaves H unchanged. Isolated graph padding and transitive tournament padding each preserve their scalar, but the next cone step does not commute. Any transfer of the graph offset proof must transport *its operation and state*, not just the oddness of its output.

A global exact value-preserving map from **all** even graphs to tournaments is impossible: even graphs realize 21 by THM-5040, while tournaments do not. A restricted map excluding that value, or a transformed scalar, remains possible. Equinumerosity alone is neither such map nor a proof of completeness.

### D. Port offsets are the promising surviving route

THM-5085, `01-canon/theorems/THM-5085-fixed-port-tails-in-every-base-and-the-ternary-bank-to-ten-million.md`, explicitly explains the distinction: graph offset moves have order-controlled increments, whereas tournament port moves depend on a response state (h,R,Delta) and deletion drift h_u. Its x=2 transport is already known from Sessions195/197/198/205 and THM-4888; only the general-x extension is new there. The finite bank up to 10^7 is not an all-order closure.

Our result suggests a concrete design rule: do not aim first to realize an arbitrary desired conflict graph. Instead, choose a carrier whose **extra-cycle sidecar closes in a small transfer state**. Port chains and bounded-span ribbons do exactly that. The new universal construction is a stress test for an alleged generic graph-to-tournament offset reduction: it must explain how it accounts for J_G and all remaining positive terms.

## Computational controls

Reproduction from repository root:

    python3 04-computation/experiments/session237_cycle_tilings/verify.py

Output: `verification.json` beside this report. Universe: all 1,099 labelled simple graphs of orders 1 through 5, no filters. For every graph the checker verifies compatible imposed arcs, all selected directed cycles, the exact induced corona, the vertex bound, existence of the extra local triangles, and the full coefficient vector of the corona polynomial identity against independent subset enumeration, with extra evaluation controls at x=1,2,3,4. Maximum constructed order is 37.

Direct Hamiltonian DP controls:
- G=K1: |T|=5, H=9; corona bound 5, strengthened bound 7.
- G=2 isolated vertices: |T|=10, H=81; bounds 25 and 29.
- G=K2: |T|=13, H=4617; bounds 21 and 25.

For K1 a second implementation enumerates all directed odd cycles, all their packings, and all permutations to check OCF and DP independently. All 16 completions of the unassigned pairs for this input are also checked; their H-values are exactly {9,11,13,15}. Hostile control: abstract K3 has I(K3,2)=7, and is deliberately **not** treated as the full conflict graph. The direct H controls use the lexicographic completion only; the all-completion assertion is proved above, not established by that finite computation.

## Provenance, corrections and remaining question

The uniform construction was derived in Session237 from the elementary edge-intersection representation of arbitrary graphs, plus private-cycle completion. We searched current canon, results, and historical reflections for induced universality, arbitrary/intersection graphs, realization, and private triangles. No identical uniform theorem was located. This is not an external novelty claim. Existing explicit induced C5/claw constructions are predecessors, not contradictions.

Historical reflection `07-reflections/the-even-graph-is-the-tournaments-cycle-half.md` uses degree-even cycle-space terminology and makes stronger “same information” assertions. Those assertions were **not imported** into this report's automorphism-even graph lane. The exact modern dictionary and THM-5040 guard the type distinction.

Question exposed: can a class of selected carriers supporting a complete offset numeration be completed so that the extra-cycle partition identity has a fixed-dimensional, nonnegative, controllable state? Port chains show a positive special case. Arbitrary induced realization does not solve this controllability problem. Chan–Pak's all-odd-except-7,21 conclusion remains open.

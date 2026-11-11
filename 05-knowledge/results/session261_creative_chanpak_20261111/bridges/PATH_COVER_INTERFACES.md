# Domination components, ordered blocks, and a literal two-apex escape

**PROVED:** the quantitative defect inequality for disjoint unions of complete bipartite components; a weighted induced-biclique bound; exponential clone amplification; and an exact two-apex alternating family (the latter with an explicit finite integer matrix certificate).
**FINITE-EXACT:** small controls and a continuation-state collision. **OPEN:** arbitrary bipartite quantitative defect and Chan–Pak coverage.

## Inheritance and live board

Accept Session260 `bridges/BIPARTITE_DEFECT.md`, `bridges/DEFECT_QUANTITATIVE.md`, `bridges/RESEARCH.md`, and `arithmetic/QUASIPERFECT_NORMS.md` as inherited dependencies. Their tools are the domination injection, positive newborn odd-cycle formula, actual Stirling substitution, and cyclic two-module identity. We do not re-audit them. The hostile is induced 2K2: one apex cannot preserve its domination count. The corrected near miss is treating a positive domination injection as if it reserved all old Hamilton paths. The neglected sidecar is ordered spanning path covers.

Board: domination components / odd-cycle packings / ordered component blocks / two-apex skeletons / operation-relative memory. Anchor: quantitative defect. Niche: multi-apex literal constructors. Wildcard: transfer an unrelated continuation-kernel theorem and a forest-complement symmetry. Sources and rejected maps are in `RANDOM_CONNECTIONS.md`.

Let U=T−z, A=N⁺(z), B=N⁻(z), and G have edge ab for a∈A,b∈B precisely when a→b. Write H for Hamilton paths, d for dominating sets, and Δ=H(T)−H(U). Empty tournament H is 1. Inherited OCF expresses Δ as the sum of weights 2^{number of cycles} of vertex-disjoint directed odd-cycle packings containing z.

## 1. Component-block theorem: an infinite non-chain class

Suppose the nonisolated components of G are complete bipartite graphs K_{p_i,q_i}, i=1,…,k. No restrictions are imposed on tournaments inside A or inside B. Isolated vertices are permitted. Set

    E_i = ((2^{p_i}−1)(2^{q_i}−1)+1)/2.

Then

    Δ ≥ 2 Σ_{∅≠J⊆[k]} |J|! ∏_{i∈J} E_i
      ≥ ∏_{i=1}^k (1+2E_i)−1 = d(G)−1.                 (1)

Consequently H(T)−d(G)≥H(U)−1 for this entire class.

### Proof

A dominating set of K_{p,q} either meets both shores, or equals one whole shore and misses the other. Thus d(K_{p,q})=(2^p−1)(2^q−1)+2=1+2E. Isolated vertices must be selected and do not change d; disconnected components multiply.

Among subsets of a p-vertex shore, the number of nonempty even subsets minus the number of odd subsets is −1. Therefore the number of pairs of nonempty subsets X⊆A_i,Y⊆B_i with |X|+|Y| even is exactly E_i.

Choose a nonempty component set J, one such pair in every component, and any linear ordering of J. In each component choose canonically a Hamilton path of U[X], followed by one of U[Y]. All A_i→B_i edges exist, so this is a directed even-vertex block beginning in A and ending in B. Between distinct components every B_i→A_j edge exists, since there are no graph edges between them. Thus the blocks concatenate in every one of the |J|! orders. Prepending z and closing the last B vertex back to z makes a directed odd cycle.

Different choices give different rooted directed cycles: their vertex sets recover the chosen subsets, and the run of component labels after z recovers the ordering. Each such cycle by itself is an OCF packing of weight 2; vertices not in the cycle are simply unused. This proves the first inequality. Finally j!≥2^{j−1} for j≥1, with equality at j=1,2 and strict inequality from j=3. Expanding the product proves the second inequality. No Hall lemma or unproved global reservation is used.

### Sharpness and stronger remainder

For k=1 and both shores transitive, U is transitive and equality holds. The inherited 2K2 lift with H(T)=13 and H(U)=5 realizes equality at k=2 as well. For k≥3, the proof alone gives the additional lower bound

    Δ−(d(G)−1) ≥ Σ_{|J|≥3}(2|J|!−2^{|J|})∏_{i∈J}E_i > 0.

In particular three independent edges force Δ≥30 whereas d−1=26. This is a genuine quantitative improvement, not merely strictness of the earlier qualitative defect. Internal cycles and unused vertices can add more positive OCF weight.

What fails for an arbitrary graph is precise: a cross-edge joining components reverses a needed B_i→A_j concatenation. Calling a collection of overlapping bicliques “components” destroys the literal block permutation. The theorem does not apply to arbitrary biclique covers.

## 2. Weighted induced-biclique lower bound

For any A₀⊆A,B₀⊆B with every A₀→B₀ edge, put W=U[A₀], V=U[B₀], R=U−(A₀∪B₀), and x(X)=Q₂(X)−H(X), where Q₂ is the ordered spanning two-path count allowing empty slots. Then

    Δ ≥ [x(W)x(V)+H(W)H(V)] H(R).                        (2)

Indeed the restriction to z,A₀,B₀ is C3[W,V,1]. The inherited exact cyclic identity gives its newborn increment x(W)x(V)+H(W)H(V). Combine every packing containing z in that restriction with any odd-cycle packing in R. Their union is a distinct packing of T containing z and their weights multiply. Cross-edges are irrelevant to this selected subset of packings.

Since Q₂(X)=Σ_{S⊆X}H(X[S])H(X\S), for nonempty X of size p,

    x(X) ≥ H(X)+2^p−2 ≥ 2^p−1.

This proves the star bound Δ≥2^q and the complete-biclique bound independently of the component proof. It also retains the extra H(R) sidecar. For a pure biclique, equality in Δ≥d−1 holds exactly when both shore tournaments are transitive: a nontransitive shore has H>1 and makes the displayed lower bound strict.

## 3. Transitive twin apexes amplify rather than repair

Replace z by a transitive m-vertex module with exactly z's external orientations. Denote its count by h(m), and let h(0)=H(U). For every m≥1,

    h(m) ≥ H(U)+2^{m−1}[H(T)−H(U)].                     (3)

Proof: take any old OCF packing containing z. In its unique z-cycle, replace z by an increasing run on any odd-cardinality subset of the m clones. The cycle remains odd, every directed edge is legal, and the number of cycles and hence the OCF weight are unchanged. There are 2^{m−1} odd subsets. Leave the other clones unused. Contracting the unique clone run recovers the old packing and selected subset, so these are disjoint objects. Packings avoiding every clone contribute H(U). All remaining packings have nonnegative weight.

Equality occurs for the cyclic triangle, where h(m)=2^m+1 and H(U)=1. Thus the exponential factor is sharp. This proof retains the full original newborn increment, rather than just one triangle witness. It gives a rigorous reason that indistinguishable extra apexes are not a small positive-defect correction.

## 4. Nevertheless two apexes escape the chain-only count image

Take U transitive on the word (AB)^k in its forward order. Add z₁→z₂ with B→z_i→A. Let F_k be the resulting actual tournament count. Then, for every integer k≥0,

    F_k = 9^{k+1}/50 + (120k+316)4^k/675 + (24k+19)/54.   (4)

The first nonempty-word values are 5,29,197,1485,12021,101949. The k=2 construction realizes 29, which the inherited all-order chain-graph/disjoint-union grammar misses. This does not make 29 a new member of the tournament spectrum; it is a concrete escape from that specific adapter.

### Literal skeleton compiler

For a single apex z, replacing z by m transitive clones gives

    h(m)=Σ_r a_r r!S(m,r).

For m=2 only r=1,2 survive: h(2)=a₁+2a₂. Here a₁ is the one-apex count. For U a transitive word, a₂ counts assignments of its vertices to three ordered bins L,M,R:

- vertices in each bin occur in their inherited transitive order;
- L is empty or ends in B;
- M is nonempty, starts in A and ends in B;
- R is empty or starts in A.

Concatenating L,z,M,z,R is exactly a two-run skeleton. The endpoints are the whole constraint; no unspecified joins remain.

Scan U from left to right. A state is (s₀,s₁,s₂)∈{0,1,2}³, recording each bin as empty, last A, or last B. On reading c∈{1,2}, put it in any one bin, except a B may not enter an empty M or R. Replace that bin state by c. All transitions have multiplicity one, with multiplicities added when state transitions coincide. The initial state is (0,0,0). Accept states with s₀≠1 and s₁=2, with no condition on s₂.

Let M=P_B P_A be the 27×27 integer matrix of a full AB block, e the initial column, and f the acceptance row. Then a₂(k)=f M^k e exactly. The entirely explicit matrix obeys

    f(M−I)²(M−4I)²(M−9I)=0.                              (5)

`path_cover_probe.py` checks every coordinate of this integer row identity, not just finitely many sequence values. Thus the sequence satisfies its order-five recurrence for all k. Its initial values are 0,1,10,83,684, and either solving that recurrence or multiplying the generating function gives

    Σ a₂(k)t^k = t(16t²−9t+1)/[(1−t)²(1−4t)²(1−9t)],
    a₂(k)=9^{k+1}/100 +(60k+8)4^k/675 −(12k+11)/108.

The one-apex word formula gives a₁(k)=(4^{k+1}+6k+5)/9. Substitution into a₁+2a₂ proves (4). The matrix certificate is a finite arithmetic proof of an all-order statement; the separate small literal Hamilton DPs are controls, not its justification.

## 5. A necessary continuation sidecar, with a literal collision

Two rooted tournaments, with rows encoding out-neighbours by binary masks, are

    T₁=(16,1,3,7,14), root 0;
    T₂=(20,1,2,7,14), root 0.

Both have (order,H(T),H(T−z),d(G_z))=(5,9,1,9), hence both have zero original domination defect. But their two-clone outputs are respectively 23 and 29; their three-clone outputs are 51 and 81. The first outside word is ABBB, the second ABAB after putting U in transitive order.

The unrelated operation-kernel criterion THM-840 therefore applies literally: this four-scalar observation is not a congruence for clone continuation. Already a₂ distinguishes the pair (7 versus 10). For a specified clone size two, (a₁,a₂) is sufficient; for arbitrary clone sizes retain the entire skeleton vector. No claim of globally minimal nonlinear memory is needed.

## 6. Controls and remaining maps

Run `python3 05-knowledge/results/session261_creative_chanpak_20261111/bridges/path_cover_probe.py`. It checks 5,405 labelled rooted tournaments of orders 1..5, actual clones of sizes 1..3, the complete-biclique subcase, the collision, the exact matrix identity, and literal alternating constructions through k=3. These are narrow new-claim controls, not a new spectrum census or an audit of inherited results.

Remaining concrete questions:

1. Replace the complete-component block order by a controlled set of admissible component permutations when cross-incidences are present. An injection needs the full incidence sidecar, not just component sizes.
2. Find a count-preserving domination source replacement whose multi-apex skeleton vector reaches a target; (4) supplies an infinite family, not a division algorithm.
3. Characterize two-apex equality in (3). Its surplus counts clone-separating odd cycles and multi-cycle packings, a positive interface potentially more useful than signed correction.
4. The arbitrary-graph inequality Δ≥d(G)−1 remains OPEN. The factorial block theorem is not a proof of Chan–Pak.

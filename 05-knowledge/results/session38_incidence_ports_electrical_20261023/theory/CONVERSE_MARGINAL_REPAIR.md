# Converse can repair an incidence signature without refining its alphabet

Status: elementary all-size sufficient theorem for explicitly reflection-constrained ordinary port counts. This is not a claim that paired-signature injectivity alone suffices.

Let C be the finite set of ordinary tournament classes at a fixed order, with converse involution c↦c*. Let S:C→Σ be a chosen signature. Form the bipartite multigraph Γ_S with left and right copies of Σ, and one edge for each class c, joining S(c) on the left to S(c*) on the right. Distinct classes remain distinct edges, including when their endpoints agree.

Fix a source ordinary class and a compatible source endpoint/color category. Its unknown outgoing counts x_c to target classes have known margins by S(c). Reflection sends this source category to its converse category, with transported endpoint scores. If exact reflection equality is required, the known signature margins of that reflected category also supply the margins of the original counts by S(c*).

Thus the same vector (x_c) is constrained by both endpoint sums on Γ_S. If some target classes are excluded by endpoint/color compatibility, restrict the edge set accordingly. This construction uses only margins already justified by reflection; it does not invent joint (S,S*) bucket counts.

## Forest criterion

If Γ_S is a forest as a multigraph, all these target counts are uniquely determined by their two families of margins. Consequently, when an authenticated actual table is feasible, every ordinary cell and hence every merged colored edge multiplicity is forced by the original signature margins together with reflection.

Proof. The difference of two solutions has zero sum at each vertex. A leaf forces its unique incident difference to vanish. Delete that edge and continue. Every edge vanishes because a finite forest is exhausted by this process. Equivalently, orient every edge from left to right, negate the right equations, and use the independence of the oriented incidence columns of a forest. Parallel edges form a cycle of length two and must not be suppressed. The proof holds over the integers and does not require positivity. QED.

The forest condition is sufficient, not necessary: capacities, nonnegativity, other margins or absent targets can force uniqueness even in a cyclic Γ_S. For a fixed feasible vector the residual-cycle criterion gives the sharper answer for this marginal subsystem.

## Singleton-converse criterion

A particularly simple sufficient condition is that each class c belongs either to a singleton S-fibre or has c* belonging to a singleton S-fibre. Then each unknown x_c is read directly from one of its two margins. This also implies that every edge of Γ_S has a degree-one endpoint, so Γ_S is a disjoint union of stars and isolated vertices.

The order-eight E signature with the neighbor triangle count reduced modulo three and exact edge incidence satisfies this condition: its four collision pairs contain eight classes, and the converse of every one has a globally singleton signature. Therefore its remaining ordinary transport ambiguities disappear when native reflection is imposed. This is a stronger and different result than merely failing to find a mirror-paired four-cycle.

## Why paired-signature injectivity is not enough

Four distinct classes may have signature pairs (a,x),(a,y),(b,x),(b,y). Every pair is distinct, yet the signed change (+1,−1,−1,+1) preserves all separate margins. If it fits the nonnegativity bounds, the class counts change. A true joint paired-signature label would distinguish the four classes, but separate signature margins do not provide those joint counts.

Therefore distinguish:

1. ordinary class separation by S;
2. class separation when the joint label (S(c),S(c*)) is explicitly supplied;
3. reconstruction from S-margins and their reflection-implied counterpart margins.

The third is the actual symmetry-repair statement proved above. It can succeed while the first fails, and the second alone does not imply it.

This theorem depends on the specified reflection rule being authenticated. It does not infer reflection compatibility from an arbitrary colored multigraph, nor assert that every feasible relaxed table is a tournament realization.

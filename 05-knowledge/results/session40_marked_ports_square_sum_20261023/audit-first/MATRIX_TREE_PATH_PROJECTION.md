# Hamiltonian edge support as a degree projection of the matrix-tree polynomial

**Status: PROVED elementary carrier, not a new claim about the classical matrix-tree theorem.** It gives an exact all-size formula. It does not claim polynomial-time Hamiltonicity or a tournament-metagraph realization theorem.

Let G=(V,E) be a finite undirected simple graph, |V|=n≥2. Give each vertex an indeterminate x_v and each edge an indeterminate y_e. Let L be the weighted Laplacian with weight w_uv=y_uv x_u x_v and let τ_G(x,y) be any principal cofactor of L. By the weighted matrix-tree theorem,

τ_G(x,y)=Σ_(T spanning tree) (∏_(e∈T)y_e)(∏_(v∈V)x_v^deg_T(v)).

For distinct s,t define

P_st(y)=[x_s x_t ∏_(v≠s,t)x_v²] τ_G(x,y).

Then P_st(y) is exactly the generating polynomial of unoriented Hamiltonian paths whose endpoint set is {s,t}, with each path weighted by the product of its edge variables. Consequently

P_G(y)=Σ_({s,t}⊂V)P_st(y)

is the Hamiltonian-path polynomial. Every path is counted once, not twice. For n=1 define P_G=1 separately.

**Proof.** A tree contributing to the coefficient has degree1 at s,t and degree2 everywhere else. A connected acyclic graph of maximum degree2 is a path, and its vertices are all V. Conversely every spanning path with those endpoints has exactly that degree monomial. Distinct trees have distinct edge monomials, so coefficients count the required paths without cancellation. ∎

For any set F⊆E, the mixed edge derivative evaluated at y=1 is

(∏_(e∈F)∂/∂y_e) P_G(y)|_(y=1)=#{Hamiltonian paths containing F}.

The polynomial is multilinear, so an edge variable can occur only once. The same identity holds endpoint-by-endpoint for P_st. Thus joint edge incidence, not merely separate marginals, is accessible in the same carrier.

Let H=P_G(1). When H>0, an edge is unused exactly when its first derivative is0, forced exactly when it equalsH, and optional otherwise. When H=0 these labels must be kept separate from vacuous universal claims. A pair of edges can each have positive participation and yet have mixed derivative0: two feasible marginals do not establish shared realization.

The incidence moment problem has a sharp finite terminal level: all mixed derivatives determine P_G by Taylor expansion about y=1. This is lossless but potentially exponential, not a small-signature theorem.

## Connectivity is a real sidecar

Degree constraints alone are insufficient. On five vertices, the disjoint union of triangle1-2-3-1 and edge4-5 has n−1=4 edges, degrees(2,2,2,1,1), and two endpoints, but it is not a spanning path. The spanning-tree cofactor eliminates it precisely because it is disconnected and cyclic. This is the exact information lost by replacing the determinant/tree condition with degree margins.

## Source and target of the analogy

Source: the edge-incidence polynomial of actual spanning trees. Target: its Hamiltonian degree sector and mixed edge moments. Preserved: exact common occurrence, endpoints, and connectedness. Lost by marginal projection: compatibility of choices. Required sidecar: mixed incidence or a valid connected-subgraph certificate.

Tournament port wiring has a parallel question about whether separately compatible endpoint observations can be glued to one actual matching. The analogy is structural, not an equality of polynomials or a proof that the same determinant enumerates tournament-metagraph wiring.

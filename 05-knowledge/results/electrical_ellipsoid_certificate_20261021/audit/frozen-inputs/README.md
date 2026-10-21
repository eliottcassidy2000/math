# Ellipsoid contacts as an exact electrical-network certificate

**Status: proof and exact controls, awaiting independent review.** This is a standard maximum-determinant duality specialized to the projected-cube carrier, with an elementary certificate proof. It is not a new general duality theorem or a Lonely Runner proof.

## 1. Carrier, coordinate map, and information retained

Let v_1,...,v_k>0, k>=3, and r_i=1/v_i. With P the orthogonal projection onto v-perp, the body D=P[-1,1]^k has the exact pair-facet description proved in the companion dominant-speed packet. Write x_i=v_i y_i and quotient by adding a common constant to all y_i. Choose y_k=0 as a coordinate gauge. The body becomes

    K={y in R^(k-1): |y_i-y_j|<=r_i+r_j for every i<j},

where the grounded coordinate y_k means zero. Explicitly, the grounded map is y_i=x_i/v_i−x_k/v_k for i<k; its inverse is x=P diag(v)(y_1,...,y_(k−1),0). This is a linear image of the projected cube. It preserves ellipsoid containment and volume comparisons; absolute determinant values depend on this coordinate choice. The arithmetic lattice must be transported separately. No lattice point, common runner time, or loneliness witness is supplied by this change of coordinates.

For each unordered pair e={i,j}, let b_e be its reduced incidence vector, with any orientation, and put c_e=r_i+r_j. A centred ellipsoid with covariance Sigma>0 lies in K iff

    b_e^T Sigma b_e <= c_e²  for every e.                    (1)

Signs of the incidence vectors cancel. There is no intrinsic tournament orientation in this graph; choosing one is only a linear-algebra gauge.

## 2. A verifiable network certificate

Choose nonnegative conductances lambda_e whose positive-edge graph is connected. Its grounded Laplacian is

    L=sum_e lambda_e b_e b_e^T >0,      d=k-1.

Define the effective resistance of any pair, including a nonedge, by

    R_e=b_e^T L^(-1) b_e.

This has its ordinary electrical meaning: injecting unit current b_e gives potential L^(-1)b_e; the voltage difference and dissipated energy both equal R_e.

**Certificate theorem.** If

    R_e<=c_e² for all pairs,
    lambda_e>0 implies R_e=c_e²,                            (2)

then Sigma=L^(-1) is the unique covariance of the maximum-volume ellipsoid in K. The ellipsoid has centre zero, even if translated competitors are allowed.

Proof. The candidate is feasible by (1). Also

    sum_e lambda_e c_e²=sum_e lambda_e R_e=tr(L L^(-1))=d.

For a centred competitor T>0, feasibility gives tr(LT)<=d. The eigenvalues of L^(1/2) T L^(1/2) are positive, so arithmetic-geometric mean gives

    det T <= (tr(LT)/d)^d / det L <= 1/det L=det Sigma.     (3)

Equality requires L^(1/2) T L^(1/2)=I, hence T=Sigma.

For a translated competitor h+T^(1/2)B, the support inequalities are

    |b_e^T h|+sqrt(b_e^T T b_e)<=c_e.

Its centred part is feasible, so the same bound holds. At equality T=Sigma. Every positive-conductance edge is tight by (2), forcing b_e^T h=0. Connectivity forces h=0. Degenerate ellipsoids have zero volume and are not maximizers. QED.

For any connected conductance choice, even one failing (2), the same calculation gives the valid upper bound

    det T <= [sum_e lambda_e c_e²/d]^d / det L             (4)

for every feasible covariance. Failure of the matching effective-resistance condition means the bound may be slack, not that the original body is empty.

## 3. Existence of such a contact network

For completeness the network certificate is not only sufficient. The feasible covariance set is bounded: grounded pair constraints bound diagonal entries, and positive semidefiniteness bounds the others. A sufficiently small multiple of the identity is strictly feasible. Maximizing log det therefore has a positive-definite optimum Sigma, unique by strict concavity.

At that optimum, the gradient Sigma^(-1) is a nonnegative combination of the active constraint matrices b_e b_e^T. Indeed, otherwise finite-dimensional conic separation supplies a symmetric direction U with tr(Sigma^(-1)U)>0 and b_e^T U b_e<=0 on every active constraint. A sufficiently small step preserves positive definiteness and all slack constraints and increases log det, a contradiction. Thus

    Sigma^(-1)=sum_e lambda_e b_e b_e^T,

with lambda_e>=0 supported on active pairs. Positive definiteness makes the positive-edge graph connected. These are exactly (2). This is the usual first-order maximum-determinant argument, not an additional hypothesis about runners.

Equivalently the grounded-Laplacian dual objective is

    sum_e lambda_e c_e² - log det L - d.

Its domain is lambda_e>=0 with L>0, equivalently a connected positive-conductance support. Its minimum equals max log det Sigma. The explicit certificate theorem proves equality whenever the resistance and contact conditions are exhibited, without requiring a numerical optimizer.

## 4. A graph-count identity, with no count-to-incidence leap

Let B be the reduced incidence matrix. Since L=B diag(lambda) B^T, Cauchy-Binet gives

    det L=sum_(spanning trees T) product_(e in T) lambda_e. (5)

A square incidence minor has determinant ±1 exactly for a spanning tree, and zero otherwise. Hence this is the weighted matrix-tree identity, with a direct proof in the present coordinates.

Together (3) and (5) link the optimal ellipsoid volume to a weighted tree count. The edge weights and the resistance inequalities are indispensable. A bare count of trees, tournament paths, or abstract cycles is not an ellipsoid or a lattice-hit certificate.

## 5. Why a cycle is necessary outside the dominant-speed cone

If a certifying contact network has exactly k−1 positive edges, connectivity makes it a tree. An edge's effective resistance in that tree is 1/lambda_e; its tightness forces lambda_e=1/c_e². Effective resistance between any two vertices is therefore

    R_ij=sum_(edges e along the tree path i to j) c_e².

Condition (2) is exactly the squared path criterion of the companion facet-tree theorem. On a length-two path a-b-c it forces

    r_b²+r_b(r_a+r_c)<=r_a r_c,

so r_b<min(r_a,r_c). A length-three path would force both r_b<r_c and r_c<r_b. Thus the only possible tree is a star, and its centre must be the dominant-speed vertex satisfying

    M²>=M(A+B)+AB,

with A,B the largest other speeds.

Consequently outside that cone EVERY certifying network must contain a cycle and have at least k positive edges. This is an exact contact-certificate obstruction, not a temporal cycle or a directed tournament cycle. The triangle-repair packet supplies the next simplest positive example: a star with one extra leaf edge, hence one triangle. Its covariance and positive contact weights satisfy (2).

## 6. Literal controls and scope

The checker uses exact rational arithmetic, independent Gaussian elimination, direct pair resistances and explicit spanning-tree expansion. It tests a dominant star at (1,2,4), and the triangle-repaired network at (1,2,3,5). Deleting the latter's extra edge leaves a connected tree but violates a resistance constraint, so connectivity alone is not the certificate.

It also tests the determinant upper bound (4) against feasible covariance samples. Normal and optimized Python runs retain every check through explicit exceptions. The code imports no repository modules.

This provides a precise information-flow interpretation: a unit source-to-sink current has a measurable resistance, active contacts saturate exact geometric budgets, and the Laplacian recovers the ellipsoid covariance. It preserves source labels and the full edge weights. It does not assert that the Gaussian covariance describes runners probabilistically, that these currents are their physical trajectories, or that a determinant improvement supplies the missing arithmetic orbit hit.

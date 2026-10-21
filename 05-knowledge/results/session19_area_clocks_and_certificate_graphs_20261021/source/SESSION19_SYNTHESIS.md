# Session 19: area clocks, physical ports, and the graph of a convex certificate

Status: synthesis in progress. The separately reviewed proofs and their scope statements are authoritative. No general Lonely Runner proof or new Hadwiger–Nelson lower bound is claimed.

## 1. What the creative analogy actually preserves

The user's proposal was to view runner counts through regular convex polygons, use Kepler's second law creatively, and retain tournament structure while seeking progress on Lonely Runner and Hadwiger–Nelson. The successful common language is not merely “a cycle.” It is a labelled configuration with a specified clock, predicate, and permitted coordinate change.

For an orientation-preserving circle homeomorphism H and 0<delta<=1/2, exact preservation of the predicate distance(x,y)>=delta for every pair is equivalent to H(x+delta)=H(x)+delta. This converts a metric question into a precise rotational-symmetry requirement. It also identifies the limit: preserving thresholds with unbounded denominators forces H to be a rotation. One fixed nonlinear clock cannot preserve all runner-count thresholds.

On a positive continuous radial curve, normalized swept area is a faithful phase coordinate. Its inverse map to polar angle preserves threshold 1/N exactly when the curve has N-fold rotational symmetry. On a tangential convex polygon, all edge support distances from the centre are equal, so swept area is proportional to perimeter arclength. A regular polygon therefore gives an exact equal-area, constant-perimeter-speed model.

Taking L=lcm(N,6) sides and then projecting radially onto the unit circle preserves three predicates simultaneously:

1. The sharp loneliness threshold 1/N.
2. The semicircle orientation used for a geometric tournament, including collision and antipodal ties.
3. Euclidean unit-distance equality on the projected unit circle, corresponding to phase difference ±1/6.

This is an exact representation change, not an orbit-hitting theorem. The physical polygon's chord lengths are not silently substituted for phase distance. A regular polygon is also an exact periodic circular-billiard orbit: straight flights and radial reflection impulses conserve angular momentum and swept-area rate. It is not an inverse-square Kepler ellipse.

## 2. Why the tempting shortcuts fail

Several explicit hostiles separate the valid map from stronger claims.

- The runners need not ever form a regular polygon. For N>=3, with velocities {0,1,...,N−2,N}, a regular N-gon would force t=j/N from speed1, then speedN collides with the stationary runner. The regular reference track is useful; requiring a regular instantaneous configuration is not equivalent to loneliness.
- An odd regular polygon can preserve its 1/N threshold while reversing a semicircle edge. An even multiple repairs this particular loss.
- A focused eccentric Kepler ellipse has no nontrivial rotational symmetry. Its equal-area-to-angle map therefore fails universal 1/N threshold preservation. At eccentric anomalies pi/2 and 3pi/2, the mean-anomaly forward difference exceeds pi while the true-anomaly difference is below pi: even a tournament edge reverses.
- Equal aggregate area sweep is too coarse. The relative-speed triples (1,2,6) and (1,3,5) have the same sum but different maximum loneliness, 2/7 and 1/2.
- A labelled geometric three-cycle can occur both when the observer is lonely and when it is not. Orientation without threshold slack does not carry the target.

The original common-start assumption is also essential. A current shifted-LRC counterexample with speeds 1,...,5 and specified nonzero starts has maximum 15/94<1/6. This fixed vector and its four maximizing times were independently verified exactly. The stronger arbitrary-start assertion cannot be imported through a geometric change of coordinates.

## 3. Physical ports give a positive colouring theorem

Every finite unit-distance graph on one circle has degree at most two and is at most three-colourable. On the radius-one circle, its components are paths or six-cycles, hence they are 2-choosable: any lists containing at least two colours can be coloured properly.

This yields a source-preserving extension theorem. Fix the actual colouring of the entire exterior X of a centre and its unit-circle neighbourhood N. For each rim vertex z, let F_z be the colours used by its actual exterior neighbours. A prescribed centre colour c extends whenever

    |F_z union {c}| <= q−2 for every z in N.

The fixed exterior colours and all geometric contacts are retained. The rim is free to recolour; prescribed rim ports must instead retain their singleton lists. Confusing those two boundary conditions would produce a false forcing claim.

At q=5, at most two exterior contacts per rim vertex guarantee extension for every centre colour. Consequently every vertex-critical six-chromatic unit-distance graph has, at every root, a neighbour with at least three neighbours outside the root's closed neighbourhood. A centre-based six-chromatic construction also needs an off-circle core of chromatic number at least four. These are necessary structural restrictions, not a construction of such a graph.

The exact clock-to-contact map further explains a limitation of the polygon picture. For common-start integer-speed runners at a reduced rational time a/q on the unit circle, unit-distance edges require 6 to divide q. The tight fourteen-runner AP configuration at time1/14 has an empty radius-one unit graph. Scaling its adjacent chord to length one produces only a cycle, not a high-chromatic graph. Conversely, preserving all unit contacts at radius two forces an angular homeomorphism to be a rotation: that radius has an irrational contact angle. Finite resonant clock symmetries cannot be promoted to plane-wide unit-distance symmetries.

A concurrent repository correction supplies another exact graph-versus-clock distinction. For the integer distance graph with forbidden differences {1,3,4,7}, vertices0,3,4,7 form a clique of size four, and an explicit period-eight four-colouring proves chromatic and fractional chromatic number four. But its best common rotation has loneliness1/5, so rotation colourings require five. Thus forgetting additive rotation equivariance changes the optimum. This is an integer distance graph, not a claimed Euclidean unit-distance embedding; the current incoming proof and its exact interval cover are retained with their original attribution.

## 4. Timed repair graphs: a real orientation, with a real obstruction

For k>=2 distinct positive integer relative speeds with common starts (N=k+1 total runners), the adaptive clock q=2k+1 has a useful exact dichotomy under the explicit all-speeds-unit-modulo-q hypothesis. Repeated signed residues leave an immediate safe grid class. Otherwise the speeds form a complete signed residue system, which forces q to be prime.

Each inverse clock then has one unsafe owner. Comparing that owner's exit time with another runner's first danger-entry time gives an intrinsic forward-corridor oriented repair graph. It is not generally complete, so it is not automatically a tournament. Reciprocal edges in this forward graph and its forward endpoint subgraph are impossible by a modular-inverse product inequality, and a sink certifies a lonely endpoint.

For a prescribed residue cycle, an explicit product deficiency is necessary and sufficient for realizing its blockers at arbitrarily large integer speed heights. The construction can keep blockers inside the first following forbidden tooth at the actual repair endpoint. Thus the cycle is not just an abstract graph artefact.

But its arrows still do not compose chronologically: arrival at a target runner's danger boundary is not arrival at that runner's own inverse clock. The timing sidecar cannot be discarded. A concrete five-speed row defeats all ten two-sided repair endpoints while t=1/24 is lonely. More strongly, a quadratic-residue family defeats every one of those proposed repairs while its actual loneliness tends to 1/2. The failed method can miss very easy rows, not just barely lonely ones.

A full signed event tape remains an exact certificate: at tied entry/exit times, retain all owners and evaluate the point itself, not only the open cells on either side. Tight arithmetic progressions demonstrate why an arbitrary ordering of tied events loses legitimate boundary-only witnesses. None of these certificates proves that all rows admit a witness without checking or controlling the remaining events.

## 5. A stronger convex-geometric result: the star certificate

For k>=3 and positive speed coordinates v=(a_1,...,a_(k−1),M), the geometric carrier is the (k−1)-dimensional projected cube D=P[-1,1]^k in the hyperplane perpendicular to v. The geometric theorem allows positive real coordinates; arithmetic lattice statements require their additional hypotheses. This is a zonotope, arithmetic when the speed direction is rational, and not generally a planar regular polygon. Its exact pair-facet inequalities are

    |v_j x_i−v_i x_j| <= v_i+v_j.

Let M be the largest speed, with A and B the two largest remaining speeds. If

    M² >= M(A+B)+AB,

an explicit projected ellipsoid with ambient diagonal entries (1+a_i/M)² and a zero entry at M is the unique maximum-volume ellipsoid in D, even among translated, non-diagonal ellipsoids. The ambient zero does not cause degeneracy after projection.

The proof identifies the pair facets incident to M with a star graph. Their coordinate map sends the enclosing parallelotope to a cube and the candidate ellipsoid to its unit ball. Hadamard's inequality bounds every ellipsoid in that cube; the displayed speed condition is exactly what makes the ball satisfy all remaining pair facets.

This also resolves the graph-shape question rather than relying on a drawing. Exactly k−1 independent pair constraints form a spanning tree. In reciprocal-speed coordinates, feasibility of its parallelotope's maximal ellipsoid requires a squared path-length inequality. Every length-two path forces its middle reciprocal weight to be smaller than both endpoint weights. A length-three path forces two contradictory inequalities. Therefore only a star can work.

For the familiar four-vertex, three-edge patterns: the star can certify the shape in its exact speed cone; the path cannot do so for any positive speeds; a triangle plus an isolated vertex is rank-deficient. This is a specific structural theorem about convex certificates, not a claim that tournaments or Collatz are isomorphic to those three graphs.

For (1,2,4), the squared-volume factor improves from 221/84 for the cited diagonal construction to75/28 for the optimal ellipsoid. This is a local volume comparison. The two ellipsoids do not contain one another, so the gain does not imply pointwise norm domination or automatically improve a global LRC speed-product bound.

### 5.1 The next contact pattern: a triangle repairs a star

When exactly one remaining leaf-pair constraint cuts the star's unit ball, an explicit two-coordinate covariance correction makes that pair tight without changing the other pair tests. A positive contact decomposition proves that the repaired ellipsoid is again the unique global maximum, including translated competitors. This completes the shape classification for every positive three-relative-speed vector: four total runners and a two-dimensional body. For k>=4 the theorem remains restricted to the one-violating-pair regime. The explicit all-dimension nonvacuity family permits repeated speed coordinates; the separate four-relative-speed example below is distinct.

For (1,2,3), the exact squared-volume factor is 22/7, compared with 64/21 for the cited construction. For (1,2,3,5), it is1488/325. In the higher-dimensional single-pair regime the optimal shape cannot be represented by ANY projected nonnegative diagonal lift: such a lift would create correlations in every leaf pair, whereas the unique optimum has only the one repaired correlation. Compactness makes the volume gap strict. This is a specific limitation of a geometric certificate class, not a disproof of the LRC method using other shapes.

### 5.2 Information flow becomes an actual electrical network

There is a precise network interpretation of the same convex problem. In grounded reciprocal-speed coordinates, pairwise width limits are r_i+r_j. Assign nonnegative conductances to labelled pairs with connected positive-conductance support, and form their positive-definite grounded Laplacian L. If every effective resistance is at most (r_i+r_j)², and every positive-conductance edge attains its bound, then L^(-1) is exactly the maximal ellipsoid covariance. The proof is a trace bound followed by arithmetic-geometric mean, including the translated-centre condition.

The converse is the usual maximum-determinant optimality condition, so this is an exact certificate framework rather than just an analogy. The matrix-tree identity expresses det L as a weighted spanning-tree sum. Both weights and resistance constraints are essential; a bare tree count is not enough.

A tree contact network must be a star by the earlier path obstruction. Therefore outside the dominant-speed cone every certifying network needs a cycle. The triangle repair supplies the simplest extra cycle. These are undirected geometric contact graphs: arbitrary incidence orientations are only a gauge, and their edges are not the timed arrows of the runner-repair graph. Electrical current is a faithful interpretation of the Laplacian certificate, not a proposed physical trajectory of the runners.

A useful cross-check joins the two mechanisms without identifying their graphs. In the quadratic-residue hostile family, all speed ratios tend to one. The right-hand side of the star condition divided by M² therefore tends to three, so the star certificate fails; the local clock repairs fail as well. Yet actual loneliness tends to1/2. Thus two certificate obstructions can coexist on increasingly easy rows. Contact-cycle complexity and repair-cycle complexity are not themselves measures of failure of loneliness.

## 6. Current literature and the exact remaining problem

The repository's older external headline requires qualification. Allikvere's September24 preprint, arXiv:2609.02604v2, claims computer-assisted proofs for fourteen and fifteen total runners. This session checked its stated claims and used a specified elementary geometric construction for comparison; it did not independently verify that entire proof or its certificate archive. The general conjecture is not settled by that paper or by this session. The mathematical statements here do not rely on promoting those claims to proved repository dependencies.

A separate current primary connection is Alon–Bucić–Davies, arXiv:2609.36220v1: for 1<=k<n, its matrix-avoidance theorem gives supremum clearance at least k/(2n) for a shifted k-parameter system under a full-spark condition, with attainment for rational rows, and supports bounds for typical norms. It is not a Euclidean-plane colouring bound. Nor can a multi-parameter clearance conclusion be restricted without proof to the original one-dimensional common-time runner orbit; the extra parameters are real freedom, not notation.

The fixed-centre zonotope hit for original LRC must not be replaced by all-coset covering. The latter corresponds to the shifted problem, where explicit counterexamples now exist. Likewise, an inner ellipsoid's exact volume is not a lattice-hit certificate. The missing arithmetic/time obligation remains visible after every faithful map.

Next useful targets are consequently specific:

1. Control the actual common-time orbit hitting the safe region, rather than its area sum or unweighted orientation.
2. Retain arrival clocks when composing local repair obligations; the cycle-product theorem alone does not supply that transport.
3. Use the exact ellipsoid and its contact graph in a valid arithmetic bound, without replacing determinant gain by norm containment.
4. Build an actual exterior context that forces five-colour saturation; low-contact free-rim gadgets are excluded by the extension theorem.

The successful theme is conservation of the right information. A useful analogy identifies both what it preserves and the coordinate whose omission would invalidate the conclusion.

## 7. Concurrent Kepler argument: closed endpoints and actual existence

Incoming THM-4737 was inspected rather than inherited automatically. Its right-continuous one-sided area-rate test misses a tight closed witness: on the triangle, relative speeds (1,2) at time1/3 are lonely, but the runner at phase2/3 is assigned to the outgoing observer-incident edge, whose rate is zero. The exact repaired test, for positive speeds and polygon count at least three, requires the maximum of the left and right area-rate limits to be positive for every runner. It retains both far incident vertices and still rejects a collision at the observer.

A second distinction concerns the least denominator N_min of a lonely rational time. Rational attainment of the maximum loneliness does not itself prove that maximum reaches the conjectured threshold. The sound invariant allows N_min=infinity when no witness exists; its finiteness is equivalent to LRC for that integer row. The maximizing-time argument does yield the conditional bound N_min<=2 max(v_i), once the threshold is reached.

For real speeds, a closed-box hit by the orbit closure describes the supremum. Density alone does not justify an exact boundary hit by the orbit itself. For distinct positive common-start real speeds, the rank-at-most-one consequence survives with its proper argument: interior hits when the relation interval is strictly longer than one, and periodic attainment in the remaining two-speed commensurable case. Haar density for incommensurable speeds is an asymptotic-time statement, not generally a one-unit-time-window measure.

These repairs preserve the useful observer-centred Kepler picture. They remove circular existence reasoning and retain precisely the endpoint and clock information required by the original question.

## 8. Publication reconciliation

While our Kepler correction was awaiting publication, concurrent commit cfcb930344470c60cf13ffc8694e990eaa074209 repaired the same closed-endpoint and existence issues. That correction was preserved. Our independent audit is being prepared alongside it with the original inspected source and a current concordance, rather than overwriting the owner's repaired theorem or duplicating its correction entry. The additional conditional denominator bound and explicit density qualification remain in the audit. A failed guarded-ref attempt did not become an extra main-branch commit.

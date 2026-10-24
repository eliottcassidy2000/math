# Balanced charts and the unavoidable distance-squared cost

**PROVED:** all-r Hamming-chart balance and the sharp uniform full-refresh path-congestion lower bound. **FINITE-EXACT:** five cross-pointed-chart connections in the complete 409-path seed atlas. **OPEN:** full-atlas all-size native connectivity and source-specific ROOT estimates.

## A complete balanced chart family inside the native cube

On the actual fixed-source cube X_r from CUBES.md, partition by Hamming weight: X_k={x:|x|=k}. These are weight charts, not the pointed partner-isomorphism charts of the full atlas. Conditional uniform measures admit the native adjacent-chart coupling

    K_(k,k+1)(x,y)=1/(r-k) if y adds one previously zero bit.

The reverse kernel has weight 1/(k+1) on deleting a one. Both marginals are exact because

    binom(r,k)(r-k)=binom(r,k+1)(k+1).

The quotient law is p_k=binom(r,k)/2^r; upward transition (r-k)/(2r), downward k/(2r), holding 1/2. It lifts to exactly the authenticated coordinate-resampling chain, with gap 1/r. Thus all levels are covered and coupled at every r, with the original pins, fixed source, distinct resources and explicit clocks. This is complete coverage of X_r, whose fraction of all source paths is (2/409)^r; it is not full-atlas coverage.

Actual nonuniform source weights cannot be replaced by those convenient conditional uniforms. There is a sharp small obstruction. At r=3 give the singleton {1} weight kappa, the pair {2,3} weight kappa, and every other state of those two levels weight one. The conditional mass of {1} is kappa/(kappa+2); its two adjacent targets together have mass 2/(kappa+2). A direct-edge balanced coupling is impossible for kappa>2. At kappa=2 equality holds and a coupling exists; at kappa<2 the symmetric remaining-edge allocations give one as well.

The obstruction is to this edge-supported conditional interface, not to connectivity or to the whole-cube Metropolis chain. A path coupling can reroute through other legal levels. In particular, independently couple the two prescribed conditional laws, assemble a reversible quotient kernel with the actual chart masses, and apply the existing `guarded_coverage_20261023/LIFTING.md` stationary-return compiler on X_r. It gives at most 2r events per return and congestion at most 4 kappa r^2 when the full weight ratio is at most kappa. This is an application of the prior stronger compiler, not a new proof of universal coverage.

## Optimal congestion theorem

Write N=2^r and P(x,x^i)=1/(2r), with uniform stationary law 1/N. For arbitrary possibly randomized routes gamma from independent uniform endpoints x,y, let ell(gamma) count native edges, including repetitions. Define length-weighted directed-edge load

    L_e = sum_(x,y) N^(-2) E[ell(gamma) times number of uses of e],
    Lambda = max_e L_e / [1/(2rN)].

Every route satisfies ell>=d_H(x,y). Summing loads over all directed edges gives E[ell^2]. Total directed nonholding conductance is 1/2. Hence

    Lambda/2 >= E[ell^2] >= E[d_H(x,y)^2] = r(r+1)/4,
    Lambda >= r(r+1)/2.

The ordered-coordinate shortest routing from CUBES.md attains equality: every directed edge is used by N/2 endpoint pairs; conditioned on use of an edge, expected route length is (r+1)/2. Its load divided by capacity is r(r+1)/2. Thus the prior upper bound is optimal among all such routings, including fractional routings, not merely an artifact of choosing a coordinate order.

Equality requires shortest routes and equal load/capacity ratios on every positive-capacity edge. This equality description concerns the displayed sum bound; it does not assert a unique routing.

More generally, for any endpoint demand pi(x)R(x,y) and any native conductances c_e, the same proof gives

    Lambda >= E_(pi R)[d_native(x,y)^2] / sum_e c_e.

For independent pi endpoints with full weight ratio at most kappa, pi(x)>=1/(kappa N). Metropolis conductance has total nonholding mass at most 1/2, so

    Lambda >= r(r+1)/(2 kappa^2).

This lower bound need not be sharp for weighted laws. It does show why constant congestion independent of depth cannot hold for broad refresh demands at bounded weight ratio. Uniformity over controllers at each r is compatible with a polynomial r-dependent estimate; depth-independent congestion for this demand is not.

## Five actual seed cross-chart bridges

`charts.py` starts with all 409 literal Hamiltonian paths of source code 29558785. Its 362 pointed charts retain actual first/last source owners and target isomorphism fixing the two corresponding marked owners. Forgetting only the target pointings gives 357 coarse fibres: 352 singleton-chart fibres and five two-chart fibres. No two-clock transposition crosses any of these five pairs.

A window move replaces a contiguous clock interval, retains the entire path outside that interval, and must remain a Hamiltonian path with the same actual source endpoints and unpointed target class. Exact bipartite transport and an independent weighted Hall test give:

| Pointed chart IDs | Sizes | Smallest balanced window |
|---|---:|---:|
|203,205|1,1|4|
|208,210|1,1|5|
|211,215|1,1|5|
|218,221|1,1|5|
|219,224|1,2|5|

At the stated window all cross edges are available, so the unequal-size case uses mass 1/2 on each target. At all smaller tested windows there are no cross edges. These are literal bounded events in this eight-owner source; the complete path records and conditional transport are in charts.normal.json.

Crossing a pointed-chart boundary loses an endpoint-fixed target transporter. Therefore these five bridges cannot automatically be tensored by the ordinal-block lemma, whose proof needs that marking to control boundary arcs. The missing sidecar is a compatible target-boundary correspondence in the larger host. This is a precise place where adding an unpointed isomorphism is too weak.

## Independent obligations

Balanced conditional endpoint mass, physical capacity along a route, and full guard preservation are separate conditions. A large spectral gap of a quotient does not create a balanced interface. A balanced interface does not bound route congestion. A path of low congestion does not authenticate missing native endpoints. The formulas above provide either exact certificates or exact obstructions at each of these steps.

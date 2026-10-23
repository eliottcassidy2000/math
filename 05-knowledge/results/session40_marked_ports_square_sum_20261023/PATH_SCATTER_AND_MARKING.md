# Endpoint conditioning, Hamiltonian Hessians, and conserved edge costs

Status: PROVED elementary finite-ensemble identities; FINITE-EXACT for the square-sum ranks below. This is not an algorithmic lower bound or a claim of new priority for covariance identities.

## 1. The exact carrier

For a finite simple graph G with at least one Hamiltonian path, count unoriented paths once. Write x(P) in {0,1}^E for edge membership and

Z(theta)=sum_P exp(sum_e theta_e x_e(P)).

Then the gradient of log Z is the vector of edge-use probabilities and its Hessian is Cov(x). Indeed differentiation weights each path by its edge indicators. For every real edge-cost vector a,

a^T Hess(log Z) a = Var(sum_e a_e x_e(P)) >= 0.

All probabilities are positive for finite real theta. Consequently the kernel consists EXACTLY of the edge costs whose sum is constant across all paths in the ensemble. Its rank is the affine dimension of {x(P)} and does not depend on theta. The same statement holds after conditioning on one actual unordered endpoint pair. Empty endpoint sectors are excluded, not assigned a covariance.

If B is the unsigned vertex-edge incidence matrix, then paths with endpoints {s,t} satisfy B x(P)=2*1-e_s-e_t. Every cost B^T alpha is therefore conserved in that endpoint sector, with value 2 sum_v alpha_v-alpha_s-alpha_t. Conservation may have additional causes: connectivity and restricted path support can impose more relations than the degree equations.

Let D contain all differences x(P)-x(P0) within each nonempty endpoint sector. Its real kernel is the intersection of the individual sector covariance kernels. Equivalently it is the kernel of any positive weighted sum of these sector covariances. Thus rank(D) is an exact measure of the remaining *linear variation after marking endpoints*. It is not the minimum number of binary edge queries needed to distinguish individual paths.

## 2. Exact square-sum controls

Use Q_n on 1..n with i+j a square and i!=j. The complete path census through32 is supplied by the square-sum packet and independently audited there. Exact rational row elimination gives:

|n|paths|endpoint pairs|total affine rank|within-endpoint rank|
|---|---:|---:|---:|---:|
|23|3|3|2|0|
|25|10|10|9|0|
|26|12|9|11|3|
|27|35|16|20|5|
|28|52|19|25|7|
|29|19|15|16|2|
|30|20|15|17|3|
|31|349|172|34|9|
|32|392|220|45|12|

At26 the within-endpoint rank is3, but the single edge-membership query {3,6} distinguishes every pair of paths sharing endpoints. There is no contradiction: one bit can distinguish two alternatives in each of three different already-marked sectors, while their three difference vectors are linearly independent in the full edge space.

At32 there are46 graph edges. The global affine rank45 is maximal under the common edge-count equation sum x_e=31. No edge is globally forced or unused. Yet conditioning on endpoints reduces the difference rank to12. The signless incidence matrix has rank32, so the common conditional conservation space has dimension34, two more than the degree-potential space.

An exact basis extension is supplied by the unit edge costs at {5,31} and {9,27}: both are individually constant within every endpoint sector, and together with the32 independent degree-potential rows they span the full34-dimensional conservation space. This is a finite32 statement, not an all-size conserved-charge law. Endpoint-conditioned constants can vary between endpoint sectors.

## 3. What transfers to tournament wiring

The reusable mechanism is conditioning on actual owners before comparing variation. A scalar histogram mixes different endpoint sectors; a shared-incidence condition restricts which changes are simultaneously feasible. On a finite candidate wiring ensemble the same covariance identity detects linear conserved costs, but it does not establish that the ensemble contains exactly the native tournament wirings. That feasibility and converse-transport sidecar must still be proved.

The electrical connection is precise but limited: degree conservation supplies a linear kernel, while the separate cotree-bit theorem uses the binary boundary kernel. Real covariance, binary cycle space, and mod8 linking/phase guards are different objects. Equal dimensions or similar-looking matrices do not identify them.

## Reproduction

python3 path_scatter.py
python3 -O path_scatter.py

The script locates the adjacent square-sum packet's frozen census, validates every path, derives all ranks by exact rational elimination, checks all endpoint degree conservation equations, and proves the stated32-dimensional basis extension by ranks. There is no fitted parameter. The paired receipt is path_scatter_receipt.json.

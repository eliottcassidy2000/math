# Shared-core exchange cubes: interaction orthogonality, exact interval selection, and capacity

Status: **PROVED scoped composition and selection criteria; FINITE-EXACT physical eight-vertex instance.** No all-target selector, unrestricted interval overlap, or Chan–Pak completeness is proved.

## 1. Inheritance and deliberate change of architecture

Closest positive predecessors are THM-4960 (`01-canon/theorems/THM-4960-twelve-bit-moving-core-consecutive-odd-intervals.md`) and THM-5000 (`01-canon/theorems/THM-5000-coprime-interface-composition-consecutive-odd-windows.md`). They provide arbitrarily long physical intervals with explicit but exponentially large baselines. Session245 `interfaces/PROOF.md`, sections 6–8, gives the h13 unit cut switch, then proves that a mixed h13/h37 ordinal chain has primitive coefficient gcd 2 but minimum spacing `50*37^(k-2)` for k>=2. That is the hostile example, not a coverage mechanism.

We change the physical architecture: several exchanges act on **one fixed deletion tournament**, rather than multiplying distinct deletion modules. The inherited quadratic ear identity is THM-4114 (`01-canon/theorems/THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature.md`), section 6, based on THM-4097. Quadraticity, submodularity, and its directed-cut interpretation are not new claims here. Full interaction rank is already known not to imply a count interval (THM-4114 section 7).

The contribution here is the exact disjoint-exchange restriction, the negative-inertia capacity bound, an iff bounded-digit interval selector retaining the physical baseline, and a literal independent-switch certificate. The unit-Lipschitz extension in section 6 admits controlled nonzero interaction. These are scoped deductions, with no literature-priority claim.

Concept board: absolute target placement; shared versus multiplied deletion core; interaction matrix; bounded digit coverage versus gcd; retained rooted response. The anchor is absolute placement, still open. The niche is shared-core exchange cancellation. The wildcard is response tomography.

After short moving-core and shared-root probes did not provide a scalable improvement, a reproducible unrelated-topic sample used seed 246071 on sorted `git ls-tree -r --name-only HEAD` paths beginning `01-canon/theorems/`, excluding filenames containing `hamilton`, `spectrum`, or `tournament`. At inherited HEAD 0e140250d07c the first three sampled files were THM-3951-affine-plane-boundary-incidence-forest-and-equianharmonic-survivor-nonentry.md, THM-3490-full-fourier-dirichlet-residue-tomography-of-periodic-polynomial-words.md, and THM-3675-russell-cylinder-critical-fold-formal-conjugacy-closure.md. The useful analogy is THM-3490's retention of a complete response array rather than its lossy scalar recurrence. Here the map is from the ear quadratic to its restriction on chosen exchange directions; it preserves **all** mixed differences and the absolute baseline, while discarding off-cube cuts. The sidecar is the literal deletion graph and cut embedding. The cheapest decisive test is a four-corner mixed difference. This is an analogy, not an imported theorem about tournaments or a claimed equivalence of the two subjects.

## 2. Literal cut polynomial and exchange cancellation

Let C be any nonempty tournament on n vertices. For a cut vector x in {0,1}^n, append a sink v, then a root r beating exactly the selected C vertices and losing to v. Write F(x) for the Hamilton-path count of this **protected closure**. Let h=H(C), and let a_i count ordered spanning two-path covers of C whose second path starts at i, allowing the first path to be empty. For i != j let W_ij be the number of such covers whose first path ends at i and second starts at j, plus the reversed endpoint version. Put W_ii=0. Then

    F(x) = h + 2 a^T x - (1/2) x^T W x.                 (1)

This is the seven-coordinate closure identity h+2A−T, with A=a^T x and T=(1/2)x^T W x. It is also the inherited quadratic ear formula applied to C plus its sink. Both paths must be nonempty in the W terms, so their endpoint vertices are distinct.

Fix a background cut x0 and k pairwise disjoint pairs (p_i,q_i), with x0(p_i)=1 and x0(q_i)=0. Put d_i=e_(q_i)−e_(p_i). Every

    x(b)=x0+sum_i b_i d_i,   b in {0,1}^k,

is a valid literal cut. Define

    g_i = F(x0+d_i)−F(x0),
    c_ij = −d_i^T W d_j.

Exact substitution in (1) proves

    F(x(b)) = F(x0)+sum_i g_i b_i+sum_(i<j)c_ij b_i b_j. (2)

In particular, the cube is affine **iff** d_i^T W d_j=0 for all i != j. Necessity follows from the four-corner difference; sufficiency follows because there are no terms of degree >=3. This is cancellation of signed endpoint incidences inside a shared deletion core. It does not need equality of individual a coordinates or a repeated module.

Every g_i and c_ij is even: each F value is odd by the tournament path-parity theorem, and c_ij is a four-corner difference. This parity assertion concerns the physical cube, not arbitrary formal W matrices.

## 3. Exact absolute interval selector

Suppose the orthogonality condition holds. Let

    B = F(x0)+sum_(g_i<0) g_i,
    w_i = |g_i|/2,
    R = sum_i w_i.

Then the exact physical count set is

    { B+2 sum_i w_i b_i : b in {0,1}^k }.               (3)

Changing the sign of a coefficient merely exchanges which physical choice is called bit zero. Zero weights may be discarded. Sort the remaining positive integers as w1<=...<=wm. The count set is **exactly all R+1 odd numbers in [B,B+2R] iff**

    w1=1 and w_j <= 1+sum_(i<j) w_i for every j>1.     (4)

For m=0 the singleton interval is covered. Sufficiency follows inductively: if the preceding subset sums fill [0,S], adjoining w joins [0,S] and [w,w+S], which touch or overlap exactly when w<=S+1. Necessity follows because if w_j>S+1, the target S+1 uses neither only earlier weights nor any later one. Thus gcd 1 is not substituted for bounded coverage.

A constructive selector processes weights from largest to smallest, subtracting w_j whenever the remaining target is at least w_j. The same inequality proves the remainder stays inside the earlier interval. Undo the sign reversals to recover the actual cut. All targets here are relative to the **computed physical baseline B**, not to an arbitrary prescribed integer.

For an orthogonal family the remaining rooted coordinates are explicit:

    deletion = root-ending count = h,
    two-path-start response U(b) = h+2 a^T x(b).

The last coordinate is generally not constant. It must be carried if these outputs are used in a subsequent rooted operation.

## 4. Capacity is negative inertia, not an isotropic-rank slogan

For every distinct i,j,

    W_ij >= 2 H(C−{i,j}) > 0.                          (5)

Here H(empty)=1. A short proof uses the inherited odd-cycle-collection expansion. In the unprotected one-root insertion, group the cycles through the new root by their first and last old vertices. Such an odd cycle is a simple **odd-length** directed path from i to j, together with the two root arcs. Its contribution is 2 H(C minus the vertices of that path) times x_i(1−x_j). The coefficient of x_i x_j is minus twice the sum over odd-length paths in either direction. The single arc between i and j contributes 2 H(C−{i,j}); all other terms are nonnegative. The quadratic coefficient is the same −W_ij as in (1): the extra sink doubles the linear A term but leaves T unchanged. This proves (5). The cycle-collection identity is inherited, not a new attribution; see CORE-PAPERS.md's Grinberg–Stanley/Irving–Omar source distinction.

For a pair exchange d_i=e_q−e_p,

    d_i^T W d_i = −2 W_pq < 0.

Consequently an interaction-orthogonal k-exchange family has negative-definite Gram matrix `diag(−2W_(p_i q_i))`. Its directions are independent and span a negative-definite subspace. Therefore

    k <= n_−(W),    and (by disjoint support) k <= floor(n/2). (6)

Here n_− is the number of negative eigenvalues with multiplicity. The directions are **not isotropic**. The inertia bound controls how many independent orthogonal exchanges fit in a given physical core; it neither produces such directions nor bounds all tournament counts.

## 5. A shared-core physical certificate

Take outgoing-neighbor masks

    C = [10,28,49,36,41,3],  h=43,
    a = [55,64,72,55,78,57],
    x0 = cut 37 = {0,2,5},
    d1 = e4−e2,  d2 = e3−e5.

The symmetric interaction matrix is

    [ 0,26,14,20,26,28]
    [26, 0,20,18,20,20]
    [14,20, 0,26,28,24]
    [20,18,26, 0,24,22]
    [26,20,28,24, 0,22]
    [28,20,24,22,22, 0].

Cancellation is literal:

    d1^T W d2 = W43−W45−W23+W25 = 24−22−26+24 = 0.

Both individual increments are 2. The four physical choices have

    bits  cut  H    U    deletion/root-end
    00    37   345  411  43
    10    49   347  423  43
    01    13   347  407  43
    11    25   349  419  43.

Every closure is a strong eight-vertex tournament. C itself is strong, every selected cut is nonempty, r reaches C through a selected vertex, and every old vertex reaches r through the sink. Thus the three consecutive odds 345,347,349 live in one fixed deletion core, with no 43-factor multiplying the switch. This is not a claim of a new small realized value or an improvement over the unbounded moving-core interval theorem. Its purpose is to certify the different cancellation architecture and its exact response loss.

## 6. Controlled nonzero interactions: the unit-Lipschitz extension

There is a useful positive relaxation of orthogonality. Normalize (2) by putting g_i=2u_i and c_ij=2v_ij. The physical normalized function is

    f(b) = (F(x(b))−F(x0))/2 = sum_i u_i b_i + sum_(i<j) v_ij b_i b_j.

It changes by at most one across **every** cube edge iff, for every i,

    −1 <= u_i+sum_(j:v_ij<0) v_ij,
     u_i+sum_(j:v_ij>0) v_ij <= 1.                    (7)

This is exact: the minimum and maximum marginal are attained by independently choosing the signs' minimizing or maximizing neighboring bits. Under (7), the range of f is every integer between its minimum and maximum. Connect minimizing and maximizing vertices by a cube path; integer values changing by at most one cannot skip an intermediate target. This proves a physical interval with exact endpoints F(x0)+2 min f and F(x0)+2 max f.

Moreover, subtracting the two inequalities gives

    sum_j |v_ij| <= 2.

Since v is integral, the graph of nonzero interactions is a disjoint union of isolated vertices, paths, cycles, and isolated edges of weight ±2. (Edges in a component with >=3 vertices have weight ±1.) This is a **sharp structural classification with the field inequalities retained**, not a statement that every degree-two interaction graph automatically satisfies (7). Fields must still pass (7).

Minimum and maximum are computable in O(k) arithmetic operations: traverse each path with a two-state dynamic program indexed by the final bit; for a cycle fix its first bit in each of two runs and include the closing edge. Sum extrema across components. Store predecessors to recover minimizing and maximizing physical cuts. A coordinate-by-coordinate path between these cuts then finds any intermediate target in at most k flips. No generic optimization oracle or mere gcd test is needed.

For example the formal chain `sum b_i−sum b_i b_(i+1)` counts selected runs and satisfies (7). This is an algebraic model illustrating the graph classification, **not** a claim that arbitrary such chains have physical tournament embeddings. The actual physical certificate above is the isolated-vertex case u=(1,1), v=0. Finding growing physical path/cycle interaction embeddings with controlled baseline is an exact next gate.

## 7. Scope, verification, and remaining constructive obligation

The standalone verifier reconstructs the core coordinates from literal Hamilton-path subset DP, checks the complete 64-cut polynomial, directly counts all four eight-vertex closures, independently enumerates all 8! permutations for each, verifies strongness, endpoint counts, response values, orthogonality, the bounded-digit selector, and the unit-Lipschitz classification/selector on explicit finite coefficient universes. Explicit raising checks survive optimized Python. `--fail` deliberately raises.

Run `python3 -I -B verify.py` and `python3 -O -I -B verify.py`; their JSON outputs agree. No discovery scan completeness is a theorem dependency.

The new physical cube removes the *local* repeated-module multiplication tax. It does not control the absolute starting count as k grows. An all-target conclusion still requires a growing sequence of literal cores and cut cubes, each passing the complete-weight criterion or (7), together with endpoint placement/overlap covering every odd target except 7,21. Neither a large negative eigenspace nor arbitrary formal coefficient data gives that embedding or baseline theorem. Chan–Pak remains open.

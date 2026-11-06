# Independent audit of the arbitrary-module common-row switch

6 November 2026. **PASS for the exact arbitrary-module profile identities and stated arithmetic scope.** The construction supplies a common padding row and two distinct endpoint values for every nonempty module. Empty modules give one endpoint value. No iterability, freely prescribed base, bounded endpoint menu, or eventual-spectrum threshold is proved.

## Inputs and physical meaning

The initial frozen proof has SHA-256 38a459ef458aeaca337d70303840fc6da571780069f2140c30c2cdb51fd3aff9. The tensor producer has SHA-256 c9ae0630f0cc48af3c282127afb12238b9ac995d0071f849a63f038de70720db, and the literal producer verifier has SHA-256 373873218674d21764eb0c63074bd05fd15a6b79d02c43c8b72f1f9209914106. The supplied six-file manifest passes. Reviewed input snapshots and complete tensor arrays are retained here. The reviewed corrected proof has SHA-256 60f72f178d670ae86a82069a0ab67dca5981131b61e83e5941bbf015abbeb2a5. It changes only one prose sentence to require a nonempty module for distinct endpoints; the original and complete diff are preserved.

A has outgoing rows (102,0,10,83,135,30,182,47), root 0 and replaceable vertex 5. C has rows (142,0,34,22,71,155,47,94), root 0 and replaceable vertex 3. Inflation by X retains all outside arcs and orients every incident bundle uniformly according to the replaced vertex. Empty inflation deletes that vertex. The root survives, and vertex 1 is a global sink in both outputs. The physical outputs remain tournaments; no strong-output assertion is justified.

There is a useful exact surgery interpretation, independently checked by the tensor script. Relabel C in the order (0,1,2,5,4,3,7,6). Its rows become (102,0,10,115,135,22,182,47), with module position 5 and root still 0. It differs from A only by reversing 5->3 to 3->5. For general X, the two physical graphs therefore differ by reversal of the whole incidence bundle between X and core vertex 3. This does not establish a selective-subset flip identity.

## Path-cover conventions and all-order certificate

For j>=1, p_j(X) counts covers by j labeled, nonempty directed paths. The slots are ordered even when their roles coincide. Reordering the slots acts freely, so j! divides p_j. The empty-support coordinate p_0 is 1 only on the empty module. It is not the constant function 1 on all modules.

Every Hamilton path in an inflated context has a unique core order and a unique set of occupied gaps. Each occupied gap contains one nonempty path of X. Conversely a permitted core order and a labeled spanning path cover reconstruct a unique Hamilton path. The gaps' empty and occupied permissions depend only on the finite core and its module incidences. This bijection works for every internal orientation of X.

The producer's cover calculations use possibly bidirected auxiliary graphs. Their use is sound: for R, merging a separator just before the root makes incoming incidences unrestricted while preserving the root's outgoing incidences. Splitting at that root recovers the arbitrary first path and the root-starting second path. The L construction is the dual. For beta, a bidirected separator after root deletion records two ordinary labeled bins, either possibly empty. These are counting devices, not tournament constructions. No tournament parity result is applied to them.

The independent tensor calculation avoids those auxiliary graphs. It partitions the actual fixed core vertices into two ordered paths. For R the second path begins at the root and its initial module gap is forbidden; for L the first path ends at the root and its final module gap is forbidden. For beta, all core cuts are allowed. A path containing no core vertex can be empty or a single nonempty module path. This directly reproduces the same labeled-cover convention.

All 11,520 core orders and 30,240 root/separator cuts are exhausted. Every coefficient of all seven profiles and both derived row coordinates agrees with the supplied full arrays, including constant terms and all vanishing higher-degree coefficients. Hence the complete tensor identity is certified for arbitrary module size, not inferred from the finite literal examples.

## Verified identities and boundaries

For C[X] minus A[X], the only nonzero profile differences are

    Delta h = Delta s = 2(p_1+p_2).

The e,L,R,d,beta differences all vanish. Consequently B=h+L-s and V=R+beta also agree. The independently verified common row is

    B = 93p_0 +253p_1 +204p_2 +48p_3,
    e = 0,
    V = 300p_0 +1022p_1 +1178p_2 +520p_3 +72p_4,
    d = 11p_0 +29p_1 +21p_2 +4p_3.

The remaining A coordinates and all C corrections exactly match the proof's displayed arrays. In particular:

- Empty X gives (h,s,e,L,R,d,beta)=(37,5,0,61,100,11,200) in both contexts, and common row (93,0,300,11). There are not two distinct endpoint choices in this boundary case.
- Singleton X gives common row (253,0,1022,29), with starts 16 and 18.
- For every nonempty tournament X, p_1=H(X)>0, so the endpoint and Hamilton gaps are strictly positive.
- For a transitive module of order m>=1, p_1=1 and p_2=2^m-2. The endpoint gap is therefore 2^(m+1)-2.

Common O/Z padding preserves the collision because both sources are nontrivial and share the actual padding row. The independently audited four-core readout 3B+2s gives an output gap of 4(p_1+p_2). Neither operation removes the dependence of the common base and both endpoints on X's path-cover moments.

## Moment kernel, not realizability

For nonempty inputs p_0 is fixed at zero. The coefficient matrix of B,V,d on (p_1,p_2,p_3,p_4) has rank three. Exact integer cofactors, followed by primitive normalization, give kernel vector

    (-1728,3420,-5427,7768).

Its endpoint changes are 1944 for A and 5328 for C. The vector is primitive over the integers. Requiring moment differences to respect j! divisibility restricts its integer multiplier to multiples of 6: the third coordinate requires an even multiplier, while the fourth requires a multiple of 3. After dividing by the factorials, the primitive unordered direction is

    (-10368,10260,-5427,1942).

These are correct arithmetic conditions. They do not provide actual tournament modules with those differences. Preserving all four moments preserves both endpoint values; preserving only the three row coordinates leaves the displayed potential direction, but its physical attainability is open. Nesting these two contexts generally changes the moments, so row equality alone is not a composition theorem.

## Independent literal verification

The independent C++ checker uses directed-prefix enumeration, rather than the producer's subset recurrence, to recover h,s,e and all induced-set path counts. Cover coordinates are formed from those literal counts. The module moments are counted independently by enumerating every module permutation and every possible cut pattern into nonempty paths, without a cover-convolution recurrence.

Optimized and debug builds agree byte-for-byte on:

- All 76 labeled module tournaments of orders 0 through 4, plus eight specified order-5 modules and four order-6 modules.
- All seven coordinates of 176 actual inflated context graphs.
- 56 literal padding profiles and 24 literal four-core readouts.
- Positive endpoint gaps for all 87 nonempty tested modules and the exact empty boundary.

The order-6 controls include a transitive module, C3->C3, and two other orientations. The checker performs 232 exact profile computations and 24,180,212 directed-prefix visits. All tested graphs have at most 13 vertices, so 64-bit counting is safe. These literal tests are additional controls; the finite core tensor bijection establishes the all-order claim.

The independent Python tensor checker passes normally and with optimization. Both producer programs are also rerun in their normal/optimized or optimized/debug modes and reproduce the frozen receipts.

## Reproduction and publication scope

Run in this directory:

    python3 check_tensors.py > TENSORS_NORMAL.json
    python3 -O check_tensors.py > TENSORS_OPTIMIZED.json
    cmp TENSORS_NORMAL.json TENSORS_OPTIMIZED.json
    g++ -std=c++17 -O3 -Wall -Wextra check_bodies.cpp -o check_bodies
    ./check_bodies > BODIES_NORMAL.json
    g++ -std=c++17 -O0 -g -D_GLIBCXX_ASSERTIONS -Wall -Wextra check_bodies.cpp -o check_bodies_debug
    ./check_bodies_debug > BODIES_DEBUG.json
    cmp BODIES_NORMAL.json BODIES_DEBUG.json

The producer replay uses PRODUCER_TENSOR.py and PRODUCER_VERIFIER.cpp; their receipts are preserved. No tournament-bank completeness assumption is required. The portable release should cite or freeze the previously audited padding and readout dependencies used by the corollaries. No selective bundle flips or later exploratory composition are included in this approval.

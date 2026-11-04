# A per-anchor criterion for even pair quotients

**PROVED:** an explicit induced-subset endpoint-parity test detects some even-order quotient converses with two-vertex children, beyond the unconditional odd-order theorem. It rejects the retained ten-name endpoint-parity hostile and supplies a physical modulo-eight separator for its twenty-name pair inflation. No claim of universal even-order separation is made.

## 1. Setup and inherited mixed identity

Let Q be a fixed named tournament of order q>=3. Replace each quotient label i by the physical pair {2i,2i+1}, directed 2i->2i+1. Let T=Q[TT2,...,TT2] and U=Q-op[TT2,...,TT2], preserving all pair interiors. The coarse attachment bits a_i,c_i are expanded uniformly over each physical pair. Ports remain x->y.

Put Delta=F_T-F_U and delta_Q=(F_Q-F_(Q-op))/2 modulo two. The division is over integers. The proved coefficientwise identity in `../odd_pair_quotient/PROOF.md`, Sections3–5, is

    [Delta]_mixed = 4[delta_Q+B_Q+B_(Q-op)]_mixed modulo8.  (1)

Here B_Q counts two vertex-disjoint directed quotient cycles, one through each port. After removing the ports it is a sum over two disjoint old directed paths P,R, each of length at least two vertices. If their endpoints are s,t and u,v, the characteristic-two contrast contributes

    a_t c_v+a_s c_u
      +a_s a_t(c_v+c_u)+c_u c_v(a_t+a_s).                 (2)

All four endpoint names are distinct. Formula(1) is an equality of complete mixed coefficient data, not a cubic truncation of the physical response.

For each quotient name j define the integer trace

    L_j = sum_(i!=j)[a_i c_j]Delta
            +sum_(j->i)[a_i a_j c_j]Delta.

## 2. The endpoint-subset criterion

For an induced tournament Q[S], let S_j(Q[S]) and E_j(Q[S]) count Hamilton paths starting and ending at the physical name j. They are zero when j is absent. Define the bit

    b_j = sum_(S subset V,
                  |S|>=2,
                  q-|S| is even and at least2)
                  (S_j(Q[S])+E_j(Q[S])) modulo2.          (3)

Then

    L_j = 4(1+b_j) modulo8.                              (4)

Consequently, if b_j=0 for at least one name j, an actual coarse query separates T and U.

**Proof.** The shared-index cubics in(2) are zero. Sum its mixed quadratic coefficients [a_i c_j] over i!=j. For fixed disjoint supports A,S of the two old paths, this sum is

    H(Q[A]) (S_j(Q[S])+E_j(Q[S])) modulo2,

where |A|,|S|>=2. The quotient paths in A and S remain disjoint; no independent endpoint-product assumption without support tracking is used. Rédei parity gives H(Q[A])=1 modulo2. If r=q-|S| names remain, the number of eligible supports A is 2^r-1-r. This is zero for r=0,1, and is odd precisely for even r>=2. Hence the quadratic column correction is exactly b_j.

The incoming [E1 theorem](https://github.com/eliottcassidy2000/math/blob/76b3b1a77d124100f03026420f6b821533cdec15/05-knowledge/results/two_port_parity_separation_20261104/PROOF.md) says the corresponding quadratic/shared-cubic trace of delta_Q is1 at every j for q>=3. Combining with(1) proves(4). A nonzero finite-difference trace implies a nonzero actual attachment contrast. QED.

The same argument explains the odd-order theorem: the total sum of the b_j is zero modulo two, because each summand in(2) has two mixed quadratic terms. If q is odd, the b_j cannot all equal1. At even q they may all be1, and(4) then gives no separator.

For example q=4 has b_j=1 for every j in every tournament: only two-element supports contribute, and there are three partners for j. This is a genuine zero-trace boundary, not a response collision. At q=5, formula(3) simplifies to b_j=outdegree_Q(j) modulo2: on a three-name support the endpoint sum at j is odd precisely when j is its source or sink.

A fixed anchor trace can be recovered from at most 2(q+1+d_j)<=4q actual queries, where d_j is the outdegree of j. For c=0 and c=e_j, use a=0, every singleton, and every pair e_i+e_j for which j->i. These sets contain every lower query needed by the finite differences. The construction does not claim optimal query counts.

## 3. A separate shared-cubic sufficient condition

The cycle-pair contrast(2) has no shared-index cubic. Therefore

    [a_i a_j c_j]Delta = 4 h_ij modulo8,
    h_ij=[a_i a_j c_j]delta_Q,      i!=j.                 (5)

The all-order shared-label identity is

    h_ij = sum_(S subset V-{j}, |S| even)
                     (S_i(Q[S])+E_i(Q[S])) modulo2.       (6)

For completeness, modulo two delta_Q is the contrast of single directed odd cycles in the quotient extension. A term contributing a_i a_j c_j must have the consecutive segment i,x,j,y,k: its attachment product is a_i(1-a_j)c_j(1-c_k). Removing x,j,y leaves a directed path on an even-size support S avoiding j, ending at i. The opposite orientation has an unavoidable c_k and cannot contribute this coefficient; a cycle using x->y cannot have two a-incidences. The length-four possibility k=i is even and excluded. Thus the coefficient from Q is the sum of E_i(Q[S]), while the quotient converse contributes the sum of S_i(Q[S]), proving(6). The empty support contributes zero. Thus any nonzero bit in(6) supplies an eight-actual-query finite-difference certificate, at odd or even quotient order. This is sufficient, not necessary. Coefficients with three distinct old labels and quartics remain part of the full observer and are not discarded.

## 4. The endpoint hostile is still visible

Use the prime ten-name quotient

    Q=(598,508,344,865,40,69,656,317,113,438).

Its all-simple-path and Hamilton-path endpoint parity vectors are balanced; the separate `../even_endpoint_boundary` packet and independent audit retain those exact boundaries. They do not imply vanishing of(3) or(6).

For this Q, the bit mask of b_j is763. Its zero positions are exactly j=2 and j=8, so both L_2 and L_8 equal4 modulo8. There are44 nonzero shared-cubic bits h_ij, including h_(3,0)=1.

More simply, its twenty-name pair inflation already has the physical query

    a={0,1},       c=empty

with

    F_T(a,c)=3 modulo8,       F_U(a,c)=7 modulo8.           (7)

These are counts of the actual22-name tournaments after adding x=20,y=21, not an abstract quotient query or relabeling. `literal_hostile.cpp` independently verifies(7) by predecessor subset DP on all22 physical names. A separate weighted module-run DP verifies all eight queries of the shared-cubic certificate; their alternating contrast is4 modulo8. The full coefficient is not inferred from endpoint imbalance.

For the weighted calculation, a visited pair used in one run has weight q1(TT2)=1, and a pair used in two runs has weight q2(TT2)=2. Three runs are impossible. Modulo8 any run word with three or more twice-used pairs contributes zero. This yields a finite dynamic program over visited pairs, at most two repeated pairs, the used ports and the last quotient label. It is independent of the odd-cycle cancellation proof. Its small controls compare against16,932 literal augmented-tournament queries.

## 5. Fixed ambient contexts and the reversal cube

The separate even-core transfer theorem in `../mod8_transfer/PROOF.md` applies to the actual augmented cores, whose common order is2q+2, hence even. Give both ports the core's exterior incidences. Its ambient multiplier modulo8 is odd, so any nonzero modulo-eight query contrast obtained here stays nonzero in every fixed ambient quotient. This consequence uses the higher-adic transfer theorem, not the weaker modulo-four transfer alone.

In particular(7), whose contrast is4 modulo8, stays4 modulo8 in every fixed ambient context. The root query is unchanged by any pair-interior converse, exactly by module run reversal. Once the quotient orientation is known, an isolated lifted query for each pair reads its converse bit independently of all other pairs. Thus this explicit even-order quotient gives2^11=2,048 distinct complete named responses in every fixed ambient context. Eleven adaptive or twenty-one fixed queries suffice: the one root query plus ten pair decoders, or decoders for both quotient orientations. This state-count conclusion is proved by the isolation argument; no exhaustive2,048-state response bank is claimed.

The unrestricted even-quotient problem remains open. In particular b_j=1 for all j, or vanishing of all shared-cubic bits, is only failure of these sufficient detectors and is not evidence of an exact full-response collision.

## 6. Reproduction

Run `python verify.py > normal.json` and `python -O verify.py > optimized.json`, then compare the files. The self-contained verifier uses Python's standard library and an available C++17 compiler; it installs nothing. Every failure is an explicit exception or nonzero process exit, including optimized execution. Both literal and weighted-run query algorithms use exact arithmetic modulo eight. The independent integer predecessor check reduces only completed path counts modulo two. All four mathematical/source files and both complete numerical data files are bound by SHA256 in the receipt. Compiled executables are temporary and are not part of the publication packet.

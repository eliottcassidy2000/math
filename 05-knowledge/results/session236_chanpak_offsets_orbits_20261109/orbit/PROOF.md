# Fifteen ternary generators: exact orbit cones, ray rigidity, and logarithmic descent

**Recovery provenance.** Reconstructed from the worker's retained exact tool-call history after the ephemeral Session236 worktree disappeared. The mathematical statements, certificates, and algorithms are preserved. This recovery proof is re-typeset, not claimed byte-identical to the lost proof. The checker has the same algorithms and success-output format; diagnostic descriptions and formatting are reconstructed. All outputs in this directory are freshly rerun, not invented historical logs.

**Status:** PROVED for the fifteen listed ternary matrices, using finite exact certificates and all-word induction. The coverage control is FINITE-EXACT. HYP-9640(A), HYP-9640(B), and tournament cofiniteness remain OPEN.

## Scope and inheritance

Number matrices M1,...,M15 in `generators_base3_15.txt` in file order. This is the inherited S18 generator list, originally at `05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/out/generators_base3_15.txt` (blob a3946ae383f63806c4e59f355f994d15665a18a8). Set s0=(1,0,1), T=M1 and v=(3,1,2). A word lists chronological application order, so (4,6) means M6 M4 s0.

By THM-5070, `01-canon/theorems/THM-5070-port-chain-transfer-for-odd-cycle-packing-counts-and-ternary-coverage-to-a-million.md`, each generator is a realizable tournament block of order at most four. This does not identify this family with all tournaments, the 57 ternary matrices, or the 56 binary matrices.

Inheritance: nonnegative port transfer is the closest mechanism; the ten restricted-family gaps are the hostile; the corrected near miss is the former bounded-length tower-direction check; E and O are the indispensable sidecars. Live concepts are inverse rays, invariant cones, offsets, descent, and the distinction between orbit coverage and tournament realization.

## 1. Elementary unimodular cone

Put a=h+E-O, b=O-E, c=E. Then (h,E,O)=(a+b,c,b+c), so the coordinate matrix is

    P=((1,1,0),(0,0,1),(0,1,1)), det(P)=-1.

Every P^{-1}M_iP is an entrywise nonnegative integer matrix, as checked and printed in `verification.out`. Since s0 has coordinates (0,1,0), every orbit state satisfies 0<=E<=O<=h+E. Each row (-1,0,3)M_i is entrywise nonnegative, so h<=3O after any step; the seed also satisfies it.

Each (8,-1,0)M_iP is entrywise nonnegative, giving E<=8h and O<=9h. The E bound is sharp: M4s0=(1,8,8). For i other than 1,5,6, (-1,6,0)M_i is entrywise nonnegative. In the exceptional cases the expressions 6E'-h' are

    i=1: (3O-h)+3(O-E),
    i=5: (3O-h)+3E,
    i=6: (3O-h)+6E+3O.

Hence h<=6E after any nonempty word. These are necessary bounds, not an orbit characterization.

## 2. Exact minimal conical hulls

Let R be the whole orbit and R+ its nonempty-word part.

**The full cone cone(R)** has the four extreme rays, in cyclic projective order:

    (1,0,1),       word (),
    (97,32,33),    word (4,6),
    (25,17,17),    word (4,2),
    (1,8,8),       word (4).

Its primitive inward facet rows are

    (-1,2,1), (-17,-824,849), (0,-1,1), (8,7,-8).

**The nonempty cone cone(R+)** has five extreme rays:

    (193,33,129),  word (4,6,1),
    (292,65,162),  word (4,6,5),
    (97,32,33),    word (4,6),
    (25,17,17),    word (4,2),
    (1,8,8),       word (4).

Its primitive inward facet rows are

    (-3039,6402,2909),
    (-1,2,1),
    (-17,-824,849),
    (0,-1,1),
    (768,1415,-1511).

**Proof.** Normalize the rays to h=1. Their cyclic lists are strictly convex polygons in coordinates (E/h,O/h): adjacent cross products give the displayed facets; each vertex lies on its two adjacent facets and strictly inside all others. These are equivalent vertex and facet descriptions. For every generator M_i, ray r_j, and facet f, verify f M_i r_j>=0. There are 240 inequalities for the full cone and 375 for the nonempty cone. The checker verifies them all using exact integers. All image h coordinates are positive, so both cones are invariant. The seed belongs to the full cone and all fifteen first-step images belong to the nonempty cone. Induction gives orbit containment. Conversely every listed ray is attained by its displayed word, giving reverse containment of the generated cones in the corresponding orbit's conical hull. Thus these are EXACT minimal cones, not sampled approximations. This finite invariant-polytope certificate proves an all-order statement.

Both cones imply the sharp bounds

    E/h<=8, O/h<=8, h/O<=97/33, (E+O)/h<=16.

Every equality is attained by a displayed ray. The nonempty cone further implies the sharp bound E/h>=33/193, attained at (193,33,129).

The integer points of the cone are not the orbit: (10,4,5) satisfies every nonempty-cone facet and h=1 mod3, while 10 is a permanent ternary tournament gap. Convexification discards essential arithmetic information. Cofiniteness is not proved.

## 3. Last-letter rigidity at the tower ray

For every i!=1, the equation M_i s=lambda v, lambda>0, has no solution with h>0 and E,O>=0. Divide by lambda. Matrices 2,3,4 always output E'=O', incompatible with v. Matrices 5,6 always output h'=2E'+O', also incompatible. The remaining exact inverse images of v are

    1: (3/2,1/2,1),
    7: (0,0,1/2),
    8: (1,-1/3,1/3),
    9: (-33/10,5/2,-7/10),
    10: (57/104,7/104,-5/13),
    11: (1,-1/3,-1/6),
    12: (0,1,-1),
    13: (0,0,1/3),
    14: (51/131,10/131,-19/131),
    15: (3/71,14/71,-5/71).

The checker independently recomputes these by rational Gaussian elimination and verifies the singular left-null separators (0,1,-1) or (1,-2,-1).

**All-word theorem.** If A is any product of these fifteen matrices and Av=lambda v for lambda>0, every letter is T. Therefore A=T^r and lambda=2^r. Indeed every generator preserves positive h and nonnegative E,O. The final letter must therefore be T by the preceding obstruction. Its inverse sends lambda v to (lambda/2)v. Induction strips the entire word.

**No reachable tower ray.** The same induction shows no orbit state from s0 is a positive multiple of v, since otherwise stripping towers would put s0 on that ray.

The h>0 condition is essential: M7(0,0,1)=2v and M13(0,0,1)=3v. This is not a full-ambient-cone nonexistence theorem.

## 4. Sharp growth and logarithmic descent

Every first row has h coefficient at least 1 and E coefficient at least 3. Thus h'>=h+3E. On the nonempty orbit, E/h>=33/193, so

    h'/h >= 292/193.

Equality is attained by T acting on (193,33,129): the image has first coordinate 292. Hence this is the SHARP universal next-step factor over all nonempty orbit states and all fifteen generators.

A word of length r>=1 ending at h has h1>=1 at its first step. Therefore

    h >= (292/193)^(r-1),
    r <= 1+floor(log_(292/193)(h)).

Every representing word satisfies this bound. The realizing tournament has order at most 1+4r, hence at most 1+4(1+floor(log_(292/193)(h))). Exact code may compare 292^(r-1)<=h*193^(r-1) rather than use floating point. The simpler bound with base 3/2 also follows from h<=6E.

This gives finite inverse descent: try a last generator, solve for integral predecessor states, restrict to the finite proved sidecar bounds and smaller h, and stop at s0. Singular matrices can have multiple predecessors; logarithmic depth is not a polynomial-time decision algorithm. Binary-encoded orbit-projection membership has certificates of O(log h) generator symbols with polynomial-time exact verification.

All steps after the first increase h by at least 3. The h=1 level is exactly (1,0,1), (1,1,1), (1,2,2), (1,4,4), (1,8,8). Higher levels are obtained from strictly smaller h levels. Forward closure with h<=H is finite and exact without saturation, since E<=8H and O<=9H. These facts describe a recursive state language, not a cofiniteness proof.

## 5. Tower offsets and finite-parallel-ray obstruction

Let ell=(1,2,1), L(s)=h+2E+O. Direct multiplication gives

    Tv=2v, ell^t T=2ell^t, ell^t v=7,
    T^3=I+v ell^t.

Consequently for k>=0,

    T^(3k)s=s+((8^k-1)/7)L(s)v.

The coefficient (8^k-1)/7 is integral. The first coordinate changes by 3((8^k-1)/7)L(s). Thus offset transport depends on L, not merely h, while the transverse class modulo Rv is fixed by T^3.

Let F be a finite union of affine lines parallel to v. If a nonempty word A maps infinitely many points of one such line into F, two distinct images lie on one target line, so subtracting gives Av parallel to v. Also Av has positive h, so its scalar is positive. The all-word theorem forces A to be a tower power. Therefore finitely many parallel-v state families cannot get additional infinite digit branches from non-tower words. This does not rule out varying directions, nonlinear families, higher-dimensional intervals, or offset moves outside these fifteen generators.

## 6. Reproduction and scope

Run `python3 verify_orbit.py` in this directory. It verifies all rational ray certificates, unimodular cone certificates, sidecar bounds, tower identities, both minimal cones, all 615 generator/ray/facet inequalities, and the sharp growth factor. The unsaturated closure through h<=3000 yields 53,732 states and exactly the gaps

    10,31,61,106,217,274,337,346,607,865.

This matches the inherited finite control; it does not extend the prior coverage record and does not imply no larger gaps. The binary56 question is not resolved here.

All verification gates use explicit `require`, never assert. `python3 -O verify_orbit.py` must have identical success output. The hostile flag `--mutate` increments M1[0][0] from 1 to 2 in memory only; normal and optimized interpreters must both reject with exit 1 at `low==(-1,-3,6)`. Fresh recovery runs and output comparisons are recorded in `robustness.out`.

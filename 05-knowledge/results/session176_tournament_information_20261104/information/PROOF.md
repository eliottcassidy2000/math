# What a tournament context can observe: optimal probes and interface information

Session 176. General claims below are proved. Numerical census claims are bounded through order five. No Hamilton-spectrum completeness, new globally attained value, or historical novelty claim is made.

## 1. The interface and its observer

A **uniform one-hole tournament context** consists of an exterior tournament R, with each exterior vertex declared either an in-neighbor or an out-neighbor of the hole. Substituting a nonempty tournament T for the hole orients every edge between a fixed exterior vertex and T the declared same way. Write C[T] for this actual tournament. No formal vector is presumed realizable.

Let q_T(k) count ordered lists of k nonempty directed paths whose vertex sets partition V(T). Paths are directed vertex words, not paths modulo reversal. Put q_T(k)=0 for k>|T|. In particular q_T(1)=H(T), and q_T(|T|)=|T|!.

Collapse each maximal consecutive T-run in a Hamilton path of C[T] to a symbol X. The remaining word uses every exterior vertex once, contains k copies of X, has no adjacent Xs, and obeys all directed adjacency rules. Let c_C(k) count these skeleton words. This gives the exact bijection

    H(C[T]) = sum_(k>=1) c_C(k) q_T(k).                 (1)

Every exterior segment between two Xs starts at an out-neighbor of T and ends at an in-neighbor of T. These are disjoint vertex classes, so the segment contains at least two vertices. Consequently, if C has m exterior vertices,

    c_C(k)=0 whenever 2(k-1)>m.                        (2)

This is independent of the exterior's internal orientation. It is a genuine interface-bandwidth restriction, not an assumption of transitive exterior order.

## 2. Optimal actual-tournament probes

For r>=0 let P_r have exterior vertices in the transitive order

    a_1,b_1,a_2,b_2,...,a_r,b_r,

with all T→a_i and b_i→T. For r=0 the context is empty. Write c_r(k)=c_(P_r)(k). Then

    H(P_r[T]) = sum_(k=1)^(r+1) c_r(k) q_T(k),
    c_r(r+1)=r!.                                      (3)

**Proof of the diagonal.** A skeleton with r+1 T-runs needs r intervening exterior segments, each containing at least two vertices. All 2r exterior vertices are therefore in these segments, exactly two per segment, with no exterior prefix or suffix. Each segment is a_i,b_j, and transitivity permits it precisely when i<=j. The segments pair all a indices with all b indices. If a permutation pi satisfies i<=pi(i) for every i, equality of the two index sums forces pi(i)=i for every i. Thus the segments are exactly (a_i,b_i), in arbitrary order. There are r! orders, all legal because they are separated by Xs. ∎

The initial rows (indexed by k=1,...,r+1) are

    r=0: 1
    r=1: 3, 1
    r=2: 9, 10, 2
    r=3: 31, 83, 50, 6
    r=4: 117, 684, 908, 336, 24.

The checker computes these by an exterior-subset skeleton DP, independently of path-cover counts.

For a tournament of at most n vertices, the n actual observations h_r=H(P_r[T]), 0<=r<n, determine its entire profile by

    q_T(r+1) = (h_r - sum_(k=1)^r c_r(k)q_T(k))/r!.     (4)

These divisions are exact for observations of actual tournaments. They are rational linear recovery on arbitrary formal input vectors.

**Optimality qualification.** To have a nonzero coefficient of the formal coordinate q(k), at least 2k-2 exterior vertices are necessary by (2), and P_(k-1) attains that bound. This does NOT say every actual tournament profile needs k separate coordinates, or that every q(k) is independent on the realizable subset. Actual profiles satisfy relations. Nor is it a lower bound on the size of an unrestricted encoding into one integer.

## 3. Contextual minimality and substitution circuits

**Theorem (exact observational equivalence).** For any two nonempty tournaments A and B,

    q_A=q_B  iff  H(C[A])=H(C[B]) for every uniform context C.

The forward implication is (1). For the reverse, set n=max(|A|,|B|) and use P_0,...,P_(n-1) and triangular inversion (4).

Thus q is the coarsest exact statistic for this observer: any state map s that determines every contextual Hamilton count must satisfy s(A)=s(B) only if q_A=q_B. This is a statement about equality classes, not bits, runtime, or arbitrary scalar encodings.

The equivalence is also a **congruence for actual tournament substitution circuits**. Replacing a module A by B with q_A=q_B leaves the complete profile of the resulting composite unchanged. To prove this, place each P_r around the composite. The resulting combined exterior is still a uniform context for the original hole, so its Hamilton count is unchanged. Recover the composite's profile by the theorem. Repeat one replacement at a time for multiple holes or a finite nested circuit. The hole is a true module at each substitution; arbitrary nonuniform attachments are outside this claim.

For fixed n, let W_n be the rational span of all q_T with |T|=n, and let d_n=dim W_n. Equation (3) shows that the Hamilton contextual response map restricted to W_n has rank d_n. Hence a linear sufficient interface state on W_n needs at least d_n scalar coordinates, and d_n coordinates suffice. By contrast, on the full formal space Q^n the rank is n. Confusing these two ranks would overstate the result.

## 4. Unbounded information rank on realized tournaments

The actual rank d_n is unbounded; no realizability assumption is needed. Define the path-color polynomial

    F_T(t)=sum_k q_T(k) binom(t,k).

For a nonnegative integer t, it counts maps assigning each vertex one of t ordered colors, together with a Hamilton path on each nonempty color class. This follows by choosing the k used colors and assigning the k ordered paths to them. If A⊕B has every edge A→B, each monochromatic path is uniquely a path in A followed by a path in B, omitting empty pieces. Consequently

    F_(A⊕B)(t)=F_A(t) F_B(t).                          (5)

Equality for nonnegative integers proves the polynomial identity. A singleton has F=t. A directed triangle has q=(3,6,6), hence F=t^3+2t.

For 0<=j<=floor(n/3), take T_j to be an ordered sum of j directed triangles and n-3j singletons. These are actual n-vertex tournaments and

    F_(T_j)(t)=t^(n-3j)(t^3+2t)^j
              =t^(n-2j)(t^2+2)^j.

Their lowest nonzero powers are the distinct exponents n-2j, with coefficients 2^j. The polynomials are linearly independent: in any dependence use the smallest exponent present, then proceed upward. The change of basis from binom(t,k) to powers of t is invertible, so their q profiles are linearly independent. Therefore

    d_n >= floor(n/3)+1.                               (6)

This proves a growing linear-interface requirement for exact uniform-context Hamilton observations. It does not preclude nonlinear or packed integer states, efficient specialized circuits, or compression for a restricted observer family.

The exhaustive labelled census through n=5 gives d_n=1,1,2,2,3, respectively, already showing why the formal rank n should not be attributed to realized tournaments. These five values are finite verified results, not a conjectured general formula.

## 5. Explicit loss of scalar information

Outgoing-neighbor masks on vertices 0,...,4:

    A=(16,9,3,5,14), q_A=(15,78,198,240,120),
    B=(24,17,3,6,12), q_B=(15,90,210,240,120).

Both have H=15. Their first nontrivial probe is the uniform cyclic wrapper T→a→b→T, for which c_1=(3,1). Thus

    H(P_1[A])=3*15+78=123,
    H(P_1[B])=3*15+90=135.

These are strong seven-vertex tournaments: each vertex can reach and be reached through the exterior directed two-step bridge. Therefore even equal-order modules with equal H can produce different strong outputs through a minimal two-vertex uniform observer. This is information relevant to spectrum synthesis, not a new global spectrum value.

Outside uniform contexts, even the full q profile is insufficient for a labelled interface. Let A be the two-vertex tournament 0→1 and B its reversal. Both have q=(1,2). Add x with x→0 and 1→x. The first becomes a directed triangle with H=3; the second is transitive with H=1. The external vertex distinguishes individual vertices, destroying the module hypothesis. Endpoint/subset data from Session 175 is therefore genuinely additional information, not redundant detail recoverable from q.

## 6. Interpretation for the spectrum program

Uniform substitution has an exact compositional information type: q. Formal circuit simplification is sound on actual modules when it respects this type, but numerical formal profiles need not be realizable. The probes demonstrate that discarded profile entries can be read back by genuine tournament contexts; the triangle family proves the required linear dimension cannot stay bounded as order grows.

This does not produce a scalable constructor for every permitted odd Hamilton count. It identifies what a proposed scalar recurrence must prove: either its restricted contexts annihilate the lost directions, its realizable states lie in a justified lower-dimensional family, or additional state is required. Nonuniform constructions require a finer type and explicit endpoint/orphan accounting. The Session 175 global exclusions 7 and 21 and unresolved prime/7p lanes remain unchanged.

## 7. Reproduction and scope

Run `python3 check.py > certificate.json` and `python3 -O check.py > certificate.optimized.json`; compare receipts. No packages are required. The checker enumerates every labelled tournament through order five; independently checks H by permutations; checks exact path-cover recurrences, skeleton transfer, probe inversion, rank computations, and explicit scalar collisions. Deep probes are checked by full subset Hamilton DP when the substituted order is at most eleven. A separate `audit.py` independently checks the explicit examples by permutations and the realized lower-rank witness family by polynomial arithmetic.

Prior mathematical context: Session 175 construction/README.md (path-cover wrapper and ordered sums), interface/PROOF.md (endpoint and orphan information), and synthesis/NOTE.md (scope of the spectrum problem). Those files are inherited results, not asserted external novelty support.

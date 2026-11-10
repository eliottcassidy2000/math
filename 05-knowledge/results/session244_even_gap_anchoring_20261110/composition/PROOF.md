# A synchronized endpoint-gap compiler and its primitive self-feeding face

Status: PROVED identities and constructions; finite exact checks PASS. Chan–Pak completeness remains OPEN.

## Inheritance and scope

Every even gap, infinitely many witnesses for each gap, the O/Z/J transfer laws, and positive endpoint universality are inherited. Sources are Session199 physical interpolation, Session221 physical pairs, THM-5001, and Session243 `constructions/GAP_UNIVERSALITY.md` at parent revision c6b359b9c0ca. Session243's pair JOJOZ/JOOOZ has lower count 9D−3E+9U, always divisible by3. The present comparison uses a different coefficient row and, more importantly, synchronizes the endpoint and deletion coordinates at the output. It is not a new proof of even-gap existence or a theorem about arbitrary absolute placement.

THM-5000's coprime moving-interface windows have exponentially growing baselines. They do not supply an additive operation on arbitrary scalar gaps. THM-5109 and HYP-9656 concern different constrained response closures and likewise do not close this placement step.

## 1. Physical transfer and exact compiler

For a rooted tournament write (D,E,U,d) for total Hamilton paths, paths ending at the root, ordered spanning two-path covers with the distinguished path starting at the root (the other path may be empty), and the total count after deleting the root.

O appends a sink and keeps the root. Z appends a vertex that beats exactly the old root and keeps the root. J makes the same graph addition but marks the new vertex. The inherited maps are

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D+U−E,d,2U−d,d),
    J(D,E,U,d)=(D+U−E,D−E,D+U,D).

Letters are executed from left to right. Direct multiplication gives

    JJJJZ: (9D−3E+11U, 3D−E+3U, 9D−3E+11U, 3D−E+3U),
    JOZJZ: (9D−E+11U,  3D−E+3U, 9D−E+11U,  3D−E+3U).

Thus the two strong outputs differ in total count by exactly 2E, while sharing both endpoint count and deletion count. Both have U=H and E=d. This is an identity of integer linear maps on the full four-coordinate state, not a numerical fit.

If the source is strong, the root reaches every vertex. Each O preserves that fact, and each J creates a new root that first reaches the old root and then everything. A final Z makes the graph strong: every old nonroot reaches the new vertex, which reaches the old root, while the old root reaches every old vertex and hence the new vertex through an old nonroot. The sources used here have at least three vertices, so this last intermediary exists.

By inherited endpoint universality, set E=N for any positive integer N, including7 and21. The compiler gives strong same-order pairs of gap2N with at most

    max(10, ceil(log2(N+1))+8)

vertices. Increasing the order of the endpoint selector while preserving its binary selection keeps E=N and makes D unbounded. Since U≥D, the lower count is at least20D−3N and is unbounded. Hence infinitely many distinct pairs arise. These global even-gap statements are inherited in strength; the new feature is the synchronized output face and coefficient row.

The lower row (9,−3,11,0) is primitive, unlike the prior common factor3. This alone does not imply primality or coprimality. Literal verified examples do demonstrate that there is no fixed composite obstruction:

- source adjacency rows [20,1,2,7,14], root2, state(9,3,11,5): counts193 and199;
- source rows [24,1,3,6,14], root0, state(11,1,13,1): counts239 and241;
- the same source rows, root3, state(11,2,20,5): counts313 and317.

Rows encode outgoing-neighbor bitmasks, with vertices numbered from0. All six output counts are prime by exact trial division, and all outputs are literal ten-vertex strong tournaments. This is finite evidence, not prime infinitude or Polignac.

## 2. Self-feeding face and exact gcd control

Restrict to physical states (h,e,h,e). The two words become

    F−(h,e)=(20h−3e,6h−e),
    F+(h,e)=(20h−e,6h−e).

The output again denotes the full state (h',e',h',e'). Both h and e are odd, since e=d is the Hamilton-path count of a deletion tournament. For any odd integers h,e, regardless of physical realizability,

    gcd(20h−3e,6h−e)=gcd(h,e),
    gcd(20h−3e,20h−e)=gcd(5h,e).

For the first identity, a common divisor of the left pair divides2h and2e by elimination. It is odd and therefore divides gcd(h,e); the reverse divisibility is immediate. For the second, the difference is2e, so oddness reduces the gcd to gcd(20h,e)=gcd(5h,e).

Consequently F− preserves primitive physical pairs. If gcd(h,e)=1 and5 does not divide e, the two compiled total counts are coprime. This is a controlled coprime family, not all-gap coprimality.

The inherited single-letter/small-word maps on this face include

    Z:   (h,e) -> (2h−e,e),
    JZ:  (h,e) -> (3h,h),
    JJZ: (h,e) -> (6h−e,2h−e).

Appending k copies of Z to both new outputs keeps their common endpoint and deletion counts and multiplies their gap by2^k. Explicitly, Z^k sends (h,e) to (2^k h−(2^k−1)e,e). This is a genuine common continuation, but not an operation that translates an arbitrary prescribed odd anchor by an arbitrary even gap.

## 3. Why these facts do not close the conjecture

Knowing E=N requires no control of D or U. The target equation for the lower output is

    M=9D−3N+11U.

The endpoint selector gives no solution of this exact equation for an arbitrary odd M. Even on the self-feeding face, the inverse of F− is

    h=(h'−3e')/2, e=3h'−10e'.

The inverse must be integral, odd, positive, and physically realizable. It does not cover the physical cone merely because its determinant is−2. F+ has determinant−14 and a further divisibility constraint. Coprimality of generated pairs supplies no missing additive gluing law.

The deletion coordinate e on the face is itself an attained total count. In particular e=7 or21 cannot occur there by the inherited exclusion theorem. Endpoint universality on unrestricted rooted tournaments cannot be silently promoted to endpoint universality on this face.

An exploratory enumeration of all face-return words of length at most8 finds536 distinct matrices. The separate `FACE_PROBE.json` records190 first-return representatives and inverse tests for **invertible matrices only**, on the broad necessary cone h≥2e+1. It is not a physical census, excludes rank-one predecessors, and asserts no global obstruction or completeness. The narrow inverse images and the uncontrolled scalar anchor remain the unresolved steps.

## 4. Reproduction and controls

From this directory:

    python3 -I -B verify.py > replay.json
    cmp replay.json RESULTS.json
    python3 -I -B -O verify.py > optimized.json
    cmp optimized.json RESULTS.json
    python3 -I -B verify.py --fail

The last command must fail nonzero. The verifier checks full matrix identities;855 literal one-letter transfers on every rooted labelled tournament of orders1 through4;24 endpoint-selected pairs N=1,...,24 with total, endpoint, deletion counts and strongness checked independently by subset DP; all four coordinates of three literal prime-pair witnesses;31375 integer gcd cases; and the exact two coefficient-row pair classes at length5 with terminal Z. `physical.py` is a credited excerpt of Session243's standard-library graph verifier. Assertions use explicit exceptions and survive Python optimization. Finite tests check the implementation, while the displayed transfer identities and arguments supply the unbounded quantifiers.

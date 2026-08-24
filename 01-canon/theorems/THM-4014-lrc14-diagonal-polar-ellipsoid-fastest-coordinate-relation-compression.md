---
id: THM-4014
title: "LRC(14) diagonal-polar ellipsoid and fastest-coordinate relation compression"
status: >
  PROVED ALGEBRA + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. A
  hypothetical primitive LRC(14) counterexample forces an H-weighted Graver
  relation for every feasible positive diagonal polar ellipsoid. The
  fastest-coordinate choice gives ||a||_1<=49 and |a_M|<=7; if a_M is
  nonzero then sum a_i^2<=193. The l1=49 boundary is rigid and requires
  M/M_2<3/2 and 4|M. A simultaneous linear bank retains every speed weight
  and reduces its top-pair support-two branch to 28 ratios. These are
  necessary reductions, not a proof of LRC(14).
source: root + lift_tower_covering / LRC14 frontier session, 2026-08-24
audit: >
  PASS. The independent audit rederived the polar vertices, exact pair
  criterion, compressed-operator polarity, strict transference map, Graver
  reduction, max-coordinate and simultaneous banks, kernel-balance cap,
  integer l1 boundary, and support-two census. It found and repaired the
  false ambient H^{-1} polarity identification using the exact hostile
  n=(1,2), H=diag(2,1). The corrected operator is A=P_V H|_V. A standard-
  library exact companion checks 156 polar vertices, all 78 pair facets,
  the hostile rank-one correction, integer optimizations, and 28 reduced
  top-pair ratios. Normal and optimized runs match the frozen LF output
  after platform-newline normalization.
depends_on:
  - THM-4009-euclidean-covering-transference-short-relation-compression
related:
  - THM-1008-lrc13-descent-floor
  - THM-1043-the-spread-ladder-and-the-refined-covering-map
  - THM-2052-finite-height-forces-high-rank-bounded-relation-code
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-4002-lrc14-signed-endpoint-cross-phase-and-fixed-scale-two-family
  - THM-4003-lrc14-scale-two-component-erosion-boundary-strip
  - THM-4004-lrc14-three-detuned-divisor-comb-profile
audit_script: 04-computation/lrc14_diagonal_ellipsoid_refinement_independent_audit_20260824.py
audit_output: 05-knowledge/results/lrc14_diagonal_ellipsoid_refinement_independent_audit_20260824.out
audit_report: 05-knowledge/results/lrc14_diagonal_ellipsoid_refinement_independent_audit_20260824.md
audit_script_sha256: dc71fc2ef24515e5eb5e6562b8e6d21f04d75b31c2f59c06d40e3545730d65c3
audit_output_sha256: bd3793f54f1ecb0a61671d01a12d2089b9dc995f4e98ad83cda5732024bf3224
audit_report_sha256: e628268f74099f922e6e04f028e86b552e5d70b095a09c70157e40247676f10e
hash_basis: raw LF bytes
---

# THM-4014 -- diagonal polarity localizes the short LRC(14) relation

**PROVED ALGEBRA + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** This
theorem strengthens THM-4009's necessary short-relation reduction. It does
not prove LRC(14).

## 1. Exact polar and feasible diagonal metrics

Let `n=(n_1,...,n_13)` be primitive with distinct positive integer entries,
let

```text
V=n^perp,                 P:R^13 -> V,
Lambda=P(Z^13),           K=P([-1,1]^13).               (1)
```

All polarity below is relative to `V`. For `u in V`,

```text
h_K(u)=sup_(z in [-1,1]^13)<u,Pz>=||u||_1,
K polar={u in V:||u||_1<=1}.                            (2)
```

The vertices of this hyperplane section of the cross-polytope are exactly

```text
v_ij^+-=+/- (n_j e_i-n_i e_j)/(n_i+n_j),   i<j.        (3)
```

There are `2*C(13,2)=156` signed vertices. Indeed no cross-polytope vertex
lies in `V`, because every `n_i` is positive, and a hyperplane-section vertex
must therefore occur where `V` crosses an edge joining `e_i` to `-e_j` or
the opposite edge.

Let `H=diag(h_1,...,h_13)` be positive. Evaluation at (3) shows that

```text
K polar subset {u in V:u^T H u<=1}                     (4)
```

holds if and only if, for every pair,

```text
(h_i-1)n_j^2+(h_j-1)n_i^2<=2n_i n_j.                  (5)
```

Necessity follows at the vertices. Sufficiency follows because the quadratic
form is convex and every point of the polar polytope is a convex combination
of its vertices. Positivity of `H` is load-bearing.

## 2. The compressed-operator polarity repair

Put

```text
A=P H|_V:V -> V.                                       (6)
```

The ellipsoid in (4) is represented on `V` by `A`, so its polar is

```text
Q_H^o={x in V:<A^(-1)x,x><=1}.                         (7)
```

It is not generally represented by ambient `H^(-1)`. Solving `Ay=x` with
`x,y in V` gives the exact identity

```text
<A^(-1)x,x>
 =x^T H^(-1)x-(n^T H^(-1)x)^2/(n^T H^(-1)n).          (8)
```

Thus the ambient-inverse ellipsoid is a generally strict subellipsoid of
(7). The equality error is not harmless: it would yield the wrong dual
quadratic form and would not imply the weighted relation below. The exact
hostile

```text
n=(1,2),       H=diag(2,1)                             (9)
```

saturates (5). On a unit vector of `V`, the restricted `H` operator is `9/5`,
so its inverse is `5/9`, whereas the ambient inverse compresses to `3/5`.
At squared radial coordinate `7/4`, the values are `35/36<1` and `21/20>1`.

## 3. Weighted transference theorem

THM-4009 writes the lonely zonotope as

```text
Z(n)=c+(3/7)K,             c=(1/2)P(1).                (10)
```

Assume `n` is a counterexample. The closed zonotope is disjoint from
`Lambda`. From (4) and polarity reversal, the closed ellipsoid

```text
c+(3/7)Q_H^o
```

is also lattice-point-free. Applying `A^(-1/2)` turns it into a closed
Euclidean ball of radius `3/7` for the lattice

```text
L=A^(-1/2)Lambda,             L*=A^(1/2)Lambda*.       (11)
```

Closed disjointness and discreteness give `mu(L)>3/7`. THM-4009's cited
twelve-dimensional Banaszczyk inequality gives

```text
lambda_1(L*)<14.                                       (12)
```

The exact dual identity in THM-4009 is

```text
Lambda*=Z^13 intersect n^perp.                         (13)
```

Therefore a counterexample forces a nonzero integer speed relation with

```text
a dot n=0,                 a^T H a<196.                (14)
```

Choose an `H`-shortest nonzero kernel vector. A nontrivial conformal
decomposition would contain a nonzero summand whose absolute coordinates are
bounded by those of `a`, strictly somewhere. Positive diagonality would make
that summand strictly `H`-shorter. Hence `a` can be chosen Graver.

The quantifier in (14) is

```text
for every feasible H, there exists an H-short Graver relation a(H).         (15)
```

Different metrics need not select the same or distinct relations.

## 4. Fastest-coordinate compression

For a chosen coordinate `i`, let `M_i=max_(j!=i)n_j` and set

```text
h_i=1+2n_i/M_i,             h_j=1 for j!=i.            (16)
```

Condition (5) is immediate and is sharp on a pair containing the largest
other speed. Let `M` and `M_2` be the largest and second-largest speeds, let
`a_M` denote the coefficient at `M`, and put `S=sum a_i^2`. Equations
(14)--(16) give one Graver relation satisfying

```text
S+2(M/M_2)a_M^2<196,
M_2 S+2M a_M^2<=196M_2-1.                             (17)
```

Consequently

```text
S<=195,
a_M!=0  =>  S<=193.                                   (18)
```

Weighted Cauchy gives

```text
||a||_1^2
 <=(12+1/(1+2M/M_2)) a^T H a
 <(37/3)196<50^2,                                     (19)
```

so

```text
||a||_1<=49.                                           (20)
```

The kernel equation improves the raw coefficient cap. Put

```text
T=sum_(j!=M)n_j^2,             r=M/M_2>1.
```

Cauchy applied to `M a_M=-sum_(j!=M)n_j a_j`, together with
`T<=12M_2^2`, gives

```text
a^T H a
 >=a_M^2(1+2r+M^2/T)
 >=a_M^2(1+2r+r^2/12)
 >a_M^2(37/12).                                       (21)
```

Since `64*(37/12)>196`,

```text
|a_M|<=7.                                              (22)
```

The exact integer boundary in (20) is rigid. At fixed `L=||a||_1` and
`t=|a_M|`, the twelve other squares are minimized by balancing `L-t` among
twelve integers. Replacing the actual maximum weight by its strict lower
limit `3`, the minima are

```text
L=49: 195, uniquely at t=1 and twelve absolute values 4;
L=50: 204.                                             (23)
```

Thus equality `||a||_1=49` requires

```text
|a_M|=1,             |a_j|=4 for every j!=M.           (24)
```

Its weighted energy is `193+2M/M_2`, so (14) requires `M/M_2<3/2`; the
kernel equation in (24) also forces `4|M`. Therefore

```text
M/M_2>=3/2 or 4 not dividing M  =>  ||a||_1<=48.       (25)
```

The weaker sufficient ratio that would force `a_M=0` using only
`T<=12M_2^2` is `M/M_2>=-12+2sqrt(621)>37`. It lies outside THM-1008's live
counterexample range `M/M_2<13` and is calibration only.

## 5. One simultaneous speed-weighted relation

Put `D=M+M_2` and choose

```text
h_i=1+2n_i/D                    for every i.            (26)
```

For every pair, the left side of (5) is

```text
2n_i n_j(n_i+n_j)/D<=2n_i n_j,                         (27)
```

with equality only on the top pair. Hence a counterexample forces one Graver
relation satisfying all weights simultaneously:

```text
D S+2 sum_i n_i a_i^2<=196D-1.                         (28)
```

It obeys `|a_M|<=9` and `||a||_1<=49`. The latter follows from
`sum_i 1/h_i<37/3`: the eleven non-top reciprocal weights are below one and
the top two contribute

```text
4(r+1)^2/((3r+1)(r+3))<4/3.                            (29)
```

This bank localizes coefficient mass but does not improve THM-4009's scalar
square cap. If the selected row is supported on the top speeds `s<t`, write
`p=s/g`, `q=t/g`, `g=gcd(s,t)`. Its coefficients are `(q,-p)` up to sign,
and (28) becomes

```text
(p+q)^2<196,                 p+q<=13.                  (30)
```

Exactly 28 coprime pairs `p<q` satisfy (30), compared with THM-4009's 47
ambient support-two ratios. This is a branch refinement; transference does
not force support two or top-pair support.

## 6. Scope and reproduction

The theorem preserves a speed-weighted metric and a fastest-coordinate
coefficient. It still forgets support, component crossing, order-two centre
character, endpoint owner, and Fourier phase. AP13 has much shorter internal
relations while remaining lonely, so no norm statement here is sufficient.
LRC(14) remains open.

The intersection with incoming THM-4003/4004 makes this boundary exact. In
THM-4003's scale-two row `2u direct-sum t(1,9)`, one has `M=9t` and
`M_2=2U`, so (25) gives `||a||_1<=48`. A selected support-two row on that top
pair would have `p+q<=13` and hence height at most twelve; because the pair
crosses the rank-eleven component cut, THM-3818 forbids this subbranch. No
parameter cell closes: the internal tail row

```text
9(t)-1(9t)=0
```

has `l1=10` and max-boost energy `82+9t/U<94`. In THM-4004's divisor lane, a
body containing `1,2,3` can likewise absorb the metric through `1+2-3=0`
while `t` is coprime to the body. To force crossing, a future metric must
control the **internal weighted minima** of both components, not only the
global maximum coefficient.

The incoming THM-4009 support-two owner audit makes the label loss explicit.
Among THM-3910's seventeen external pair types, fifteen are compatible with
the old square cap; `(8,21)` and `(9,11)` are incompatible only when the
selected support is extra/extra. A body/extra pair involving `tp` necessarily
satisfies `tp<=13U`, while AP11 has a norm-five body/body absorber for every
external type. Consequently neither the 47 old ratios nor the 28 top-pair
ratios may be intersected with the seventeen types without carrying support
owner and component incidence.

The independent exact report is
[the diagonal-ellipsoid audit](../../05-knowledge/results/lrc14_diagonal_ellipsoid_refinement_independent_audit_20260824.md).
Reproduce its frozen output with

```text
python3 -B 04-computation/lrc14_diagonal_ellipsoid_refinement_independent_audit_20260824.py
python3 -B -O 04-computation/lrc14_diagonal_ellipsoid_refinement_independent_audit_20260824.py
```

Both streams match the recorded LF output after platform-newline
normalization. **QED.**

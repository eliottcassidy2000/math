# One-pair triangle repair of the projected-cube ellipsoid

**Status: PROVED by the explicit certificate below, with an exact finite
companion. Session 19, 2026-10-21. Separate addendum to the frozen
`session19-ellipsoid-optimum` packet; that packet is unchanged. No
literature-priority or global LRC-bound claim.**

## 1. Scope and result

Let `v=(a_1,...,a_m,M)>0`, `M=max v_i`, `m=k-1>=2`, `H=v^perp`, and
`D=P[-1,1]^k`. Use the base packet's invertible star coordinates

```
z_i=(M x_i-a_i x_k)/(M+a_i),     i=1,...,m.              (1)
```

The exact description of `L(D)` consists of axis slabs `|z_i|<=1` and,
for every leaf pair `i,j`,

```
|A_ij z_i-B_ij z_j| <= 1,
A_ij=a_j(M+a_i)/[M(a_i+a_j)],
B_ij=a_i(M+a_j)/[M(a_i+a_j)].                           (2)
```

These coefficients are positive and at most one: for example `A_ij<=1`
is equivalent to `a_j<=M`. Suppose **exactly one** leaf pair, labelled
`a,b`, violates the unit-ball test:

```
A^2+B^2>1,
A_ij^2+B_ij^2<=1 for every other pair.                 (3)
```

Set

```
D_ab=-M^2+M(a+b)+ab,
c=(A^2+B^2-1)/(2AB)=D_ab/[(M+a)(M+b)] > 0,
S=I_m+c(e_a e_b^T+e_b e_a^T).                          (4)
```

Then `L^(-1)(S^(1/2) B_2^m)` is the **unique maximum-volume ellipsoid
inside D, including translated ellipsoids**. Its squared volume divided
by `vol(B_2^m)^2` is

```
Delta_triangle
 = Delta_star (1-c^2)
 = 4[M(a+b)+ab]/||v||_2^2
     * product_(i notin {a,b}) (1+a_i/M)^2.             (5)
```

Here `Delta_star=M^2 product_i(1+a_i/M)^2/||v||_2^2` is a universal
star-parallelotope volume bound, not the volume of a feasible ellipsoid
under (3). The strict `c>0` distinguishes this repair from the base cone.

If the leaves are ordered `a>=b>=d>=...`, then (3) is equivalent to

```
M^2 < M(a+b)+ab,
M^2 >= M(a+d)+ad                                      (6)
```

when `m>=3`; there is no second condition when `m=2`. Indeed,
`M(s+t)+st` is increasing in both arguments, and `(a,d)` is the largest
remaining pair after `(a,b)`. Ties are counted with multiplicity; a tied
second leaf may make the single-violation regime empty, as (6) correctly
records.

## 2. Feasibility and positive definiteness

Since `0<A,B<=1`, one has `|A-B|<1`, and therefore

```
0 < c < 1.
```

Thus `S` is positive definite, with eigenvalues `1+c,1-c` on the affected
two-dimensional block and eigenvalue one elsewhere. Its axis supports are
all one. The exceptional pair's squared support is

```
A^2+B^2-2ABc=1.
```

Every other pair functional uses at most one of the affected coordinates,
so its squared support is unchanged from the unit ball and is at most
one by (3). This proves containment. The support-two nature of the facet
normals is essential: an unrecorded higher-support constraint could see
the changed off-diagonal entry.

## 3. Explicit global optimality certificate

Let `n=A e_a-B e_b` be the exceptional slab normal. Define

```
lambda_pair=c/[AB(1-c^2)],
lambda_a=(1-c A/B)/(1-c^2),
lambda_b=(1-c B/A)/(1-c^2),
lambda_i=1 for i notin {a,b}.                          (7)
```

Every weight is strictly positive. For example

```
1-c A/B=(1+B^2-A^2)/(2B^2)>0
```

because `A<=1` and `B>0`; the other formula is symmetric. Direct block
multiplication gives

```
S^(-1)=sum_i lambda_i e_i e_i^T + lambda_pair n n^T.    (8)
```

All normals in (8) are contact normals: `e_i^T S e_i=1` and
`n^T S n=1`. Taking the trace after multiplying by `S` gives

```
sum_i lambda_i+lambda_pair=m.                         (9)
```

Now let `h+T^(1/2)B_2^m` be any ellipsoid in `L(D)`, with covariance
`T>=0`. Every contact normal `w` in (8) satisfies

```
|w dot h|+sqrt(w^T T w)<=1,
so w^T T w <= (1-|w dot h|)^2 <= 1.                   (10)
```

Multiply (10) by the positive weights and sum. Equations (8)--(9) yield

```
tr(S^(-1)T)<=m.
```

Apply arithmetic--geometric mean to the nonnegative eigenvalues of
`S^(-1/2) T S^(-1/2)`:

```
det(T)/det(S) <= [tr(S^(-1)T)/m]^m <= 1.               (11)
```

If equality holds, all eigenvalues are one, hence `T=S`; equality in the
strictly positive axis terms in (10) also forces `h_i=0` for every `i`.
Thus the maximizer is unique, including its center. Singular competitors
have zero volume. This elementary certificate needs no external
optimality theorem or numerical optimizer.

## 4. Determinant simplification

The base packet proves the Jacobian factor `Delta_star=1/det(LL^T)`.
The repaired covariance has determinant `1-c^2`. To simplify, let

```
K=(M+a)(M+b)=M^2+M(a+b)+ab.
```

Then

```
K^2-D_ab^2=4M^2[M(a+b)+ab].                            (12)
```

Cancel the two corresponding star factors to obtain (5). All expressions
are homogeneous of degree zero in the speeds, as required because
scaling `v` does not change `H`, `P`, or `D`.

## 5. Complete three-speed classification

When `k=3`, there is only one leaf pair, so the base star case and this
repair cover **every positive triple**, after choosing a largest speed
as `M`. Write the other two speeds as `a,b`.

If `D_ab<=0`, the unique maximum ellipsoid is the base star ellipsoid,
with diagonal lift

```
(q_a,q_b,q_M)=((1+a/M)^2,(1+b/M)^2,0).
```

If `D_ab>0`, it is represented by the strictly positive diagonal lift

```
q_M=D_ab/(ab),
q_a=(a+b)/b+a(b-a)/(bM),
q_b=(a+b)/a+b(a-b)/(aM).                               (13)
```

Indeed, in star coordinates a projected diagonal lift has covariance

```
S_ii=[M^2 q_i+a_i^2 q_M]/(M+a_i)^2,
S_ij=a_i a_j q_M/[(M+a_i)(M+a_j)]  (i!=j),              (14)
```

and (13) gives diagonal entries one and off-diagonal entry `c`. Positivity
is transparent from

```
q_a=1+a/M+(a/b)(1-a/M)>0,
q_b=1+b/M+(b/a)(1-b/M)>0.
```

At `D_ab=0`, (13) agrees continuously with the star lift. The exact
maximum determinant is therefore

```
if M^2 >= M(a+b)+ab:
    Delta_max=(M+a)^2(M+b)^2/[M^2(M^2+a^2+b^2)],
otherwise:
    Delta_max=4[M(a+b)+ab]/(M^2+a^2+b^2).               (15)
```

For `(a,b,M)=(1,2,3)`, one gets `c=1/10`, the lift
`(5/3,7/3,1)`, and `Delta_max=22/7`. The construction used in Allikvere's
[arXiv:2609.02604v2, Section 3.1](https://arxiv.org/html/2609.02604v2#S3.SS1)
has determinant `64/21` here, so the squared volume ratio is `33/32`.
This is only a fixed-shape comparison. For tied largest speeds that
source construction can already be optimal, so no universal strict
comparison is asserted for every positive triple in this region.

## 6. A precise failure of the entire diagonal class in higher dimension

When `k>=4`, there are at least three leaf coordinates. Formula (14) shows
that all off-diagonal entries of a diagonal lift are a common scalar
`q_M` times strictly positive factors. They are therefore all zero, or all
nonzero. The repaired `S` in (4) instead has exactly one nonzero unordered
off-diagonal pair. It cannot be represented by **any real diagonal Q**,
and in particular not by a nonnegative diagonal lift.

The nonnegative diagonal feasible set is compact: every coordinate has a
finite upper bound from a pair constraint. Its volume maximum is attained.
By uniqueness in (11), attaining the global optimum would require the
same covariance `S`, which has just been ruled out. Consequently

```
maximum volume over all ellipsoids in D
  > maximum volume over nonnegative projected diagonal lifts             (16)
```

for every `k>=4` in the single-violating-pair regime. This does not quantify
the gap or identify the best diagonal lift there.

A distinct primitive four-speed example is `(1,2,3,5)`. The pair `(2,3)`
is the sole violation (`25<31`); the next pair passes (`25>=23`). The
repaired star covariance has only `S_23=S_32=3/28` off the diagonal, and

```
Delta_max=1488/325.
```

Thus restricting to diagonal ambient lifts genuinely misses the global
maximum even in this small, explicit projected-cube family. The
overcomplete contact set is a star plus one leaf edge, hence contains a
triangle; this is consistent with the base packet's impossibility of a
nonstar **basis-tree** certificate.

The regime is nonempty in every `k>=4`, not just four dimensions. Take
`v=(1,...,1,2,2,4)`, with `k-3` entries equal to one. Only the pair `(2,2)`
violates: it has `A=B=3/4`, `c=1/9`. The mixed `(1,2)` pairs have squared
normal length `17/18`, and `(1,1)` pairs have `25/32`. The repaired chart
determinant is always `80/81`. This repeated-speed family establishes
nonvacuity; the separate four-speed example above has distinct speeds.

## 7. Reproduction and boundaries

Run `python3 -B check_exact.py` and `python3 -B -O check_exact.py` from any
directory using the script's absolute path if necessary. Both outputs
must equal `check_exact.out`. The standard-library companion enumerates
nondecreasing integer speeds in `1..8` for dimensions `k=3..6`, records
how many have zero, one, or more violating pairs, and tests (4)--(14) only
on the explicitly declared single-pair region. All three-speed rows are
covered by the complete piecewise formula. Positive weights, exact
inverse certificates, normalized covariance, determinant identities,
strictly positive triple lifts, diagonal nonrepresentation in higher
dimension, threshold equality and the two concrete comparisons are
checked with rational arithmetic. These finite tests do not replace the
universal proof above.

The dimension-independent family is additionally replayed for `k=4..10`.
Other tight nonviolating pair facets may occur at the boundary of the
single-pair regime; they need no positive weight in (8), and do not affect
the certificate or its uniqueness conclusion.

Multiple violating leaf pairs are outside the theorem: repairs then
interact through positive definiteness and the inverse certificate.
The diagonal obstruction is not an LRC obstruction. No lattice-point
claim, lonely-time witness, computational-prime audit, literature-priority
claim, or improved global LRC speed-product constant follows from this
packet. Larger volume also does not imply containment or uniform
pointwise norm improvement, as the base packet's hostile demonstrates.

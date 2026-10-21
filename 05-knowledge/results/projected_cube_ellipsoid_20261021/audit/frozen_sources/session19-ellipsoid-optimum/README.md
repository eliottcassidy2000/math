# Exact ellipsoid of a dominant-speed projected cube

**Status: PROVED by the elementary arguments below; finite exact companion
provided. Session 19, 2026-10-21. Staged outside the repository for independent
review. No literature-priority claim.**

## Result and scope

Let `k >= 3`, let `v=(a_1,...,a_(k-1),M)` have positive real coordinates,
with `M=max v_i`, and let `A >= B` be the two largest coordinates other than
`M`, counting multiplicity. Write `H=v^perp`, let `P` be orthogonal projection
onto `H`, and put `D=P[-1,1]^k`. If and only if

```
M^2 >= M(A+B)+AB,                                      (1)
```

the following ellipsoid is contained in `D`:

```
q_i^*=(1+a_i/M)^2  (i<k),       q_k^*=0,
Q^*=diag(q_i^*),                 E^*=P(Q^*)^(1/2) B_2^k. (2)
```

When (1) holds, `E^*` is the **unique maximum-volume ellipsoid among all
ellipsoids contained in D, including translated ellipsoids**. Its squared
volume divided by `vol(B_2^(k-1))^2` is

```
Delta^* = M^2 product_(i<k) (1+a_i/M)^2 / ||v||_2^2.      (3)
```

The ambient zero in `Q^*` does not make the projected ellipsoid degenerate.
The induced covariance on the `(k-1)`-dimensional space `H` is positive
definite. For `k>=3`, (2) is also the unique maximizing nonnegative diagonal
lift. If (1) fails, every maximizing nonnegative diagonal lift is strictly
positive; the present result does not identify that interior optimizer.

Moreover, a parallelotope made from exactly `k-1` independent pair facets of
`D` can certify its maximum ellipsoid only when those facets form a star.
No nonstar spanning tree works for any positive speed vector. Thus changing
the facet-basis tree does not extend this certificate beyond (1).

This is a fixed-shape convex-geometric theorem. It is not an LRC proof, does
not improve a uniform speed-product constant, and does not certify any
external computational proof. Its source comparison is local to the stated
ellipsoid construction.

## 1. Inheritance and the corrected carrier

The closest repository mechanism is
`THM-4014-lrc14-diagonal-polar-ellipsoid-fastest-coordinate-relation-compression`
in `01-canon/theorems/`: the exact polar vertices, pair feasibility test, and
the compressed covariance rather than an ambient inverse. Its canonical
hostile uses `v=(1,2), Q=diag(2,1)`: the inverse of the compressed covariance
is `5/9`, whereas compression of the ambient inverse gives `3/5`.

`THM-4009-euclidean-covering-transference-short-relation-compression` supplies
the speed-dependent inradius viewpoint, and
`THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction`
supplies the projected-lattice/dual interpretation. None is needed as an
unproved premise here; the finite-dimensional identities are rederived.

The new operation is degeneration of the ambient diagonal lift while
retaining a nondegenerate projected ellipsoid. The least-used sidecar is the
full labelled facet basis, not its determinant alone. The corrected near
miss is to confuse `Q^(-1)|H` with `(PQP|H)^(-1)` or to discard `q_M=0`
because an ambient inverse is unavailable.

The concept board is: (i) pair polar vertices, (ii) compressed covariance,
(iii) volume/covolume, (iv) facet-tree paths, (v) LRC semantic loss. The
source-to-target map preserves containment and exact Euclidean volume;
it loses lattice-point incidence, translation phase and lonely-time
witnesses. Those require additional arithmetic data and arguments.

## 2. Exact pair criterion, also at zero diagonal entries

All polar bodies in this note are taken inside `H`. For `y in H`,

```
h_D(y)=||y||_1,
h_(E_Q)(y)^2=y^T Q y,
D^o={y in H: ||y||_1 <= 1}.                            (4)
```

The vertices of `D^o` are exactly

```
 +/- r_ij,       r_ij=(v_j e_i-v_i e_j)/(v_i+v_j), i<j. (5)
```

Indeed, no cross-polytope vertex lies on `H`, because all speeds are
positive. Every section vertex lies on an edge; the edges meeting `H`
join coordinate vertices of opposite signs and distinct indices, and their
intersections are (5). Conversely each such intersection is a vertex: it
is the unique solution on its edge, whose intersection with the section is
an exposed face. Equivalently, on a fixed sign orthant the two affine
constraints leave section vertices with support two, and positive-speed
balance forces opposite signs.

The quadratic form in (4) is convex for every nonnegative diagonal `Q`,
including singular `Q`. A convex function on a polytope is bounded above by
its largest vertex value. Therefore

```
E_Q subset D
  iff q_i v_j^2 + q_j v_i^2 <= (v_i+v_j)^2 for all i<j. (6)
```

By polarity and the bipolar theorem, the same vertex description gives the
exact primal facet presentation

```
D={x in H: |v_j x_i-v_i x_j| <= v_i+v_j for all i<j}.   (7)
```

No integrality, distinctness, primitive normalization, or LRC assumption is
used in (4)--(7).

An independent direct proof of (7) avoids polarity: `x in H` belongs to
`D` exactly when some `s` has `x+s v in [-1,1]^k`. This means that the
closed intervals `[(-1-x_i)/v_i,(1-x_i)/v_i]` share a point. Finitely many
intervals have a common point exactly when each pair does, and their
pairwise compatibility is precisely (7).

## 3. The dominant-speed condition is exact

At (2), every pair involving `M` saturates (6). For `i,j<k`, direct expansion
gives the slack

```
(a_i+a_j)^2 - q_i^* a_j^2 - q_j^* a_i^2
 = (2 a_i a_j/M^2) [M^2-M(a_i+a_j)-a_i a_j].            (8)
```

Since `M(a+b)+ab` is strictly increasing in each positive argument, all
these slacks are nonnegative exactly when the largest pair `A,B` satisfies
(1). An equivalent normalization is

```
(1+A/M)(1+B/M) <= 2.                                  (9)
```

Equality is allowed and just adds contact facets. The condition implies
`M>A+B`, hence the maximum speed is unique. Repeated smaller speeds cause
no difficulty.

## 4. Global maximum, by a parallelotope certificate

For `i<k`, set

```
r_i=(M e_i-a_i e_k)/(M+a_i),
L:H -> R^(k-1),        (Lx)_i=r_i dot x,
C={x in H: |r_i dot x|<=1 for all i<k}.                 (10)
```

Each `r_i` belongs to `H`, and these `k-1` rows are linearly independent
(look at coordinates `1,...,k-1`). Thus `L` is an isomorphism, `L(C)` is the
unit cube, and (7) gives `D subset C`. Also `LP=L`. At (2),

```
L P (Q^*)^(1/2) = [I_(k-1) | 0],
L(E^*)=B_2^(k-1).                                    (11)
```

This identity simultaneously proves nondegeneracy on `H` and identifies
the candidate as the maximal ellipsoid of `C`.

For completeness, let `c+T B_2^d` be any ellipsoid inside `[-1,1]^d`,
where `d=k-1`. Coordinate support gives

```
|c_i| + ||row_i(T)||_2 <= 1.
```

Hadamard's determinant inequality yields

```
|det T| <= product_i ||row_i(T)||_2
        <= product_i (1-|c_i|) <= 1.                  (12)
```

Equality requires `c=0` and orthonormal rows, so the ellipsoid is exactly
`B_2^d`. Singular ellipsoids have zero volume and cannot be maximizers.
Under (1), (8) proves `E^* subset D subset C`; transport (12) by `L` to
obtain the asserted global maximum and uniqueness. This proves the result
without invoking a characterization theorem for John's ellipsoid.

## 5. Determinants and the independent diagonal proof

For any nonnegative diagonal `Q`, choose an orthonormal matrix `U` with
image `H`. Cauchy--Binet gives, including zero entries,

```
det(U^T Q U)
 = sum_i v_i^2 product_(j != i) q_j / ||v||_2^2.        (13)
```

The squared `(k-1)`-minors of `U` equal `v_i^2/||v||_2^2`, since its unit
normal is `v/||v||_2`. Unlike an inverse formula, (13) is a polynomial
identity and remains valid on the boundary. Substitution of (2) gives (3).

There is a second proof restricted to nonnegative diagonal lifts. Fix
`t=q_k>=0` and define

```
A_i=(1+a_i/M)^2,   b_i=a_i^2/M^2,   x_i=b_i/A_i.
```

The `M`-pair constraints imply `0<=q_i<=A_i-b_i t` and
`0<=t<=T=min_i 1/x_i`. The numerator in (13) is coordinatewise
nondecreasing on the nonnegative orthant, so replacing every `q_i`, `i<k`,
by its cap gives an upper bound even if the cap vector violates another
pair constraint. With `f(t)=product_(i<k)(1-x_i t)`, this bound is

```
N_cap(t)=M^2 product_i A_i [f(t)-t f'(t)],
N_cap'(t)=-M^2 product_i A_i t f''(t),
f''(t)=2 sum_(i<j) x_i x_j product_(h != i,j)(1-x_h t). (14)
```

On `[0,T]`, every term of `f''` is nonnegative. Since `k-1>=2`, it is
strictly positive on `[0,T)`, including the empty-product case `k=3`.
Consequently `N_cap(t)<N_cap(0)` for every `t>0`, including `t=T` by
continuity and integration on an initial positive interval. When `t=0`,
the numerator is `M^2 product_i q_i`, uniquely maximized at all the caps.
This proves uniqueness within the diagonal cone and checks the singular
endpoint separately from any denominator cancellation.

For `k=2`, the feasible set is `M^2 q_1+a_1^2 q_2<=(M+a_1)^2`.
Every point on its nonnegative boundary segment maximizes the determinant;
all yield the same interval `D`. The ellipsoid is unique, but its diagonal
lift is not. For `k>=3`, diagonal covariance restriction is injective:
if a diagonal difference `d` vanishes on `H`, the pair vectors imply
`d_i/v_i^2+d_j/v_j^2=0` for every pair, and any triangle forces every entry
to vanish.

## 6. Exact boundary classification of diagonal maximizers

The diagonal feasible set in (6) is closed and bounded. It contains `Q=I`,
whose compressed determinant is one. A maximizing diagonal lift therefore
has positive determinant and at most one zero coordinate.

Suppose a maximum has `q_p=0`; all other `q_j` are positive. If a pivot
constraint for the pair `p,i` is slack, consider the direction

```
delta q_p=1,
delta q_i= v_i^2/v_p^2,
delta q_j=-v_j^2/v_p^2  for j notin {p,i}.              (15)
```

For sufficiently small positive step all coordinates remain nonnegative.
Other pivot constraints and all nonpivot pairs involving `i` are unchanged;
pairs involving neither `p` nor `i` decrease. The single increasing pivot
constraint starts with slack. Thus this is a feasible direction. The
directional derivative of the numerator of (13) is

```
2 v_i^2 product_(j notin {p,i}) q_j > 0,                (16)
```

a contradiction. Every pivot constraint must therefore be tight, forcing
`q_j=(1+v_j/v_p)^2`. Its feasibility is precisely (8), with `p` as pivot;
this forces `v_p` to be the unique largest speed and gives (1).
Conversely, (1) supplies the zero-coordinate maximum already proved.
Thus the dominant-speed condition is equivalent to the existence of any
boundary maximizing nonnegative diagonal lift. Outside it, every such
maximum is strictly interior. This is not a description of the full
maximum ellipsoid outside the cone.

## 7. Why no other facet-basis tree can extend the certificate

Put `r_i=1/v_i` and `y_i=x_i/v_i`. Equation (7) becomes

```
|y_i-y_j| <= r_i+r_j.                                  (17)
```

For a graph on the `k` labels, its selected pair functionals have rank
`k-c`, where `c` is its number of connected components. This is the usual
rank of the incidence matrix, unchanged by the invertible diagonal change
of coordinates. Restricting to `H` loses no rank: the common null direction
before restriction is the speed direction, transverse to `H`. More
generally its kernel has dimension `c`, and intersection with `H` has
dimension `c-1`, giving the same rank. In particular, `k-1` selected pair
functionals form a basis exactly when the graph is a spanning tree.

Given a spanning tree `T`, orient its edges arbitrarily and use coordinates

```
z_uv=(y_u-y_v)/(r_u+r_v),   uv in T.
```

The tree's enclosing parallelotope is the cube in these coordinates, and
its unique maximum ellipsoid is the unit ball. For any pair `i,j`, the
path telescoping identity expresses `(y_i-y_j)/(r_i+r_j)` as a signed sum
of edge coordinates with coefficient magnitudes `(r_u+r_v)/(r_i+r_j)`.
Hence this ellipsoid lies in `D` exactly when

```
sum_(uv in path_T(i,j)) (r_u+r_v)^2 <= (r_i+r_j)^2
for every pair i,j.                                   (18)
```

On a length-two path `a-b-c`, (18) reduces to

```
r_b^2+r_b(r_a+r_c) <= r_a r_c,
```

which implies `r_b<min(r_a,r_c)`: if, for example, `r_b>=r_a`, its left
side is strictly greater than `r_a r_c`. A length-three path `a-b-c-d`
would require both `r_b<r_c` and `r_c<r_b`, impossible. Every nonstar tree
has a length-three path. Therefore only stars can satisfy (18). For a star
with centre `p`, (18) is exactly the dominant-speed condition with pivot
`v_p`. This classifies all pair-facet basis certificates, not all possible
John contact configurations. Outside the cone a certificate based on pair
facet contacts cannot consist of only `k-1` independent facet pairs whose
enclosing parallelotope has the same ellipsoid.

## 8. Exact source comparison and normalization

Jaan Allikvere, *Fourteen and fifteen lonely runners*,
[arXiv:2609.02604v2](https://arxiv.org/html/2609.02604v2#S3.SS1),
Section 3.1, Lemmas 3.1--3.2, uses the same `D,E,H,P` normalization, with
`q_i=1+2u_i-u_i^2`, `u_i=v_i/M`, and gives its covariance/covolume formula.
The source's family is the cap family in (14) at `t=2`, whereas (2) is
`t=0`. Thus (14) proves strict volume improvement over that construction
for every `k>=3` in the dominant-speed cone. This does not contradict its
inclusion lemma, which makes no optimality claim.

For `v=(1,2,M)`, condition (1) is `M^2>=3M+2`; among integer `M>=2`, it
holds exactly for `M>=4`. Formula (3) becomes

```
Delta^*=(M+1)^2(M+2)^2/[M^2(M^2+5)].                   (19)
```

For the distinct primitive example `(1,2,4)`:

```
Q^* diagonal          = (25/16, 9/4, 0),
source Q diagonal     = (23/16, 7/4, 2),
Delta^*               = 75/28,
Delta_source          = 221/84,
Delta^*-Delta_source  = 1/21,
volume ratio          = 15/sqrt(221) > 1.              (20)
```

For primitive integer `v`, `covol(P Z^k)=1/||v||_2`. In the norm with
unit ball `E^*`, its covolume is therefore

```
1/[M product_(i<k)(1+a_i/M)].                          (21)
```

At `(1,2,4)` this is `2/15`, against `2/sqrt(221)` for the source's
ellipsoid. The source defines its `LRC(k)` with `k` nonzero speeds, hence
`k+1` runners; that convention must not be silently conflated with a
runner-count index. The current preprint's much broader LRC claims are not
inputs to this proof and are not audited here. No claim is made that a
pointwise shape improvement changes its optimized global constants.

The gain is in volume, **not containment or uniform norm dominance**. In
the star coordinates (10), the source covariance for `(1,2,4)` is
`[[1,2/15],[2/15,1]]`, whereas the optimum covariance is the identity.
The directions `(1,1)` and `(1,-1)` give squared supports respectively
`34/15` and `26/15`, compared with `2` for the optimum. Neither ellipsoid
contains the other. In particular, one cannot claim that every weighted
short-relation inequality improves merely from the determinant gain.

## 9. Hostiles, exact companion, and remaining frontier

The example `(1,2,3)` violates (1): its small-pair polar vertex evaluates
to `89/81>1` at (2). A repeated-speed example `(1,1,1)` also fails. Equality
is nonvacuous over positive rationals: `(3,6,12)` is strictly inside, while
`(2,3,6)` satisfies `36=30+6` and has an extra contact pair. The diagonal
lift remains unique there. The `k=2` continuum is an explicit hostile to
unqualified uniqueness of the lift. All tree orientations give the same
squared path criterion. A triangle plus an isolated vertex has rank two,
not three, and is not a facet basis.

Run from any directory, using Python 3 and its standard library:

```
python3 -B /path/to/check_exact.py
python3 -B -O /path/to/check_exact.py
```

Both streams must equal `check_exact.out`. The companion checks rational
determinant identities through an independent nonorthonormal-basis formula,
the exact containment threshold, cube coordinates and contact identities,
the cap-polynomial identity and derivative, translated cube controls,
boundary improvement directions, the source comparison, and all labelled
four-vertex spanning trees on a declared positive-integer box. Finite
checks support the written proof; they are not substitutes for its
universal quantifiers. `MANIFEST.json` freezes owned artifact hashes and
reproduction commands. No repository or remote state is changed.

The remaining geometric question is the maximum ellipsoid outside the
dominant-speed cone, where this facet-basis certificate is impossible.
The independent LRC obligations are unchanged: volume is not lattice
incidence, contact geometry is not a lonely-time witness, and a local
ellipsoid improvement is not a uniform speed-product improvement.

Method cards used: `Attack a proposed bound before extending it`,
`Test structured adversaries, not only random samples`, and `Verify the
consequence, not the model's own assumptions`, in
`00-navigation/META-PATTERNS.md`. The singular projected shape and the
nonstar-tree obstruction are the decisive structured tests. The present
packet does not promote a new universal method card.

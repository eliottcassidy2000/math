---
status: PROVED ALGEBRA + VERIFIED-EXACT independent audit; not canonized; LRC(14) OPEN
source: lift_tower_covering independent audit / root LRC(14) session, 2026-08-24
tags:
  - lonely-runner
  - lrc14
  - geometry-of-numbers
  - polarity
  - ellipsoid
  - weighted-relation
  - transference
---

# Independent audit: diagonal ellipsoids refine the short-relation reduction

## Decisive verdict

The proposed refinement survives, but one displayed polarity step must be
repaired.

Let `V=n^perp`, let `P` be orthogonal projection to `V`, and put

```text
K=P([-1,1]^13).
```

For a positive diagonal matrix `H=diag(h_i)`, define the compressed operator

```text
A=P H|_V : V -> V.                                      (1)
```

The polar of `{u in V : u^T H u<=1}` is governed by `A^(-1)`, not by the
ambient matrix `H^(-1)`.  Generically `H` does not preserve `V`.  Therefore
the correct inellipsoid is

```text
Q_H^o={x in V : <A^(-1)x,x><=1}.                        (2)
```

The smaller ambient set `{x in V:x^T H^(-1)x<=1}` is still contained in `K`,
but using only that smaller set does not imply the desired `H`-weighted dual
bound.  Retaining `(2)` repairs the proof completely: a hypothetical LRC(14)
counterexample forces a nonzero integer Graver relation

```text
a dot n=0,                   a^T H a<196.                (3)
```

This is a genuine refinement of
[THM-4009](../../01-canon/theorems/THM-4009-euclidean-covering-transference-short-relation-compression.md),
but this audit does not alter or promote the theorem.  LRC(14) remains open.

The companion exact audit is
[`04-computation/lrc14_diagonal_ellipsoid_refinement_independent_audit_20260824.py`](../../04-computation/lrc14_diagonal_ellipsoid_refinement_independent_audit_20260824.py),
with frozen output
[`05-knowledge/results/lrc14_diagonal_ellipsoid_refinement_independent_audit_20260824.out`](lrc14_diagonal_ellipsoid_refinement_independent_audit_20260824.out).

## 1. Polar vertices and the exact pair criterion

For `u in V`,

```text
h_K(u)=sup_(z in [-1,1]^13) <u,Pz>
      =sup_z <u,z>
      =||u||_1.
```

Hence, with polarity taken inside `V`,

```text
K^o={u in V:||u||_1<=1}.                                (4)
```

The ambient `l1` ball is the cross-polytope with vertices `+-e_i`.  None of
those vertices lies in `V`, because every `n_i` is positive.  A vertex of its
hyperplane section must therefore be the crossing of `V` with an edge joining
`e_i` to `-e_j`, or with the opposite edge.  Thus the vertices of `(4)` are
exactly

```text
v_ij^+-=(+/-)(n_j e_i-n_i e_j)/(n_i+n_j),    i<j.       (5)
```

There are `2*C(13,2)=156` of them.  Each has `l1` norm one and is orthogonal
to `n`.

Assume `H` is positive diagonal and write `delta_i=h_i-1`.  At `(5)`,

```text
(v_ij)^T H v_ij=(h_i n_j^2+h_j n_i^2)/(n_i+n_j)^2.
```

Consequently this value is at most one exactly when

```text
delta_i n_j^2+delta_j n_i^2<=2n_i n_j.                  (6)
```

The quadratic form is convex.  Every point of the polytope `(4)` is a convex
combination of its vertices, so bounds at all vertices bound the whole
polytope.  Conversely, containment bounds every vertex.  Therefore `(6)` for
every pair is **necessary and sufficient** for

```text
K^o subset {u in V:u^T H u<=1}.                          (7)
```

Positivity of `H` is load-bearing; pair inequalities alone should not be used
to define an ellipsoid if arbitrary negative diagonal entries are allowed.

## 2. The ambient-inverse trap and its exact repair

On `V`, the quadratic form in `(7)` is represented by `A=P H|_V`.  Its polar
is `(2)`.  One can see the difference without choosing a basis.  For `x in V`,
solve `Ay=x`.  Since `P(Hy)=x`, there is a scalar `lambda` with

```text
Hy=x+lambda n.
```

The constraint `y in V` gives

```text
lambda=-(n^T H^(-1)x)/(n^T H^(-1)n),
```

and hence

```text
<A^(-1)x,x>
 =x^T H^(-1)x-(n^T H^(-1)x)^2/(n^T H^(-1)n).           (8)
```

The second term is nonnegative.  Thus the originally displayed ambient
ellipsoid is a valid but generally strict subellipsoid:

```text
{x in V:x^T H^(-1)x<=1} subset Q_H^o subset K.          (9)
```

It is not the polar ellipsoid.  An exact two-coordinate hostile makes the
loss visible.  Take

```text
n=(1,2),       V=span((2,-1)),       H=diag(2,1).
```

The sole pair inequality is an equality.  On a unit vector in `V`, the
restricted `H` operator is `9/5`, so its inverse is `5/9`; the ambient
`H^(-1)` compression is instead `3/5`.  At squared radial coordinate `7/4`,

```text
(7/4)(5/9)=35/36<1,       but       (7/4)(3/5)=21/20>1. (10)
```

The point is in the true polar ellipsoid and outside the ambient-inverse one.
This refutes equality, not the surviving inclusion `(9)`.

Moreover, starting from only the smaller ellipsoid in `(9)` would produce a
dual norm associated with `(P H^(-1)|_V)^(-1)`.  Since

```text
(P H^(-1)|_V)^(-1) <= A,
```

a bound for that smaller quadratic form does **not** imply `a^T H a<196`.
The compressed operator cannot be discarded and recovered later.

## 3. Correct transformation, polarity, and strictness

Let

```text
Q_H={u in V:<Au,u><=1}.
```

Equations `(7)` and polarity reversal give

```text
Q_H^o={x:<A^(-1)x,x><=1} subset K.                      (11)
```

As in THM-4009, write the lonely zonotope as

```text
Z(n)=c+(3/7)K,              c=(1/2)P(1).
```

For a counterexample, the **closed** zonotope is disjoint from
`Lambda=P(Z^13)`.  Therefore the closed ellipsoid

```text
c+(3/7)Q_H^o
```

is lattice-point-free.  Apply `A^(-1/2)` on `V`.  It sends this ellipsoid to
the closed Euclidean ball of radius `3/7` about `A^(-1/2)c`, and sends the
lattice to

```text
L=A^(-1/2)Lambda.
```

Closed-ball disjointness and discreteness give the strict covering-radius
inequality

```text
mu(L)>3/7.                                               (12)
```

The transformed dual is

```text
L*=A^(1/2)Lambda*.                                      (13)
```

Using Banaszczyk's inherited twelve-dimensional Euclidean inequality
`2 mu(L)lambda_1(L*)<=12`, equations `(12)`--`(13)` give some nonzero
`a in Lambda*` with

```text
||A^(1/2)a||_2<14.
```

THM-4009's exact dual identity says
`Lambda*=Z^13 intersect n^perp`.  Finally, since `a in V`,

```text
||A^(1/2)a||_2^2=<a,P(Ha)>=a^T H a<196,                 (14)
```

which proves `(3)`.  Every strict sign is accounted for: the closed
lattice-free ellipsoid gives `mu>3/7`, which gives `lambda_1<14`, which gives
the strict weighted square threshold `<196`.

Choose `a` shortest for the `H`-norm.  If it had a nontrivial conformal
decomposition `a=u+v` by integer kernel vectors, then
`|u_i|<=|a_i|` in every coordinate, strictly somewhere.  Positivity and
diagonality of `H` would give `u^T H u<a^T H a`, a contradiction.  Thus the
weighted-shortest row can be chosen Graver.

## 4. The max-coordinate boost

For a chosen coordinate `i`, put

```text
M_i=max_(j!=i)n_j,
h_i=1+2n_i/M_i,                 h_j=1 for j!=i.          (15)
```

For a pair containing `i`, condition `(6)` reduces to `n_j<=M_i`; for every
other pair it is trivial.  Hence `(15)` is feasible, with equality at a
largest other speed.

Let `M` be the largest speed and `M_2` the second-largest.  Applying `(14)` to
the max-coordinate boost gives a Graver relation with ordinary square sum
`S=sum a_i^2` satisfying

```text
S+2(M/M_2)a_M^2<196,
M_2 S+2M a_M^2<=196M_2-1.                              (16)
```

Immediate exact consequences are:

- `S<=195`;
- if `a_M!=0`, then `S<194`, hence `S<=193`;
- `(M_2+2M)a_M^2<=196M_2-1`, hence the raw bound
  `|a_M|<=8`;
- `||a||_1<=49`.

For the last line, weighted Cauchy gives

```text
||a||_1^2
 <=(12+1/(1+2M/M_2)) a^T H a
 <(37/3)196<50^2.                                      (17)
```

The integer optimization is sharper about the boundary.  At fixed
`L=||a||_1` and `t=|a_M|`, the other twelve squares are minimized by balancing
`L-t` among twelve integers.  Replacing the actual max weight by its strict
lower limit `3`, the exact minima are

```text
L=49:  195, uniquely at t=1 and twelve absolute coefficients 4;
L=50:  204.                                             (18)
```

Thus `L<=49` with room to spare at `50`.  Equality `L=49` can occur under the
energy bound only in the rigid shape

```text
|a_M|=1,              |a_j|=4 for all j!=M.             (19)
```

It additionally requires `M/M_2<3/2`, because its energy is
`192+1+2M/M_2`, and the relation equation forces `M=0 mod 4`.  Therefore

```text
M/M_2>=3/2 or 4 not|M   =>   ||a||_1<=48               (20)
```

for the relation selected by the max-coordinate boost.

The kernel equation improves the raw coefficient cap.  Put

```text
T=sum_(j!=M)n_j^2,                 r=M/M_2.
```

If `a_M!=0`, Cauchy applied to
`M a_M=-sum_(j!=M)n_j a_j` gives

```text
a^T H a
 >=a_M^2(1+2r+M^2/T)
 >=a_M^2(1+2r+r^2/12)
 >a_M^2(37/12).                                      (21)
```

Since `64*(37/12)>196`, the genuinely relation-aware cap is

```text
|a_M|<=7.                                              (22)
```

This is stronger than the requested immediate `<=8` bound and uses the
kernel balance, not ellipsoid containment alone.

Equation `(21)` also says `a_M=0` whenever

```text
1+2M/M_2+M^2/T>=196.                                  (23)
```

Using only `T<=12M_2^2` gives the simpler sufficient threshold

```text
M/M_2>=-12+2sqrt(621)=37.8397431775... .               (24)
```

This is calibration only.  By
[THM-1008](../../01-canon/theorems/THM-1008-lrc13-descent-floor.md),
a live hypothetical counterexample already has `M/M_2<13`; the regime in
`(24)` is closed independently.  More generally `(23)` cannot fire below
`13`, because `T>=M_2^2` makes its left side at most `(1+M/M_2)^2<196` there.
No LRC progress should be claimed from the avoidance threshold.  The live
consequences are the weighted relation, `(16)`, `(17)`, `(20)`, and `(22)`.

The boosts for different indices have the quantifiers

```text
for every i, there exists a weighted-short Graver relation a^(i).
```

They do not supply one relation satisfying all thirteen boosted bounds, nor
thirteen distinct or independent relations.  Indeed, simultaneously spending
the full pair budget on two boosted coordinates violates `(6)` on their pair.

## 5. The simultaneous linear bank

Let

```text
D=M+M_2,                   h_i=1+2n_i/D for every i.    (25)
```

For every pair,

```text
(h_i-1)n_j^2+(h_j-1)n_i^2
 =2n_i n_j(n_i+n_j)/D
 <=2n_i n_j,                                           (26)
```

because `D` is the maximum pair sum.  Thus `(25)` is feasible, with equality
only on the top pair when speeds are distinct.  Unlike the coordinate boosts,
it yields one relation carrying all weights simultaneously:

```text
D S+2 sum_i n_i a_i^2<=196D-1.                         (27)
```

In particular,

```text
|a_i|<=floor(sqrt((196D-1)/(D+2n_i))),
|a_M|<=9.                                              (28)
```

Weighted Cauchy again gives `||a||_1<=49`.  To see the universal constant,
the eleven non-top reciprocal weights are less than one, while, for
`r=M/M_2>1`, the top two sum to

```text
(r+1)/(3r+1)+(r+1)/(r+3)
 =4(r+1)^2/((3r+1)(r+3))<4/3.                          (29)
```

Thus `sum_i 1/h_i<11+4/3=37/3`, as in `(17)`.

This bank does **not** improve the scalar square cap already in THM-4009.  If
`m` is the smallest speed, `(27)` only gives

```text
S<196D/(D+2m).                                         (30)
```

Writing `r=M/m` and `s=M_2/m`, the difference between the multiplier in
`(30)` and THM-4009's exact-inradius multiplier is

```text
(r+s)/(r+s+2)-(r^2+1)/(r+1)^2
 =2(rs-1)/((r+s+2)(r+1)^2)>0.                          (31)
```

So `(30)` is strictly weaker.  The new information is the localization term
`sum n_i a_i^2`, the max-coordinate cap, and the improved `l1` cap.

As a consistency check, the best uniform diagonal increment is

```text
delta=2mM/(m^2+M^2).
```

It saturates `(6)` on the extreme pair and gives
`H=((m+M)^2/(m^2+M^2))I`, exactly recovering THM-4009's speed-dependent
Euclidean inradius.  The diagonal method extends that geometry rather than
contradicting it.

There is also a sharp support-two refinement.  Suppose the selected relation
is supported on speeds `s<t`, with `g=gcd(s,t)`, `p=s/g`, and `q=t/g`.  Its
coefficients are `(q,-p)` up to sign, and `(27)` becomes

```text
p^2+q^2+2((s+t)/D)pq<196.                              (32)
```

If the support is the top pair, `s+t=D`, so `(32)` is

```text
(p+q)^2<196,                  p+q<=13.                  (33)
```

There are exactly `28` coprime pairs `p<q` satisfying `(33)`, versus the `47`
ambient support-two ratios under THM-4009's `p^2+q^2<=195`.  This is a branch
refinement, not a claim that the transference-selected relation must use the
top pair or have support two.

## 6. Typed map and loss ledger

```text
source:     closed lattice-free lonely zonotope in V=n^perp
target:     H-weighted shortest integer relation
map:        polar pair facets -> compressed ellipsoid -> A^(-1/2) lattice
            -> Euclidean covering transference -> A^(1/2) dual vector
preserved:  counterexample implies one Graver a with a^T H a<196
destroyed:  endpoint owner, phase, component, odd centre character, support
sidecars:   H choice; compressed operator A; THM-1008 ratio range;
            quantifier identifying which boost selected which relation
hostile:    n=(1,2), H=(2,1) separates A^(-1) from ambient H^(-1)
boundary:   coordinate boosts cannot be combined at full strength;
            linear weights recover simultaneity but weaken scalar norm
```

The weighted relation remains necessary only.  The AP speed row has many very
short relations and is lonely, so no diagonal norm threshold by itself proves
LRC(14).

## 7. Exact controls and reproduction

The standard-library companion verifies:

- all `156` signed polar vertices for `n=(1,...,13)`;
- all `78` pair identities for the max boost and linear bank, including their
  unique equality facets;
- the compressed inverse equation and rank-one loss exactly over `Q`;
- the exact hostile `(10)`;
- the integer optimizations `(18)` and the relation-aware coefficient check;
- the `28` top-pair ratios; and
- the scalar comparison identity `(31)`.

Run in normal and optimized modes:

```bash
python3 -B 04-computation/lrc14_diagonal_ellipsoid_refinement_independent_audit_20260824.py
python3 -B -O 04-computation/lrc14_diagonal_ellipsoid_refinement_independent_audit_20260824.py
```

Both runs must byte-match the frozen output.  The only external mathematical
input is the Banaszczyk inequality already cited and scoped by THM-4009; this
audit makes no new literature claim.

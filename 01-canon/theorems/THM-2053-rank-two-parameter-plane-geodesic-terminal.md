---
id: THM-2053
title: "Rank-two parameter planes have an explicit large-direction LRC terminal"
status: >
  PROVED from the standing cited LRC through thirteen total runners and
  Ungar's 1982 planar direction theorem. The torus-geodesic estimate is
  elementary and unconditional. Every rational two-plane containing a
  positive row has a full-support repeat projection, hence torus margin at
  least 1/13; every primitive parameter direction of norm at least 91 times
  the plane Lipschitz constant is LRC14-safe. This makes every THM-2052 plane
  finite in parameter space, but does not enumerate the disks or prove
  LRC(14).
source: codex-2026-07-21-DC2-LRC14-termination
depends_on:
  - THM-2052
  - LRCUpTo13
  - "P. Ungar, 2N noncollinear points determine at least 2N directions, JCTA 33 (1982), 343-347"
related:
  - HYP-4342
  - HYP-4346
  - HYP-8846
---

# THM-2053 -- large directions in every relation plane

## 1. Parameter planes and their torus

Let `U` be a two-dimensional rational subspace of `Q^n`, and choose a
`Z`-basis `u,z` of the saturated lattice `U intersect Z^n`. Write

```text
c_i=(u_i,z_i) in Z^2,
F_U(x,y)=min_i ||u_i x+z_i y||,
M_T(U)=max_((x,y) in (R/Z)^2) F_U(x,y),                 (1)
L(U)=max_i sqrt(u_i^2+z_i^2).
```

For a primitive nonzero parameter direction `d=(a,b) in Z^2`, put

```text
v(d)_i=a u_i+b z_i,
M(d)=max_(t in R/Z) min_i ||v(d)_i t||.                 (2)
```

Thus `M(d)` is the ordinary one-dimensional lonely-runner maximin of the row
selected by `d`, while `M_T(U)` permits the two parameters to move
independently.

## 2. General primitive-geodesic rate

For every primitive `d=(a,b)`,

```text
M_T(U)-L(U)/(2 sqrt(a^2+b^2)) <= M(d) <= M_T(U).        (3)
```

This is the all-direction version of HYP-4342's `(1,N)` rate.

### Proof

Let

```text
C_d={(at,bt) mod Z^2:t in R/Z}.
```

The right inequality in (3) is immediate because `C_d` is a subset of the
two-torus and `F_U(at,bt)` is exactly the objective in (2).

For the left inequality, put `n=(b,-a)`. Since `gcd(a,b)=1`, the subcircle
`C_d` is the kernel of the primitive character

```text
(x,y) |-> bx-ay mod Z.                                  (4)
```

Given any torus point `X`, choose a lift to `R^2` and an integer `m` with

```text
|n dot X-m|<=1/2.
```

Moving `X` in the normal direction by

```text
-(n dot X-m)n/||n||_2^2
```

lands in the kernel (4). Hence every torus point is within Euclidean distance

```text
1/(2||n||_2)=1/(2 sqrt(a^2+b^2))                       (5)
```

of `C_d`.

The circle-distance function is one-Lipschitz, so each character
`(x,y)|->||u_i x+z_i y||` is `||c_i||_2`-Lipschitz. Their pointwise minimum
`F_U` is therefore `L(U)`-Lipschitz. Apply (5) to a maximizer of `F_U`; its
nearby point on `C_d` loses at most the error in (3). QED.

The proof uses only elementary torus geometry. It does not invoke LRC and it
does not assume that the coordinates `v(d)_i` are positive or distinct.

## 3. A checkable generic floor

Call `U` **full-support-repeat exposed** if some primitive direction
`d_0=(a_0,b_0)` has

```text
p_i=a_0u_i+b_0z_i !=0 for every i,
|p_i|=|p_j| for some i!=j.                              (6)
```

For `n=13`, the standing cited LRC through thirteen total runners implies

```text
M_T(U)>=1/13.                                           (7)
```

Indeed, the thirteen nonzero integers in (6) have at most twelve distinct
absolute values. Apply the cited theorem to that set of at most twelve speeds.
It supplies a time `t_0` at which every distinct speed has distance at least
`1/(k+1)>=1/13`. Repetitions and sign changes do not change circle distance,
so

```text
F_U(a_0t_0,b_0t_0)>=1/13,
```

which proves (7).

Combining (3) and (7) gives the explicit terminal

```text
sqrt(a^2+b^2)>=91 L(U)  ==>  M(d)>=1/14,               (8)
```

because

```text
1/13-1/14=1/182
```

and `L/(2||d||)<=1/182` under (8). Equality in (8) is allowed: LRC needs a
weak `1/14` witness.

## 4. Every genuine two-plane has the repeat direction

The word "generic" can be removed. Assume `U` contains a row with every
coordinate positive. Then no column `c_i` is zero. Since `u,z` are independent,
the columns span `R^2`.

Form the finite centrally symmetric point set

```text
P={+c_i,-c_i:1<=i<=n} subset R^2,                      (9)
```

discarding repetitions, and let `N=|P|`. It is non-collinear, contains no
zero, and `N` is even. Let `D` be the set of radial directions spanned by the
columns. Every radial direction contributes at least the two distinct points
`p,-p`, so

```text
|D|<=N/2.                                              (10)
```

Ungar's planar direction theorem says that `N=2m` non-collinear points
determine at least `N` distinct secant directions.

Suppose, toward a contradiction, that (6) fails. Take distinct `p,q in P`. If
their secant direction `span(p-q)` is not in `D`, choose a primitive integer
direction `d_0` perpendicular to `p-q`. (All points are integral, so such a
direction is integral.) Then

```text
d_0 dot p=d_0 dot q,
```

and this common equality gives two equal absolute projected speeds. Moreover
`d_0 dot c_k` is nonzero for every `k`, because a zero would make `c_k`
parallel to `p-q`, putting the secant direction in `D`. Thus `d_0` would
satisfy (6), a contradiction.

Therefore failure of (6) would force **every** secant direction of `P` to
belong to `D`. The number of directions determined by `P` would be at most
`|D|<=N/2`, contradicting Ungar's lower bound `N`. Hence (6) always holds.

Reference for the exact external input: P. Ungar, *2N noncollinear points
determine at least 2N directions*, Journal of Combinatorial Theory, Series A
33 (1982), 343--347, DOI `10.1016/0097-3165(82)90045-0`.

## 5. Consequence for the THM-2052 atlas

THM-2052 places every hypothetical primitive counterexample in the kernel of
an eleven-dimensional bounded relation code. Each such kernel has dimension
at most two, and only finitely many kernels occur.

For a one-dimensional kernel, its saturated integer lattice has only one
primitive generator up to sign, so that atlas cell is already a single exact
row. For a two-dimensional kernel, Section 4 supplies the repeat direction;
choose a saturated basis. A primitive positive speed row has a primitive parameter
direction: if `(a,b)=g(a',b')` with `g>1`, every speed is divisible by `g`.
Equation (8) leaves only the finite disk

```text
a^2+b^2 < (91L(U))^2.                                  (11)
```

Every row in (11) can be decided exactly by the pair-sum ruler theorem.
Therefore the entire THM-2052 atlas has an explicit finite parameter-space
terminal **without harvesting a twelfth bounded relation**.

What remains is computational/structural rather than infinitary: construct or
compress the finite atlas, choose controlled saturated bases, and discharge
the finite disks (11). This theorem does not claim those disks have been
enumerated and does not prove LRC(14).

## 6. Quantifier repair and assumption challenge

HYP-4346 says that if infinitely many independent directions share finitely
many scale-only templates, then the plane collapses. Applied naively, it only
produces *some* escaping direction and has the wrong quantifier for a specified
target row. Equations (3) and (8) are the repair: the repeat direction plus
the settled lower-dimensional theorem raises the two-torus to `1/13`, and
then every sufficiently long **specified** target direction inherits the
`1/14` bound with an explicit loss.

This also challenges the earlier assumption that the only terminal above
rank eleven is a twelfth bounded relation. Rank twelve is still a valid
maximal-minor terminal, but the torus-geodesic disk is already a terminal at
rank eleven. The next decisive problem is to make (11) small enough to check:
use lattice-reduced saturated bases and exploit the bounded triple-code shape,
rather than search for another local scalar invariant.

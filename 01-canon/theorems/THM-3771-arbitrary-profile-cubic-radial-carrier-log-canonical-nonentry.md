---
id: THM-3771
title: "Arbitrary-profile cubic radial-carrier log-canonical nonentry"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT HOSTILE AUDIT.
  For every nonzero radial profile u and every nonzero polynomial phi of
  degree at most two, the family U=Xu(XT), W=U+3XT+r, Q=U phi(W) has a
  birational log-canonical chart and an explicit complete rational-mate
  torsor.  The candidate below proves an exact smoothness criterion and
  claims that a polynomial mate exists precisely on the boundary where
  both u and phi are constant.  It has no canonical force until the proof,
  boundary atlas, exact companion, and metadata are independently audited.
source: root + jc_sparse_direct_search / 2026-08-23
audit: PENDING INDEPENDENT HOSTILE AUDIT.
depends_on: []
related:
  - THM-3551-one-ray-planar-jacobian-mate-no-go
  - THM-3598-danielewski-rational-exact-polar-graph-family-and-classification
  - THM-3755-composite-monomial-generic-fibre-residue-obstruction
  - THM-3758-quadratic-radial-carrier-rational-exact-split-fibre-nonentry
  - THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate
script: 04-computation/jc2_arbitrary_profile_cubic_radial_carrier_thm3771.py
output: 05-knowledge/results/jc2_arbitrary_profile_cubic_radial_carrier_thm3771.out
script_sha256: 80ca106c48dac4dbc25d0fa881c8afaee6c846f9f4cf7856208bd8902f29093f
output_sha256: 16ed5478381242cd7ce76edaf735f7eed5b57de75a0e0d16665f17a71df17aa6
semantic_sha256: fa1f0e9b4ffef79f16e86a6aff42eb5bf7a779e211b9d820ee212a9aa82fc6ac
hash_basis: raw LF bytes
---

# THM-3771 -- arbitrary radial dressing reaches a vertical-address wall

**RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT HOSTILE
AUDIT.**  This file is self-contained.  In particular, it does not use the
presently reserved THM-3770.

Let `k` be an algebraically closed field of characteristic zero.  Choose

```text
u in k[z]\{0},                  r in k,
phi(W)=a+bW+cW^2 in k[W]\{0},                         (1)
```

where `phi` is allowed to have degree zero or one.  In `k[X,T]` put

```text
z=XT,              U=Xu(z),              W=U+3z+r,
Q=U phi(W).                                                (2)
```

Then the following are equivalent:

```text
Q is smooth on A^2;

u(0)phi(r)!=0,  u and phi are squarefree,
and gcd(u(z),phi(3z+r))=1.                                (3)
```

Here a nonzero constant is squarefree.  Under `(3)`, for every
`lambda in k*`:

1. the rational function

   ```text
   P0=-lambda W/(3Q)                                      (4)
   ```

   satisfies `J(P0,Q)=lambda`;

2. every rational solution of `J(P,Q)=lambda` is uniquely of the form

   ```text
   P=P0+H(Q),                    H in k(Q);                (5)
   ```

3. a polynomial solution exists if and only if both `u` and `phi` are
   constant.  On that boundary `Q` is a nonzero multiple of `X`.

If at least one profile is nonconstant, `Q` is therefore a smooth reducible
noncoordinate with an explicitly integrated rational Hamiltonian time form,
but it is not a planar Keller component.

The promised cubic radial-carrier form is literal.  With `v=3z+r`, expansion
of the quadratic in `(1)` gives

```text
Q=XA(z)+X^2B(z)+X^3C(z),
A=u phi(v),            B=u^2 phi'(v),            C=c u^3. (6)
```

Thus `c!=0` is the genuine cubic sector, while `c=0` records its sharp
quadratic and linear boundaries.  A normalized integral hostile is

```text
u=1+z,                 r=3,                 phi(W)=W^2-5,

Q=X^6T^3+6X^5T^3+9X^4T^3+3X^5T^2+18X^4T^2+27X^3T^2
  +3X^4T+18X^3T+22X^2T+X^3+6X^2+4X.                    (7)
```

It is smooth of total degree nine, has the rational mate `-W/(3Q)`, and has
no polynomial mate.

## 1. The birational chart is log-canonical

Work first in the independent chart `(X,z)`.  Since

```text
J_(X,T)(X,z)=X,                                         (8)
```

direct differentiation gives

```text
J_(X,z)(U,W)=3u(z),
J_(X,T)(U,W)=3Xu(z)=3U.                                 (9)
```

No information has been discarded: `(U,W)` is a birational coordinate pair.
Indeed

```text
z=(W-U-r)/3,          X=U/u(z),          T=z/X.         (10)
```

Consequently

```text
k(X,T)=k(U,W),
J(W,Q)=-3U phi(W)=-3Q.                                 (11)
```

Equation `(11)` proves `(4)` immediately.  It also proves completeness.
Writing `L=Q`, one has

```text
k(U,W)=k(L,W),                  U=L/phi(W),
J(-,Q)=-3L partial_W.                                      (12)
```

The constant field of this derivation is exactly `k(L)`.  Hence
`J(P,Q)=lambda` is equivalent to

```text
partial_W P=-lambda/(3L),                               (13)
```

whose complete rational solution is `(5)`.  This proves both rational
existence and the full primitive torsor, rather than a bounded response jet.

## 2. Exact smoothness, including every boundary

The chart determinant in `(9)` is nonzero wherever `U!=0`.  In that locus,

```text
dQ=phi(W)dU+U phi'(W)dW                                (14)
```

vanishes exactly when `phi` has a repeated root.  Conversely, if `rho` is a
repeated root of `phi`, choose `U!=0` so that

```text
z=(rho-U-r)/3
```

is not a root of `u`.  This is possible because `k` is infinite.  Formula
`(10)` then supplies a source point with `W=rho`, and `(14)` vanishes there.
Thus squarefreeness of `phi` is necessary and sufficient off `U=0`.

The zero set of `U=Xu(z)` has one axis component and one hyperbola for every
root `beta` of `u`:

```text
D_0: X=0, z=0, W=r;
D_beta: z=beta, W=r+3beta.                              (15)
```

At the generic point of these components,

```text
D_0:       dU=u(0)dX;
D_beta:    dU=Xu'(beta)dz.                              (16)
```

On `U=0`, equation `(14)` reduces to `dQ=phi(W)dU`.  Therefore the axis is
safe exactly when `u(0)phi(r)!=0`; every root component is safe exactly when
the root is simple and `phi(r+3beta)!=0`.  These conditions are precisely
`(3)`.

This also proves necessity without hidden genericity:

- `u(0)=0` or `phi(r)=0` makes the origin critical;
- a multiple root of `u` makes `(16)` vanish on its component;
- a common root of `u(z)` and `phi(3z+r)` kills `(14)` there;
- a multiple root of `phi` gives the off-axis critical point constructed
  above.

There are no other source points.  This proves the smoothness equivalence.

## 3. The zero fibre remembers all vertical addresses

Assume `(3)`.  Since `Q` is smooth, every irreducible component of `Q=0` is
reduced.  Its value of `W` belongs to the exact finite set

```text
S={r}
  union {r+3beta : u(beta)=0}
  union {rho : phi(rho)=0}.                            (17)
```

Every displayed value occurs: the first two groups come from `(15)`, and
each divisor `W-rho=0` above a root of `phi` is nonempty.  Moreover all
values in `S` are pairwise distinct.  Squarefreeness separates roots within
each profile; `u(0)!=0` separates `r` from the `u`-addresses;
`phi(r)!=0` separates `r` from the `phi`-addresses; and the gcd condition in
`(3)` separates the two non-axis groups.  Hence

```text
|S|=1+deg(u)+deg(phi).                                  (18)
```

At the generic point of a component with address `w_D in S`, `(4)` has
simple `1/Q` principal coefficient

```text
-lambda w_D/3.                                         (19)
```

In the complete torsor `(5)`, a pole of `H` at zero of order at least two
creates an uncancellable higher pole on every component.  A regular `H`
does not alter `(19)`.  A simple pole `h_-1/Q` adds the same coefficient
`h_-1` on every component, so it regularizes all components only if all
addresses in `S` agree.  By `(18)`, this happens exactly when
`deg(u)=deg(phi)=0`.

If either profile is nonconstant, no rational mate is regular on `A^2`, and
therefore no polynomial mate exists.  The factorization

```text
Q=Xu(z)phi(W)                                           (20)
```

is then nontrivial, so `Q` is also a noncoordinate.

If `u=u0` and `phi=phi0` are nonzero constants, then

```text
Q=u0 phi0 X,                    P=-lambda T/(u0 phi0)   (21)
```

is a polynomial Keller pair.  Equivalently, the target shear
`H(Q)=lambda r/(3Q)` removes the only nonconstant denominator in `(4)`.
This proves the iff and shows that both the constant-`u` and constant-`phi`
faces are sharp: either face remains obstructed until the other profile is
constant too.

## 4. Exact controls and scope

The companion named in the metadata verifies the universal chart,
quadratic expansion, inverse, Jacobian sign, and primitive; direct Groebner
smoothness for constant, linear, and quadratic profile controls; all five
failure mechanisms in Section 2; the complete address sets in Section 3;
the positive coordinate boundary `(21)`; the expanded degree-nine hostile
`(7)`; and full polynomial-mate linear systems for that hostile through
total mate degree eight.  Normal and optimized executions are required to
byte-match the frozen transcript before promotion.

The theorem classifies only the log-canonical dressing `(2)`.  It does not
classify arbitrary triples `A,B,C` in `(6)`, and it does not prove the planar
Jacobian conjecture.  Its constructive message is instead exact: arbitrary
radial freedom in `u` preserves rational integrability, but it creates one
new zero-fibre address per root, moving a polynomial mate farther away rather
than closer.  **QED, conditional on independent hostile audit.**

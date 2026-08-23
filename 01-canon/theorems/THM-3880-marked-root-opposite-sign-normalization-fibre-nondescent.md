---
id: THM-3880
title: "Marked-root descent is exactly sign matching on nodal-cuspidal normalizations"
status: >
  PROVISIONAL STRENGTHENING CANDIDATE + VERIFIED-EXACT; awaiting independent
  hostile audit of the full iff.  Universally, opposite nonzero marked-root
  signs over one normalization fibre obstruct every polynomial carrier.
  More sharply, if the forced rational carrier extends regularly and the
  plane curve has only ordinary nodes and A2 cusps, it descends if and only
  if no node fibre has opposite nonzero signs.  Same-sign node values and
  every A2 derivative condition are automatic, including points on A=0.
  The regular-extension pole gate and higher singularity jets remain open.
source: jc_sparse_direct_search / post-THM-3876 intrinsic sign-monodromy abstraction, 2026-08-23
audit: >
  STRENGTHENING SELF-AUDITED; independent hostile audit pending.  The
  assertion-free symbolic companion checks the cusp identity, both
  global-sign carrier formulas and their exact jumps, same-sign node values,
  the A=0 node value, both A2 derivative seams, the THM-3876 specialization,
  and positive z=0/A=0 controls.  Normal and optimized runs byte-match the
  frozen 19-gate transcript.  The local-to-global iff uses the standard
  completed local rings of an ordinary node and an A2 cusp; no finite census
  substitutes for that proof.
related:
  - THM-3859-marked-root-polynomial-graph-companion-puncture-obstruction
  - THM-3864-integrated-three-cusp-conductor-seminormal-three-direction-gate
  - THM-3876-monomial-marked-root-normalization-nondescent
script: 04-computation/jc2_marked_root_opposite_sign_nondescent_thm3880.py
output: 05-knowledge/results/jc2_marked_root_opposite_sign_nondescent_thm3880.out
script_sha256: cca8a8d35d9e00ff3af6ffef08e8410d124fc84b025d74e8f260476e7c5cf506
output_sha256: 44f8bd1022260ee5c8ae760be91e3a4e07be43df5f4cac69d14cc3bb4f02e28f
semantic_sha256: 47182169c96e351ac2cbf7f476be52b554ab1abe596639afb59121af10112006
hash_basis: raw LF bytes
---

# THM-3880 -- the marked-root sign and conductor theorem

**PROVISIONAL STRENGTHENING CANDIDATE + VERIFIED-EXACT; awaiting independent
hostile audit.**  Work over an algebraically closed field `k` of
characteristic zero.  For `b in k[A,C]`, put

```text
Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b,                       (1)

P=1+(2/3)AC,                         u=1+AC+A^2b.               (2)
```

Let `Gamma subset A2_(A,C)` be an irreducible reduced affine curve and

```text
nu: X -> Gamma                                                    (3)
```

its normalization.  Use the same letters `A,C` for their pullbacks to the
integral normal curve `X`, and assume that the marked square has a regular
root

```text
z in O(X),                         z^2=P.                       (4)
```

There are two conclusions.

### Universal opposite-sign obstruction

If two distinct geometric points `p,q in X(k)` satisfy

```text
nu(p)=nu(q),                 z(q)=-z(p),                 z(p)!=0, (5)
```

then no polynomial `b(A,C)` with `Delta_b|_Gamma=0` exists.

### Complete ordinary-node/A2 criterion under regular extension

Assume additionally that `A` is not the zero function, every singularity of
`Gamma` is an ordinary node or an `A2` cusp, and fix
`epsilon in {+1,-1}`.  On `A!=0`, define the forced rational carrier

```text
B_epsilon=(epsilon z^3-1-AC)/A^2 in k(X).                      (6)
```

Suppose `(6)` extends to a regular function `B_epsilon in O(X)`.  Then

```text
B_epsilon in O(Gamma)
 iff
no ordinary-node fibre {p,q} has z(q)=-z(p)!=0.                (7)
```

Because `Gamma` is affine in the `(A,C)`-plane, membership on the left means
exactly that `B_epsilon` is the restriction of some polynomial `b(A,C)`;
that lift automatically satisfies `Delta_b|_Gamma=0`.

The hypotheses in the second statement are load-bearing.  The theorem does
not prove that `(6)` has no poles, and singularities beyond nodes and `A2`
cusps can impose higher conductor jets.  Those are now the two precise open
directions.

## 1. One global sign on an integral normalization

The universal cusp identity is

```text
A^2 Delta_b=27(P^3-u^2).                                      (8)
```

Suppose first that a polynomial carrier exists and write `B=nu^*b`.  From
`(4)` and `(8)`, one has in the domain `O(X)`

```text
u^2=z^6.
```

Therefore `(u-z^3)(u+z^3)=0`, and integrality forces one sign globally:

```text
u=epsilon z^3,                         epsilon in {+1,-1}.      (9)
```

This is not a sign chosen independently on each branch.

Conversely, under the regular-extension hypothesis `(6)`, equation `(9)`
holds on the dense open set `A!=0` and hence on all of `X`.  Substitution in
`(8)` gives `A^2 Delta_B=0`.  Since `A` is a nonzero function in the domain,

```text
Delta_B=0 in O(X).                                             (10)
```

Thus every regular extension in the second statement is already a genuine
normalization-level discriminant solution.

There is also an intrinsic recovery of the square root from a putative
carrier.  If `P` is not the zero function, then `z_b=u/P in k(X)` satisfies
`z_b^2=P`.  At a zero of `P` of order `d`, one has

```text
2 ord(u)=3d,
ord(z_b)=ord(u)-ord(P)=d/2>=0.                                 (11)
```

Hence `d` is even and `z_b` is regular.  Any prescribed root `(4)` differs
from `z_b` by one global sign.  This makes `(9)` intrinsic.

## 2. Opposite signs contradict descent without dividing by A

If `p,q` lie over one point of `Gamma`, every descended function has equal
values there.  For a polynomial carrier,

```text
A(p)=A(q),       C(p)=C(q),       B(p)=B(q),       u(p)=u(q). (12)
```

But `(5)` and `(9)` give

```text
u(q)=epsilon z(q)^3=-epsilon z(p)^3=-u(p).                     (13)
```

Equations `(12)--(13)` imply `z(p)=0`, contrary to `(5)`.  This proves the
universal obstruction.

No division occurs here.  If the common point lies on `A=0`, then `P=1`
and the opposite roots are `+1,-1`; equation `(13)` still gives opposite
nonzero `u`-values.  A denominator-based proof would miss this seam.

If `Gamma` itself is the vertical line `A=0`, then `z^2=1` in the domain
`O(X)`, so `z` is globally `+1` or globally `-1`; `(5)` cannot occur.  This
trivial smooth boundary is outside `(6)` only because division by `A^2` is
unnecessary there.

## 3. Exact forced values away from A=0

Since

```text
AC=(3/2)(z^2-1),                                               (14)
```

the two signs in `(6)` have the factored forms

```text
A^2B_+(z)= (z-1)^2(2z+1)/2,                                   (15)

A^2B_-(z)=-(z+1)^2(2z-1)/2.                                  (16)
```

At two opposite addresses, ordered as `z -> -z`,

```text
B_epsilon(-z)-B_epsilon(z)
 =-2 epsilon z^3/A^2.                                         (17)
```

This is the explicit nonzero jump when `Az!=0`; Section 2 is its regular
extension across `A=0`.

## 4. Local conductor conditions for nodes and A2 cusps

Let `R=O(Gamma)` and `S=O(X)`.  The finite quotient `S/R` is supported at
the singular points.  Its completed local descriptions are standard:

```text
ordinary node:
  R_hat=k[[x,y]]/(xy) subset k[[x]] direct_sum k[[y]],
  membership iff the two constant terms agree;                 (18)

A2 cusp:
  R_hat=k[[tau^2,tau^3]] subset k[[tau]],
  membership iff the coefficient of tau is zero.               (19)
```

Equivalently, `(19)` asks for zero first derivative in a normalization
parameter.  Because completion is faithfully flat and `S/R` has finite
length, a regular function in `S` belongs to `R` exactly when it passes
`(18)--(19)` at every singularity.

## 5. Same-sign node values are automatic

Consider an ordinary node with normalization addresses `p,q`.  Since
`z(p)^2=z(q)^2=P(nu(p))`, their root values are equal or opposite.  If the
forbidden case in `(7)` does not occur, then

```text
z(p)=z(q)                                                       (20)
```

(including the common value zero).

If the node has `A!=0`, formula `(6)` and equality of `A,C,z` immediately
give

```text
B_epsilon(p)=B_epsilon(q).                                    (21)
```

If the node lies on `A=0`, use `(10)` rather than `(6)`.  Specializing the
original discriminant gives

```text
0=Delta_B|_(A=0)=9C^2-54B,
B=C^2/6.                                                       (22)
```

Both branches have the same `C`, so `(21)` holds there as well.  Thus every
allowed node passes exactly the conductor condition `(18)`.

## 6. Every A2 first-jet condition is automatic

Let `p` be an `A2` cusp and use a local normalization parameter.  Write a
prime for its derivative.  Both plane coordinates have zero derivative:

```text
A'(p)=C'(p)=0.                                                 (23)
```

First suppose `A(p)!=0`.  Differentiate `z^2=P`:

```text
2zz'=(2/3)(A'C+AC')=0.                                        (24)
```

Differentiate the numerator identity in `(6)`:

```text
2AA'B+A^2B'=3 epsilon z^2z'-A'C-AC'.                          (25)
```

At the cusp, the right side is zero: `A'=C'=0`, and
`z^2z'=z(zz')=0` by `(24)`.  Since `A(p)!=0`, equation `(25)` gives

```text
B'(p)=0.                                                       (26)
```

This argument includes `z(p)=0`; it never divides by `z`.

Now suppose `A(p)=0`.  Differentiate the polynomial identity `(10)`
directly.  On substituting `A=A'=C'=0`, every term vanishes except

```text
-54B'(p)=0.                                                    (27)
```

Thus `(26)` also holds on the `A=0` seam.  Every cusp passes `(19)`.

Sections 5--6 prove the sufficiency direction of `(7)`.  The necessity
direction is Section 2, applied to a polynomial representative of any
element of `R`.  This completes the iff.

## 7. Sharp z=0 and A=0 controls

At `z=0`, one has

```text
P=0,                         2AC+3=0.                          (28)
```

Both signs give `u=0`, and the jump `(17)` vanishes.  This is a genuine
positive boundary.  On the irreducible hyperbola in `(28)`, take

```text
b_0=2C^2/9.
```

Then

```text
Delta_(b_0)=-(C^2/3)(2AC+3)^2.                                (29)
```

Likewise `A=0` is not forbidden when its sign is constant.  The profile
`b_A=C^2/6` satisfies

```text
Delta_(b_A)=-(AC^3/4)(3AC+4),                                 (30)
```

and contains the vertical component with `z=1`.  The obstruction is the
opposite-sign gluing, not either divisor by itself.

## 8. THM-3876 is exactly the sign-flip specialization

After gcd reduction, THM-3876 has

```text
A=r^M,
C=6r^N(1+r^(M+N)),
z=1+2r^(M+N).                                                  (31)
```

For `M>=3`, choose primitive `zeta^M=1`, put

```text
eta=zeta^N,
U=r^(M+N)=-1/(eta+1).                                         (32)
```

At the colliding addresses `r,zeta r`,

```text
z(r)=(eta-1)/(eta+1),               z(zeta r)=-z(r)!=0.       (33)
```

Thus the universal obstruction applies.  Using
`U^2=A^2r^(2N)`, the jump `(17)` becomes exactly

```text
B(zeta r)-B(r)=-2r^(2N)(eta-1)^3/(eta+1),                     (34)
```

the THM-3876 primitive-root difference.

## 9. Exact replay and remaining frontier

Run

```text
python3 04-computation/jc2_marked_root_opposite_sign_nondescent_thm3880.py
python3 -O 04-computation/jc2_marked_root_opposite_sign_nondescent_thm3880.py
```

and compare both streams byte-for-byte with the frozen output.  The exact
companion checks the universal identities and jumps, both node seams, both
cusp derivative seams, `(29)--(30)`, and the THM-3876 specialization.

The theorem converts the former broad same-sign uncertainty into two exact
frontiers:

1. **regular-extension/pole gate:** decide when the rational function `(6)`
   has no pole at any finite normalization point, especially above `A=0`;
2. **higher conductor jets:** for singularities worse than `A2`, determine
   whether the marked identities force all missing Taylor coefficients.

For ordinary nodes and `A2` cusps, once regularity is known, there is no
additional hidden conductor obstruction beyond the sign packet `(7)`.

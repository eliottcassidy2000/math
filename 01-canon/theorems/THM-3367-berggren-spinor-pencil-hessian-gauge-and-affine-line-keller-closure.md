---
id: THM-3367
title: "Berggren spinor pencil Hessian gauge and affine-line Keller closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  After the fixed
  middle-branch target gauge, the full three-weight Berggren spinor pencil is
  the universal symmetric 2 by 2 pencil.  Polynomial row-integrable
  pullbacks are exactly fixed target transforms of binary gradient maps, and
  their determinant is minus the Hessian determinant.  No affine rank-two
  coefficient pullback has constant determinant.  Every integrable
  constant-nonzero pullback whose coefficient image lies in an affine line
  is an explicit tame shear, in every degree.  The top homogeneous Hessian
  channel of every remaining Keller pullback is one isotropic pure-power
  ray.  This is a binary symmetric subproblem and not a reduction of JC(2).
source: root/factorial-jacobian-lrc-threebranch-2026-08-14
depends_on:
  - THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit
related:
  - THM-2063-one-fiber-linear-planar-keller-pairs
  - THM-2811-special-image-linear-intertwiner-rigidity-and-jc-degree-wall
  - THM-3010-ballot-column-newton-ratios-and-metallic-alternation
  - THM-3362-odd-berggren-determinant-ray-factorial-closure
  - HYP-8905-binary-symmetric-jc2-subcase-and-separate-descent-programs
  - MISTAKE-237
script: 04-computation/berggren_spinor_hessian_gauge_thm3367.py
output: 05-knowledge/results/berggren_spinor_hessian_gauge_thm3367.out
script_sha256: cc56744421ed39e47ca0d53d9f71c763292ecd4af91bda0f78f0ea01f8fab30a
output_sha256: 4d3ee069caf0ae34ebe88410ec52f1613000566831645339d6e2bd0131af8e17
hash_basis: LF-normalized bytes
---

# THM-3367 -- Berggren spinor pencil Hessian gauge and affine-line Keller closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem identifies the exact planar content of the three Berggren
spinor branches.  It yields a useful representation and two sharp closure
gates, but it does not turn the ternary tree into a proof of the planar
Jacobian Conjecture.

## 1. Inheritance and conventions

The closest proved source is
[THM-3357](THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md),
with its Euclid-parameter branch matrices

```text
L=[ 0 1;-1 2],       M=[0 1;1 2],       R=[1 0;2 1].               (1)
```

The closest planar closure is
[THM-2063](THM-2063-one-fiber-linear-planar-keller-pairs.md): a planar
Keller pair with one affine source fibre in one output-pencil direction is
tame.  The corrected near miss is MISTAKE-237: a binary symmetric-Hessian
calculation is not a bridge from NC2/GMC to full `JC(2)`, and the stable
symmetric reduction of a general planar pair raises dimension.  The
least-used sidecar activated here is **row integrability** of a matrix-valued
polynomial.  A determinant identity alone does not supply it.

Work over `C`.  Source coordinates are `r=(s,t)^T`, and Jacobians have output
polynomials as rows:

```text
Jac(P,Q)=[P_s P_t;Q_s Q_t].                              (2)
```

For coefficient polynomials `a,b,c in C[s,t]`, put

```text
B(a,b,c)=aL+bM+cR.                                      (3)
```

These coefficients are now arbitrary source polynomials.  This change of
type forgets positivity, primitivity, branch words, and ternary-tree
ancestry.  Nothing below claims that `(a,b,c)` is a branch-count vector.

## 2. The middle branch gauges the pencil to the full Hessian frame

Let

```text
D=[-1 0;0 1],                    J=[0 1;1 0].             (4)
```

Direct multiplication gives

```text
L=MD,                            R=MJ,                    (5)
D^2=J^2=I,                       DJ+JD=0.                 (6)
```

Consequently

```text
B(a,b,c)=M S(a,b,c),                                      (7)
S(a,b,c)=aD+bI+cJ=[b-a c;c a+b].                         (8)
```

Thus the fixed left gauge `M^(-1)` sends the three Berggren matrices to a
basis of the complete three-dimensional space `Sym_2(C)`.  Since
`det(M)=-1`, the exact determinant signs are

```text
det S=b^2-a^2-c^2,
det B=a^2-b^2+c^2=:q(a,b,c).                             (9)
```

The quadratic form `q` is nondegenerate.  Equations (5)--(9), rather than a
shared word such as "determinant", are the lawful algebraic connection.

## 3. Row integrability is exactly the binary Hessian condition

**Theorem (integrability dictionary).**  For `a,b,c in C[s,t]`, the following
are equivalent.

1. There are `P,Q in C[s,t]` with `B(a,b,c)=Jac(P,Q)`.
2. There is `H in C[s,t]` such that

```text
a=(H_tt-H_ss)/2,       b=(H_tt+H_ss)/2,       c=H_st.   (10)
```

Under this correspondence one may take

```text
(P,Q)=M grad(H)=(H_t,H_s+2H_t).                         (11)
```

Moreover

```text
det Jac(P,Q)=q(a,b,c)=-det Hess(H).                     (12)
```

### Proof

Suppose first that `B=Jac(P,Q)`, and put

```text
G=M^(-1)(P,Q)^T.                                       (13)
```

Because `M` is constant, `Jac(G)=M^(-1)B=S`.  Matrix `S` is symmetric, so

```text
partial_t G_1=S_12=S_21=partial_s G_2.                 (14)
```

The polynomial one-form `G_1 ds+G_2 dt` is closed.  Integrating `G_1` in
`s` and using (14) to choose the remaining polynomial in `t` gives a
polynomial `H` with `G=grad(H)`.  Equating `Hess(H)` with (8) gives (10),
and multiplying by `M` gives (11).

Conversely, (10) makes `S=Hess(H)`, so

```text
B=M Hess(H)=Jac(M grad(H)).                             (15)
```

Taking determinants in (15) proves (12).  Linear terms of `H` absorb all
constant output translations, so there is no omitted constant case.  QED.

Therefore a row-integrable pullback of (3) is Keller exactly when

```text
q(a,b,c)=kappa in C*,       det Hess(H)=-kappa.         (16)
```

The fixed target map `M` preserves polynomial invertibility.  This is an
exact representation of the **binary symmetric** Keller subproblem, not of
all planar Keller pairs.

## 4. No affine rank-two coefficient slice is Keller

Let an affine coefficient pullback have the form

```text
w(r)=(a,b,c)^T=w_0+A r,       A in Mat_(3 by 2)(C).     (17)
```

**Theorem (affine-rank gate).**  If `q(w(r))` is constant, then

```text
rank(A)<=1.                                               (18)
```

In particular, no affine rank-two pullback has constant determinant, even
before row integrability is imposed.

### Proof

Let `W=im(A)`, and polarize

```text
q(x_1,x_2,x_3)=x_1^2-x_2^2+x_3^2                      (19)
```

to its nondegenerate symmetric bilinear form.  The homogeneous quadratic
part of `q(w_0+Ar)` is `q(Ar)`.  Constancy therefore gives `q|_W=0`, and
polarization gives `W subset W^perp`.  Nondegeneracy in dimension three
gives

```text
dim W<=dim W^perp=3-dim W,
```

so `2 dim W<=3` and `dim W<=1`.  QED.

The exact companion also checks the generic algebraic version.  If `u,v`
are the two columns of `A`, the equations

```text
q(u)=q(v)=<u,v>_q=0                                    (20)
```

make the square of every `2 by 2` minor of `[u v]` reduce to zero in the
generic Groebner ideal.  Over `C`, every minor itself therefore vanishes.

Rank zero in (17) gives a constant Jacobian matrix and, when (16) holds, an
affine automorphism.  The rank-one case is part of the all-degree closure
below.

## 5. Every integrable affine-line image is an explicit tame shear

The next theorem is not restricted to affine weights or bounded degree.

**Theorem (affine-line image closure).**  Suppose the coefficient image is
contained in one affine line.  In the nonconstant presentation write

```text
w(r)=w_0+v h(r),       v in C^3 nonzero,
h in C[s,t] nonconstant.                                (21)
```

Assume that `q(w)=kappa` is constant and that `B(w)` is row-integrable.  Put

```text
C=S(w_0),                         N=S(v).                (22)
```

Then there are `g in C^2` nonzero, `rho in C*`, a one-variable polynomial
`F`, a vector `d`, and a scalar `e` such that

```text
N=rho g g^T,                    ell=g^T r,               (23)
H(r)=(1/2)r^T C r+F(ell)+d^T r+e,                       (24)
g^T adj(C)g=0.                                           (25)
```

Here `H` is the potential in Section 3.  If `kappa!=0`, then `C` is
invertible,

```text
g^T C^(-1)g=0,                                           (26)
```

and the resulting Keller map `Phi=M grad(H)` is tame.  Given a target
`y=Phi(r)`, its inverse is

```text
u=M^(-1)y-d,
sigma=g^T C^(-1)u,
r=C^(-1)(u-F'(sigma)g).                                  (27)
```

### Proof

Let the polar form of `q` be

```text
<x,y>_q=x_1y_1-x_2y_2+x_3y_3.                           (28)
```

Since a nonconstant polynomial in `C[s,t]` is transcendental over `C`, the
identity

```text
q(w_0+v h)=q(w_0)+2<w_0,v>_q h+q(v)h^2=kappa           (29)
```

forces

```text
q(v)=0,              <w_0,v>_q=0,              q(w_0)=kappa. (30)
```

The map `S:C^3 -> Sym_2(C)` in (8) is a linear isomorphism.  By (9),
`det N=-q(v)=0`; since `v!=0`, `N` is a nonzero symmetric rank-one matrix.
Over `C` it has the form in (23).

Row integrability and Section 3 give

```text
Hess(H)=C+N h.                                           (31)
```

Third derivatives of a polynomial potential are fully symmetric.  Writing
`h_j=partial_j h`, (31) therefore gives

```text
N_ij h_k=N_ik h_j.                                      (32)
```

Choose `i` with `g_i!=0` in (23).  Equation (32) reduces to

```text
g_j h_k=g_k h_j,                                        (33)
```

so `grad(h)` is everywhere parallel to `g`.  Complete `ell=g^T r` to a
linear coordinate system `(ell,m)`.  Equation (33) says `partial_m h=0`,
and hence

```text
h=f(ell)                         for some f in C[T].     (34)
```

Choose `F` with `F''=rho f`.  The Hessian of the difference between `H` and
the first two terms in (24) is zero, so that difference is affine.  This
proves (24).

For every scalar `z`, the rank-one determinant identity gives

```text
det(C+z rho g g^T)=det C+z rho g^T adj(C)g.             (35)
```

Constancy and nonconstancy of `f` prove (25).  If `kappa!=0`, then by (9)

```text
det C=-q(w_0)=-kappa!=0,
```

and (25) is equivalent to (26).

Finally

```text
grad(H)=C r+F'(ell)g+d.                                 (36)
```

For `u` as in (27),

```text
C^(-1)u=r+F'(ell)C^(-1)g,
g^T C^(-1)u=ell+F'(ell)g^T C^(-1)g=ell.                (37)
```

Thus `sigma=ell` and (27) follows.  More structurally, put
`h_0=C^(-1)g`.  The map before the invertible linear factor `C` is

```text
r -> r+F'(g^T r)h_0,              g^T h_0=0.           (38)
```

It fixes `ell` and shears one complementary coordinate by a polynomial in
`ell`, so it is elementary triangular after a linear change.  Multiplication
by `C` and `M`, and translation by `d`, prove tameness.  QED.

If `v=0` or `h` is constant, the coefficient map is constant; a nonzero
constant determinant gives the affine boundary already noted.  If
`kappa=0`, equations (23)--(25) and the normal form (24) still hold, but `C`
is singular and no inverse is claimed.  These are the exact degenerate
boundaries.

Combining Sections 4 and 5 proves:

```text
every row-integrable affine-weight Keller pullback is tame;              (39)
every row-integrable constant-q affine-line-image pullback is tame.       (40)
```

Equation (40) covers arbitrary degree.  Its map always has a linear
output-pencil member after affine changes, so it is consistent with the
broader closure in THM-2063; (23)--(38) supply the transfer-native normal
form and inverse.

## 6. The top Hessian channel collapses to one isotropic ray

The affine-line theorem does not say that every polynomial solution of
`det Hess(H) in C*` has line image.  Lower Hessian layers may introduce new
directions.  The top layer, however, is rigid.

**Binary Hessian-zero lemma.**  Let `G(s,t)` be a nonzero homogeneous binary
form of degree `d>=2`.  Then

```text
det Hess(G)=0        iff        G=lambda ell^d          (41)
```

for `lambda in C*` and a nonzero linear form `ell`.

### Proof of the lemma

Write, on the chart `t!=0`,

```text
G(s,t)=t^d f(z),                    z=s/t.              (42)
```

Direct differentiation gives

```text
det Hess(G)
=(d-1)t^(2d-4)[d f f''-(d-1)(f')^2].                  (43)
```

If `f` is constant, (42) is a pure power.  Otherwise choose a finite root
`z_0` of `f`, of multiplicity `m>=1`, and write its first coefficient as
`A!=0`.  The coefficient of `(z-z_0)^(2m-2)` in the bracket in (43) is

```text
A^2[d m(m-1)-(d-1)m^2]=A^2 m(m-d).                    (44)
```

It must vanish, so `m=d`.  Since `deg f<=d`, this forces
`f=A(z-z_0)^d`.  The converse follows because the Hessian of a pure power
has rank at most one.  QED.

Now let `H in C[s,t]` satisfy `det Hess(H) in C*`, and put `d=deg H`.

- If `d<=1`, its Hessian determinant is zero, so this case is impossible.
- If `d=2`, its Hessian is constant and nonsingular, so `grad(H)` and
  `M grad(H)` are affine automorphisms.
- If `d>=3`, the unique total-degree `2d-4` part of `det Hess(H)` is
  `det Hess(H_d)`, where `H_d` is the top homogeneous form.  It must vanish.

The lemma therefore gives

```text
H_d=lambda ell^d.                                       (45)
```

Writing `ell=g^T r`,

```text
Hess(H_d)=lambda d(d-1)ell^(d-2) g g^T.                (46)
```

Since `S` is a linear isomorphism, the degree-`d-2` coefficient triple has
the exact form

```text
(a,b,c)_(d-2)=v ell^(d-2),        v!=0,        q(v)=0. (47)
```

Thus a second coefficient direction can enter only in a lower Hessian
layer.  A candidate outside the affine-line closure must be nonlinear and
must retain the first lower layer at which the coefficient image leaves
the top isotropic ray.

This also fixes the scope of the THM-3010 contact.  The middle matrix `M`
has the silver norm-minus-one recurrence, but (7) shows that it is a fixed
legal target gauge and becomes the identity under `M^(-1)`.  Its iterated
metallic alternation is therefore not a planar-Keller invariant.  Independently,
(45)--(47) give only one leading direction, so the two-direction metallic
maximal-alternation circuit cannot occur in the top Hessian channel.  A
lower-layer, target-shear-compatible circuit would be a necessary sidecar;
none is claimed here.

## 7. Hostiles, connection contract, and the separate cubic determinant

### Constant determinant without integrability

Let

```text
p=s,                 r=-t(2+st),                 c_0=1+st,
a=(p+r)/2,           b=(r-p)/2.                         (48)
```

Then

```text
q(a,b,c_0)=pr+c_0^2=1.                                 (49)
```

Nevertheless

```text
S=[-s 1+st;1+st -t(2+st)],                             (50)
partial_t S_11-partial_s S_12=-t,
partial_t S_12-partial_s S_22=s+t^2.                   (51)
```

Both Hessian curl equations fail.  This is a two-parameter polynomial map
into the determinant quadric `q=1`, but it supplies no planar Jacobian pair.
It is the canonical hostile to dropping the row-integrability sidecar.

### Positive affine-line control

For

```text
H=st+s^6/30                                             (52)
```

the weights are `(-s^4/2,s^4/2,1)`, `q=1`, and

```text
Phi=(s,t+s^5/5+2s),
Phi^(-1)(u,v)=(u,v-u^5/5-2u).                          (53)
```

This tests a non-affine, arbitrary-degree line image rather than only the
constant or quadratic boundary.

### The FC determinant is a different pencil

THM-3357 Section 7 and THM-3362 use

```text
P(a,b,c)=det(aT_L+bT_M+cT_R)
        =(a-b-c)(a+b-c)(a+b+c),                        (54)
```

where the `T` matrices are `3 by 3` nodewise Pythagorean transfers.  The
quadratic lift agrees with each branch separately; it is not additive under
linear mixing.  The cheapest exact hostile is

```text
det(L+M+R)=q(1,1,1)=1,
det(T_L+T_M+T_R)=P(1,1,1)=-3.                          (55)
```

Therefore (54)'s FC(3) factorial-moment law does not transfer through (7)
to a planar Keller determinant.  Conversely, the Hessian dictionary for the
`2 by 2` spinor pencil says nothing about factorial moments of (54).

The full connection contract is

```text
source:     three-weight Berggren spinor pencil B
target:     binary symmetric matrix/Hessian pencil S=M^(-1)B
map:        fixed left target gauge, followed by row integration
preserves:  determinant up to the sign det(M)=-1; Keller and invertibility
loses:      positivity, tree words, ancestry, Pell iteration, branch owners
sidecar:    both Hessian curl equations and lower-layer coefficient directions
hostiles:   (48)--(51) for integrability; (55) for false 2x2/3x3 mixing. (56)
```

MISTAKE-237 remains load-bearing.  This theorem does not supply a
two-variable symmetric reduction of an arbitrary planar Keller pair, does
not reverse any Special-Image implication, and does not connect NC2 or GMC
to `JC(2)`.  It makes no literature-priority or novelty claim.

## 8. Exact verification

The companion uses exact arithmetic only.  Its explicit verifier universe is:

1. generic matrix identities over `Z[a,b,c]`;
2. a generic six-symbol affine-rank ideal, where all three squared minors
   reduce to zero;
3. one dense quartic potential over `Q[s,t]` for (10)--(12);
4. the rational control (52)--(53) and a second control over `Q(i)` with
   `C=I`, `g=(1,i)`;
5. the constant and zero-determinant boundary controls;
6. the exact hostile (48)--(51);
7. formula (43) with generic coefficients for every `2<=d<=8`, with a pure
   power and a non-pure hostile at each degree; and
8. the cubic factorization (54) and determinant-mixing hostile (55).

The finite degree range in item 7 audits the implementation; the universal
quantifier in the lemma follows from (42)--(44), not from that range.
Ordinary and optimized runs must byte-match the stored transcript:

```text
python3 04-computation/berggren_spinor_hessian_gauge_thm3367.py
python3 -O 04-computation/berggren_spinor_hessian_gauge_thm3367.py
```

Both end in `PASS`.

**End of proof.**

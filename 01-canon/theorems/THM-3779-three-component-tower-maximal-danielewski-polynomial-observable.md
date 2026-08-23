---
id: THM-3779
title: "Three-component tower maximal Danielewski polynomial observable"
status: >
  PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.  For every
  m>=2 in the THM-3774 rational Keller tower, V=UP and the newly recognized
  E=P(V-1) are polynomial, and the complete intersection of the rational
  target field with the source polynomial ring is the smooth exponent-one
  Danielewski ring k[U,V,E]/(UE-V(V-1)).  The induced unramified map has
  exact image the surface minus one point.  Over C the nonzero logarithmic
  symplectic class excludes every Darboux pair in this maximal observable,
  hence excludes every polynomial Keller pair obtained by any rational
  target-field change.  The m=1 boundary is false: it omits a divisor and
  x^2 is an additional polynomial target-field observable.  At m=3 the
  cover has S4 monodromy and its cubic resolvent realizes the V4 quotient.
source: root + jc_quartic_c3_construct / target-field observable search, 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The exact companion verifies E's source
  formula, the relation and all three Poisson minors for m=1,...,8, the
  two arm parametrizations, finite-root recovery identities, critical
  double-root residual degrees, the m=1 inverse-discriminant hostile, and
  the m=2/m=3 discriminants, branch factorizations, cubic resolvent, scaled
  polynomial resolvent, and logarithmic chart transition.  Normal and
  optimized replay and an independent audit of the image/DVR argument are
  still required before status promotion.
depends_on:
  - THM-3774-three-component-rational-keller-cover-tower
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
related:
  - THM-3777-three-component-tower-reciprocal-payment-horizontal-divisor-transfer
  - THM-3561-rational-keller-danielewski-polynomial-completion
script: 04-computation/jc2_three_component_maximal_danielewski_observable_thm3779.py
output: 05-knowledge/results/jc2_three_component_maximal_danielewski_observable_thm3779.out
script_sha256: 2d028f2a2e30b266c4ab2192e6f7bb1d34b0ffd365c10ce271e12925e771fc65
output_sha256: 5e1ae1396b85b6e89c3d2caf1fea08ea24ca5a8c22e734a7c2d684e4492522f0
semantic_sha256: 9d7beea76b64741acea7c87e547b03ef4f672dfef7295272b3605de287c639b8
hash_basis: raw LF bytes
---

# THM-3779 -- the full polynomial target-field observable is Danielewski

**PROVED + VERIFIED-EXACT + PENDING INDEPENDENT HOSTILE AUDIT.**  This
theorem closes the complete rational target-field route for the first
non-Galois members of THM-3774.  It does not construct a planar Jacobian
counterexample and says nothing about polynomial functions outside the
rational target field.

Let `k` be an algebraically closed field of characteristic zero, let `m>=2`,
and put

```text
A=1+xy,
B=1+x^(2m+1)A^m,
U=xAB,
P=((m+1)B-1)/(mx).                                    (1)
```

THM-3774 proves

```text
J(U,P)=1,
k(x,y)=k(U,P,t),
t^(m+1)-mPt+mU=0,             [k(x,y):k(U,P)]=m+1,    (2)
```

where `t=x^2A`.  The theorem below determines exactly which elements of
`k(U,P)` are polynomials in `x,y`.

## 1. The missing observable

Set

```text
C_m=(m+1)B-1,                   V=UP.                  (3)
```

The blowdown `V` is polynomial, but it is not the maximal polynomial
observable.  Define

```text
N_m = y
    + ((2m+1)/m)x^(2m)A^(m+1)
    + ((m+1)/m)x^(4m+1)A^(2m+1).                      (4)
```

Direct expansion gives the two exact divisibilities

```text
V-1=xN_m,
E:=P(V-1)=C_m N_m/m in k[x,y].                        (5)
```

Thus `E` is a polynomial source function lying in the rational target field.
It satisfies

```text
UE=V(V-1).                                             (6)
```

In particular `E` does not belong to `k[U,V]`: otherwise `(6)` would make
`U` divide `V(V-1)` in the polynomial ring `k[U,V]`.  Algebraic
independence of `U,V` follows from `(2)` and `P=V/U`.

Let

```text
D=k[u,v,e]/(ue-v(v-1)),             Y=Spec D.          (7)
```

The assignment `(u,v,e)=(U,V,E)` is injective.  Indeed, the ring in `(7)`
is a domain of dimension two, while its image has fraction field
`k(U,V)=k(U,P)` of transcendence degree two.  Hence `(6)` is the complete
relation among the three observables.

## 2. A smooth unramified target

The three source Jacobian minors are

```text
J(U,V)=U,
J(U,E)=2V-1,
J(V,E)=E.                                              (8)
```

They satisfy the explicit unit identity

```text
(2V-1)^2-4UE=1.                                       (9)
```

The same identity proves that `Y` is smooth: its gradient
`(e,u,1-2v)` has no common zero.  Equations `(8),(9)` show that

```text
phi_m:A2_(x,y) -> Y,          (x,y) |-> (U,V,E),       (10)
```

is everywhere unramified, hence locally quasi-finite.  It is dominant and
has generic degree `m+1` by `(2)`.

## 3. The exact image for m>=2

Let `(u,v,e)` be a geometric point of `Y`.

### 3.1 The open set u!=0

Put `p=v/u` and consider

```text
f_m(T)=T^(m+1)-mpT+mu.                                (11)
```

For a root `t` of `(11)`, set

```text
s=u-t^(m+1).                                          (12)
```

A root with `s=0` must obey

```text
u=t^(m+1),                 p=((m+1)/m)t^m.             (13)
```

Because `u!=0`, such a root is nonzero.  For fixed `(u,p)` it is the unique
distinct root satisfying `(13)`: it equals `(m+1)u/(mp)`.  Moreover

```text
f_m'(t)=0,
f_m''(t)=m(m+1)t^(m-1)!=0,                             (14)
```

so it has multiplicity exactly two.  Since `deg f_m=m+1>=3`, at least one
other root remains.  That root has `s!=0`.  If `(13)` has no solution, every
root is already good.

For any good root the exact reconstruction is

```text
x=t/s,
A=s^2/t,
y=s(s^2-t)/t^2,
B=u/s.                                                 (15)
```

Here `t!=0` because `f_m(0)=mu!=0`.  Formulas `(12),(15)` give

```text
xA=s,               x^2A=t,              xAB=u,
P=1/x+((m+1)/m)t^m=u/t+t^m/m=p.                       (16)
```

Thus every geometric point with `u!=0` has a source preimage.  Notice the
load-bearing use of `m>=2`: the bad root consumes two sheets, and a third
sheet is needed.

### 3.2 The two arms

When `u=0`, relation `(7)` gives `v=0` or `v=1`.  Every point of the
`v=1` arm is hit by

```text
x=0,                      y=e,                         (17)
```

because `(U,V,E)|_(x=0)=(0,1,y)`.  Every point of the `v=0` arm with
`e!=0` is hit by

```text
x=-1/e,                   y=e,                         (18)
```

for which `A=0` and `(U,V,E)=(0,0,e)`.  The other zero-fibre component
`B=0` gives the compatible parametrization `E=1/(mx)`.

The point

```text
o=(0,0,0) in Y                                             (19)
```

has no preimage.  Indeed `U=xAB`, and the three factors are pairwise
coprime.  On their zero components the exact restrictions are

```text
x=0:  (V,E)=(1,y),
A=0:  (V,E)=(0,-1/x),
B=0:  (V,E)=(0, 1/(mx)).                              (20)
```

The latter two parameters never vanish.  Sections 3.1 and 3.2 prove

```text
phi_m(A2)=Y minus {o}                       for every m>=2. (21)
```

This is a set-theoretic statement on geometric points, and hence a
scheme-theoretic surjectivity statement over the displayed open target.

## 4. Exact polynomial intersection

Work inside `k(x,y)` and put

```text
K=k(U,P)=Frac D,
R_m=k[x,y] intersection K.                            (22)
```

The inclusion `D subset R_m` follows from `(5)--(7)`.  For the converse,
take `F in R_m` and a height-one prime `q` of the normal ring `D`.
Its generic point is not the omitted codimension-two point `o`, so `(21)`
and local quasi-finiteness supply a height-one prime `r` of `k[x,y]`
mapping dominantly to `q`.  If `pi_q` is a uniformizer, then

```text
ord_r(F)=e(r/q) ord_q(F),             e(r/q)>0.        (23)
```

The left side is nonnegative because `F` is a source polynomial.  Hence
`ord_q(F)>=0`.  This holds for every height-one prime of `D`.  Smoothness
makes `D` normal, and the Krull intersection of its height-one DVRs gives

```text
k[x,y] intersection k(U,P)
   =k[U,V,E]
   ~=k[u,v,e]/(ue-v(v-1))                 (m>=2).      (24)
```

This is the maximal polynomial observable in the rational target field,
not merely another convenient intermediate ring.

## 5. The sharp m=1 failure

The restriction `m>=2` cannot be removed.  At `m=1`, take `(u,p)=(1,2)`.
Then

```text
f_1(T)=T^2-2T+1=(T-1)^2,                              (25)
```

and its only root has `s=1-T^2=0`.  Thus the point
`(u,v,e)=(1,2,2)` of `Y` is omitted.  In fact the whole quadratic
discriminant branch is omitted, and the source identity

```text
P^2-4U=1/x^2,
x^2=1/(P^2-4U) in k[x,y] intersection k(U,P)           (26)
```

exhibits an observable outside `D`.  So codimension-one coverage, not just
the formulas `(5)--(9)`, is load-bearing in `(24)`.

## 6. The exponent-one obstruction closes all target-field words

Now take `k=C`.  Change coordinates on `Y` by

```text
b=1-v,                    c=u,                    e=e. (27)
```

Then `(7),(8)` become the standard exponent-one, two-arm Danielewski model

```text
ce=b(b-1),
{b,c}=c,             {c,e}=-(2b-1),             {b,e}=-e. (28)
```

On `c!=0` its symplectic form is

```text
omega=db wedge dc/c.                                     (29)
```

Its pullback by `(10)` is `dx wedge dy`, because
`db wedge dc=dU wedge dV=U dx wedge dy`.

THM-3600, specifically its exponent-one boundary `(23),(24)`, proves that
`[omega]` is nonzero in algebraic de Rham cohomology whenever the squarefree
polynomial has at least two roots.  Here the roots are `0,1`.  Concretely,
the two arm-plane coordinates are glued by

```text
a_1=a_0-c^(-1),                                        (30)
```

and their local Liouville potentials differ by `-dc/c`, whose residue is
nonzero.  Therefore `omega` is not exact.

If `F,G in D` satisfied `{F,G}=kappa!=0`, then

```text
omega=kappa^(-1)dF wedge dG
     =kappa^(-1)d(F dG),                               (31)
```

contradicting `(29),(30)`.  Thus `D` has no polynomial Darboux pair.
Together with `(24)` this gives the exact all-word consequence:

```text
F,G in k(U,P),   F,G in k[x,y]
        => J_(x,y)(F,G) is not a nonzero constant.      (32)
```

In particular no rational symplectic target change of `(U,P)` can make both
outputs polynomial.  This includes target words of arbitrary length, not
only the three-shear orientation of THM-3777.  It does not constrain a
construction that leaves `k(U,P)`.

## 7. Cubic and quartic cover anatomy

The finite normalization of `Y` in `k(x,y)` adds the ramified sheets that
the unramified nonproper map `(10)` leaves at infinity.

### 7.1 The cubic m=2 member

For `m=2`, `(11)` and its discriminant are

```text
f_2(T)=T^3-2PT+2U,
disc(f_2)=4(8P^3-27U^2).                              (33)
```

THM-3774 proves irreducibility.  The cuspidal factor is irreducible and not
a square, so the Galois closure has group `S3`.  On the normalization of the
branch curve,

```text
P=(3/2)r^2,                  U=r^3,
f_2(T)=(T-r)^2(T+2r).                                 (34)
```

The double root has `s=U-r^3=0` and lies on the missing finite-normalization
boundary.  The simple root `-2r` supplies the source point when `r!=0`.
At `r=0` the target is the omitted point `o`.  This explains how the field
extension is ramified while `(10)` is everywhere unramified.

### 7.2 The quartic m=3 member and its V4 quotient

For `m=3`,

```text
f_3(T)=T^4-3PT+3U,
disc(f_3)=27(256U^3-81P^4).                           (35)
```

The standard cubic resolvent is

```text
R_3(Z)=Z^3-12UZ-9P^2,                                 (36)
```

and has the same discriminant as `(35)`.  It is irreducible over `k(U,P)`:
solving `(36)` for `U` gives

```text
U=(Z^3-9P^2)/(12Z),                                   (37)
```

a coprime rational map of degree three over `k(P)`.  The irreducible quartic
from THM-3774 and irreducible resolvent leave `A4` or `S4` as the generic
group; the nonsquare irreducible discriminant in `(35)` excludes `A4`.
Hence

```text
Gal(f_3^split/k(U,P))=S4.                              (38)
```

The action of `S4` on the three pair-partitions has kernel `V4`.  Thus the
splitting field of `(36)` is exactly the fixed field of `V4`, with quotient

```text
S4/V4 ~= S3.                                           (39)
```

Along the branch normalization,

```text
P=(4/3)r^3,                  U=r^4,
f_3(T)=(T-r)^2(T^2+2rT+3r^2),
R_3(Z)=(Z+2r^2)^2(Z-4r^2).                            (40)
```

Again the quartic double root is the missing `s=0` boundary sheet, while the
two residual roots remain visible on the source.

There is a polynomial presentation of the resolvent over `D`: if
`Z_0=UZ`, then

```text
Z_0^3-12U^3Z_0-9UV^2=0,                               (41)
```

whose discriminant is

```text
27U^2(256U^7-81V^4).                                  (42)
```

This is an integral cubic-resolvent model whose Galois closure is the
`V4`-quotient field, not another source polynomial.  Indeed a
resolvent-root field is not contained in the quartic root field inside the
generic `S4` closure.  The genuinely new target-field polynomial is `E`,
and `(24)` proves that `U,V,E` exhaust all of them.

## 8. Scope and exact controls

The companion named in the metadata checks `(4)--(9),(15)--(20)` for
`m=1,...,8`, including the critical double-root quotient of degree `m-1`.
It checks the hostile identity `(26)`, all formulas `(33)--(42)`, the
standard two-arm coordinate change, and the logarithmic transition.  The
normal and optimized runs must byte-match the frozen transcript before
promotion.

The theorem proves:

1. the complete polynomial target-field intersection for every `m>=2`;
2. the complete no-Darboux/no-rational-target-word consequence over `C`;
3. the exact `S3` cubic and `S4 -> S3` quartic-resolvent anatomy; and
4. the codimension-one `m=1` failure and its first extra observable.

It proves no statement about polynomial pairs using functions outside
`k(U,P)`, no planar Jacobian counterexample, and no positive-characteristic
analogue.  **QED, conditional only on independent hostile audit.**

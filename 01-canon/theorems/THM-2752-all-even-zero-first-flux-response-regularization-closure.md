---
id: THM-2752
title: "All-even zero-first-flux response regularization closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the full
  chosen-sheet split degree-22 exact-square-prefix
  response family, every all-even member with zero first flux is physically
  empty, including reducible and nonreduced members.  The five transverse
  infinity points are excluded by the exact-prefix lift resultant.  At the
  remaining point P_infty every normalization branch has
  min(v(q),v(s))>=4v(h), while the exact homogeneous cancellation
  R_25+(d/2)F_23 has order at least 28v(h); hence the affine response is
  regular there.  A physical component either has a forbidden smooth pole or
  makes its nonconstant response globally regular on a projective curve.
  Combined with THM-2725 and THM-2745, this closes the entire chosen-
  sheet split degree-22 response family, not the broader split branch,
  JC(2), or DC(2).
source: jc-even-zero-flux-next-2026-07-28
audit: jc-even-zero-flux-audit-2026-07-28 (independent Faber reconstruction, valuation and unequal-slope audit, weighted-index-cover typing, componentwise projective argument, normal/optimized replay, and hash audit: ACCEPT)
addendum_audit: root-2026-07-28-all-degree-even-response (independent recurrence review and symbolic control through m=50 of the post-audit all-M=4L+2 local-order, face-gcd, and q^3 addendum: ACCEPT)
depends_on:
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
related:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2725-split-even-faber-nonzero-first-flux-unified-closure
  - THM-2745-highest-odd-faber-componentwise-exact-prefix-closure
  - THM-2747-highest-odd-reduced-boundary-divisor-and-one-ended-factorization-atlas
  - THM-2755-all-even-zero-flux-componentwise-global-regularity-closure
script: 04-computation/jc2_degree22_all_even_zero_flux_response_regularization_20260728.py
output: 05-knowledge/results/jc2_degree22_all_even_zero_flux_response_regularization_20260728.out
script_sha256: 7e92b7f241965f9f8cdc4619488237101fc0b5b897c4493897d2ba9ad20a07af
output_sha256: 8f8a2389bb407ef3555f3bbb222ec218d5928d54851bfbb336756ba8af9874cc
hash_basis: LF-normalized bytes
---

# THM-2752 -- the all-even zero-flux response becomes regular at its last end

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The odd response argument fails at the all-even boundary for a structural
reason: all five even rows have one common slope-four Newton face, and their
third responses cancel the first flux to one additional local order.  That
failure is itself the obstruction.  It makes the response regular at the
only infinity point not already forbidden by the polynomial exact prefix.

## 1. Exact statement

Work over `C` with the chosen-sheet split polynomial exact-square-prefix
degree-22 response family in

```text
P(1,2,3,4)_[h,d,q,s].
```

Suppose all odd Faber coefficients vanish and the chosen-sheet first flux is
zero.  Thus

```text
F_23=Phi_22+a_14 h^8 Phi_14+a_10 h^12 Phi_10
             +a_6 h^16 Phi_6+a_2 h^20 Phi_2=0,

G_24=Psi_22+a_14 h^8 Psi_14+a_10 h^12 Psi_10
             +a_6 h^16 Psi_6+a_2 h^20 Psi_2-W h^24=0.       (1)
```

The coefficients `a_2,a_6,a_10,a_14,W` are arbitrary elements of `C`.  Let

```text
R_25=R_22+a_14 h^8 R_14+a_10 h^12 R_10
             +a_6 h^16 R_6+a_2 h^20 R_2.                   (2)
```

Then no reduced irreducible component of `(1)` carries the response map of a
polynomial split exact-square-prefix Keller trajectory.  Since the source is
reduced, the same conclusion holds for reducible or nonreduced members.

## 2. The coefficient-independent top boundary

At `h=0`, every lower row and both flux constants disappear.  After removing
nonzero rational scalars, the two top equations are

```text
Phi_22=unit*q*P,                 Psi_22=unit*Q,

P=-84d^2q^4s+3dq^6+280dq^2s^3-21q^4s^2-84s^5,

Q=-224d^3q^6+3360d^2q^4s^2-336dq^6s-3360dq^2s^4
  +3q^8+560q^4s^3+224s^6.                              (3)
```

They are coprime.  Hence `(1)` has no surface component and no curve
component contained in `h=0`; every projective curve component meets `h=0`
in the finite top intersection.

When `q=0`, `(3)` gives `s=0`, leaving only

```text
P_infty=[h:d:q:s]=[0:1:0:0].                           (4)
```

On `q!=0`, pass to the `q=1` index chart and the coarse invariants

```text
z=s^3,                         rho=d/s^2.
```

The two equations become

```text
Pbar=3rho-21+z(-84rho^2+280rho-84),

Qbar=3+z(-336rho+560)
       +z^2(-224rho^3+3360rho^2-3360rho+224).           (5)
```

Their lexicographic basis is one linear equation in `rho` and

```text
p_5(z)=20141047808z^5-14386462720z^4+1089822720z^3
       -21288960z^2-35910z+81.                         (6)
```

The exact companion verifies that `(6)` is squarefree, that `z`, `rho`, and
the Jacobian of `(5)` are nonzero in the quotient, and that the stripped top
response

```text
Tbar=3+z(9rho^2-105rho+70)
       +z^2(-84rho^3+280rho^2-84rho)                   (7)
```

is a unit there.  Thus the rest of the top boundary consists of exactly five
transverse coarse points.  At each, `h` has order one and

```text
R_aff=R_25/h^25
```

has pole order `25`.  This calculation uses only the top degree-22 row; it
does not import a conclusion stated under a nonzero-odd-seed hypothesis.

## 3. Every branch over `P_infty` has slope at least four

Pass to the finite `d=1` local index cover at `(4)`.  Let `v` be the discrete
valuation on any normalization branch, and put

```text
e=v(h)>0,             u=v(q),             v_0=v(s),
a=min(u,v_0),                                           (8)
```

allowing value `+infinity` for an identically zero coordinate.  For

```text
m=2,6,10,14,22
```

the local `(q,s)` orders of both `Phi_m` and `Psi_m` are respectively

```text
k_m=1,2,3,4,6,                                         (9)
```

and

```text
(22-m)+4k_m=24.                                        (10)
```

Suppose first that `a<4e`.  Every lower even row has valuation at least

```text
(22-m)e+k_m a
 =6a+(6-k_m)(4e-a)>6a,                                 (11)
```

and `W h^24` also has valuation strictly greater than `6a`.  Thus only the
degree-six local faces of the top row can occur initially.

If `u<v_0`, the pure `q^6` monomial of the `Psi_22` face is its unique
lowest term and has nonzero coefficient.  If `v_0<u`, the pure `s^6`
monomial is uniquely lowest.  Both cases contradict `G_24=0`.  If
`u=v_0=a`, the two initial equations are

```text
q s(q^2-3s^2)(3q^2-s^2)=0,
(q^2-s^2)(q^4-14q^2s^2+s^4)=0.                        (12)
```

The two homogeneous forms in `(12)` have gcd one, so they have no common
nonzero projective direction.  This is the final contradiction.  The same
unequal-valuation argument covers exactly one of `q,s` being identically
zero.  Therefore every nontrivial branch satisfies

```text
min(v(q),v(s))>=4v(h).                                 (13)
```

If `q=s=0` identically, treat the branch separately below; no minimum needs
to be assigned.

The residual `mu_2` on the `d=1` index cover may exchange two lifted
branches.  Passing to that finite cover multiplies all displayed orders by
the same ramification index.  Inequality `(13)` and regularity of an invariant
rational function therefore descend to the coarse normalization.

## 4. The even response gains four orders

The homogeneous, correctly typed cancellation uses `d`, of weight two, not
`h^2`.  Exact Faber recurrence gives

```text
D_m=R_m+(d/2)Phi_m,

D_2 =0,
D_6 =-q^3/4,
D_10=5q^3s/8,
D_14=7q^3(dq^2-5s^2)/32,
D_22=-33q^3(6d^2q^4-84dq^2s^2+3q^4s+70s^4)/1024.     (14)
```

On the `d=1` chart, every nonzero `D_m` has `(q,s)` order at least
`k_m+1`.  Because the first flux is zero, summing `(14)` gives the exact
homogeneous identity

```text
R_25+(d/2)F_23
  =D_22+a_14h^8D_14+a_10h^12D_10
         +a_6h^16D_6+a_2h^20D_2.                      (15)
```

On a branch of `(1)`, `F_23=0`.  Equations `(10)`, `(13)`, and `(15)` give

```text
v(R_25)>=(22-m)e+4(k_m+1)e=28e                        (16)
```

term by term.  Consequently

```text
v(R_aff)=v(R_25/h^25)>=3e>0.                           (17)
```

Thus the response is regular and actually vanishes at every normalization
point over `P_infty`.  If `q=s=0` identically, every expression in `(14)`
vanishes identically, so `(17)` holds directly.

This is exactly why the highest-odd argument disappears.  An odd seed
produces a transverse response pole at `(4)`; the all-even bank instead has
the cancellation `(15)`.  A nonzero first flux would add
`-(lambda d/2)h^23` to the right side and destroy `(16)`.

## 4A. The `q^3` cancellation is an all-degree parity syzygy

The common factor in `(14)` is not a degree-22 accident.  For every positive
integer

```text
m=4ell+2,                                                (17a)
```

define the Faber coefficients by the formal series

```text
(1+2dz^2+qz^3+(d^2-s)z^4)^(m/4)=sum_(n>=0)c_n z^n.      (17b)
```

The coefficient recurrence used throughout the exact companions gives

```text
Phi_m=4c_(m+1),
R_m=4c_(m+3)+2d c_(m+1),
D_m=R_m+(d/2)Phi_m=4(c_(m+3)+d c_(m+1)).                (17c)
```

Put

```text
N=m/2,             alpha=m/4,             beta=alpha-1,
A(y)=1+2dy+(d^2-s)y^2.                                  (17d)
```

Expanding

```text
(A(z^2)+qz^3)^alpha
 =sum_(j>=0) binom(alpha,j)q^j z^(3j)A(z^2)^(alpha-j) (17e)
```

shows first that both coefficients in `(17c)` contain only odd powers of
`q`: their indices are odd, whereas `A(z^2)` is even in `z`.  Write

```text
A(y)^beta=sum_(n>=0)b_n y^n.
```

The coefficient of `q` in `D_m/4` is

```text
alpha(b_N+d b_(N-1)).                                  (17f)
```

Differentiate `A^beta`.  The identity

```text
A(A^beta)'=beta A' A^beta
```

gives the coefficient recurrence

```text
(n+1)b_(n+1)+2d(n-beta)b_n
 +(d^2-s)(n-1-2beta)b_(n-1)=0.                         (17g)
```

At `n=N-1`, the last coefficient is zero because
`2beta=N-2`, while the middle coefficient is `dN`.  Hence

```text
N(b_N+d b_(N-1))=0.                                    (17h)
```

Characteristic zero and `(17f)` kill the linear term.  Parity already
killed the constant and quadratic terms, so

```text
R_m+(d/2)Phi_m is divisible by q^3
in Q[d,q,s] for every m=4ell+2.                         (17i)
```

This lemma explains the common `q^3` visible in `(14)`.  It does **not** by
itself give the extra local-order invoice `(9)`--`(10)`, align another
degree's coefficient gaps, or close any broader split family.  The
degree-22 proof still uses the five exact rows and the slope-four valuation
argument above.

There is also a stronger all-degree explanation of that local-order invoice.
On the `d=1` chart write

```text
1+2z^2+qz^3+(1-s)z^4=(1+z^2)^2+(qz^3-sz^4).          (17j)
```

For `m=4ell+2`, put `k=(m+2)/4=ell+1`.  Expanding the `m/4` power by
perturbation order `r` gives

```text
sum_(r>=0) binom(m/4,r)(qz^3-sz^4)^r
                         (1+z^2)^(m/2-2r).             (17k)
```

For `r<k`, the last exponent is a positive integer and the whole summand has
`z`-degree at most `m`.  It contributes nothing to `c_(m+1)`, `c_(m+2)`, or
`c_(m+3)`.  At `r=k`, the last factor is `(1+z^2)^(-1)`.  The numerator
`(qz^3-sz^4)^k` has degree at most `4k=m+2`; hence, if its quotient has
coefficients `b_n`, the identity

```text
(1+z^2)sum b_nz^n=(qz^3-sz^4)^k
```

gives

```text
b_(m+3)+b_(m+1)=0.                                    (17l)
```

By `(17c)`, the entire order-`k` contribution to `D_m` cancels.  Therefore

```text
ord_(q,s)(Phi_m)=ord_(q,s)(Psi_m)=k,
ord_(q,s)(R_m+(1/2)Phi_m)>=k+1                         (17m)
```

on `d=1`.  Up to one common nonzero scalar, the two leading flux faces are

```text
Phi_face=sum_(a odd) binom(k,a)(-1)^((a-1)/2)
                              q^a(-s)^(k-a),

Psi_face=sum_(a even) binom(k,a)(-1)^(a/2)
                              q^a(-s)^(k-a).            (17n)
```

With `zeta=-s+iq` and `zetabar=-s-iq`, these are the odd and even parts of
`zeta^k`: their sum and difference recover `zeta^k` and `zetabar^k`.
The two linear forms are coprime, so the faces in `(17n)` have gcd one.
Moreover a pure `s^k` term occurs in `Psi_face`, while a pure `q^k` term
occurs in one of the two faces.  This proves the unequal- and equal-slope
obstruction uniformly for every such `m`.

Consequently, for a top degree `M=4L+2` and lower rows of the same congruence
class,

```text
(M-m)+4k_m=M+2,
(M-m)+4(k_m+1)=M+6.                                   (17o)
```

Thus the zero-first-flux response is locally regular at `P_infty` in this
whole formal all-even family: its numerator has order at least `M+6`, versus
denominator order `M+3`.  This is only a local regularization theorem.  It
does not classify the other `h=0` points, provide their exact-prefix
resultants, derive the response chart, or close any degree other than the
audited degree `22` application above.

## 5. A smooth response pole cannot lift

Let `X` be the normalization of the reduced irreducible component containing
the generic physical image.  THM-2723 gives on the source line

```text
U R_source'=kappa,                         kappa!=0,   (18)
```

and classifies `R_source` as either a nonconstant affine-linear function
with `U` constant or

```text
U=u_0(x-a)^M,
R_source=r_0+r_1(x-a)^(1-M),               M>=2.      (19)
```

The nonconstant source map extends to a finite surjective morphism

```text
gamma:P1_x -> X.                                        (20)
```

Suppose `X` contains one of the five smooth boundary points.  Its response
pole of order `25` pulls back with order `25e_gamma>1`.  The constant-`U`
alternative in `(18)` has only a simple pole, so `(19)` holds and the unique
source pole is the finite point `a`.

Write the polynomial exact prefix as

```text
P_source=H^2+L,
H=U^2z_source^2+Bz_source+C,       L=Az_source+E.
```

On the split sheet put

```text
h w=(hU)z_source+beta.
```

Coefficient comparison gives the two exact identities

```text
beta^2+d=h^2C,
q beta-s=h^4E.                                         (21)
```

All polynomial coefficients are regular at the finite point `a`.  The first
identity in `(21)` says explicitly that
`beta^2=h^2C-d` belongs to the source DVR.  Thus `v(beta)>=0`, so `beta` is
regular.  Taking residues at a smooth boundary point gives, for
`omega=beta(a)`,

```text
omega^2+d_0=0,                    q_0 omega-s_0=0.     (22)
```

Here `q_0d_0s_0!=0`.  Substitution `d=-omega^2`, `s=q omega` in the two top
forms yields, after removing nonzero factors,

```text
56omega^3+3q,
7168omega^6+896omega^3q+3q^2.                         (23)
```

Their exact `q`-resultant is

```text
-76608 omega^6.                                       (24)
```

It cannot vanish because `d_0=-omega^2!=0`.  Therefore no physical component
can contain a smooth boundary point.  This rederives the lift obstruction in
the present all-even scope rather than citing an odd-member conclusion.

## 6. Projective regularity closes the last component

The rational primitive in `(18)` has exactly one pole on the projective
source.  Thus even before `(24)`, surjectivity of `(20)` prevents `X` from
containing two distinct smooth response-pole points: their inverse fibres
would be disjoint nonempty pole fibres.  Section 5 excludes the remaining
possibility of exactly one.

Every boundary point of `X` must consequently lie over `P_infty`.  By
Section 4, `R_aff` is regular at every such normalization point.  It is
already regular on the affine chart `h!=0`.  Hence `R_aff` is a global
regular function on the projective integral curve `X`, so it is constant.
Its pullback `R_source` is constant, contradicting `(18)`.  This proves the
candidate theorem.

If the complete intersection is reducible, the argument applies only to the
reduced irreducible generic-image component.  If it is nonreduced, the map
from the reduced source kills nilpotents and factors through the target
reduction.  Embedded or nilpotent structure cannot change `(18)`, the
boundary valuation, or the exact identities `(21)`.

## 7. Synthesis and exact boundary

This result combines with:

1. THM-2725, which closes the all-even nonzero-first-flux bank; and
2. THM-2745, which closes every member with a nonzero odd coefficient,

to empty the entire chosen-sheet split degree-22 response family `(1)` with
the odd columns restored.

This does not derive that response family for every planar Keller pair,
close another reduced degree, perform integral two-adic raising, treat a
nonpolynomial or non-exact prefix, prove the broader split branch, or prove
`JC(2)` or `DC(2)`.

## 8. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree22_all_even_zero_flux_response_regularization_20260728.py
python3 -O 04-computation/jc2_degree22_all_even_zero_flux_response_regularization_20260728.py
```

Both executions byte-match

```text
05-knowledge/results/jc2_degree22_all_even_zero_flux_response_regularization_20260728.out.
```

The companion independently reconstructs every Faber row, the universal
five-point top boundary, the coprime `P_infty` faces, all five homogeneous
cancellations `(14)`, the slope-four divisibility of `(15)`, and the smooth
resultant `(24)`.  It uses explicit `require` checks rather than
optimization-sensitive `assert` nodes.

Current LF-normalized SHA-256 values are

```text
script 7e92b7f241965f9f8cdc4619488237101fc0b5b897c4493897d2ba9ad20a07af
output 8f8a2389bb407ef3555f3bbb222ec218d5928d54851bfbb336756ba8af9874cc
```

QED.

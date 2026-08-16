---
id: THM-3506
title: "Fixed Keller five-face norm transform and the 271/99 boundary"
status: >
  PROVED (conditional Newton/resultant face transform) + VERIFIED-EXACT
  (fixed-map packet and finite-sheet gate).  For the fixed sporadic Keller
  tower, the first untested cleared norm has exposed pair (271,99), not
  (259,87): G=L^43 N(J) has top face C x^271(3xz-2y)^99.  Moreover
  v_L(N(G))=-271, so R_5=L^271 N(G) is polynomial and coprime to L, and its
  exposed pair is (1699,615).  This theorem itself does not prove the two
  renewal faces for G; THM-3513 subsequently proves both for that fixed
  polynomial.  THM-3522 subsequently proves that complete packets renew in
  this fixed inverse chart whenever the next cleared norm is polynomial,
  closing the packets of R_5 and R_6.  THM-3523 subsequently closes the next
  finite-sheet/polynomiality gate and gives the packet of R_7.  THM-3528
  subsequently proves raw polynomial complete packets at all levels; later
  finite-sheet units, image primes, and every general Jacobian claim remain
  open.
source: codex/tropical-keller-norm/2026-08-16
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness
  - THM-3498-level-four-old-boundary-cancellation-and-degree81-discriminant-gate
  - THM-3504-level-four-sporadic-keller-image-prime-and-four-component-nonproperness
related:
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
  - THM-3513-fixed-G-hybrid-newton-renewal-faces
  - THM-3521-fixed-R5-finite-sheet-unit-and-next-old-L-clearing
  - THM-3522-fixed-keller-five-face-renewal-propagation
  - THM-3523-fixed-R6-finite-sheet-unit-and-next-old-L-clearing
script: 04-computation/keller_tropical_norm_face_recurrence_probe_20260816.py
output: 05-knowledge/results/keller_tropical_norm_face_recurrence_probe_20260816.out
script_sha256: fe1b03de9061c997a1abba6b88753e589c0a01f9f92762e49fe7ea0504ce9797
output_sha256: 0dd07d1af0621a9f767e9c803e805de61ee428fb6a980b004cd3f06625082b52
hash_basis: raw LF bytes; exact rational values use the ASCII numerator/denominator convention
---

# THM-3506 -- the fixed Keller five-face norm transform

**PROVED (conditional transform) + VERIFIED-EXACT (fixed-map instance).**

Let `F:C^3->C^3` be the fixed sporadic Keller map of THM-2473.  On the
target write

```text
L=27a^2c^2-18abc+16a+b^3c-b^2,
T=4-3bc,
S=27ac^2-9bc+8,
E(w)=Lw^3+Tw-2c.                                      (1)
```

Write `N(P)` for the cubic function-field norm of a source polynomial
`P(x,y,z)` through the inverse chart of THM-3495.  Thus `N(P)` is the
product of `P(q(w))` over the three roots of `E`.

## 1. Five faces, not one scalar exponent

For a monomial `x^i y^j z^k`, put

```text
lambda=i-k,       beta=i-j-2k,       gamma=i-j-5k.      (2)
```

Let `e,m` be nonnegative integers with `3|m`, and suppose all exponents
below are nonnegative.  A polynomial `P` has the five-face packet
`A(e,m)` when, up to independently nonzero rational scalars, its indicated
complete faces are

```text
max lambda=e:    x^e(3xz-2y)^m,                         (3)
min lambda=-e:   y^(3e-2m) z^e,                         (4)

min beta=-5e+2m:
  y^(3e-2m) z^(e-m)
  (y^2+27z)^(2m/3)(y^2+108z)^(m/3),                   (5)

max k=2e-2m/3:
  x^(2e-4m/3) z^(2e-2m/3),                             (6)

min gamma=-8e+2m:
  z^e(27x^2z+y^3)^(e-2m/3).                            (7)
```

The scalars in (3)--(7) need not agree.  The phrase *complete face* is
load-bearing: an extra equal-weight term can change a resultant even when
the displayed monomial remains present.

Let

```text
Q=L^e N(P),
e'=7e-2m,              m'=3e-2m.                       (8)
```

Assume `Q` is polynomial.  Then the first three faces of `Q` are, again up
to nonzero rational scalars,

```text
max lambda=e':    x^e'(3xz-2y)^m',                      (9)
min lambda=-e':   y^(3e'-2m') z^e',                    (10)

min beta=-5e'+2m':
  y^(3e'-2m') z^(e'-m')
  (y^2+27z)^(2m'/3)(y^2+108z)^(m'/3).                 (11)
```

In particular, the norm-face transform on the retained state is

```text
 [e']   [ 7 -2 ] [e]
 [m'] = [ 3 -2 ] [m].                                  (12)
```

Equations (9)--(11) do **not** by themselves prove that `Q` has the two
renewal faces (6)--(7).  If those two faces hold at every successive rung,
then (12) iterates rigorously.  Since the matrix has trace `5` and
determinant `-8`, each coordinate satisfies

```text
u_(n+2)=5u_(n+1)+8u_n.                                 (13)
```

This is the precise conditional induction inside this theorem.  THM-3528
subsequently discharges its polynomiality hypothesis at every raw rung.

### Conditional Cassini and projective sidecar

Start the formal matrix orbit at

```text
v_0=(e_0,m_0)=(1,0),       v_(n+1)=M v_n,
M=[[7,-2],[3,-2]].                                      (13a)
```

Whenever the renewal hypotheses make this an actual face-packet orbit,
the determinant is an exact Cassini invariant:

```text
det(v_n,v_(n+1))
 =det(M)^n det(v_0,v_1)=3(-8)^n.                       (13b)
```

Induction modulo `6` in (12) gives

```text
e_n=1 (mod 6),       m_n=3 (mod 6),       n>=1.         (13c)
```

For `n>=1`, any common divisor of `e_n,m_n` divides
`det(v_(n-1),v_n)=3(-8)^(n-1)`.  Equation (13c) excludes both prime
divisors `2` and `3`; hence

```text
gcd(e_n,m_n)=1.                                         (13d)
```

Thus `r_n=m_n/e_n` is reduced.  However (13b) has magnitude
`3*8^n`, never `1`, so consecutive fractions are not Farey neighbors and
the orbit is not a Stern--Brocot edge path.

Projectivizing (12) gives

```text
r |-> phi(r)=(3-2r)/(7-2r).                             (13e)
```

On `[0,1]`, `phi` maps into `[1/5,3/7]` and
`phi'(r)=-8/(7-2r)^2`, so it is an orientation-reversing contraction.
Consequently the ratios alternate around and converge to the unique fixed
point in this interval,

```text
r_*=(9-sqrt(57))/4.                                     (13f)
```

Since `e_n,m_n` are coprime odd integers for `n>=1`, Euclid's parity
division gives a primitive Pythagorean triple

```text
((e_n^2-m_n^2)/2, e_n m_n, (e_n^2+m_n^2)/2).
```

The first three are

```text
(20,21,29),       (812,645,1037),
(31820,26829,41621).
```

This is a Fibonacci-style/Cassini sidecar to the conditional matrix orbit,
not an identification with the Fibonacci sequence or a Farey-tree path.

## 2. The top face from the inverse chart

Use the target scaling

```text
(a,b,c)=(A/t,B,Ct),        t -> 0,
h=3AC-2B.                                                   (14)
```

Then `L=(16A)t^-1+O(1)`.  The Newton polygon of `E` has one linear edge
and one quadratic edge.  Exact substitution in the inverse formulas of
THM-3495 gives

```text
linear root:
  w=(C/2)t+O(t^2),
  (q_x,q_y,q_z)=((C/2)t,-h/2,A/t)+higher,               (15)

quadratic roots:
  w=W t^(1/2)+higher,       4AW^2+1=0,
  (q_x,q_y,q_z)=(W t^(1/2),
                  6AW t^(-1/2),
                 -26A t^-1)+higher.                    (16)
```

On the linear root, (4) contributes exactly

```text
constant * A^e h^(3e-2m) t^-e.                         (17)
```

On the two quadratic roots, the product of the two copies of (5) has
order `t^(-5e+2m)`.  Its leading coefficient cannot vanish: after putting
`A=1` and `W^2=-1/4`,

```text
y^2=-9,
y^2+27z=-711=-9*79,
y^2+108z=-2817=-9*313.                                 (18)
```

Thus the quadratic resultant is a nonzero Gaussian norm.  Multiplication
by `L^e` adds order `t^-e`; the total order is

```text
-e + (-5e+2m) - e = -(7e-2m)=-e'.                     (19)
```

The only `h`-dependence comes from (17), with exponent
`3e-2m=m'`.  Weighted homogeneity supplies `A^e'`.  This proves (9), and
it also proves the weaker top-face conclusion from just (4)--(5), without
the renewal faces (6)--(7).

With the product-of-roots norm convention, the raw coefficient in (9) is

```text
16^e c_- (-1)^(m')/2^(m')
 * Norm_(W^2=-1/4)( in_beta(P)(W,6W,-26) ),             (20)
```

where `c_-` is the coefficient in (4).  Formula (20) is the coefficient
used by the exact companion.

## 3. The opposite face from the Chebyshev end

Keep (14), but now let `t -> infinity`.  After division by the dominant
coefficient, the inverse cubic tends to

```text
(Bw)^3-3(Bw)-2=(Bw-2)(Bw+1)^2.                         (21)
```

Writing `q=Bw`, the three Puiseux roots have leading labels
`q=2,-1,-1`.  At a generic point of this toric divisor the inverse chart
has

```text
q=2:     q_y~-B/2,       q_z~-B^3 C t/8,
q=-1:    q_y~ B,         q_z~ B^3 C t.                 (22)
```

Put

```text
d=2e-2m/3,       p=2e-4m/3.                            (23)
```

The monomial face (6) is nonzero on all three labels in (21), so each
copy of `P` has exact order `t^d`.  Their product and `L^e` have order

```text
3d+e=7e-2m=e'.                                         (24)
```

The powers of `B,C` are

```text
B^(3e+9d-3p) C^(e+3d)
 =B^(15e-2m) C^e'
 =B^(3e'-2m') C^e'.                                    (25)
```

This is exactly the monomial (10).  The double limiting root in (21)
causes no cancellation because the norm multiplies the two Puiseux
factors and the face monomial is nonzero at their common leading label.

## 4. The beta face and the two exceptional factors

Use instead

```text
(a,b,c)=(At,B/t,C/t^2),       t -> 0.                   (26)
```

Now `L~B^3C t^-5`, and the same limiting cubic (21) appears after writing
`w=(q/B)t`.  Exact inverse-chart leading terms are

```text
q=2:
  (q_x,q_y,q_z)~(2t/B,-B/(2t),-B^3C/(8t^5)),

q=-1:
  (q_x,q_y,q_z)~(-t/B,B/t,B^3C/t^5).                   (27)
```

For `K=27x^2z+y^3`, these give

```text
q=2:      K~-B(B^2+108C)/(8t^3),
q=-1:     K~ B(B^2+27C)/t^3.                           (28)
```

Put `r=e-2m/3=m'/3`.  Substitution of the gamma face (7) into the one
simple and two double labels gives, after multiplication by `L^e`,

```text
constant *
B^(15e-2m) C^(4e)
(B^2+27C)^(2r)(B^2+108C)^r.                            (29)
```

Since

```text
15e-2m=3e'-2m',       4e=e'-m',       r=m'/3,          (30)
```

equation (29) is precisely (11).  Its weight is
`-29e+6m=-5e'+2m'`, proving both the face and its extremality.

Sections 2--4 derive (12) from the inverse chart and the two resultant
geometries.  No numerical interpolation of `1,7,43` is used.

## 5. Exact packet for L, H, and J

Retain the canonical normalizations

```text
L N(L)=H/64,
L^7 N(H)=J/2^35,
G=L^43 N(J).                                            (31)
```

Exact extraction from `L`, the frozen `H` polynomial, and THM-3495's
`66,146`-term coefficient ledger gives the following five-face data.  All
omitted coefficients are nonzero.

| polynomial | `(e,m)` | beta face after its scalar is removed | z-top face | gamma face after its scalar is removed |
|---|---:|---|---|---|
| `L` | `(1,0)` | `y^3 z` | `x^2 z^2` | `z(27x^2z+y^3)` |
| `H` | `(7,3)` | `y^15 z^4(y^2+27z)^2(y^2+108z)` | `x^10 z^12` | `z^7(27x^2z+y^3)^5` |
| `J` | `(43,15)` | `y^99 z^28(y^2+27z)^10(y^2+108z)^5` | `x^66 z^76` | `z^43(27x^2z+y^3)^33` |

The opposite `lambda` faces are respectively

```text
y^3z,       y^15z^7,       y^99z^43.                   (32)
```

The companion independently recovers the two already-proved transforms in
(31), coefficient-for-coefficient at the exposed face.  This calibrates all
normalizations before the first untested step.

Applying (12) to `J` gives

```text
(e_G,m_G)=(7*43-2*15, 3*43-2*15)=(271,99),             (33)
```

and therefore

```text
in_max-lambda(G)=C_G x^271(3xz-2y)^99,       C_G!=0.    (34)
```

For the exact normalization `G=L^43N(J)`, `C_G` is a negative integer with
numerator bit length `1179`, and the SHA256 of the ASCII string
`numerator/denominator` is

```text
fa3ed2569db09620551a553a22c8392e193ba12ec8783ad66456c3af9fe6f49f. (35)
```

The transform also proves two previously unexpanded faces of `G`:

```text
in_min-lambda(G)=c_- y^615 z^271,                       (36)

in_min-beta(G)=c_beta y^615 z^172
 (y^2+27z)^66 (y^2+108z)^33,                            (37)
```

with `c_-,c_beta` nonzero.  In particular `min beta(G)=-1157`.

Thus the proposed scalar laws

```text
e_next=6e+1,       m_next=2e+1                          (38)
```

are refuted at their first untested point.  They held for the first two
steps because the hidden second coordinate happened to take the values
`m=0,3`; the actual state transform is (12).

## 6. The next old-boundary valuation

THM-3498 gives the generic divergent inverse sheets over `(L)`:

```text
x=u^-1,
y=D/S+O(u),
z=-3(D/S)u+O(u^2),
3xz-2y=-11D/S+O(u),                                    (39)
```

where `D/S` is a unit and `v_L(u)=1/2`.  Equation (34) therefore gives
exact valuation `-271/2` on each divergent sheet.

The finite sheet must still be checked: a zero there would add a positive
valuation and spoil exactness.  At the canonical point

```text
(a,b,c)=(2/27,1,1) on L=0,
q=(2,5/6,-7/8),        F(q)=(2/27,1,1),                 (40)
```

the companion evaluates `G(q)=L(q)^43N(J)(q)` in three independent good
reductions.  The exact ledgers are

```text
prime     N(J)(q)     L(q)     G(q)     cubic discriminant
 101         39        16       27              56
 103         89        12       98              91
 107         51        38       27              14.     (41)
```

Every denominator and cubic discriminant is a unit.  Nonzero reduction at
even one good prime proves that the rational value is nonzero; all three are
retained as hostile controls.  Hence the finite sheet is a generic unit and

```text
v_L(N(G))=-271.                                         (42)
```

Since the inverse cover is finite etale over `Q[a,b,c,L^-1]`, the same UFD
localization argument as THM-3498 gives

```text
R_5:=L^271N(G) in Q[a,b,c],       gcd(R_5,L)=1.          (43)
```

Finally, (36)--(37) and the top-face part of Section 2 require no renewal
faces.  They give one further exact exposed face:

```text
in_max-lambda(R_5)
 =C_5 x^1699(3xz-2y)^615,       C_5!=0,                 (44)

(1699,615)=(7*271-2*99,3*271-2*99).                    (45)
```

Equation (44) does not by itself prove `v_L(N(R_5))=-1699`; that requires a
new finite-sheet unit test.  THM-3521 subsequently supplies that test and
proves the valuation.  Neither (43) nor THM-3521 proves that `R_5` is
irreducible or an image equation.

## 7. Failure modes and induction boundary

The exact transform identifies five distinct ways an extrapolation can
fail.

1. **Linear-edge cancellation.**  The complete `min lambda` face can vanish
   after (15), changing both the top weight and the power of `3xz-2y`.
2. **Quadratic-resultant cancellation.**  The `min beta` face can vanish in
   `Q(A)[W]/(4AW^2+1)`.  A visible face monomial is insufficient; its
   Gaussian norm must be nonzero.
3. **Chebyshev-end cancellation.**  A nonmonomial `z`-top face can vanish at
   either label `Bw=2` or `Bw=-1`, changing (10).
4. **Gamma-face cancellation.**  The complete `min gamma` face can vanish
   after (27), destroying one of the factors in (11).
5. **Finite-sheet vanishing.**  Even a correct exposed face controls only
   the two divergent old-`L` sheets.  A zero on the finite sheet changes the
   exact norm valuation and denominator clearing.

There are also two algebraic bookkeeping hazards.  A raw resultant for the
nonmonic cubic differs from the function-field norm by a leading-coefficient
power, and polynomiality of `L^eN(P)` must come from the finite-cover
localization plus an exact old-boundary valuation.  Neither may be inferred
from a face picture alone.

This theorem proves the three transported faces (34), (36), (37), but not
the two renewal faces of `G`.  THM-3513 subsequently derives both by two
hybrid Newton limits, so the fixed `G` has the complete packet `A(271,99)`.
THM-3522 subsequently proves that a complete packet renews through this
fixed inverse chart whenever the next cleared norm is polynomial.  Together
with THM-3521's following finite-sheet gate, it gives complete packets
`A(1699,615)` for `R_5` and `A(10663,3867)` for `R_6`.  THM-3523 subsequently
closes `L^10663N(R_6)` and gives `A(66907,24255)` for `R_7`.  THM-3528 later
proves polynomiality at every raw rung; a separate finite-sheet unit remains
necessary for `L`-coprimality and image-prime arguments.

## 8. Monoid comparison and exact scope

Composition degree and this Newton state use different operations.  For the
fixed Keller map, fibre degrees multiply under composition, giving `3^n`.
By contrast, `(e,m)` is an affine boundary packet transported by a cubic
function-field norm; its update is the signed matrix (12), and its validity
depends on face resultants and finite-sheet units.  Norms are multiplicative
on functions, but admissible Newton packets are not known to be closed under
norm or product.

Thus (12) is compatible grammar for one fixed composition orbit, not a
degree law for the Keller monoid.  It proves no statement about arbitrary
Keller maps, atom degrees, `KDeg(n)`, `JC(2)`, `DC(2)`, or LRC.  It also does
not prove a fifth nonproperness component, a fifth image prime, degree-`243`
separability, or a fifth discriminant square class.

Reproduce the exact packet and finite-sheet gate with

```text
python 04-computation/keller_tropical_norm_face_recurrence_probe_20260816.py
python -O 04-computation/keller_tropical_norm_face_recurrence_probe_20260816.py
```

**QED.**

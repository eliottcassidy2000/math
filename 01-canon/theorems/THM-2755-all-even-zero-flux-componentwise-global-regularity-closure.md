---
id: THM-2755
title: "All-even zero-flux componentwise global-regularity closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The all-even,
  zero-first-flux residual of the chosen-sheet split polynomial exact-square-
  prefix degree-22 Faber response family is physically empty, including
  vertical, reducible, and nonreduced members.  The five smooth infinity
  points fail the polynomial exact-prefix resultant.  At P_infty, coprime
  residual top faces and the common slope-four weight force
  min(v(q),v(s))>=4v(h), so q/h^3 extends regularly and vanishes on every
  boundary branch of a physical component.  Projective global regularity
  forces q=0, while every even third response is q-divisible, contradicting
  U R_Q'=kappa.  With THM-2719, THM-2725, and THM-2745 this closes the full
  chosen-sheet split polynomial exact-prefix reduced-degree-22 family, not
  the broader split branch, JC(2), or DC(2).
source: root/all-even-zero-flux-componentwise-regularity-closure-2026-07-28
audit: >
  componentwise-global-regularity-DVR-2026-07-28 (independent multinomial
  reconstruction, degenerate/vertical/component audit, false-shortcut
  correction, and normal/-O replay: ACCEPT); all-even-zero-flux-second-
  hostile-2026-07-28 (independent generalized-multinomial, parity,
  slope-four DVR/index-quotient, smooth-prefix, normalization, and
  reducible/nonreduced audit: ACCEPT); degree22-frontier-scope-audit-2026-07-28
  (independent coefficient-bank exhaustion, chosen-sheet typing, and
  upstream/downstream frontier audit: ACCEPT)
depends_on:
  - THM-2719-full-split-odd-faber-generic-normalization-genus-four-hundred-nineteen
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
  - THM-2725-split-even-faber-nonzero-first-flux-unified-closure
  - THM-2745-highest-odd-faber-componentwise-exact-prefix-closure
related:
  - THM-2747-highest-odd-reduced-boundary-divisor-and-one-ended-factorization-atlas
  - THM-2752-all-even-zero-first-flux-response-regularization-closure
script: 04-computation/jc2_degree22_all_even_zero_flux_componentwise_closure_thm2755.py
output: 05-knowledge/results/jc2_degree22_all_even_zero_flux_componentwise_closure_thm2755.out
script_sha256: e4c47380efd336a4cb054c499e230c0dbf869a21321b3b48234c817608787b95
output_sha256: 3671979838572dcde99f0cd122f5e88cdc530f4345e71bc532939ed68dff4190
hash_basis: LF-normalized bytes
---

# THM-2755 -- the last all-even degree-22 edge is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

**Scope:** the chosen-sheet, split, polynomial exact-square-prefix,
constant-coefficient reduced Faber family of degree `22`, with every odd
Faber coefficient zero and chosen-sheet first flux `lambda=0`.  Normalize
`a_22=1`.  The exact target-shear quotient removes the coefficients
`a_(4k)`, and the legal translation `P -> P+c` kills `a_18`.
THM-2719's remaining bank is therefore precisely the four even coefficients
`B,C,D,E` and the eleven odd coefficients.  This is the last edge left open by
`THM-2725-split-even-faber-nonzero-first-flux-unified-closure` and
`THM-2745-highest-odd-faber-componentwise-exact-prefix-closure`.  It is not a
derivation of splitness or reduced degree `22`: THM-2202 derives the
polynomial exact prefix only after entering the quartic, twice-odd terminal
reduction, and no theorem places an arbitrary Keller pair in a quartic
source-fibre chart or forces this reduced degree.  Once splitness is given,
choosing `U` chooses a sheet presentation and replacing `U` by `-U` flips
the odd bank; uniformity in all odd coefficients closes both presentations.

The broader split branch remains open in reduced degrees `6,10,14,18` and
all degrees at least `26` (degree `2` is tame), and this theorem supplies no
degree-raising or degree-descent map.  Genuine-nonsplit degrees at least
`26`, other source-fibre/Newton/Jelonek branches, the refuted canonical
first-grade `2`-adic raising route, `JC(2)`, and `DC(2)` also remain open.

The proof is deliberately componentwise.  It does not assume that the whole
response complete intersection is irreducible or reduced.

## 1. Family, parity, and the two component types

Work over `C` in

```text
P(1,2,3,4)_[h,d,q,s].
```

Put `c_22=1`, `c_14=B`, `c_10=C`, `c_6=D`, and `c_2=E`.  The all-even
zero-first-flux member is

```text
F_23 = sum_(j=22,14,10,6,2) c_j h^(22-j) Phi_j,

G_24 = sum_(j=22,14,10,6,2) c_j h^(22-j) Psi_j - W h^24,

R_25 = sum_(j=22,14,10,6,2) c_j h^(22-j) R_j.          (1)
```

The quartic Faber recurrence has the exact parity

```text
Phi_j(d,-q,s) = -Phi_j(d,q,s),
Psi_j(d,-q,s) =  Psi_j(d,q,s),
R_j(d,-q,s)   = -R_j(d,q,s)                            (2)
```

for every even `j`.  Consequently there are homogeneous polynomials `A_20`
and `T_22` such that

```text
F_23=q A_20,                         R_25=q T_22.       (3)
```

This factorization is global, not merely a tangent-cone statement.

At `h=0`, after harmless nonzero scalings, the residual top forms are

```text
A_infty = -84 d^2 q^4 s + 3 d q^6 + 280 d q^2 s^3
          -21 q^4 s^2 -84 s^5,

G_infty = -224 d^3 q^6 +3360 d^2 q^4 s^2 -336 d q^6 s
          -3360 d q^2 s^4 +3 q^8 +560 q^4 s^3+224 s^6.
                                                                  (4)
```

Exact polynomial gcds give

```text
gcd(A_infty,G_infty)=1,              gcd(q,G_infty)=1. (5)
```

Indeed, in the stripped normalization of `(4)`, on `q=0`,

```text
A_infty=-84s^5,                       G_infty=224s^6.  (6)
```

The companion also checks the equivalent recurrence-normalized values
`-(693/128)s^5` and `(231/256)s^6`.  Therefore neither
`(A_20,G_24)` nor `(q,G_24)` has a common hypersurface factor for any values
of `B,C,D,E,W`: a positive-weight common factor would restrict at `h=0` to a
common factor in `(5)`, while a factor whose restriction is zero would be
divisible by `h`, impossible because both top forms are nonzero.  Thus the
full member `(qA_20,G_24)` has no surface component.  Both possible generic
component types are pure projective curves.

Let a physical polynomial Keller trajectory be given.  The third Hamiltonian
identity from
`THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity` is

```text
U (R_source)' = kappa,                         kappa != 0.          (7)
```

It makes the coefficient map nonconstant.  The reduced source kills all
target nilpotents.  The closure of its generic image is consequently a
reduced integral curve component, and its function field is a domain.  From
`qA_20=0`, exactly one of the following applies generically:

```text
q=0,                                             (vertical type)
A_20=0 and q!=0.                                 (residual type)  (8)
```

The vertical type is immediately impossible: `(3)` gives `R_source=0`, in
contradiction with `(7)`.  It remains to exclude a residual component.

## 2. Complete residual support at infinity

The residual boundary is

```text
A_infty=G_infty=0 in P(2,3,4)_[d,q,s].                 (9)
```

Equation `(6)` shows that its `q=0` support consists only of

```text
P_infty=[d:q:s]=[1:0:0].                              (10)
```

On the `q=1` index chart, use the coarse `mu_3` invariants

```text
rho=d/s^2,                         z=s^3.              (11)
```

After removing nonzero monomials, `(9)` becomes

```text
Pbar=3rho-21+z(-84rho^2+280rho-84),

Qbar=3+z(-336rho+560)
       +z^2(-224rho^3+3360rho^2-3360rho+224).          (12)
```

Its lexicographic basis consists of one linear equation in `rho` and

```text
p5(z)=20141047808 z^5-14386462720 z^4+1089822720 z^3
      -21288960 z^2-35910 z+81.                       (13)
```

The exact ideals obtained by adjoining `z`, `rho`, or the Jacobian of `(12)`
are units, and `gcd(p5,p5')=1`.  Hence `(13)` gives exactly five distinct
coarse points, all with

```text
q != 0,                       s != 0,                 d != 0,       (14)
```

and the residual intersection is transverse there.  In particular `h` is a
local parameter on every residual curve germ through one of these points.
Together with `(10)`, this is the complete set-theoretic support of `(9)`.

The stripped top response on the chart is

```text
Tbar=3+z(9rho^2-105rho+70)
       +z^2(-84rho^3+280rho^2-84rho).                 (15)
```

The ideal `(Pbar,Qbar,Tbar)` is the unit ideal.  Thus `R_22` is nonzero at
all five smooth points, and `R_aff=R_25/h^25` has pole order `25` there.

## 3. Universal exclusion of every smooth infinity point

Let `Y` be the normalization of a residual component carrying the physical
image.  The rational source map extends, by properness, to a nonconstant
morphism

```text
gamma:P1_x -> Y;                                      (16)
```

it is finite and surjective.  Suppose `Y` contains one of the five smooth
points.  At every source point above it, the pullback of `R_aff` has pole
order `25e>1`, where `e` is the ramification index.  The rational-primitive
classification in THM-2723 therefore excludes constant `U`, because that
case makes `R_source` affine-linear with pole order one.  The unique response
pole on the source is instead a finite zero `x=a` of

```text
U=u_0(x-a)^m,                          m>=2.           (17)
```

The polynomial exact-square prefix supplies more information than the
abstract response curve.  Write

```text
P_source=H^2+L,
H=U^2 z_source^2+B_src z_source+C_src,
L=A_src z_source+E_src,

w=U z_source+B_src/(2U),
H=w^2+d_aff,                    L=q_aff w-s_aff.       (18)
```

Pull the homogeneous response coordinates to the source DVR at `a`, passing
to a finite weighted-index cover if needed, and put

```text
beta=h B_src/(2U).                                    (19)
```

Direct coefficient comparison in `(18)` gives the two exact identities

```text
beta^2+d=h^2 C_src,
q beta-s=h^4 E_src.                                   (20)
```

All source coefficients on the right are regular at `a`.  Since `d` is
regular, the first identity makes `beta` regular in the DVR.  If `omega` is
its residue, then at a smooth infinity point `(20)` reduces to

```text
d=-omega^2,                       s=q omega.           (21)
```

Here `q,d,s` are all nonzero, so `omega!=0`.  Substitution in the two top
forms factors them as

```text
A_infty=-8 omega^2 q^5(56 omega^3+3q),

G_infty=q^6(7168 omega^6+896 omega^3 q+3q^2).          (22)
```

After removing the nonzero factors, their resultant in `q` is

```text
Res_q(56 omega^3+3q,
      7168 omega^6+896 omega^3 q+3q^2)
    =-76608 omega^6 != 0.                              (23)
```

This contradiction is coefficient-independent.  Therefore a physical
residual component contains none of the five smooth boundary points.

## 4. The uniform `P_infty` valuation lemma

It remains to understand every normalization branch above `(10)`, including
all exceptional parameter values.  Work on a finite `d=1` index cover and
let `v` be the branch valuation.  Set

```text
a=v(h)>0,
b=min(v(q),v(s))>0,                                   (24)
```

with `v(0)=infinity`.  Since this is the residual type, `q` is not identically
zero and `b` is finite.

For a row of Faber degree `j`, let `gap=22-j`.  The exact minimum ordinary
degrees in `(q,s)` at `d=1` are:

```text
 j    gap    ord(Phi_j/q)    ord(Psi_j)    ord(R_j)
22      0          5              6             6
14      8          3              4             4
10     12          2              3             3
 6     16          1              2             2
 2     20          0              1             1.     (25)
```

For the first two equations, every lower row ties the degree-22 row exactly
at the slope `b=4a`:

```text
gap*a+ord(Phi_j/q)*4a = 20a,
gap*a+ord(Psi_j)*4a   = 24a.                           (26)
```

The term `W h^24` also ties the second equation at `24a`.  Therefore, if
`b<4a`, the degree-22 lowest faces strictly precede every lower row, every
higher `(q,s)`-degree term in the degree-22 row, and `W h^24`.  The two
coefficients at orders `5b` and `6b` would have to vanish simultaneously at
the leading projective direction `[q_*:s_*]`.  They are

```text
A_5=-(231/128)s(q^2-3s^2)(3q^2-s^2),

G_6=-(231/256)(q^2-s^2)(q^4-14q^2s^2+s^4).            (27)
```

But

```text
gcd(A_5,G_6)=1.                                       (28)
```

They have no common point of `P1_[q:s]`, including the directions where one
coordinate vanishes.  This contradiction proves the uniform inequality

```text
min(v(q),v(s)) >= 4v(h).                              (29)
```

No genericity, reducedness of the total member, slope-four chart, or
assumption `ord(h)=1` entered the proof.  In particular, on every residual
normalization branch over `P_infty`,

```text
v(q_aff)=v(q/h^3)=v(q)-3v(h) >= v(h)>0.                (30)
```

The calculation was made on a finite index cover; finite ramification scales
all valuations by the same positive integer, so regularity and strict
positivity descend to the coarse normalization.

## 5. Global regularity closes every residual component

Every projective curve component meets `h=0`.  Indeed, a component disjoint
from `h=0` would be both projective and contained in the affine chart
`h!=0`, hence would be zero-dimensional.  Equations `(5)` also exclude a
curve component contained inside `h=0`.

By Section 3, every boundary point of a physical residual component lies
over `P_infty`.  On the affine chart `h!=0`, the function

```text
q_aff=q/h^3                                             (31)
```

is a regular coordinate.  Its pullback to the normalization remains regular
there.  Equation `(30)` proves that it extends regularly, and in fact
vanishes, at every omitted boundary point.  It is therefore a global regular
function on the complete integral curve `Y`.

A global regular function on a complete integral complex curve is constant.
Since `(31)` vanishes at the nonempty boundary, it is identically zero.  This
contradicts the residual assumption `q!=0`; equivalently, `(3)` would give
`R_source=0`, contradicting the third flux `(7)`.  No residual physical
component exists.

Together with the vertical contradiction in Section 1, this closes every
all-even `lambda=0` member in the stated chosen-sheet family.

## 6. Reducible, nonreduced, and vertical boundaries

The proof never assumes that `(F_23,G_24)` is integral.

* A map from the reduced source factors through the target reduction.
* Its nonconstant generic image has a reduced irreducible curve closure;
  because the ambient complete intersection is pure of dimension one, that
  closure is an irreducible component.
* In that component's function field, `qA_20=0` gives the exhaustive
  vertical/residual dichotomy `(8)`.
* Nilpotents and embedded structure cannot alter the source identity `(7)`,
  the polynomial coefficient identities `(20)`, or a valuation on the
  normalized reduced component.
* The gcds `(5)` prevent a hidden common surface or a curve trapped in the
  boundary.  Every projective component therefore acquires at least one of
  the boundary points already audited.
* A component on which `q=0`, including the extreme locus `q=s=0`, has
  `R_25=0` by global parity and is killed directly by `(7)`.

Thus reducibility, generic nonreducedness, embedded points, and vertical
components supply no escape.

## 7. Corrections and sharp failure boundaries

Two tempting shortcuts are false and are not used above.

### 7.1 The exceptional slope-four response cancels

At equality in `(29)`, write `q=h^4 Q` and `s=h^4 S` on the `d=1` cover.
All five even rows enter the same leading face.  If `A_ex` is the resulting
leading face of `A_20`, the exact response face is

```text
R_ex=-(Q/2) A_ex.                                     (32)
```

Hence `R_ex` vanishes identically on the residual leading equation
`A_ex=0`.  The apparent order-24 numerator, and therefore the apparent simple
pole of `R_25/h^25`, is not a uniform response-pole certificate.  Degenerate
parameters can also produce the literal vertical locus `q=s=0`, where every
even response vanishes.  Section 5 uses `q_aff`, not this cancelled response
face.  Independently, THM-2752 starts from this same cancellation, proves the
next exact Faber order, and shows that `R_25/h^25` is regular at `P_infty`;
thus the false shortcut is the claimed leading simple pole, not the
response-regularization mechanism.

### 7.2 A simple finite response pole does not force constant `U`

The rational-primitive equation itself has the hostile example

```text
U=(x-a)^2,                  R=-1/(x-a),
U R'=1.                                                (33)
```

Thus a simple finite pole is compatible with nonconstant `U`.  Only a pole
of order greater than one excludes the constant-`U` alternative; Section 3
uses the exact smooth-point order `25e`, while the `P_infty` closure avoids a
response-pole assertion entirely.

## 8. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree22_all_even_zero_flux_componentwise_closure_thm2755.py
python3 -O 04-computation/jc2_degree22_all_even_zero_flux_componentwise_closure_thm2755.py
```

The two transcripts byte-match.  The companion reconstructs Faber rows in
two independent ways (Newton recurrence and generalized multinomial
extraction), checks every parity and `q` factor, the two top gcds, the exact
five-point quotient, smoothness and response nonvanishing, table `(25)`, the
slope-four margins, gcd `(28)`, resultant `(23)`, cancellation `(32)`, and
hostile `(33)`.

```text
script_sha256 = e4c47380efd336a4cb054c499e230c0dbf869a21321b3b48234c817608787b95
output_sha256 = 3671979838572dcde99f0cd122f5e88cdc530f4345e71bc532939ed68dff4190
hash_basis    = LF-normalized bytes
```

Two independent hostile audits reconstructed the Faber rows by generalized
multinomial extraction rather than trusting the recurrence, checked the
weighted index quotient and every component boundary, and reproduced the
normal, optimized, and stored transcript bytes.  They found and repaired the
false simple-pole shortcut recorded in Section 7 and found no remaining
hypothesis gap.  QED.

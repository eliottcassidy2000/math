---
id: THM-3522
title: "Fixed Keller five-face renewal propagation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For the fixed sporadic Keller inverse chart, if
  P has the complete five-face packet A(e,m) and Q=L^e N(P) is polynomial,
  then Q automatically has the complete packet
  A(7e-2m,3e-2m).  In particular, the two renewal faces propagate; no
  separate finite-sheet condition enters that implication.  Hence the fixed
  R_5 and R_6 have full packets A(1699,615) and A(10663,3867), with explicit
  nonzero renewal scalars below.  THM-3523 and THM-3527 subsequently clear
  the next two rungs and apply this theorem to give A(66907,24255) for R_7
  and A(419839,152211) for R_8; THM-3525/3526 close the distinct degree-243/
  729 separability gates.  THM-3528 subsequently discharges polynomiality
  for every complete packet and proves the fixed raw packet orbit at all
  levels; THM-3529 subsequently proves all complete packets are finite-sheet
  units, so every positive raw rung is L-coprime.  Image status, arbitrary
  Keller maps, and every general Jacobian claim remain open at this theorem
  state.
source: codex/five-face-renewal-propagation/2026-08-16
depends_on:
  - THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness
  - THM-3506-fixed-keller-five-face-norm-transform-and-271-99-boundary
  - THM-3513-fixed-G-hybrid-newton-renewal-faces
  - THM-3521-fixed-R5-finite-sheet-unit-and-next-old-L-clearing
related:
  - MISTAKE-415
  - THM-3523-fixed-R6-finite-sheet-unit-and-next-old-L-clearing
  - THM-3525-level-five-degree243-separability-and-discriminant-square-class
  - THM-3526-level-six-degree729-separability-and-discriminant-square-class
  - THM-3527-fixed-R7-finite-sheet-unit-and-next-old-L-clearing
  - THM-3528-fixed-keller-all-level-cleared-norm-polynomiality-and-finite-sheet-defect
  - THM-3529-fixed-keller-complete-packet-finite-sheet-unit
scripts:
  - 04-computation/keller_five_face_renewal_propagation_probe_20260816.py
  - 04-computation/keller_five_face_renewal_propagation_independent_audit_20260816.py
  - 04-computation/keller_five_face_renewal_vieta_independent_audit_20260816.py
outputs:
  - 05-knowledge/results/keller_five_face_renewal_propagation_probe_20260816.out
  - 05-knowledge/results/keller_five_face_renewal_propagation_independent_audit_20260816.out
  - 05-knowledge/results/keller_five_face_renewal_vieta_independent_audit_20260816.out
script_sha256:
  - 5fd3c27bf49f8fab5e96b3d3fe608b91b86151ca265d35878225f8ce6aa2f05e
  - 2fafef6bca64174b71fac22aa89e5bb713847ccc47318635c61b9684c5b063cd
  - e2651175dbe4ef21553c5d4a7949ac5ba19ddb50d31183cce322d9cc74f34770
output_sha256:
  - 9a3eba81a00e35a5c99c2285e8ea7a10b6b1764b6c8d56b4d5276c449bc1b8b8
  - 9485d14a9eb764374ca6abb25fdecafc75a210ce720bdcf4c0f2e7d0536d8dc8
  - 91c05274cf307ddd6d341b55b8d94616aaa4b8278792b2e5a85fb97dba38a1cd
semantic_sha256:
  - 8b6a447c98e4e7f6bfc493818696d4a9193b4da47ab7b2f9e0368e9155940a91
  - 9de2b0a149105263ee1b3a1fba01424f9c7ff274368c689cdc5737fb340bf804
  - 7468b48f87fab27aacf77af954ea5b63189eb8ba22ec1402b9df57d6d1378375
hash_basis: raw LF bytes; exact large scalars are pinned by their prime exponents
---

# THM-3522 -- the two renewal faces propagate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

Retain the fixed sporadic Keller map, target polynomial `L`, inverse cubic,
and cubic function-field norm `N` of THM-3495.  Thus, temporarily writing the
target coordinates as `(a,b,c)`,

```text
L=27a^2c^2-18abc+16a+b^3c-b^2,
E(w)=Lw^3+(4-3bc)w-2c.                                (1)
```

Let `e,m` satisfy the packet hypotheses of THM-3506: they are nonnegative,
`3|m`, and all displayed packet exponents are nonnegative.  For a source
monomial `x^i y^j z^k`, put

```text
gamma=i-j-5k.                                          (2)
```

It is enough for the renewal argument to assume that the following two
complete faces of `P` hold:

```text
r=e-2m/3,       p=2e-4m/3=2r,
d=2e-2m/3=e+r,

in_max-k(P)=kappa x^p z^d,                             (3)
in_min-gamma(P)=eta z^e(27x^2z+y^3)^r,                (4)
```

where `kappa,eta` are nonzero rationals.  Put

```text
Q=L^e N(P),
e'=7e-2m,              m'=3e-2m,                      (5)
r'=e'-2m'/3,           p'=2e'-4m'/3,
d'=2e'-2m'/3,          gamma'=-8e'+2m'.               (6)
```

If `Q` is polynomial, then its two complete renewal faces are

```text
in_max-k(Q)=kappa' x^p' z^d',                          (7)

in_min-gamma(Q)
 =eta' z^e'(27x^2z+y^3)^r',                           (8)
```

with the explicit nonzero scalars

```text
rho=p-3d=-4e+2m/3,
eta'   =kappa^3 2^rho,
kappa' =kappa^3 2^rho 27^r'.                          (9)
```

Negative powers in (9) are rational powers, not Laurent variables.  In
particular, if `P` has the full packet `A(e,m)`, THM-3506 supplies the other
three complete faces of `Q`, and (7)--(8) give

```text
A(e,m) --L^e N(-), polynomial--> A(7e-2m,3e-2m).      (10)
```

This is a one-step closure theorem for the fixed inverse chart.  It is not a
claim that the polynomiality hypothesis always holds.

## 1. The common endpoint is forced and scalar-compatible

Expanding (4), its terms are indexed by `0<=ell<=r`:

```text
x^(2ell) y^(3(r-ell)) z^(e+ell).                       (11)
```

The unique term of largest `z`-degree occurs at `ell=r`; it is

```text
x^(2r)z^(e+r)=x^p z^d.                                (12)
```

Thus the faces in (3)--(4) meet exactly at one monomial.  Since they are
faces of the same polynomial, their apparently independent normalizing
scalars obey the load-bearing compatibility

```text
kappa=eta 27^r.                                        (13)
```

This resolves the equal-weight issue: the scalars in the definition of
`A(e,m)` need not be numerically equal, but they cannot be chosen
inconsistently at an overlap of complete faces.

Now define

```text
delta_6=i-j-6k=gamma-k,
delta_8=i-j-8k=gamma-3k.                               (14)
```

Every monomial of `P` has `gamma>=-8e+2m` and `k<=d`.  If

```text
u=gamma-(-8e+2m)>=0,       v=d-k>=0,                  (15)
```

then its gaps above the two lower bounds are respectively

```text
delta_6-[(-8e+2m)-d]=u+v,
delta_8-[(-8e+2m)-3d]=u+3v.                           (16)
```

Equality in either bound therefore forces `u=v=0`.  The complete `k=d`
face is the singleton (3), so both hybrid initial forms consist of exactly
`kappa x^p z^d`.  Explicitly,

```text
min delta_6=-10e+8m/3,
min delta_8=-14e+4m.                                  (17)
```

No unlisted equal-weight source term can enter either norm limit.

## 2. The maximum-`z` target face

Use the first THM-3513 hybrid scaling

```text
(a,b,c)=(A,B,Cs^3),       w=q/s,       s -> infinity. (18)
```

Exact weighted extraction from THM-3495's inverse numerators gives the
generic residual cubic and inverse coordinates

```text
27A^2Cq^3-2=0,
(q_x,q_y,q_z)~(q/s,-s/q,-Cs^6/q^3).                  (19)
```

The residual cubic is separable on the generic open `AC!=0`, and none of its
roots is zero.  By (3), (16), and (19), each inverse branch has the exact
leading term

```text
P(q_x,q_y,q_z)
 ~kappa C^d q^rho s^(-p+6d).                          (20)
```

There is no sign in (20), because

```text
d=2(e-m/3)
```

is even.  Branchwise root-independence is generally false: for `P=L` and
`P=H`, the exponents `rho` are `-4` and `-26`, neither divisible by three.
The norm only needs the product of the three labels.  The cubic in (19) is
nonmonic, so Vieta gives

```text
product(q)=2/(27A^2C).                                 (21)
```

Together with

```text
L^e~(27A^2C^2)^e s^(6e),                              (22)
```

equations (20)--(22) give

```text
Q(A,B,Cs^3)
 ~kappa^3 2^rho 27^(e-rho)
   A^(2e-2rho) C^(2e+3d-rho)
   s^[6e+3(-p+6d)].                                   (23)
```

Direct simplification using (5)--(6) gives

```text
e-rho=r',       2e-2rho=p',
2e+3d-rho=d',   6e+3(-p+6d)=3d'.                      (24)
```

Because `Q` is polynomial, substitution (18) grades its monomials by three
times their `c`-degree.  Equation (23) is its complete generic leading
coefficient in `A,B,C`; it is a single monomial and is nonzero.  Equations
(23)--(24) prove (7) and the scalar `kappa'` in (9), including completeness.

## 3. The minimum-`gamma` target face and nonmonic Vieta factor

Use the second hybrid scaling

```text
(a,b,c)=(At,B/t,C/t^5),       w=qt,       t -> 0,
D=27A^2C+B^3.                                         (25)
```

Exact weighted extraction gives

```text
Dq^3-3Bq-2=0,
(q_x,q_y,q_z)~(qt,-1/(qt),-C/(q^3t^8)).               (26)
```

The discriminant of the residual cubic is

```text
-2916 A^2 C D,                                         (27)
```

so it is generically separable and again has no zero root.  The unique
`delta_8` initial monomial gives, on each branch,

```text
P(q_x,q_y,q_z)
 ~kappa C^d q^rho t^(p-8d).                           (28)
```

The leading coefficient `D` in (26) cannot be dropped.  Its raw resultant
with `q` is `2`, while the function-field norm is

```text
product(q)=2/D.                                        (29)
```

Equivalently, (29) is the determinant of multiplication by `q` in the monic
three-dimensional quotient.  Since

```text
L^e~(CD)^e t^(-8e),                                   (30)
```

taking the norm in (28) yields

```text
Q(At,B/t,C/t^5)
 ~kappa^3 2^rho C^(e+3d) D^(e-rho)
   t^[-8e+3(p-8d)].                                   (31)
```

The exponent identities are

```text
e+3d=e',       e-rho=r',
-8e+3(p-8d)=-8e'+2m'=gamma'.                          (32)
```

For polynomial `Q`, the power of `t` is exactly the target weight
`i-j-5k`.  Therefore (31)--(32) are the complete initial form (8), with
scalar `eta'` in (9).  The common endpoint coefficient of (8) is

```text
eta' 27^r'=kappa',                                    (33)
```

exactly matching (7).  Thus the two output faces are not merely separately
nonzero; their overlap is automatically coefficient-compatible.

The nonmonic factor is decisive.  Treating (26) as monic without dividing by
`D` would leave exponent `e` instead of `e-rho=r'`.  At the calibrated
`J -> G` step this would give `43` rather than `205`, immediately
contradicting THM-3513.

## 4. Polynomiality, denominators, and the finite sheet

The hybrid calculation takes place over the generic residual fields in
(19) and (26).  Their constant term is `-2`, so every occurrence of
`q^rho`, including negative `rho`, is a unit.  Their leading coefficients
are retained by (21) and (29), and the apparent denominators cancel into
the nonnegative output exponents in (24) and (32).  Packet admissibility
implies `m<=e`; hence `e',m',r',p',d'` are all nonnegative and `3|m'`.

Polynomiality of `Q` is nevertheless a genuine global hypothesis.  It is
what turns the generic Laurent asymptotics (23) and (31) into complete faces
of a polynomial in the target coordinates.  This theorem does not derive
polynomiality from a face picture.  THM-3528 subsequently derives it from
finite-etale norm regularity and the nonnegative finite-sheet defect.

No additional finite-sheet unit condition is used in (18)--(33).  A
finite-sheet test at the old divisor `(L)` determines the exact `L`-adic
valuation and coprimality after clearing; it is a different gate.  Once
polynomiality is already known, renewal propagation needs neither a finite
specialization nor a good-reduction witness.  For the fixed tower, THM-3521
supplies precisely the finite-sheet input needed at the `R_5` rung; its files
are not altered here.

## 5. Exact fixed-tower consequences

THM-3513 proves that the fixed

```text
G=L^43N(J)
```

has packet `A(271,99)` and maximum-`z` coefficient

```text
kappa_G=3^1128/2^117.                                  (34)
```

THM-3506 proves that `R_5=L^271N(G)` is polynomial.  Applying (7)--(10)
therefore gives the full packet

```text
R_5 has A(1699,615),                                   (35)
```

with its two newly proved complete faces

```text
in_max-k(R_5)
 =(3^7251/2^1369)x^2578z^2988,                        (36)

in_min-gamma(R_5)
 =(3^3384/2^1369)
   z^1699(27x^2z+y^3)^1289,                           (37)

min gamma(R_5)=-12362.                                 (38)
```

THM-3521 proves that

```text
R_6=L^1699N(R_5)
```

is polynomial.  Applying the closure theorem again gives

```text
R_6 has A(10663,3867),                                 (39)
```

and the complete renewal faces

```text
in_max-k(R_6)
 =(3^46008/2^10493)x^16170z^18748,                    (40)

in_min-gamma(R_6)
 =(3^21753/2^10493)
   z^10663(27x^2z+y^3)^8085,                          (41)

min gamma(R_6)=-77570.                                 (42)
```

The first three faces in (35) and (39) come from THM-3506; equations
(36)--(37) and (40)--(41) supply the two renewal faces.  Thus renewal is no
longer an independent open gate at these two rungs.

## 6. Exact audit and boundary

The primary exact companion reconstructs both hybrid initial systems directly from
THM-3495's inverse numerators.  It verifies the residual coordinate formulas
modulo both cubics, their generic discriminants, the raw resultants, and the
multiplication-matrix norms.  It checks all symbolic exponent identities,
the boundary packet `(e,m)=(3,3)`, and the exact scalar calibrations

```text
L -> H/2^6,       H -> J/2^35,       J -> G.           (43)
```

The first two controls have `rho=-4,-26` and therefore hostilely rule out
the unnecessary branchwise root-independence assumption.  Ordinary and
optimized replays agree line-for-line with the stored transcript.

A disjoint implementation starts again from THM-3495's literal inverse
numerators, uses separate weighted-initial and quotient-matrix code, checks
both hybrid gap identities on a `17 by 17` hostile grid, and expands the
literal gamma face for all `2,500` admissible packets with `e<=120`.  It
independently obtains (36)--(42), the normalization controls in (43), and
the two successive Cassini determinants `-1536,12288`.  Its verdict is
`SOUND; promotion recommended at fixed conditional scope`.

A second clean-room Vieta audit checks all `15,250` admissible packets with
`1<=e<=300`; `10,167` of them have `rho` nonzero modulo three.  Split-prime
controls at `109,127,163` verify the product-of-powers identity separately
for the `J`, `G`, and `R_5` input states.  Dropping the nonmonic leading
coefficient fails the required `A,C,D` exponents in every nonzero bank row.

The theorem closes the renewal implication at every one-step application of
the fixed inverse chart for which polynomiality is supplied.  Standing alone,
it does **not** prove:

- polynomiality, an exact old-`L` denominator, or a finite-sheet unit for the
  next norm `L^10663N(R_6)`;
- irreducibility, squarefreeness, or image-equation status of `R_5` or `R_6`;
- a fifth nonproperness component or the degree-`243` separability gate;
- its polynomiality hypothesis or an unconditional all-level packet recurrence;
- a classification of arbitrary Keller maps, `JC(2)`, `DC(2)`, LRC, or any
  general Jacobian-conjecture statement.

The first bullet is an internal boundary of this theorem.  THM-3523
subsequently supplies that finite-sheet and polynomiality input by a separate
exact computation, and then uses the implication proved here to give the
complete packet `A(66907,24255)` for `R_7`.  THM-3525/3526 subsequently close
the degree-243/729 gates, and THM-3527 clears `R_8` with packet
`A(419839,152211)`.  THM-3528 then supplies polynomiality and complete raw
packets at every later rung.  THM-3529 subsequently proves the separate
finite-sheet-unit and all-level old-`L`-coprimality statements.  Image status,
later separability, and general-map bullets stay open.

Reproduce the exact certificate with

```text
python 04-computation/keller_five_face_renewal_propagation_probe_20260816.py
python -O 04-computation/keller_five_face_renewal_propagation_probe_20260816.py
python -B 04-computation/keller_five_face_renewal_propagation_independent_audit_20260816.py
python -B -O 04-computation/keller_five_face_renewal_propagation_independent_audit_20260816.py
python -B 04-computation/keller_five_face_renewal_vieta_independent_audit_20260816.py
python -B -O 04-computation/keller_five_face_renewal_vieta_independent_audit_20260816.py
```

**QED.**

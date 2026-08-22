# Independent hostile audit of fixed-chart five-face renewal propagation

**Verdict: SOUND; SUBSEQUENTLY PROMOTED AS THM-3522 AT THE STATED
CONDITIONAL SCOPE.**

The proposed renewal lemma survives a clean derivation from the literal
inverse-chart numerators.  There is no hidden denominator, missing branch,
equal-weight cancellation, or nonmonic-resultant error.  Once a source
polynomial has a complete packet `A(e,m)` and

```text
Q=L^e Norm(P)
```

is independently known to be polynomial, its two renewal faces propagate.
Together with the first-three-face transform in
[THM-3506](../01-canon/theorems/THM-3506-fixed-keller-five-face-norm-transform-and-271-99-boundary.md),
this gives a complete packet for `Q`.  The fixed consequences are complete
packets `A(1699,615)` for `R_5` and `A(10663,3867)` for `R_6`.

This audit did not itself promote the then-reserved
`THM-3522-fixed-keller-five-face-renewal-propagation` namespace or edit the
candidate files.  Its conclusions were subsequently combined with the
independent Vieta audit and promoted in THM-3522.

## 1. Inheritance and hostile posture

The closest proved mechanisms are:

- [THM-3495](../01-canon/theorems/THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness.md),
  which freezes the cubic inverse chart and its exact numerators;
- [THM-3506](../01-canon/theorems/THM-3506-fixed-keller-five-face-norm-transform-and-271-99-boundary.md),
  which transports the first three faces but explicitly leaves renewal
  conditional;
- [THM-3513](../01-canon/theorems/THM-3513-fixed-G-hybrid-newton-renewal-faces.md),
  which proves the two hybrid limits for the fixed input `J`; and
- [THM-3521](../01-canon/theorems/THM-3521-fixed-R5-finite-sheet-unit-and-next-old-L-clearing.md),
  which supplies polynomiality of `R_6` but explicitly does not claim the
  renewal faces of `R_5`.

The canonical hostile is MISTAKE-415 in
[`MISTAKES.md`](../01-canon/MISTAKES.md): the exposed pair alone is not a
closed scalar state.  The proposed proof does not repeat that extrapolation.
It uses both complete renewal faces, their marked common endpoint, both exact
inverse degenerations, and the separate hypothesis that the cleared norm is
polynomial.

The least-used but load-bearing sidecars are the leading coefficient of the
nonmonic residual cubic and the distinction between renewal and finite-sheet
polynomiality.  Both were attacked directly.

## 2. The common endpoint is forced and unique

For a packet `A(e,m)`, put

```text
a=2e-4m/3,
b=2e-2m/3,
h=e-2m/3,
gamma_0=-8e+2m.
```

The two renewal faces are

```text
in_max-z(P)=c x^a z^b,
in_min-gamma(P)=d z^e(27x^2z+y^3)^h,
```

with `c,d` nonzero.  Since `a=2h` and `b=e+h`, the `ell=h` term of the
second face is

```text
d 27^h x^a z^b.
```

It is the unique term of that binomial face with `z`-degree `b`.  Because the
top-`z` face is the singleton `c x^a z^b`, the intersection of the two
complete faces is exactly this monomial.  In particular, the coefficients
are not independent at the intersection:

```text
c=d 27^h.                                               (1)
```

For any monomial `x^i y^j z^k` of `P`, define

```text
gamma=i-j-5k,
delta_6=gamma-k,
delta_8=gamma-3k.
```

Completeness gives `gamma>=gamma_0` and `k<=b`, hence the exact gap formulas

```text
delta_6-(gamma_0-b)
  =(gamma-gamma_0)+(b-k),

delta_8-(gamma_0-3b)
  =(gamma-gamma_0)+3(b-k).                              (2)
```

Both summands are nonnegative.  Equality in either line forces
`gamma=gamma_0` and `k=b`, hence the unique common endpoint.  This proves
completeness of both hybrid source initial forms and rules out equal-weight
cancellation before any norm is taken.

The independent companion also expands the literal minimum-`gamma` face for
all `2,500` admissible packets with `e<=120`; every intersection is the same
singleton.  That finite sweep is a hostile control, not the proof of (2).

## 3. Both inverse hybrid scalings reconstruct from THM-3495

The audit starts from THM-3495's exact identities

```text
2S q_y=Y0+6Lw-3(9ac-b)Lw^2,
8S q_z=Z0+6LA1w-9LA2w^2.
```

It independently extracts weighted initial forms of all eight polynomials
`L,T,S,9ac-b,Y0,A1,A2,Z0`.

### 3.1 Top-`z` chart

Under

```text
(a,b,c)=(A,B,Cs^3),       w=q/s,       s -> infinity,
```

the inverse cubic has initial equation

```text
C(27A^2Cq^3-2)=0.
```

The exact numerator reductions modulo `27A^2Cq^3-2` give

```text
(q_x,q_y,q_z)
  ~(q/s,-s/q,-Cs^6/q^3),
L~27A^2C^2s^6.                                        (3)
```

The leading `S` is `27AC^2s^6`, so it is a generic unit in the residual
coefficient field.  The apparent inverse-chart denominator therefore does
not add a hidden factor; after reduction, the ratios in (3) are exact.

This scaling assigns `s`-weight `-delta_6` to a source monomial.  Equation
(2) therefore makes `x^a z^b` its unique leading source term on every
residual branch.

### 3.2 Minimum-`gamma` chart

Under

```text
(a,b,c)=(At,B/t,C/t^5),   w=qt,        t -> 0,
D=27A^2C+B^3,
```

the inverse cubic and coordinates reduce to

```text
Dq^3-3Bq-2=0,
(q_x,q_y,q_z)
  ~(qt,-1/(qt),-C/(q^3t^8)),
L~CDt^-8.                                               (4)
```

Here `S~27AC^2t^-9`, again a generic residual unit.  A source monomial has
`t`-weight `delta_8`, so the second line of (2) again leaves the common
endpoint as the unique leading term.

The audit reduces both claimed `q_y` and `q_z` identities modulo their
residual cubics.  It does not import the candidate implementation or the
fixed-`J` face ledger.

## 4. The nonmonic Vieta factor is correct

For the cubic in (4), the raw resultant with `q` is

```text
Res_q(Dq^3-3Bq-2,q)=2.
```

That raw resultant is not the function-field norm of `q`: it retains one
factor of the leading coefficient `D`.  In the monic quotient, the
multiplication-by-`q` matrix on the basis `(1,q,q^2)` has determinant

```text
Norm(q)=2/D.                                            (5)
```

An independent matrix determinant and resultant computation both reproduce
(5).  The top-`z` cubic similarly gives

```text
Norm(q)=2/(27A^2C).                                    (6)
```

The hostile monic substitution `Norm(q)=2` fails immediately: it would give
`D`-exponent `e` rather than `e-rho` below.  On `J -> G`, that would be `43`
rather than the independently proved exponent `205` in THM-3513.

Negative and non-multiple-of-three powers of `q` cause no root-labelling
problem.  The endpoint exponent is an integer and

```text
product_i q_i^rho=(product_i q_i)^rho
```

because every residual root is nonzero.  The earlier exact transitions
`L -> H/64` and `H -> J/2^35` have endpoint exponents `-4` and `-26`, so
they are direct hostile controls against an illicit branchwise cube-root
simplification.

## 5. General exponents and scalars

Set

```text
rho=a-3b=-4e+2m/3,
e'=7e-2m,
m'=3e-2m.
```

Substituting the unique endpoint into (3), multiplying over the three roots,
and multiplying by `L^e` gives

```text
Q(A,B,Cs^3)
 ~c^3 27^e(2/27)^rho
   A^(2e-2rho) C^(2e+3b-rho)
   s^(6e-3a+18b).                                      (7)
```

The packet identities are

```text
2e-2rho       =2e'-4m'/3,
2e+3b-rho     =2e'-2m'/3,
6e-3a+18b     =3(2e'-2m'/3).                           (8)
```

Similarly, (4)--(5) give

```text
Q(At,B/t,C/t^5)
 ~c^3 2^rho C^(e+3b)D^(e-rho)
   t^(-8e+3a-24b),                                    (9)
```

where

```text
e+3b          =e',
e-rho         =e'-2m'/3,
-8e+3a-24b    =-8e'+2m'.                              (10)
```

The sign from `q_z~ -C/q^3` is `(-1)^(3b)`.  It is always `+1` because
`b=2(e-m/3)` is even in every admissible packet.  Thus the complete renewal
faces are

```text
in_max-z(Q)
 =c^3 27^e(2/27)^rho
   x^(2e'-4m'/3)z^(2e'-2m'/3),                        (11)

in_min-gamma(Q)
 =c^3 2^rho
   z^e'(27x^2z+y^3)^(e'-2m'/3).                       (12)
```

Their common endpoint coefficients agree because

```text
e-rho=e'-2m'/3,
c_gamma 27^(e-rho)=c_top.                              (13)
```

Equations (2), (7), and (9) show that these are complete faces, not isolated
surviving terms.

## 6. Is polynomiality alone enough?

Yes, with the source packet and fixed inverse chart understood as hypotheses.
The two generic one-parameter substitutions read exactly the target
`z`-degree and target `gamma=i-j-5k` weight.  If `Q` is a polynomial, their
generic leading expressions are its complete initial forms.  No finite-sheet
unit is additionally needed for renewal.

Polynomiality is not a consequence of this face argument.  Without it, the
same asymptotics describe only a localized rational function and cannot be
promoted to polynomial Newton faces.  In the fixed tower, old-`L` valuation,
finite-sheet nonvanishing, and the UFD localization argument are the separate
mechanism establishing polynomiality.  This distinction is exactly where a
future rung can still fail.

Thus the lawful implication is

```text
complete A(e,m) packet
+ polynomial Q=L^e Norm(P)
=> complete A(7e-2m,3e-2m) packet.                    (14)
```

It is not an unconditional recurrence for arbitrary localized norms.

## 7. Exact normalization controls

The scalar law in (11)--(12) reproduces three pre-existing exact transitions:

```text
L -> H/2^6:
  H top coefficient   =2^2 3^24,
  H gamma coefficient =2^2 3^9;

H -> J/2^35:
  J top coefficient   =2^15 3^171,
  J gamma coefficient =2^15 3^72;

J -> G:
  G top coefficient   =3^1128/2^117,
  G gamma coefficient =3^513/2^117.
```

The last row exactly matches THM-3513.  The first two rows are especially
useful controls because their negative endpoint powers are not multiples of
three.

## 8. Fixed `R_5` and `R_6` consequences

THM-3513 proves `G` has `A(271,99)`, while THM-3506 independently proves

```text
R_5=L^271 Norm(G)
```

is polynomial.  Applying (14) gives the full packet `A(1699,615)` and the
exact renewal faces

```text
in_max-z(R_5)
 =(3^7251/2^1369)x^2578z^2988,

in_min-gamma(R_5)
 =(3^3384/2^1369)
   z^1699(27x^2z+y^3)^1289.                            (15)
```

THM-3521 independently proves

```text
R_6=L^1699 Norm(R_5)
```

is polynomial and `L`-coprime.  The new `R_5` packet therefore propagates to
the full packet `A(10663,3867)` with

```text
in_max-z(R_6)
 =(3^46008/2^10493)x^16170z^18748,

in_min-gamma(R_6)
 =(3^21753/2^10493)
   z^10663(27x^2z+y^3)^8085.                           (16)
```

The state determinants are `-1536` and `12288=-8(-1536)`, agreeing with the
determinant of the THM-3506 face matrix.

## 9. Why THM-3506 honestly left renewal open

No concealed denominator explains the old boundary.  THM-3506 proved three
output faces from five input faces, but did not prove that the two input
renewal faces themselves generate the two output renewal faces.  At that
point the hybrid singleton mechanism had not been abstracted and audited;
THM-3513 subsequently established it only for the fixed `J -> G` transition.

The new argument identifies the missing state: not the exposed pair alone,
but the marked intersection of the complete max-`z` and min-`gamma` faces.
Once that intersection is retained, the two hybrid valuations are positive
linear combinations of its two endpoint gaps, as in (2).  This is a genuine
closure proof and not a fit to the visible sequence.

THM-3506's warning remains essential in two places:

1. the full packet, not merely `(e,m)` or one displayed monomial, is required;
2. polynomiality at the next rung still needs a separate finite-sheet gate.

## 10. Promotion recommendation and exact boundary

Promotion of the reserved THM-3522 namespace is recommended with two layers
kept explicit:

1. **Conditional fixed-chart lemma:** (14), with exact scalars (11)--(12).
2. **Fixed-tower corollary:** complete packets and faces (15)--(16), using
   THM-3506, THM-3513, and THM-3521 for the needed polynomiality statements.

The promotion must continue to exclude:

- polynomiality or finite-sheet nonvanishing of
  `L^10663 Norm(R_6)`;
- an unconditional all-level face orbit;
- irreducibility, squarefreeness, or image-equation status of `R_5` or `R_6`;
- a fifth image component or degree-`243` separability;
- arbitrary Keller maps, `JC(2)`, `DC(2)`, or any general Jacobian claim.

This was the next tower gate at the time of the audit.  THM-3523 subsequently
closes that finite old-`L` sheet, proves polynomial and `L`-coprime
`R_7=L^10663N(R_6)`, and gives `R_7` the complete packet
`A(66907,24255)`.  The next analogous polynomiality gate is now the
finite-sheet unit for `L^66907N(R_7)`, not another renewal-face extraction.

## 11. Independent artifact and replay

The disjoint companion is

```text
04-computation/keller_five_face_renewal_propagation_independent_audit_20260816.py
```

with frozen output

```text
05-knowledge/results/keller_five_face_renewal_propagation_independent_audit_20260816.out
```

and semantic digest

```text
9de2b0a149105263ee1b3a1fba01424f9c7ff274368c689cdc5737fb340bf804.
```

The raw LF-byte SHA256 values are

```text
script  2fafef6bca64174b71fac22aa89e5bb713847ccc47318635c61b9684c5b063cd
output  9485d14a9eb764374ca6abb25fdecafc75a210ce720bdcf4c0f2e7d0536d8dc8
```

Reproduce with

```text
python -B 04-computation/keller_five_face_renewal_propagation_independent_audit_20260816.py
python -B -O 04-computation/keller_five_face_renewal_propagation_independent_audit_20260816.py
```

Normal and optimized runs are required to agree byte-for-byte.

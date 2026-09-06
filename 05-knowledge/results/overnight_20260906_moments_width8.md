# Two-rung trinomial certificates for smaller endpoint degree at most eight

**Status: PROVED, with independently audited FINITE-EXACT symbolic
certificates.** This establishes an unbounded
family statement, not a larger finite coefficient/height census.

Let a complex Laurent polynomial have at most three nonzero terms, exact
extreme exponents `-M<0<N`, and `min(M,N)<=8`. Then

```text
CT(f^m)!=0 for some 1<=m<=M+N.                       (1)
```

More precisely, for every collided first support return in this class,
the first and doubled support-return polynomials have no common torus
zero. An infinite family with smaller endpoint degree three attains
`m=M+N` by cancelling its first support return. Arbitrary-support
`min(M,N)>=3` and
unrestricted trinomial two-rung coprimality remain **OPEN**.

Here the smaller endpoint degree is `min(M,N)`; the total exponent width
`M+N` is unbounded. The artifact filenames retain `width8` as shorthand for
the smaller-endpoint condition only.

## Inheritance, board, and scope of the advance

The closest mechanisms are [THM-4417's width-two parabolic bound](../../01-canon/theorems/THM-4417-width-two-laurent-first-return-parabolic-critical-bound.md)
and [the exact trinomial semigroup classification](synthesis_20260905_moments_trinomial.md).
Both actual proofs were read. The latter's two residue carries are the
sidecar used here; THM-4417 gives context, not the proof of the new result.

The canonical hostile `(-13,1,8)` has two first-level channels but five
second-level channels. The corrected near miss is MISTAKE-544 / THM-1755:
tunability does not force a bounded short relation. The recovered older
technique is [THM-1680's coefficient gcd](../../01-canon/theorems/THM-1680-tnc-trinomial-gcd-decision.md),
combined with [THM-2639's two-rung certificate](../../01-canon/theorems/THM-2639-gmc-equal-mass-two-rung-persistent-collision-certificate.md).
The latter's free-semigroup assumption is kept explicit.

The five-object board was: small-root product; local root correspondence;
integer channel line; residue carries; the resultant as a polynomial in
the unbounded return parameter. The positive signal was a carry-sensitive
family with two first channels and four second channels. A symbolic
remainder closes it, and shifted coefficient positivity extends that
mechanism to a finite list of families with unbounded parameter.

The theorem includes some previously proved boundary families. Those are
controls, not thirty separate new discoveries. Its added scope is the
complete trinomial `min(M,N)=3,...,8` strip and the new carry-sensitive
proof mechanism. No literature-priority claim is made.

## 1. Normalize the support, not the scalar equations

Reflect the variable if necessary so the singleton sign is negative, and
divide all charges by their common content. The exact three-term support
then has the form

```text
(-a,b,c), 1<=b<c, gcd(a,b,c)=1.
```

Reflection and common charge scaling preserve every constant term. The
original width bound follows from a bound on the primitive width because
the latter is no larger. A neutral term is detected at `m=1`; fewer than
three terms are the standard binomial/one-term boundary.

Nonzero scalar and variable dilations normalize the two extreme
coefficients to one. Indeed, for original coefficients `alpha,beta,gamma`,
choose `lambda^(a+c)=alpha/gamma` and multiply `f(lambda*u)` by
`lambda^a/alpha`. Thus it is enough to study

```text
f(u)=u^(-a)+t*u^b+u^c, t!=0.                         (2)
```

The normalization multiplies moments by nonzero scalars and preserves
their zero sets. It does not replace a scalar sum by termwise equations.

Write

```text
g=gcd(a+b,a+c), A=(a+b)/g, B=(a+c)/g.
```

Then `1<=A<B`, `gcd(A,B)=gcd(a,g)=1`. By the proved classification,
collision at the first support return occurs exactly when

```text
a-A*B belongs to <A,B>.                              (3)
```

In that case its mass is `g`, `a>=AB`, `g>B`, and

```text
2g<=a+c=gB.                                         (4)
```

If the first row is singleton, its multinomial and nonzero coefficient
monomial already detect the polynomial. The two-charge subset `{-a,b}`
has a support return at mass `(a+b)/gcd(a,b)`, so the first support-return
mass is at most `(a+b)/gcd(a,b)<=a+c`. Consequently only (3) needs further
work.

## 2. The complete symbolic row retains both carries

At mass `k*g`, balance is the exact equation

```text
A*y+B*z=k*a, x=k*g-y-z.
```

The coefficient polynomial is therefore

```text
T_k(g,t)=sum_(Ay+Bz=ka) (k*g)_(y+z) t^y/(y!*z!),     (5)
```

where `(q)_j=q(q-1)...(q-j+1)` is the falling factorial. The sum runs over
all nonnegative integer pairs, including the later channels supplied by
both residue carries. The actual hypotheses imply `g>a/A`, so every
multinomial in (5) is a positive rational number at admissible integral `g`.

Let `y_k` be the smallest middle exponent in that row. All other middle
exponents differ from it by multiples of `B`. After removing the torus
monomial and setting `v=t^B`, write

```text
T_k(g,t)=t^y_k C_k(g) P_k(g,t^B),                    (6)
```

where `C_k` is the polynomial content in `Q[g]` and `P_k in Q[g,v]` is
primitive as a polynomial in `v`. These contents have no zero at an
admissible `g`: otherwise all coefficients of the nonempty positive row
would vanish. Thus the common torus-zero question is exactly the common
nonzero-root question for `P_1(g,v),P_2(g,v)`.

This is a polynomial family over `Q[g]`, with `g` an indeterminate.
No bound is imposed on `g` in computing its resultant.

## 3. A readable one-carry proof and sharp equality family

Take

```text
a=3, A=1, B=2,
g>=4, gcd(g,3)=1,
charges=(-3,g-3,2g-3).
```

The first row has two channels, but the second has four. The first-level
semigroup is therefore insufficient for the second polynomial. Formula
(5), with the extra channel retained, gives

```text
T_1=(g)_3 t^3/6+(g)_2 t,
T_2=(2g)_6 t^6/720+(2g)_5 t^4/24
    +(2g)_4 t^2/4+(2g)_3/6.
```

At a torus zero of the first row, `t^2=-6/(g-2)`, and direct division gives

```text
T_2=2g(g-1)(2g-1)(23g^2-47g+20)/(15(g-2)^2).          (7)
```

This is positive for every `g>=4`, since writing `g=s+4` gives
`23g^2-47g+20=23s^2+137s+200`. Thus the two rows never vanish together.

All support returns have masses divisible by `g`. The first row can be
cancelled by either nonzero square root of `-6/(g-2)`, so the first nonzero
constant term is then exactly `2g=a+c`. This gives an unbounded sharp
three-term family with smaller endpoint degree three. The smallest
control is `(-3,1,5)`: the first row is `4t^3+12t`, and the second row
equals `560` when `t^2=-3`.

## 4. Thirty symbolic families close every sole-negative `a<=8`

Under (3), `AB<=a<=8`. The complete possibilities are

```text
(a,A,B)=(a,1,B), 2<=a<=8, 2<=B<=a;                  28 families,
(a,A,B)=(6,2,3),(8,2,3).                            2 families.
```

Indeed, `A>=3` would give `AB>=12`; if `A=2`, coprimality leaves only
`B=3`, and the semigroup condition permits `a=6,8`, but not `a=7`.
There are exactly thirty families. The least admissible integer parameter is

```text
g0=floor(a/A)+1.                                    (8)
```

For each family form the resultant of the two primitive polynomials in
(6), clear constant rational denominators, remove constant integer content,
and choose positive leading coefficient. Denote that integer polynomial
by `R_(a,A,B)(g)`. The frozen symbolic certificate proves

```text
R_(a,A,B)(s+g0) has every coefficient strictly positive. (9)
```

All `352` coefficients in the complete thirty-family bank are positive.
The largest resultant degree is `32`, at `(a,A,B)=(8,1,2)`.
Consequently the resultant is nonzero at every real `g>=g0`, in particular
at every admissible integer parameter. Equations (6) and the resultant
identity exclude every common complex torus zero of the first two rows.

This is a finite exact polynomial proof obligation with an unbounded
consequence. It does not extrapolate a sequence of specialized gcds.
The manifest stores both primitive moment polynomials, both contents, the
full integer resultant coefficients, its factorization, and the complete
shifted positive coefficient list. Every shifted identity is reconstructed
exactly before acceptance.

| Sole-negative charge `a` | Symbolic families | Resultant degrees |
|---|---:|---|
| 2 | 1 | 2 |
| 3 | 2 | 3, 4 |
| 4 | 3 | 8, 4, 6 |
| 5 | 4 | 10, 6, 6, 8 |
| 6 | 6 | 18, 16, 9, 8, 10, 2 |
| 7 | 6 | 21, 16, 9, 8, 10, 12 |
| 8 | 8 | 32, 20, 24, 12, 10, 12, 14, 2 |

## 5. The opposite-endpoint reduction is finite and complete

It remains possible that the singleton negative charge is large, while the
positive endpoint is small. For `a>8` and `c<=8`, collision implies

```text
g divides c-b, so g<=c-1,
B<g, so a+c=gB<=g(g-1),
9<=a<=(c-1)(c-2)-c.                                (10)
```

Thus the remaining supports lie in an explicitly bounded rectangle derived
from necessary inequalities. The script exhausts all primitive supports in
this rectangle: `255` triples. Exactly five satisfy (3):

```text
(-9,1,6), (-11,1,7), (-13,1,8), (-20,1,8), (-12,3,8).
```

Each is tested by independent literal Laurent convolution at masses `g`
and `2g`, yielding a gcd which is only a power of the nonzero middle
coefficient. The manifest freezes the two actual coefficient polynomials,
all channel vectors, and a nonzero compressed resultant for each support.
This includes the two-extra-channel hostile `(-13,1,8)`.

The old finite height-sixty bank is not used to fill this proof step.
The five cases are rederived from (10), computed directly, and checked
against formula (5). Together with Section 4, this covers both possible
orientations of `min(a,c)<=8`. Equations (3)--(4) then prove (1).

## 6. Exact map and the finite-line inspiration boundary

The channel set at a fixed mass is the intersection of the two affine
equations `charge=0` and `mass=m` with the nonnegative integer octant.
For three terms it lies on one integer line. The explicit map to its
integer parameter preserves every channel and its multinomial coefficient;
the map to `v=t^B` factors out exactly a torus monomial. These are the
target-preserving maps used in the proof.

A finite-line incidence or no-three-in-line analogy does not itself
control the target: all the channels here are collinear, and rows can have
arbitrarily many points. The cheap hostile is the previously proved family
`(-2h,1,2h+2)`, whose first row has `h+1` channels. The needed information
is the arithmetic line parameter, its two carry residues, and the actual
multinomial weights. No theorem from the Guy--Kelly paper is imported into
this result, and no claim about its content or priority is made.

The connection to the width-two dynamics is a scope intersection. A
three-small-root product divided by the selected root is a product of two
other branches, so THM-4417's root-swap map does not extend by the same
formula. The present route instead bounds the complete two-rung ideal on
the finite normalized-support families. A root correspondence for arbitrary
many terms with smaller endpoint degree three remains a different open problem.

For a horizontal Wick face, choose `H>=c` and send a charge `q` to the
monomial `Z^H W^(H-q)`. Its moment is `(Hm)! CT(f^m)`. This exact map
transfers the new first-return bound to three-term horizontal faces and
to their lowest-face seed in THM-2022. The good prime for general Gaussian
moments remains coefficient dependent; a uniform final Gaussian detection
bound is not supplied.

## 7. Reproduction and audit ledger

```powershell
python 04-computation/overnight_20260906_moments_width8.py
python -O 04-computation/overnight_20260906_moments_width8.py
```

- [Source](../../04-computation/overnight_20260906_moments_width8.py)
- [Exact output](overnight_20260906_moments_width8.out)
- [Byte-matching optimized output](overnight_20260906_moments_width8_optimized.out)
- [Full symbolic certificates](overnight_20260906_moments_width8_certificates.json)
- [Independent proof and complete certificate audit](overnight_20260906_smith_moments_audit.md)

Frozen SHA-256 values, over LF bytes:

```text
source:       8fd1f99f8907b18a454b5524bf0aab08c68406ec9a8c7c1734ce1e18bc3b9df4
both outputs: dd52d856395f080314ec0bd1fc1e59ea54fbec2f20f19f0a9fa8443eadecc405
certificates: 0e0006d444f38326ad0c6ee9fe30a6bd65c4457b9d29a96c77370c384cac659f
```

The certificate bank has thirty symbolic resultants and `352` positive
shift coefficients, five exhaustively derived opposite-endpoint cases,
thirty separate raw-convolution rows across five symbolic families, ten
exceptional raw rows, and two sharp controls. All checks use explicit
failure gates that remain live under optimization. Polynomial identities
are compared by coefficients/expressions, avoiding an irrelevant `QQ`
versus `ZZ` polynomial-object equality mismatch caught in the first run.

The independent audit reconstructs every symbolic resultant using
degree-bounded integer Sylvester determinants: 674 evaluations certify
thirty polynomial identities in the unbounded parameter, and an independent
shift reconstructs all 352 positive coefficients. Its 63,843 exact checks
also cover the opposite-endpoint reduction and torus normalization. A
separate root-agent proof review accepts the normalization, complete case
coverage, and the sharp one-carry formula.

The proof closes the trinomial class with smaller endpoint degree at most
eight, with no bound on the opposite endpoint or coefficients. Uniform trinomial coprimality for larger
smaller endpoint degree, coefficient-uniform general Gaussian detection, and
arbitrary-support Laurent `min(M,N)>=3` remain **OPEN**. The research stopping point is
the proved finite family reduction and audited positive symbolic bank,
before increasing the endpoint-degree range or conjecturing universal positivity.

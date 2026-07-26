---
id: THM-2386
title: "Degree-eighteen H4 common-root elimination"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Every
  genuine four-transposition degree-eighteen H_4 survivor satisfies
  Res_y(P,Q)!=0. A multiple P-root lies on one of the complete walls
  closed by THM-2345 or THM-2359. At a simple P-root, a simple Q-root is
  a hidden totally ramified three-cycle and contradicts the genus-zero
  H_4 signature. The only remaining Q-multiple incidence is an exact
  one-parameter family; after removing its forced cubic discriminant
  root, both a complete specialized-subresultant atlas and an
  independent coefficient-comparison basis [1] exclude the required
  H_3 S_3^2 residual. This proves only uniform P,Q coprimality on genuine
  H_4 survivors. It does not empty the coprime H_4 locus, close degree
  eighteen, or imply JC(2) or DC(2).
source: codex-2026-07-26-h4-common-root-elimination
depends_on:
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - THM-2359-degree-eighteen-perfect-quartic-wall-closure
related:
  - THM-2341-degree-eighteen-deep-wall-local-genus-split
script: 04-computation/jc2_degree18_h4_common_root_elimination_thm2386.py
output: 05-knowledge/results/jc2_degree18_h4_common_root_elimination_thm2386.out
script_sha256: 398112a04c9cb766a259c301ffb4e214224b7d82545b33642b60d8edf7a7a787
output_sha256: 1fa1f203fbb73071f1e02f1fcc5354d49c36b80538461e87e0f629f5ac308bd7
hash_basis: working-tree bytes (LF)
---

# THM-2386 -- every genuine H4 survivor has coprime P and Q

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Retain THM-2332's structured covariants

```text
P
 =245y^4+1890By^2-24300B^2+122472D,

Q
 =539y^6+11340By^4+183708Cy^3
  +(72900B^2-367416D)y^2
  +(2361960BC+2480058W)y.                          (1)
```

Write

```text
F=4P^3+49Q^2.                                      (2)
```

A **genuine H4 survivor** here means a putative degree-eighteen Keller
trajectory in THM-2332's connected higher-support branch whose
normalization has genus zero and genuine ramification signature
`(2,2,2,2)`. In particular,

```text
F=H4*S4^2,

H4 squarefree,                    deg(H4,S4)=(4,4), (3)
```

where the four roots of `H4` are the four transposition branch values of
the normalized cubic cover. This qualification matters: a naked
polynomial identity of the shape (3) does not determine local
normalization. THM-2341's excluded rational ratio is the canonical
hostile example; its even discriminant root hides a genuine
three-cycle.

This theorem proves

```text
Res_y(P,Q)!=0                                      (4)
```

uniformly on every genuine H4 survivor.

## 1. A multiple P-root is already on a closed wall

Suppose `r` is a common root of `P` and `Q`. Direct differentiation
gives

```text
P'(y)=140y(7y^2+27B).                              (5)
```

If `P'(r)=0` and `r=0`, then `P(0)=0` is exactly

```text
126D=25B^2.                                        (6)
```

THM-2345 proves that the complete degree-eighteen wall (6), including
all its edges, carries no Keller trajectory.

If `P'(r)=0` and `r!=0`, equations (5) and `P(r)=0` give

```text
7r^2+27B=0,               504D=115B^2.             (7)
```

THM-2359 proves that the complete perfect-quartic wall in (7) carries no
Keller trajectory. These invocations use the complete-wall conclusions
of THM-2345 and THM-2359, not merely their generic normalized charts.

Consequently a common root on a genuine survivor must satisfy

```text
P(r)=Q(r)=0,                       P'(r)!=0.         (8)
```

In particular `r!=0`, because `Q(0)=0` identically and the `r=0`
incidence is (6).

## 2. A simple P,Q overlap is a hidden three-cycle

Use `s=y-r` and THM-2332's depressed cubic coordinate

```text
v^3+p(y)v+q(y)=0,

p=(16/964467)P,            q=(64/703096443)Q.       (9)
```

Assume first that `Q'(r)!=0`. Then both `p` and `q` have valuation one
in `C[[s]]`. Equation (9) is Eisenstein at `s`: every nonleading
coefficient is divisible by `s`, while its constant coefficient is not
divisible by `s^2`. Thus the normalization has one point over `r` with

```text
e=3,                       ramification contribution 2. (10)
```

The same obstruction has a useful Kummer description. Since `P,Q` are
both simple at `r`, equation (2) has order two there. In (3), `H4` is a
unit and `S4` has order one. On the unramified quadratic resolvent
`t^2=H4`, put

```text
A_+=7Q+S4*t,                A_-=7Q-S4*t.
```

Then

```text
A_+ A_-=-4P^3.                                    (11)
```

At each of the two resolvent points over `r`, both factors vanish; their
valuation sum is three, while their sum and difference each have
valuation one. Hence at each point their valuations are exactly

```text
{ord(A_+),ord(A_-)}={1,2}.                         (12)
```

The resolvent involution swaps the labels `A_+` and `A_-`. The nonzero
valuations modulo three are the Kummer shadow of the same totally
ramified cubic point.

A genuine H4 cover already has four transposition values, contributing
four to Riemann--Hurwitz. The additional contribution two in (10)
would give total ramification at least six and therefore normalization
genus at least one. This contradicts genuine-survivor genus zero.
Therefore every remaining common root satisfies

```text
Q'(r)=0.                                           (13)
```

This is the local hidden-C3 exclusion. It is stronger than a parity
check on `F`: at a simple P,Q overlap, `ord_r(F)=2`, so the
three-cycle is invisible in the squarefree part of the polynomial
discriminant.

## 3. The Q-multiple incidence is one parameter

The covariants are weighted homogeneous. For `r!=0`, replace `y` by
`r*y` and divide `P,Q` by `r^4,r^6`. Equivalently replace

```text
(B,C,D,W) by (B/r^2,C/r^3,D/r^4,W/r^5).           (14)
```

Thus the common root may be normalized to `r=1`. The three equations

```text
P(1)=Q(1)=Q'(1)=0                                  (15)
```

are linear in `(C,D,W)` and have coefficient determinant

```text
111598077326956416 !=0.                            (16)
```

Their unique solution is

```text
C=-35(81B+7)/26244,

D=5(4860B^2-378B-49)/122472,

W=(4050B^2+395B+7)/39366.                          (17)
```

Substitution factors the two covariants exactly:

```text
P=35(y-1)(y+1)(54B+7y^2+7),

Q=7y(y-1)^2
   (1620By+405B+77y^3+154y^2+231y+63).             (18)
```

Write

```text
p3=P/(y-1),                  q4=Q/(y-1)^2.
```

The two endpoint values are

```text
p3(1)=140(27B+7),            q4(1)=525(27B+7).     (19)
```

The excluded multiple-P case is exactly `27B+7=0`. Off that wall,
equations (18)--(19) give

```text
ord_1(F)=3,

F=(y-1)^3 R9,

R9:=4p3^3+49(y-1)q4^2,            R9(1)!=0.        (20)
```

Squarefreeness in (3) now forces one copy of `y-1` into `H4` and one
copy into `S4`. Hence

```text
H4=(y-1)H3,                  S4=(y-1)S3,

R9=H3*S3^2,                 deg(H3,S3)=(3,3).      (21)
```

In particular,

```text
deg gcd(R9,R9')>=3.                                (22)
```

The rest of the proof excludes (22) uniformly in `B`.

## 4. Complete specialized-subresultant atlas

The residual has constant degree and leading coefficient:

```text
deg_y(R9)=9,              lc_y(R9)=73060029,

R9=343*T9,                lc_y(T9)=213003.          (23)
```

Exact factorization over `Q[B]`, up to a nonzero rational constant,
gives

```text
Disc_y(T9)
 =(27B+7)^6
  (54B+7)^3
  (1215B+91)^3
  (2105352B^2+518049B+31997)^3
  Phi13(B).                                        (24)
```

The quadratic and `Phi13` are irreducible over `Q`, and

```text
Phi13(B)
 =2581517833820205830323568640000000B^13
 +27433606764111716665224806400000000B^12
 +22675744145697996379882215628800000B^11
 +7555894329856417179696508849920000B^10
 +1039840519524574652314595234244000B^9
 -40352748361733406218104457728875B^8
 -35168809165204401406565197022160B^7
 -4785349181326208520124913840124B^6
 -143144590616746631993379033504B^5
 +35604040498012586373759791190B^4
 +5212992075701472327189319920B^3
 +330046225412419196564414580B^2
 +10674830277647551598569056B
 +143991306409917519520541.                        (25)
```

Because the leading coefficients of `T9,T9'` are nonzero constants,
their subresultants specialize without a degree-drop exception. Let
`Sres_j` be the degree-`j` subresultant. Reduction of every coefficient
in `Q[B]` modulo each irreducible factor of (24) gives the complete
atlas

| parameter component | exponent in (24) | `(Sres_0,Sres_1,Sres_2)` | specialized gcd degree |
|---|---:|---|---:|
| `27B+7` | 6 | `(0,0,nonzero)` | 2 |
| `54B+7` | 3 | `(0,nonzero,nonzero)` | 1 |
| `1215B+91` | 3 | `(0,nonzero,nonzero)` | 1 |
| `2105352B^2+518049B+31997` | 3 | `(0,nonzero,nonzero)` | 1 |
| `Phi13(B)` | 1 | `(0,nonzero,nonzero)` | 1 |

At the three rational components the monic gcds are, respectively,

```text
(y-1)^2,                    y,                    y+1. (26)
```

Every complex parameter at which `R9` has a repeated root lies on
exactly one of the irreducible components in (24). The table therefore
proves the uniform bound

```text
deg gcd(R9,R9')<=2.                                (27)
```

This contradicts (22).

An independent audit reconstructed the complete specialization
argument. It additionally certified the quadratic factor through its
nonzero discriminant and `Phi13` through an irreducible
degree-thirteen reduction modulo `61`, independently of the
companion's rational factorization.

## 5. Independent coefficient-comparison certificate

There is a second certificate which does not use the discriminant
factorization (24). If (21) held, rescale `S3` to be monic and write

```text
S3=y^3+a2*y^2+a1*y+a0,

H3=73060029y^3+h2*y^2+h1*y+h0.                    (28)
```

The coefficients of `y^8,y^7,y^6` in `H3*S3^2-R9` solve
triangularly:

```text
h2=-73060029(2a2-3),

h1=-453789(
      -4320B+322a1-483a2^2+966a2-966),

h0=-16807(
      233280Ba2-287550B+8694a0
      -26082a1a2+26082a1+17388a2^3
      -39123a2^2+52164a2-38080).                   (29)
```

After substitution, let `E5,...,E0` be the primitive integer
coefficients of `y^5,...,y^0`. Their
`(number of terms,total degree)` signatures are

```text
(16,4), (20,5), (23,5), (25,5), (17,5), (13,5).   (30)
```

An exact grevlex computation over `Q`, in variable order
`(a0,a1,a2,B)`, gives

```text
Groebner(E5,E4,E3,E2,E1,E0)=[1].                  (31)
```

Thus the coefficient equations have no solution even over the
algebraic closure. This independently excludes every factorization
(21), without using squarefreeness of `H3`.

## 6. Scope and reproduction

The alternatives for a common root are now exhausted:

```text
P multiple:       complete wall closed by THM-2345 or THM-2359;

P,Q both simple:  hidden e=3 point, incompatible with genuine H4;

P simple,
Q multiple:       one-parameter residual excluded by (27) and (31).
```

Therefore every genuine four-transposition degree-eighteen H4 survivor
satisfies `Res_y(P,Q)!=0`.

This is a coprimality theorem only. The locus with `Res_y(P,Q)!=0` is
not proved empty here. No degree-eighteen closure, Jacobian-conjecture
consequence, or Dixmier-conjecture consequence follows.

Run

```bash
python3 04-computation/jc2_degree18_h4_common_root_elimination_thm2386.py
python3 -O 04-computation/jc2_degree18_h4_common_root_elimination_thm2386.py
```

Both transcripts are byte-identical to the stored output. The companion
reconstructs (15)--(20), factors the full discriminant in (24), proves
the irreducibility assertions, reduces `Sres_0,Sres_1,Sres_2` modulo
all five factors, checks the three rational gcds, derives (29), and
computes the exact unit basis (31). The Eisenstein/Kummer local
normalization and the complete-wall invocations are the mathematical
proof above, not computer assumptions. An independent audit also
reconstructed the wall split, local ramification, weighted
normalization, subresultant atlas, and coefficient certificate without
finding a missing parameter or quantifier. QED.

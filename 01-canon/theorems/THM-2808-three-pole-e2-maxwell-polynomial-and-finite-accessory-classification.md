---
id: THM-2808
title: "Three-pole e=2 Maxwell polynomial, complete Nielsen atlas, and closure of the two-critical-point response layer"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For every ordered positive three-pole partition (a,b,c), the
  balanced e=2 response chamber is cut out by an explicit degree-(N-3)
  Maxwell polynomial.  Its two false confluent critical points form an
  exact quadratic factor; every remaining root is simple and admissible.
  A four-gap Nielsen atlas independently has N-3 ordered classes.  Together
  with THM-2799 and THM-2805 this closes the abstract balanced e=2 response
  layer, not Keller-chart entry, JC(2), or DC(2).
source: root/jc-e2-three-pole-maxwell-2026-07-28
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2799-one-pole-two-double-zero-chebyshev-response-classification
  - THM-2805-general-two-pole-e2-maxwell-eliminant-and-nielsen-classification
related:
  - THM-2800-two-pole-two-double-zero-stieltjes-recurrence-and-first-nielsen-pair
  - THM-2768-modular-c2-c3-a4-s4-bass-serre-quotient
script: 04-computation/jc2_three_pole_e2_maxwell_thm2808.py
output: 05-knowledge/results/jc2_three_pole_e2_maxwell_thm2808.out
script_sha256: 053958bfbc85baf06e40eb6d88d56ad23a2a03a8fbb04d9f8cf39cf9b5574f99
output_sha256: f8bb83a8acef3e3d541a94219506bc7dd6534e65672607715faa851a5967522d
independent_script: 04-computation/jc2_three_pole_e2_maxwell_nielsen_thm2808.py
independent_output: 05-knowledge/results/jc2_three_pole_e2_maxwell_nielsen_thm2808.out
independent_script_sha256: 0634fece763f1bb4a3e5c908c8f4a9d69e7fefbedc35f009fe20320bdb605c41
independent_output_sha256: aea811ce152ef9c9bc70df56e27239686c17e2ff348008c6628be414257e638f
hash_basis: LF-normalized bytes
---

# THM-2808 -- three-pole e=2 Maxwell polynomial, complete Nielsen atlas, and closure of the two-critical-point response layer

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

The last possible pole count in the balanced two-double-zero response layer
has one apparent accessory parameter.  It is not a positive-dimensional
modulus: it is a one-variable equal-critical-value, or Maxwell, coordinate.
Two independent exact companions arrived concurrently.  The first checks
unordered passport representatives deeper in degree; the second checks
every ordered passport in a smaller range and reconstructs the labelled
Nielsen orbits directly.

## 1. Three poles reduce to one Maxwell condition

Work over an algebraically closed field of characteristic zero.  Fix

```text
a,b,c>=1,                 N=a+b+c>=4,
e=2,                      h=3,
pole partition=(a,b,c).                                  (1)
```

THM-2796 gives `s=N-4` simple zeros and third partition `(N)`.  Label and
normalize the poles to `0,1,lambda`, the totally ramified third-fibre point
to infinity, and `F(infinity)=1`.  Put

```text
D=x^a(x-1)^b(x-lambda)^c,
T=x(x-1)(x-lambda).                                    (2)
```

If the monic numerator is `B=S E^2`, total ramification at infinity says

```text
F=B/D,                    B-D=-v,              v!=0.   (3)
```

Thus the two roots of `E` are distinct non-pole critical points of `D`
with common critical value `v`.  The logarithmic derivative is

```text
D'/D=K/T,

K=Nx^2-Ux+a lambda,
U=(a+c)+(a+b)lambda.                                  (4)
```

Divide `D` by the monic quadratic `E=K/N`:

```text
D=H E+v_0(lambda)+R(lambda)x.                         (5)
```

If `gamma,delta` are the roots of `E`, then

```text
R(lambda)
 =[D(gamma)-D(delta)]/(gamma-delta).                  (6)
```

Away from a collision, the two critical values agree exactly when `R=0`.
An explicit all-degree recurrence is

```text
sigma=U/N,                   pi=a lambda/N,
h_(-1)=0, h_0=1, h_1=sigma,
h_j=sigma h_(j-1)-pi h_(j-2),                         (7)

x^m=h_(m-1)x-pi h_(m-2)                 modulo E.
```

Consequently

```text
R =
 sum_(0<=j<=b,0<=k<=c)
 (-1)^(b-j+c-k) binom(b,j)binom(c,k)
 lambda^(c-k) h_(a+j+k-1).                            (8)
```

## 2. Remove exactly the two confluent critical points

The collision polynomial is

```text
Delta(lambda)=Disc_x(K)
 =U^2-4Na lambda.                                      (9)
```

It has two simple roots, neither `0` nor `1`, because

```text
Disc_lambda(Delta)=-16Nabc!=0,
Delta(0)=(a+c)^2,              Delta(1)=(b+c)^2.      (10)
```

Moreover,

```text
Delta divides R.                                      (11)
```

At a root of `Delta`, write `E=(x-u)^2`.  The point `u` is not a pole by
`(4)` and `(10)`.  Differentiating `(5)` at `u` gives `R=D'(u)=0`, since
`D'=NDE/T`.

The factor occurs only once.  Write locally

```text
E=(x-u)^2-w^2,             w^2=Delta/(4N^2),
P=x^(a-1)(x-1)^(b-1)(x-lambda)^(c-1).
```

Since `D'=NP E`, termwise antidifferentiation gives

```text
R
 =N/(2w) integral_(-w)^w P(u+t)(t^2-w^2)dt,

lim_(Delta->0) R/Delta=-P(u)/(6N)!=0.                 (12)
```

The integral is only compact notation for a polynomial identity.  Define
the genuine Maxwell polynomial

```text
Q_(a,b,c)(lambda):=R(lambda)/Delta(lambda).           (13)
```

Then `gcd(Q,Delta)=1`, and

```text
deg Q_(a,b,c)=N-3.                                    (14)
```

Indeed, as `lambda` tends to infinity,

```text
gamma~(a+b)lambda/N,             delta~a/(a+b),

[lambda^(N-1)]R
 =(-1)^c (a+b)^(a+b-1)c^c/N^(N-1),

[lambda^(N-3)]Q
 =(-1)^c (a+b)^(a+b-3)c^c/N^(N-1),                  (15)
```

so neither leading coefficient vanishes.

## 3. Every Maxwell root is simple and admissible

Neither `0` nor `1` is a root of `Q`.  At those parameters one critical
root merges with a pole, while the other has nonzero value; `(6)` and
`(10)` give `Q(0)Q(1)!=0`.

Continue the two simple critical points locally and set

```text
J(lambda)=D_lambda(gamma(lambda))-D_lambda(delta(lambda)).
```

At `J=0`, the critical-point derivative terms vanish, and

```text
J'
 =cv(gamma-delta)/
   [(gamma-lambda)(delta-lambda)] !=0.                (16)
```

Here `c,v`, the critical-point difference, and both denominators are
nonzero.  Therefore every zero of `Q` is simple.

For such a root, `(5)` has `R=0` and `v_0=v`.  Put

```text
B=D-v,                    E=K/N,
S=B/E^2,                  C=Nv.                       (17)
```

Criticality gives `E^2|B`.  Any repeated root of `B` is a root of `D'`.
The pole roots have value `B=-v`, while the only remaining critical points
are the two roots already removed by `E^2`.  Thus `S` is squarefree and
disjoint from `ET`.  There are no hidden passport gates.

Finally set

```text
F=B/D,
G=C E/(2DT),
V=4 S D T^2/C^2.                                     (18)
```

Because `D'=NDE/T`, one obtains

```text
F'=2G,                    F=VG^2,
2VG'+V'G=2.                                         (19)
```

For general prescribed `kappa!=0`, rescale `G,V` as in THM-2796.
Conversely, `(3)--(6)` show that every normalized response in `(1)` gives
a root of `Q`.  Hence `Q` is already a complete ordered classification.

## 4. Four cyclic gaps give the same `N-3`

Fix the totally ramified inertia

```text
rho=(0 1 ... N-1).
```

The zero inertia is a product of two disjoint transpositions, represented
by two chords.  Its product with `rho` has three cycles exactly when the
chords are noncrossing.  If the four positive cyclic gaps are
`g_1,g_2,g_3,g_4`, the two noncrossing pairings have cycle lengths

```text
(g_1,g_3,g_2+g_4),          (g_2,g_4,g_1+g_3).       (20)
```

Label the pole cycles by the ordered parts `(a,b,c)`.  Choosing one part
as the sum of two opposite gaps gives `part-1` positive splits.  Equivalently,
the rooted encodings number `2[(a-1)+(b-1)+(c-1)]`, and forgetting which
opposite single-cycle gap was written first divides by two.  Therefore

```text
# ordered Nielsen classes for (a,b,c)=N-3.            (21)
```

Their defects sum to `2N-2`, so Riemann existence realizes all of them as
genus-zero covers.  The ordered poles and totally ramified point make the
normalization unique.  Thus the Nielsen count independently saturates
`deg Q=N-3` and supplies a second proof that all roots are admissible.

## 5. Unmarked count and closure of the `e=2` response layer

For repeated pole multiplicities, unmarked maps are the orbits of the
`N-3` roots under the **actual** multiplicity stabilizer in the anharmonic
`S_3` action.  Dividing naively by its order is false.  For `(2,2,2)`,

```text
Q proportional to (lambda-2)(lambda+1)(2lambda-1),
```

and `S_3` is transitive, so three ordered roots give one unmarked map.

Summed over all unordered three-pole partitions, Burnside gives

```text
H_3(N)
 =[2 binom(N,4)+1_(2|N) 2 binom(N/2,2)]/N
 =(N-1)(N-2)(N-3)/12 + 1_(2|N)(N-2)/4
 =sum_(m=4)^N floor((m-2)/2)^2.                      (22)
```

The identity fixes two noncrossing pairings on every four-subset.  The
only other contribution is the half-turn for even `N`, fixing two pairings
on each pair of antipodal vertex pairs; order-four rotations swap the two
pairings.

The passport equation `N-r=e-h+1>=0` forces `1<=h<=3` for `e=2`.
THM-2799, THM-2805, and the present theorem therefore exhaust the entire
abstract balanced `e=2` response layer:

```text
H_1(N)=floor((N-2)/2),
H_2(N)=binom(N-2,2),
H_3(N) given by (22).                                (23)
```

Every present `h=3` map is nonsplit.  For `N>=5`, `deg S=N-4>0`; at
`N=4`, the pole multiset `(2,1,1)` has odd parts.  This is precisely
THM-2796's squareclass criterion.

## 6. The quartic `C2` and `C3` factors are separate

For `(a,b,c)=(1,1,2)`,

```text
lambda=1/2,
E=x^2-x+1/8,              v=-1/64,
F=E^2/[x(x-1)(x-1/2)^2].                             (24)
```

The branch generators may be taken as a four-cycle and
`(0 1)(2 3)`.  They generate the order-eight square group `D_4`, with

```text
D_4/V_4 isomorphic to C_2.                            (25)
```

THM-2800's adjacent quartic `h=2` response instead has monodromy `A_4`,
with

```text
A_4/V_4 isomorphic to C_3.                            (26)
```

This validates the modular-group heuristic and identifies its missing
sidecar: the `C_2` and `C_3` quotients live on different accessory strata,
not in one `S_3=S_4/V_4` action.  The complete quartic `e=2` layer has the
`h=1` function of `x^2`, the `h=2` `A_4` map, and the `h=3` `D_4` map;
none has `S_4` monodromy.  A quartic Keller-resolvent program needs a
genuine cross-stratum gluing mechanism or another response layer.

## 7. Exact controls

The primary companion checks every unordered positive pole partition
through `N=10`, the response identities and gates, pole-swap covariance,
and the unmarked Burnside formula through `N=20`.  The independently
written companion checks all 55 ordered passports through `N=8`, directly
enumerates every labelled Nielsen passport through `N=12`, and verifies
the quartic `D_4/V_4` quotient.

Both use exact rational polynomial and permutation arithmetic, contain no
Python `assert` node, and have identical normal, optimized, and stored
transcripts.  Run

```text
python 04-computation/jc2_three_pole_e2_maxwell_thm2808.py
python -O 04-computation/jc2_three_pole_e2_maxwell_thm2808.py
python 04-computation/jc2_three_pole_e2_maxwell_nielsen_thm2808.py
python -O 04-computation/jc2_three_pole_e2_maxwell_nielsen_thm2808.py
```

The finite controls support but do not replace the all-degree proofs.

## 8. Scope and failure boundaries

This theorem closes an **abstract balanced response layer**.  It does not
show that any response enters a polynomial Keller chart, satisfies the
inherited Faber flux equations, or comes from a Weyl-algebra endomorphism.
It proves neither `JC(2)` nor `DC(2)`.

Three boundaries are load-bearing:

1. using `R` before dividing by `Delta` includes two false confluent
   critical-point solutions;
2. allowing one of `a,b,c` to vanish degenerates `(10)` and returns to the
   two-pole chamber; and
3. quotienting repeated pole parts by stabilizer size rather than its
   actual action miscounts symmetric maps such as `(2,2,2)`.

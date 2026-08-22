---
id: THM-3263
title: "Degree-seven retuned quartic infinity wall and odd critical-resultant inertia"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In each of
  THM-3257's two cubic accessory fields, enlarge the clutch to
  B_(r,s,t)=1+r*x^7+s*x^8+t*x^9. Keeping the old degree-eight tuning is
  transverse: the degree-seven term reopens the preceding x^97 layer.
  Retuning s as a function of r gives one strict-transform line; one unique
  r=r_* cancels its next coefficient and exposes a nonzero quartic carry.
  The saturated critical resultant then has degree 51, with 51 squarefree
  boundary-disjoint critical points, while four reciprocal roots escape in
  one odd local 4-cycle. This is critical-resultant inertia, not polynomial
  inverse-cover or Jelonek inertia, and supplies no Keller map or JC(2)
  consequence.
source: root/creative-synthesis-cont/2026-08-03
audit: >
  The exact companion pins fourteen artifacts from THM-3237, THM-3257,
  THM-3057, THM-3059, THM-3064, THM-3066 and THM-3068; independently expands
  the universal three-parameter reciprocal jet; reconstructs both accessory
  fields; proves every tuning denominator and quartic carry nonzero by exact
  norms; computes the degree-51 characteristic-zero quotients; and supplies
  two good reductions certifying squarefreeness and boundary-disjointness.
  Normal and optimized executions byte-match the frozen transcript, and the
  source has no assertion node or floating literal. A separate audit
  rederived the reciprocal coefficients from the response-square-root jet,
  independently recomputed all eight accessory norms, replayed the finite
  reductions, and verified the Newton edge and 4-cycle without importing the
  discovery formulas.
depends_on:
  - THM-3257-degree-eight-tuned-cubic-infinity-wall-and-three-root-critical-escape
related:
  - THM-3057-tame-quartic-inertia-clutch-index-resonance
  - THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor
  - THM-3068-c3-escape-inverse-pole-ledger-and-reciprocal-cofactor-shift
script: 04-computation/jc_degree7_retuned_quartic_infinity_wall_wildcard_20260803.py
output: 05-knowledge/results/jc_degree7_retuned_quartic_infinity_wall_wildcard_20260803.out
script_sha256: aabdb06854ef9045d31c4fdbe76e2933282f44b49230ce574a77591af9a6df56
output_sha256: 2ad12c392ff4c3228a7a242cffa58449a734339e102628d3e73386883a9bd0b2
hash_basis: LF-normalized bytes
---

# THM-3263 -- degree seven buys odd inertia, not a Keller map

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3257 shows that a degree-eight term, dormant on the generic degree-nine
face, becomes live after restriction to the primitive infinity wall. It tunes
a cubic reciprocal carry and makes three critical-resultant roots escape in
one even 3-cycle. The next lower coefficient reveals both a stronger
strict-transform staircase and a sharper stopping boundary.

The direct extension fails first: adding degree seven while retaining the old
degree-eight tuning reopens the preceding wall. After the required retuning,
however, one unique degree-seven value exposes a quartic carry and an odd
4-cycle.

## 1. Three-parameter clutch and normalized jets

Retain either of THM-3257's characteristic-zero cubic accessory fields,
response pair (V,A), monic degree-44 passport boundary, and saturated critical
resultant. Write

~~~text
V=a^2*x^16*(1+r1/x+r2/x^2+r3/x^3+r4/x^4+...),
B_(r,s,t)=1+r*x^7+s*x^8+t*x^9.                         (1)
~~~

As in THM-3257, a=2/Gamma. Define

~~~text
Kappa =4*r2-3*r1^2,
Lambda=4*r2-r1^2,
Rho   =r1^4+16*r1*r3-16*r2^2-8*r1,                    (2)

Sigma =
 128-128*r3+256*r2*r4-256*r3^2-64*r1^2*r4
 +128*r1*r2*r3-64*r2^3-32*r1^3*r3
 +48*r1^2*r2^2-12*r1^4*r2+r1^6.                       (3)
~~~

Let K_(r,s,t) be the universal critical resultant before division by the
passport boundary and H_(r,s,t) its saturated quotient. Because the boundary
is monic, after a string of top cancellations the first surviving coefficient
of K is the leading coefficient of H.

## 2. The fixed-s near miss is false

On the primitive degree-nine wall

~~~text
t=t_*=a,                                                (4)
~~~

the coefficients of x^99 and x^98 vanish. The next universal layer is

~~~text
[x^97]K=-2*a^11*(a*Kappa+4*r1*s-8*r).                  (5)
~~~

THM-3257 used

~~~text
s_old=-a*Kappa/(4*r1).
~~~

Substituting this value into (5) gives

~~~text
[x^97]K=16*a^11*r.                                     (6)
~~~

Thus the tempting claim that degree seven first changes the old cubic carry
is false. It is transverse to the old wall and reopens the preceding layer.
Equation (6) is the minimal witness and identifies the missing coordinate:
s must move with r.

## 3. Retuned strict transform and unique quartic point

The unique repair of (5) is

~~~text
s=s(r)=(8*r-a*Kappa)/(4*r1).                           (7)
~~~

Along (7), the next coefficient is

~~~text
[x^96]K=-a^12*(Rho+8*r*Lambda/a)/r1.                   (8)
~~~

When Lambda is nonzero, it has the unique zero

~~~text
r_*=-a*Rho/(8*Lambda),
s_*=(8*r_*-a*Kappa)/(4*r1).                            (9)
~~~

At this point the next carry is

~~~text
[x^95]K=-3*a^12*Sigma/(8*Lambda).                      (10)
~~~

In both accessory fields r1, Lambda, Rho and Sigma are nonzero. Consequently
the exact degree staircase on the primitive wall is

~~~text
deg H=53  when s!=s(r),
deg H=52  when s=s(r) and r!=r_*,
deg H=51  at (r,s)=(r_*,s_*).                          (11)
~~~

As in THM-3257, the generic off-wall degree is 55; special lower
cancellations away from the displayed top face are not classified here.

## 4. Exact accessory controls

The nonzero gates and tuned parameters have the following exact norms.

| passport | Norm(Lambda) | Norm(Sigma) | Norm(r_*) | Norm(s_*) |
|---|---:|---:|---:|---:|
| (4,1,1,1) | 12448544/30625 | -157881139404422797328384/28722900390625 | 2796986463061/2489708800000 | -15118672709/6224272000 |
| (3,2,1,1) | -42471424/275625 | 64559323245631953174528/28722900390625 | 8958651810434375/82564448256 | -1836115611739609375/14090999169024 |

The two exact degree-51 residual coefficient digests are respectively

~~~text
3c2ac8e2cbd13ee2568f2565bd780e2b0b8bb11facad583cc4b626e1f2265626,
2808b49be3b3ac590efdd39e1015da59e620d45b11bd103271e21e82bdf6d499. (12)
~~~

Good reductions give:

| passport | (p,u) | (t_*,s_*,r_*) | deg H | deg gcd(B,g) | deg gcd(H,g) | deg gcd(H,H') | [x^51]H |
|---|---:|---:|---:|---:|---:|---:|---:|
| (4,1,1,1) | (113,85) | (85,103,91) | 51 | 0 | 0 | 0 | 41 |
| (3,2,1,1) | (101,64) | (89,87,55) | 51 | 0 | 0 | 0 | 9 |

Here g=S*T is the inherited owner boundary. The reductions preserve degree
and have unit denominators, so the characteristic-zero residuals are
squarefree and boundary-disjoint. The constant leading y-coefficient in the
first gradient equation therefore leaves 51 reduced transverse critical
points in each accessory field.

## 5. Quartic escape and odd local inertia

Fix (r,s)=(r_*,s_*) and put

~~~text
delta=t-a,              epsilon=Gamma*t-2=Gamma*delta,
w=1/x,                  Hhat=w^55*H(1/w).              (13)
~~~

The transverse derivative of the leading coefficient is -16*a^11. Together
with (10), the completed initial form is

~~~text
Hhat=
 -16*a^11*delta
 -3*a^12*Sigma*w^4/(8*Lambda)
 +O(delta^2,delta*w,delta*w^2,delta*w^3,w^5).          (14)
~~~

Its single compact Newton edge gives

~~~text
w^4 ~ -128*Lambda*delta/(3*a*Sigma)
    = -64*Lambda*epsilon/(3*Sigma).                    (15)
~~~

Thus exactly four units of critical-resultant intersection multiplicity move
to x=infinity. A complex meridian around delta=0 cycles the four Puiseux
branches, so the escaping critical roots carry one 4-cycle. Its permutation
sign is odd.

## 6. Typed consequence and stopping boundary

The strict-transform mechanism is exact:

~~~text
degree eight tunes the first restricted wall
  -> degree seven reopens the preceding layer
  -> retuning degree eight makes degree seven live one layer later
  -> one degree-seven value exposes a quartic edge.     (16)
~~~

This proves that parity of local critical-resultant inertia is tunable in the
inherited family. It does not prove that odd infinity inertia is available
for a generically finite polynomial map. The typed connection is

~~~text
source:      reciprocal jet of THM-3257's critical resultant;
operation:   add r*x^7, restrict t=a, retune s=s(r), divide and recompute;
target:      first nonzero strict-transform coefficient and Newton edge;
preserved:   vanishing order, escaping intersection multiplicity, local
             permutation of critical-resultant roots;
destroyed:   inverse-sheet labels, affine polynomial-map realization,
             branchwise Keller cofactors and global chart;
sidecar:     a lawful second polynomial coordinate and pointed inverse-
             different/cofactor units on every inverse branch.           (17)
~~~

THM-3066 proves that a symmetric product quotient cannot replace the missing
sheetwise cofactor ratios, and THM-3068 gives a disconnected Laurent hostile.
Most directly, 51 Morse critical points remain at the tuned wall, so the
displayed scalar polynomial has no constant-Jacobian mate in the current
chart. The 4-cycle is not Jelonek inertia and yields no result about JC(2) or
DC(2).

## 7. Exact reproduction

Run

~~~text
python3 04-computation/jc_degree7_retuned_quartic_infinity_wall_wildcard_20260803.py
python3 -O 04-computation/jc_degree7_retuned_quartic_infinity_wall_wildcard_20260803.py
~~~

and compare LF-normalized bytes with the declared transcript. The companion
uses exact symbolic, rational and finite-field arithmetic only and performs
all dependency, source-AST, degree, gcd, norm and local-law checks without
randomness or optimization-sensitive assertions.

QED.

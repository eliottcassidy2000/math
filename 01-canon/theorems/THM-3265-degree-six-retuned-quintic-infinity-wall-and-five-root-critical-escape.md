---
id: THM-3265
title: "Degree-six retuned quintic infinity wall and five-root critical escape"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In each of
  THM-3263's two accessory fields, enlarge the clutch to
  B_(q,r,s,t)=1+q*x^6+r*x^7+s*x^8+t*x^9. On the primitive wall, q first
  reopens the x^96 layer after s is retuned; r must then be retuned as a
  function of q. One unique q=q_* cancels the next coefficient and exposes a
  nonzero quintic x^94 carry. The saturated critical resultant has degree 50,
  with 50 squarefree boundary-disjoint critical points, while five reciprocal
  roots escape in one even local 5-cycle. This is critical-resultant inertia,
  not inverse-cover or Jelonek inertia, and supplies no Keller map or JC(2)
  consequence.
source: root/creative-synthesis-cont/2026-08-03
audit: >
  The exact companion pins sixteen artifacts from THM-3237, THM-3257,
  THM-3263 and the cofactor/inertia boundary chain; expands the universal
  reciprocal jet through x^94; extracts a canonical 58-term quintic-edge
  invariant; reconstructs both accessory fields; proves every denominator
  and leading carry nonzero by exact norms; computes the degree-50 quotients;
  and supplies two good squarefree boundary-disjoint reductions. Normal and
  optimized executions byte-match the frozen transcript and the source has
  no assertion node or floating literal. An independent SymPy expansion
  rederived the complete tuning chain and quintic invariant without importing
  the candidate; an independent characteristic-zero reconstruction matched
  all norms and quotient digests, and coefficientwise reduction of those
  actual quotients reproduced both finite controls and the Newton law.
depends_on:
  - THM-3263-degree-seven-retuned-quartic-infinity-wall-and-odd-critical-resultant-inertia
related:
  - THM-3057-tame-quartic-inertia-clutch-index-resonance
  - THM-3068-c3-escape-inverse-pole-ledger-and-reciprocal-cofactor-shift
  - THM-3233-fraction-free-exceptional-quotient-tail-renormalization
script: 04-computation/jc_degree6_retuned_quintic_infinity_wall_wildcard_20260803.py
output: 05-knowledge/results/jc_degree6_retuned_quintic_infinity_wall_wildcard_20260803.out
script_sha256: 8b0c633366869aa97c207e9ba3535f9b393d5b8441b6cbda5419638401b7a1a0
output_sha256: c7b08d089213a6b5737880a537af8c72c4febd9a061172ddf0b1989a285fd578
hash_basis: LF-normalized bytes
---

# THM-3265 -- degree six retunes the quartic wall to a quintic escape

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3237, THM-3257 and THM-3263 expose a nested reciprocal wall: degree nine
creates a square-root face, degree eight tunes a cubic edge, and degree seven
must first reopen and then retune that wall before exposing a quartic edge.
This theorem performs the next step exactly.

The degree-six coefficient again enters one layer earlier than the naïve
carry. After two successive retunings it has one unique value that exposes a
quintic edge and five-root local escape.

## 1. Four-parameter clutch and invariants

Write the leading response jet as

~~~text
V=a^2*x^16*(1+r1/x+r2/x^2+r3/x^3+r4/x^4+r5/x^5+...),
B=1+q*x^6+r*x^7+s*x^8+t*x^9.                           (1)
~~~

Retain the invariants from THM-3263,

~~~text
Kappa =4*r2-3*r1^2,
Lambda=4*r2-r1^2,
Rho   =r1^4+16*r1*r3-16*r2^2-8*r1,                    (2)

Sigma =
 128-128*r3+256*r2*r4-256*r3^2-64*r1^2*r4
 +128*r1*r2*r3-64*r2^3-32*r1^3*r3
 +48*r1^2*r2^2-12*r1^4*r2+r1^6,                       (3)
~~~

and define the next normalized invariant

~~~text
Theta=r1^3-4*r1*r2+8*r3+8.                             (4)
~~~

As before, a=2/Gamma. Let K be the critical resultant before division by the
monic degree-44 passport boundary and H its saturated quotient.

## 2. Three nested retunings

On the primitive wall t=a, the x^99 and x^98 coefficients vanish and

~~~text
c97=[x^97]K=-2*a^11*(a*Kappa+4*r1*s-8*r).              (5)
~~~

The first strict transform is therefore

~~~text
s=(8*r-a*Kappa)/(4*r1).                                (6)
~~~

Along (6), the new coefficient q appears one layer earlier than the old
quartic carry:

~~~text
c96=-a^11*(a*Rho+8*r*Lambda-32*q*r1)/r1.               (7)
~~~

When Lambda is nonzero, the second retuning is

~~~text
r(q)=(32*q*r1-a*Rho)/(8*Lambda),
s(q)=(8*r(q)-a*Kappa)/(4*r1).                          (8)
~~~

Along (8), the next layer is

~~~text
c95=-3*a^11*(a*Sigma+64*q*Theta)/(8*Lambda).            (9)
~~~

When Theta is nonzero, it has the unique zero

~~~text
q_*=-a*Sigma/(64*Theta),
r_*=r(q_*),               s_*=s(q_*).                  (10)
~~~

At (10), define the exact integral polynomial Psi by

~~~text
Psi=-16*Theta^2*a^(-12)*c94,        c94=[x^94]K.        (11)
~~~

After the substitutions (6), (8) and (10), the companion proves that Psi is
a 58-term polynomial in r1,...,r5 of total degree 11, with canonical
expression digest

~~~text
6b1814756894c97cfe6a31a3bd4ae51bc1dc591879211ec5f67d65f2a960e860.
~~~

The first surviving coefficient is

~~~text
c94=-a^12*Psi/(16*Theta^2).                             (12)
~~~

Setting q=0 in (8)--(9) recovers THM-3263's degree-seven tuning and quartic
carry exactly. The full on-wall degree staircase is

~~~text
deg H=53  before (6),
deg H=52  after (6) but before (8),
deg H=51  after (8) with q!=q_*,
deg H=50  at q=q_*.                                    (13)
~~~

## 3. Exact accessory gates

The three load-bearing invariant norms are:

| passport | Norm(Theta) | Norm(Sigma) | Norm(Psi) |
|---|---:|---:|---:|
| (4,1,1,1) | -37200002176/5359375 | -157881139404422797328384/28722900390625 | -572298646716642082095739982919499776/301716116790771484375 |
| (3,2,1,1) | 249238609408/48234375 | 64559323245631953174528/28722900390625 | -4402683590765317259172542063793800401125376/1697153156948089599609375 |

They are all nonzero. The tuned-parameter norms and exact degree-50 quotient
digests are:

| passport | Norm(q_*) | Norm(r_*) | Norm(s_*) | H digest |
|---|---:|---:|---:|---|
| (4,1,1,1) | 2107940627729338801/186000010880000000 | -3018470836755969/1860000108800000 | 1568249062443/1162500068000 | 7926bec47cfbfad38a1de91dfeebf84b824da3c780f478240bf66bc6e9981d26 |
| (3,2,1,1) | -8978749567656192625/9187932097216512 | 8285100359018853125/484519856689152 | -4924346007190859375/161506618896384 | 631e99ab7332025f3a338fd15919afa90cd9e542e2ed3c314dc0a4b4867015ec |

Good reductions give:

| passport | (p,u) | (t_*,s_*,r_*,q_*) | deg H | deg gcd(B,g) | deg gcd(H,g) | deg gcd(H,H') | [x^50]H |
|---|---:|---:|---:|---:|---:|---:|---:|
| (4,1,1,1) | (113,85) | (85,1,43,54) | 50 | 0 | 0 | 0 | 97 |
| (3,2,1,1) | (101,64) | (89,75,80,100) | 50 | 0 | 0 | 0 | 56 |

Here g=S*T is the inherited owner boundary. The reductions preserve degree
and denominators, so the characteristic-zero quotients are squarefree and
boundary-disjoint. Fifty reduced transverse critical points remain in each
accessory field.

## 4. Quintic escape and even local inertia

Fix (q,r,s)=(q_*,r_*,s_*), and put

~~~text
delta=t-a,              epsilon=Gamma*t-2=2*delta/a,
w=1/x,                  Hhat=w^55*H(1/w).              (14)
~~~

The transverse leading derivative is -16*a^11. Together with (12), the
completed initial form is

~~~text
Hhat=
 -16*a^11*delta
 -a^12*Psi*w^5/(16*Theta^2)
 +O(delta^2,delta*w,delta*w^2,delta*w^3,delta*w^4,w^6). (15)
~~~

Its one compact Newton edge gives

~~~text
w^5 ~ -256*Theta^2*delta/(a*Psi)
    = -128*Theta^2*epsilon/Psi.                         (16)
~~~

Thus exactly five units of critical-resultant intersection multiplicity move
to x=infinity. A complex meridian cycles the five Puiseux branches, producing
one 5-cycle. Its sign is even.

## 5. Mechanism and stopping boundary

The successive strict transforms now read

~~~text
degree 9: primitive square-root wall;
degree 8: one tuning exposes a cubic edge;
degree 7: reopen, retune, expose a quartic edge;
degree 6: reopen, retune twice, expose a quintic edge.   (17)
~~~

This proves one more step of the coefficient staircase in the two inherited
accessory fields. It does not prove an all-degrees pattern, and it exposes why
local inertia parity is not the load-bearing planar-JC obstruction: the signs
alternate while the inverse-cover and cofactor debts do not move.

The five escaping objects are roots of a critical resultant, not normalization
sheets of a generically finite polynomial map. No branchwise Keller cofactor,
pointed inverse-different unit, affine second coordinate, polynomial inverse
cover or Jelonek component is supplied. The 50 surviving critical points
already prevent a constant-Jacobian mate in the current chart. Nothing here
proves or refutes JC(2) or DC(2).

## 6. Exact reproduction

Run

~~~text
python3 04-computation/jc_degree6_retuned_quintic_infinity_wall_wildcard_20260803.py
python3 -O 04-computation/jc_degree6_retuned_quintic_infinity_wall_wildcard_20260803.py
~~~

and compare LF-normalized bytes with the declared transcript. All symbolic,
algebraic-field and finite-field calculations are exact and dependency-pinned.

QED.

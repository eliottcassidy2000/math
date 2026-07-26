---
id: THM-2350
title: "Owner-pivot dual dipoles and magnetic target Dirichlet energy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For every
  explicit owner-aligned pivot packet of THM-2309, the two target
  characters have unique omitted-unit gauge representatives
  e_a-e_(k_a) and e_b-e_(k_b). Equivalently the target quotient map is
  pi(r)=(r_a-r_(k_a),r_b-r_(k_b)). Thus canonical target twists are two
  disjoint balanced target-unit dipoles, and an atomic marked relation
  lands by the difference of its left and right dipole residues plus the
  deepest-comb translation; transported-word modes vanish from this
  formula. The inverse-character boundary of THM-2343/2344 is the
  one-dimensional space of covariantly constant sections of a magnetic
  C_13 square C_13
  connection. Its exact magnetic Dirichlet energy equals the ordinary
  Dirichlet energy of the full current H and a Laplacian-weighted
  nonzero-target energy. A single nonzero dipole edge defect therefore
  proves nonzero target survival. No such defect is proved nonzero for a
  canonical row; no word-matching component, all-unit visible aggregate,
  scalar exclusion, terminal phase, or LRC(14) closure is asserted.
source: codex-2026-07-25-owner-pivot-dual-dipole
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2343-deep-comb-affine-target-catalyst
  - THM-2344-correlation-inverse-rigidity-and-aligned-tooth-twist-hostile
  - THM-2349-first-depth-one-delayed-shallow-restart
related:
  - THM-2307-dual-rank-six-reconstruction-spectrum-and-selector-no-go
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2340-owner-word-anova-target-landing
script: 04-computation/lrc14_owner_pivot_dual_dipole_thm2350.py
output: 05-knowledge/results/lrc14_owner_pivot_dual_dipole_thm2350.out
script_sha256: 59748a8bb55fe87f336e9e2ea5b0722a326dfe8d6f367b46019aaa0a416328fe
output_sha256: b9c101db48f7f11381111adb81cd179d597cb7955b2e6153b4c2f4fb81c63175
hash_basis: working-tree bytes (LF)
---

# THM-2350 -- the target plane is a pair of balanced dipoles

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2344 proves that an arbitrary transported word on one aligned active
axis can preserve the bad inverse character.  The actual owner-pivot
quotient does not translate one coordinate at a time.  Its two character
axes have a canonical, sparse, and different physical meaning:

```text
target a minus its graft unit k_a,
target b minus its graft unit k_b.                  (1)
```

This theorem derives that normal form, identifies the target gain of
every atomic relation as a pair of dipole residues, and rewrites the
remaining inverse-character problem as flatness for a magnetic
connection on the `13 x 13` twist torus.

The reformulation is lossless.  It replaces a generic `169`-variance
question by local nearest-neighbour defects.  It does not prove that a
canonical defect is nonzero.

## 1. The owner-pivot matrix modulo thirteen

Use THM-2309's labelled coordinates

```text
U={six guard/unit labels},
B={j,a,b},
```

where `j` is the selected source owner and `a,b` are the two target
blockers.  Modulo thirteen,

```text
w_u!=0       for u in U,
w_j=w_a=w_b=0.                                     (2)
```

Fix an omitted unit `u_0`, choose distinct graft units

```text
k_a,k_b in U minus {u_0},
```

and let `L` be the row space of THM-2309's six-row packet
`R_(j,u_0)`.  Put

```text
q=w_(u_0)!=0,
alpha_k=w_k/q.                                     (3)
```

After dividing rows by `q`, the packet has the following exact form:

```text
-e_j,

alpha_k e_(u_0)-e_k
    for the three ungrafted units k,

alpha_(k_a)e_(u_0)-e_(k_a)-e_a,

alpha_(k_b)e_(u_0)-e_(k_b)-e_b.                    (4)
```

The six rows are independent.  With

```text
K=w^perp subset F_13^9,
```

THM-2309 gives

```text
K=L direct_sum span(e_a,e_b).                      (5)
```

## 2. Exact dual normal form

Characters of `K/L` are represented by

```text
L^perp/<w>.                                        (6)
```

Let `ell in L^perp`.  Pairing it with (4) gives

```text
ell_j=0,

ell_k=alpha_k ell_(u_0)
    for each ungrafted unit k,

ell_(k_a)=alpha_(k_a)ell_(u_0)-ell_a,

ell_(k_b)=alpha_(k_b)ell_(u_0)-ell_b.              (7)
```

Subtract

```text
(ell_(u_0)/q)w
```

from `ell`.  This is exactly the quotient gauge in (6), and it makes the
`u_0` coordinate zero.  Equations (7) then kill every ungrafted unit
coordinate and leave

```text
ell(s,t)
 =s(e_a-e_(k_a))+t(e_b-e_(k_b)),

(s,t)=(ell_a,ell_b) in F_13^2.                     (8)
```

The representative (8) is unique in the gauge `ell_(u_0)=0`.
Consequently the two target character axes are the disjoint dipoles

```text
D_a=e_a-e_(k_a),
D_b=e_b-e_(k_b).                                   (9)
```

They are also minimum-support representatives of their character axes.
Adding a nonzero multiple of `w` turns on every ungrafted unit
coordinate, while adding zero leaves the two entries in (9).

This is the precise physical action forgotten by the abstract label
`G^ isomorphic to F_13^2`.

## 3. The quotient map is dipole residue

For `r in K`, write its THM-2309 decomposition

```text
r=l+x e_a+y e_b,                 l in L.            (10)
```

Pairing with (8) gives both

```text
ell(s,t).r=sx+ty
```

and

```text
ell(s,t).r
 =s(r_a-r_(k_a))+t(r_b-r_(k_b)).                  (11)
```

Therefore the target quotient map itself is

```text
pi(r)
 =(r_a-r_(k_a), r_b-r_(k_b)).                     (12)
```

Define the dipole residues of any integer mode vector `z` by

```text
delta_a(z)=z_a-z_(k_a) mod 13,
delta_b(z)=z_b-z_(k_b) mod 13.                    (13)
```

For a full atomic relation address from THM-2334,

```text
r_full=u+R beta+m e_d-v,
R=13^(lambda_j+1),
d in {a,b}                                        (14)
```

is the deepest target label.  Since `13|R`, equations (12)--(14) give

```text
pi(r_full)
 =(
    delta_a(u)-delta_a(v)+m 1_(d=a),
    delta_b(u)-delta_b(v)+m 1_(d=b)
   ).                                             (15)
```

Thus the word harmonic `beta` is target-neutral, but the present and bare
endpoint modes enter through **paired differences**, not raw target
coordinates.  Equation (15) is the atomic version of THM-2343's affine
translation law.

## 4. What a target twist does to the physical products

Under the representative (8), the present endpoint factors are translated
in opposite pairs:

```text
a:       +s/13,             k_a:       -s/13,
b:       +t/13,             k_b:       -t/13,       (16)
```

and every other present factor is unshifted.  The transported-word shifts
are

```text
R s/13, -R s/13, R t/13, -R t/13,
```

all integral, so the word factors remain physically unshifted.  The same
dipoles act on the bare endpoint.

Every canonical factor at the four labels in (16) is a nonconstant
centered interval or complement.  Hence THM-2344's one-active-coordinate
hostile is not a literal canonical owner-pivot twist unless its paired
factor is artificially made inert.  This still does not prove
nonconstancy: a multi-factor Fourier coefficient may remain in one
dipole-residue class or two signed residue arrays may be convolution
inverses.

The exact next analytic object is therefore the **dipole-residue width**
of the endpoint coefficient arrays, together with enough sign or phase
information to prevent a shifted convolution inverse.

## 5. Magnetic derivatives of the phase-free response

Relabel the target axes as

```text
(o,d)=(other target, deepest target),
p=(0,m),
zeta=exp(2*pi*i/13).                               (17)
```

Write the phase-free response as `K(s,t)` in the dual dipole coordinates
from (8).  Define two covariant nearest-neighbour derivatives:

```text
(nabla_o K)(s,t)
 =K(s+1,t)-K(s,t),

(nabla_d^m K)(s,t)
 =K(s,t+1)-zeta^(-m)K(s,t).                       (18)
```

The magnetic Dirichlet energy is

```text
E_mag(K)
 =1/169 sum_(s,t)
   (|nabla_o K(s,t)|^2+|nabla_d^m K(s,t)|^2).      (19)
```

It is nonnegative, and

```text
E_mag(K)=0

iff K(s+1,t)=K(s,t)
    and K(s,t+1)=zeta^(-m)K(s,t) for every s,t

iff K(s,t)=c zeta^(-m t) for some c.               (20)
```

By THM-2343, (20) is exactly the zero-only full-current boundary.  The
bad inverse character is therefore the one-dimensional space of
covariantly constant sections of (18).

This is a flat connection: going thirteen steps around the deep cycle
has holonomy

```text
(zeta^(-m))^13=1.                                  (21)
```

The deepest-comb modulation is its gauge trivialization.  Indeed, put

```text
H(s,t)=zeta^(m t)K(s,t).                           (22)
```

Then

```text
H(s+1,t)-H(s,t)
 =zeta^(m t)nabla_o K(s,t),

H(s,t+1)-H(s,t)
 =zeta^(m(t+1))nabla_d^m K(s,t).                  (23)
```

Consequently

```text
E_mag(K)
 =1/169 sum_(s,t)
   (|H(s+1,t)-H(s,t)|^2
    +|H(s,t+1)-H(s,t)|^2).                         (24)
```

The inverse-character obstruction is ordinary constancy after the lawful
deep-phase gauge.

## 6. Exact spectral invoice

Use THM-2343's normalized inverse transforms

```text
A_K(x,y)
 =1/169 sum_(s,t)zeta^(-sx-ty)K(s,t),

A_H(q)=A_K(q-p).                                   (25)
```

Finite Parseval applied to (18) gives

```text
E_mag(K)
 =sum_(x,y)
   [4 sin^2(pi*x/13)
    +4 sin^2(pi*(y+m)/13)]
   |A_K(x,y)|^2                                    (26)

 =sum_q
   [4 sin^2(pi*q_o/13)
    +4 sin^2(pi*q_d/13)]
   |A_H(q)|^2.                                     (27)
```

The bracket in (27) vanishes only at `q=0`.  Thus (27) is a
Laplacian-weighted version of THM-2334's nonzero-target variance, with
the exact comparison

```text
4 sin^2(pi/13) sum_(q!=0)|A_H(q)|^2
 <=E_mag(K)

 <=8 sin^2(6*pi/13) sum_(q!=0)|A_H(q)|^2.          (28)
```

No cancellation is hidden in (26)--(28): every summand is nonnegative.

## 7. The cheapest canonical escape observables

Equations (18)--(20) make any single violated torus edge a complete
certificate of nonzero target survival.  At the trivial twist, the two
cheapest dipole-aligned tests are

```text
D_o=K(1,0)-K(0,0),

D_d=K(0,1)-zeta^(-m)K(0,0).                       (29)
```

If either value is nonzero, then

```text
E_mag(K)>0
```

and some nonzero full target coefficient survives.  Computing (29)
requires only the untwisted current and one physical translation of each
balanced dipole (16), rather than all `169` twists.

The converse is deliberately not asserted: both displayed edges may
vanish while another torus edge detects the obstruction.  A one-point
alternative from THM-2344 is also available:

```text
K(0,1) real  =>  escape,                           (30)
```

because `(0,1)` detects `p` and a nonzero real scalar times
`zeta^(-m)` is not real.

The carrier here is the Cayley square torus with a flat character
connection, not a tournament.  Its vertices are twist states, its edges
are the two owner-pivot dipoles, its preserved predicate is target
nonconstancy, and its lost data under a tournament shadow would include
translation, phase, amplitude, and word-support masks.

No edge in (29) is proved nonzero.  THM-2349 supplies the delayed marked
current on all `165` rows, including repeated-first and strict-resonance
arms, so this magnetic reformulation is now uniform on the whole ledger.
It excludes no row, and LRC(14) remains open.

## 8. Exact companion

The companion reconstructs THM-2309's six-row packet over `F_13`,
verifies the complete `13^3` dual normal form and unique `169`-element
gauge section, proves the quotient kernel equals the row space, and checks
the dipole projection, transported-word invisibility, and atomic
deep-translation formula.  Over an exact finite field containing a
primitive thirteenth root, it builds the
`338 x 169` magnetic incidence matrix, verifies rank `168`, checks its
inverse-character kernel, explicit gauge equivalence with ordinary torus
incidence, and the unique zero spectral character.  Every load-bearing
check raises under ordinary and optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_owner_pivot_dual_dipole_thm2350.py
python3 -O 04-computation/lrc14_owner_pivot_dual_dipole_thm2350.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_owner_pivot_dual_dipole_thm2350.out
```

byte-for-byte after LF normalization.

## 9. Independent audit

An independent proof and exact-computation audit checked the six-row
owner-pivot rank calculation, the quotient/gauge typing, the two balanced
dipole representatives, and the atomic residue formula.  It separately
verified the magnetic gauge transformation, incidence rank `168`, the
one-dimensional covariantly constant space, the spectral weights and
both inequalities in (28).  Ordinary and optimized executions are
byte-identical to the stored transcript, and the source/output hashes
match the metadata.

QED.

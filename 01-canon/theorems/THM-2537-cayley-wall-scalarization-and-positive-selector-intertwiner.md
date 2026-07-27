---
id: THM-2537
title: "Cayley wall scalarization and positive selector intertwiner"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The cyclic-tournament Cayley current of
  THM-2532 has an exact signless-incidence boundary: (I+P_tau)C_tau=P_tau-I.
  Hence the two current values at any occupied-to-empty tau-edge sum to
  minus its root-independent owner weight.  Applied fibrewise to the
  THM-2531 canonical wall on THM-2527/2533's actual owner-marked Boolean
  predecessor, this turns every positive selector/threshold cylinder into
  a strictly negative two-branch Radon correlation, and their finite sum is
  exactly the positive odd scalar O_tau(-4tau).  Equivalently the selected
  boundary rho=delta_s-delta_t has centred Cayley dual
  2/13*1-delta_s-delta_t=-B_tau rho; denominator thirteen is sharp on the
  augmentation lattice, although the two endpoint packets give an integral
  representative modulo constants.  Globally, the positive occupied-to-
  empty and empty-to-occupied boundary transforms have equal mass and their
  difference is both (P_tau-I)(ge) and (I+P_tau)C_tau(ge).  On the anchored
  deep comb this becomes (P_1-I)H=P_1 gamma^+-gamma^- and commutes with the
  lawful target transform.  The canonical selected-source packet is a
  positive rational Boolean subset of the original carrier, retains its
  terminal word and late owner, has all twelve root colours, and therefore
  re-enters THM-2349's carrier lemma to give twelve owner-marked deepest-
  comb triangles.  Untwisted averaging, either endpoint alone, and a fixed
  edge without the selector all have exact hostiles.  The theorem closes
  the same-predecessor phase-to-positive-boundary intertwiner, not the
  chronological step: the selected head is still an empty predecessor,
  not a semantic arrival or target-active terminal event.  No scalar row is
  excluded and LRC(14) remains open.
source: codex-2026-07-27-cayley-wall-scalarization
depends_on:
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2527-owner-weighted-all-mode-odd-bank-and-boolean-cut-coordinate
  - THM-2531-prime-necklace-guard-boundary-selector
  - THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient
  - THM-2533-owner-weighted-phase-and-mixed-gain-radon-ladders
  - THM-2536-deep-comb-selector-flow-target-drift-and-centered-toothpick-boundary
related:
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2534-positive-kakeya-boundary-transform-and-crofton-reconstruction
  - THM-2535-boundary-tooth-clock-intertwiner-and-neutral-collapse
script: 04-computation/lrc14_cayley_wall_scalarization_thm2537.py
output: 05-knowledge/results/lrc14_cayley_wall_scalarization_thm2537.out
script_sha256: 9654d07b1ed4a735201e24239dbae49e4ba7370763eca430a2af44df7906dd7b
output_sha256: 99a84966fc9df2c64697fed415937995a26b92bf31b6e2d5d7c3baf59633dae6
hash_basis: working-tree bytes (LF)
---

# THM-2537 -- the selected wall scalarizes the actual Cayley phase

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2533 retains the uncontracted owner-attached predecessor phase, while
THM-2531 selects a literal occupied-to-empty wall on the same full root
mask.  Their missing composition is the shortest identity in the Cayley
algebra.  If `C_tau` is THM-2532's cyclic-tournament current, then

```text
signless edge incidence after C_tau
  = ordinary cyclic gradient.                                      (1)
```

At a selected `1 -> 0` wall, the right side is exactly `-1`.  Thus the
signed phase already has a positive scalar consumer: take the negative sum
of its two selected endpoint values.  This is pointwise, before Fourier
transform, integration, or autocorrelation.

The result closes a genuine same-ancestry seam.  It deliberately does not
rename the empty endpoint as an arrival.  Both endpoints are inverse
branches at one predecessor horizon.

Throughout, roots lie in `F_13`, vectors use the unnormalised real pairing

```text
<x,y>=sum_(v in F_13)x_v y_v,                                  (2)
```

and all statements hold coefficientwise for rational, real, complex, or
Hilbert-valued profiles whenever the displayed pairing makes sense.

## 1. The Cayley current is a signless-incidence inverse

Use THM-2532's convention

```text
(P_tau p)_v=p_(v+tau),

C_tau=P_tau-P_tau^2+...+P_tau^11-P_tau^12.                    (3)
```

Its Cayley law is

```text
(I+P_tau)C_tau=P_tau-I=C_tau(I+P_tau).                        (4)
```

Define the signless incidence and ordinary gradient of the oriented
thirteen-cycle by

```text
S_tau=I+P_tau,                  D_tau=P_tau-I.                 (5)
```

Then (4) is the exact factorisation

```text
S_tau C_tau=D_tau.                                            (6)
```

Because the cycle has odd length, `S_tau` is invertible.  Thus `C_tau` is
the unique cyclic-tournament current whose signless edge incidence is the
ordinary gradient, modulo the constant kernel of `D_tau`.  This is a typed
discrete Hodge statement, not merely a spectral analogy.

In coordinates, if `t=s+tau`, (6) says for every profile `p`

```text
(C_tau p)_s+(C_tau p)_t=p_t-p_s.                              (7)
```

This local two-coordinate identity is the load-bearing form below.

## 2. Positive global boundary transform

Let

```text
e:F_13->{0,1},
g>=0                                                               (8)
```

where `g` is constant on the root fibre.  In the LRC application `g=g(z)`
is the common late owner or owner--word bit at base point `z`.  Put `q=ge`
and define the two nonnegative directed boundary vectors

```text
B^+_tau(v)=g e_v(1-e_(v+tau)),

B^-_tau(v)=g(1-e_v)e_(v+tau).                                 (9)
```

The first is the packet of occupied-to-empty tails; the second indexes an
empty-to-occupied boundary by its empty tail.  Directly at each root,

```text
B^-_tau-B^+_tau
 =g(e_(v+tau)-e_v)
 =D_tau(ge)
 =S_tau C_tau(ge).                                           (10)
```

Thus the owner-attached signed Cayley phase and the difference of two
positive Boolean boundary packets are the same object after the signless
incidence map.  No inverse, averaging, or target choice is involved.

The two positive packets have equal total mass:

```text
sum_v B^+_tau(v)=sum_v B^-_tau(v).                            (11)
```

This is cyclic telescoping, or equality of the numbers of `1 -> 0` and
`0 -> 1` walls.  If `g` is Boolean, as on the live owner event, then

```text
||D_tau(ge)||_2^2
 =2 sum_v B^+_tau(v)
 =2 sum_v B^-_tau(v).                                        (12)
```

For a general nonnegative weight, (12) holds with `g^2` in place of `g`
inside both boundary packets.  This distinction is the only use of
Booleanity in the energy formula.

Equations (10)--(12) are the global phase-to-positive-boundary
intertwiner.  They retain every boundary, not yet the canonical one.

They also lose no nontrivial root phase.  With
`zeta=exp(2 pi i/13)`, the multiplier of `D_tau` on root colour `a` is

```text
zeta^(a tau)-1,                                               (12a)
```

which is nonzero for every `a!=0`.  Thus all twelve owner-attached modes
of THM-2533 survive in the signed difference of the two disjoint positive
boundary packets.  Their equal untwisted masses explain why the zero root
still cancels.

## 3. One selected wall gives an exact positive scalar

Let `e` be nonconstant and let THM-2531 select, in the oriented slope
`tau`,

```text
alpha=alpha_tau(e),             ell=q_tau(e),

s=alpha+(ell-1)tau,             t=s+tau,                      (13)

e_s=1,                          e_t=0.
```

Centre the owner-weighted profile:

```text
p^g=ge-(g/13)(sum_v e_v)1.                                   (14)
```

Since `C_tau` kills constants, `C_tau p^g=C_tau(ge)`.  Evaluating
(7) on the selected wall gives the pointwise identity

```text
(C_tau p^g)_s+(C_tau p^g)_t=-g,                              (15)

-[(C_tau p^g)_s+(C_tau p^g)_t]=g>=0.                         (16)
```

When `g>0`, (16) is strict.  The selector has converted a signed,
zero-sum, owner-attached Radon current into its exact positive owner mass.
It uses only two current coordinates and the same root mask which defined
that current.

There is a useful dual form.  Put

```text
rho_(s,t)=delta_s-delta_t,

u_(s,t)=2/13 1-delta_s-delta_t.                              (17)
```

Let `B_tau=C_tau^(-1)` on the centred module, with THM-2532's sawtooth
normalisation.  Then

```text
B_tau rho_(s,t)=delta_s+delta_t-2/13 1,

u_(s,t)=-B_tau rho_(s,t),                                   (18)

C_tau u_(s,t)=-rho_(s,t).                                   (19)
```

Indeed `delta_t=P_tau^(-1)delta_s`, and (4) also gives

```text
C_tau(I+P_tau^(-1))=I-P_tau^(-1).                           (20)
```

Applying (20) to `delta_s` proves (18)--(19).  Skew-adjointness now gives

```text
<C_tau p,u_(s,t)>
 =<p,rho_(s,t)>
 =p_s-p_t.                                                   (21)
```

At the selected wall, (21) is one.  Since `C_tau p` has zero total, its
pairing with `u` is also the negative of the sum at the two endpoints,
recovering (16).

### The denominator thirteen is sharp, but only after centring

Let

```text
Lambda={x in Z^13:sum_v x_v=0}.                              (22)
```

THM-2532 proves

```text
C_tau Lambda
 ={y in Lambda:sum_v v y_v=0 mod 13}.                        (23)
```

For the wall boundary,

```text
sum_v v rho_(s,t)(v)=s-t=-tau!=0 mod 13.                    (24)
```

Therefore `rho_(s,t)` is not in `C_tau Lambda`.  Its unique centred
rational inverse is the vector in (18), whose two endpoint coordinates are
`11/13` and whose other eleven coordinates are `-2/13`.  Its exact
denominator is thirteen; multiplying by thirteen is necessary and
sufficient for an augmentation-lattice lift.

There is no contradiction with the positive two-endpoint consumer.
The noncentred integral vector `delta_s+delta_t` already satisfies

```text
C_tau(delta_s+delta_t)=rho_(s,t).                            (25)
```

It has total two and cannot be made integral and centred by adding a
constant vector.  Scalar evaluation only sees the class modulo constants,
so the Boolean union of the two endpoint packets is an integral dual
representative; a centred source potential pays the sharp denominator
thirteen.

## 4. Exact selector-cylinder scalarisation of the positive odd bank

Return to the live THM-2527/2533 event.  At base point `z`, write

```text
e_r(z)=F((z+r)/13),

g(z) in {0,1},

W((z+r)/13)=g(z)e_r(z),                                     (26)

p^W_r(z)=g(z)[e_r(z)-1/13 sum_u e_u(z)].                    (27)
```

The equality in (26) is exactly the root independence of the sufficiently
late owner factor in THM-2533.  In its notation,

```text
(mathcal C_tau W)((z+r)/13)=(C_tau p^W(z))_r.                (28)
```

Let `Lambda_(a,ell)^tau` be THM-2531's selector cylinder and let
`Psi_tau(e)` be THM-2527's integer cut score.  On that cylinder put

```text
s_(a,ell)=a+(ell-1)tau,

t_(a,ell)=a+ell tau.                                        (29)
```

For `1<=j<=98`, define the finite Boolean layer

```text
A_(j,a,ell)
 ={z:g(z)=1,
      Psi_tau(e(z))>=j,
      z in Lambda_(a,ell)^tau}.                              (30)
```

Equation (15) gives, with no cancellation,

```text
mu(A_(j,a,ell))
 =-integral 1_(A_(j,a,ell))(z)
    [(C_tau p^W(z))_(s_(a,ell))
     +(C_tau p^W(z))_(t_(a,ell))] dz.                       (31)
```

The `156` selector cylinders partition the mixed-mask locus and

```text
Psi_tau(e)=sum_(j=1)^98 1_(Psi_tau(e)>=j).                   (32)
```

Consequently THM-2527's fixed positive odd coordinate has the exact
phase-selector expansion

```text
O_tau(-4tau)
 =integral g(z)Psi_tau(e(z))dz
 =sum_(j,a,ell)mu(A_(j,a,ell))                               (33)

 =-sum_(j,a,ell)integral 1_(A_(j,a,ell))
    [(C_tau p^W)_s+(C_tau p^W)_t]dz.                         (34)
```

Every nonnull summand in (34) is strictly negative before the leading
minus sign.  This is stronger than nonvanishing of a Fourier norm: it is a
finite positive scalarisation of the actual charged phase on the exact
THM-2531 Boolean atoms.  On every one of the `165` live rows, THM-2527's
explicit late-owner threshold makes (33) positive.

There is also a literal circle-level Boolean consumer.  For a fixed layer
`A=A_(j,a,ell)`, let

```text
iota_r(z)=(z+r)/13,

A^s=iota_s(A),                 A^t=iota_t(A).                 (35)
```

The two branch packets lie in disjoint thirteenth intervals.  Equations
(28) and (31), including the Jacobian `1/13`, give

```text
integral_(A^s union A^t) mathcal C_tau W(x)dx
 =-mu(A)/13.                                                  (36)
```

Thus an honest Boolean union of two selected branch packets pairs
nontrivially with the owner-attached Radon current.  Untwisted Haar
averaging still gives zero; the selector is the necessary twist.

For a rational owner weight rather than a Boolean owner, THM-2531's finite
threshold expansion of `g` applies first.  Equations (31)--(36) then hold
on each Boolean owner layer and sum with its common rational factor.

## 5. Covariance, reflection, and complement

The construction uses no preferred root.  Under the affine relabelling

```text
r'=U r+H,                       tau'=U tau,

p'_(Ur+H)=p_r,                  U!=0,                          (37)
```

THM-2531 sends `(s,t)` to `(Us+H,Ut+H)`, while direct substitution in (3)
gives

```text
(C_(U tau)p')_(Uv+H)=(C_tau p)_v.                            (38)
```

Therefore (7), (15), the dual packet, and every cylinder mass are affine
covariant.  Translation fixes the slope; multiplication transports it.

Reflection is the case `U=-1`.  It is lawful only with

```text
tau -> -tau,

(s,t)->(H-s,H-t).                                            (39)
```

Keeping `tau` fixed while reflecting the roots is not a covariance law.

Complement has a separate exact local law.  The reverse edge `t -> s` is
an occupied-to-empty `-tau` wall of `1-e`.  Its centred profile is `-p`,
and

```text
C_(-tau)(-p)=C_tau p.                                       (40)
```

Hence the same two current entries again sum to `-1` in the complement
chart.  THM-2531's **canonical** complement wall need not be this reversed
edge; it is the old word tournament's minimum boundary.  The theorem
applies independently to that wall.  The local reversal law must not be
misstated as equality of the two canonical selectors.

## 6. Deep-comb target line: the correct shifted selector identity

In one THM-2530/2536 target-anchored deep cell, let

```text
H_a=integral F e_a,

gamma^+_a=integral F e_a(1-e_(a-1)),

gamma^-_a=integral F e_a(1-e_(a+1)).                         (41)
```

The global boundary packets (9) at `tau=1` are exactly

```text
B^+_1=gamma^-,

B^-_1=P_1 gamma^+.                                          (42)
```

Integrating (10) therefore gives

```text
(P_1-I)H
 =P_1 gamma^+-gamma^-
 =(I+P_1)C_1H.                                               (43)
```

The shift on `gamma^+` is load-bearing.  THM-2536's opposite-selector
divergence is instead

```text
delta=gamma^+-gamma^-.                                      (44)
```

It detects adjacent-pair flow.  Equations (43) and (44) are different
shadows and cannot be substituted for one another.

All operators in (43) act only in the relative root coordinate, so they
commute with THM-2365/2536's lawful Fourier transform across the `169`
target cells.  On every target character `(b,q)`, therefore,

```text
(P_1-I)Hhat(b,q)
 =P_1 gamma^+hat(b,q)-gamma^-hat(b,q)
 =(I+P_1)C_1 Hhat(b,q).                                     (45)
```

This is an exact target-line Radon intertwiner.  It does not force a
nonzero target character: THM-2536's full-packet non-cover hostile makes
all three profiles target-cell invariant.  A scalar-cover/nonlacunarity
input must still force target variation or supply a semantic target-active
endpoint.

For the target-anchored selector chord `a -> 0`, the same local identity
uses the edge slope `tau_edge=-a`:

```text
(C_(-a)p)_a+(C_(-a)p)_0=p_0-p_a=-1.                          (46)
```

Equation (46) scalarizes every selected star edge, but its slope varies
with the marker class and need not be the physical guard slope.  This is
the all-slope diagonal behind THM-2535's clock-chart construction.

## 7. The boundary-selected source packet retains the old carrier

The pointwise bridge has a useful analytic consequence on every live row.
Let

```text
F=1_(E_j) 1_(Q_(j,sigma))(13^k x)
```

be THM-2349's positive rational Boolean carrier, and let `g` be the same
sufficiently late Boolean owner--word block used in THM-2527/2533.  Thus

```text
W=Fg.                                                        (47)
```

Use the unique branch representation `x=iota_r(z)` and define the
canonical selected-source and selected-head packets

```text
S^sel_tau(iota_r(z))
 =g(z)1_(e(z) nonconstant)1_(r=s_tau(e(z))),

T^sel_tau(iota_r(z))
 =g(z)1_(e(z) nonconstant)1_(r=t_tau(e(z))).                 (48)
```

They are rational Boolean step functions.  On every active fibre each has
exactly one root, and

```text
0<=S^sel_tau<=W<=F,

T^sel_tau W=0,

integral S^sel_tau=integral T^sel_tau
 =1/13 integral g 1_(e nonconstant)>0.                       (49)
```

Positivity is THM-2527's positive first threshold layer.  The singleton
root fibre of `S^sel_tau` has nonzero Fourier coefficient in every nontrivial
root colour, so

```text
E_a S^sel_tau!=0                         for all a!=0.        (50)
```

Moreover `S^sel_tau<=W` means that the source packet literally retains

```text
the THM-2349 terminal word sigma,
the same late owner block g,
the old shallow carrier,
and avoidance of the old deepest comb.                       (51)
```

Apply THM-2349's abstract shallow-carrier lemma with its old carrier
`e=F` and the new positive subset `f=S^sel_tau`.  For every prescribed
`kappa in F_13^*`, it gives an owner-marked `91`-unit deepest-comb triangle
whose marked coefficient is a coefficient of `S^sel_tau`.  Hence the boundary-
selected event, not merely the unselected owner event, carries all twelve
old shallow colours and retains the terminal word and late owner.

This is a source-side strengthening.  The packet `T^sel_tau` is disjoint from
`W`; it does not inherit the terminal word.  Reapplying the carrier lemma
to `S^sel_tau` therefore does not turn `T^sel_tau` into a chronological arrival.

If the late owner is refined by THM-2517's seven scheduler cells, intersect
`S^sel_tau` and its marker classes with each cell.  These are exactly the
positive selector/scheduler intersections which feed THM-2535's
clock-holonomy cycle.  Equation (15) holds in every clock cell.  The
remaining THM-2535 warning is unchanged: the scheduler clock is not the old
THM-2449 source/deep residue, and forgetting the chart kills its signed
incidence table.

## 8. Sharp hostiles and what scalar cover must now force

Four elementary controls locate the exact boundary.

### Untwisted averaging is identically blind

Every row of `C_tau` sums to zero, so

```text
sum_v(C_tau p)_v=0                                           (52)
```

for every profile.  A scalar Haar average with no branch twist cannot see
the current even when all twelve root modes survive.  Equation (36) works
because the selected two-branch packet is nonconstant.

### Neither endpoint alone is uniform

At `tau=1`, take the singleton mask `e=delta_0`.  Its selected wall is
`0 -> 1`, but

```text
((C_1e)_0,(C_1e)_1)=(0,-1).                                 (53)
```

Thus the selected source coordinate alone can vanish.  For the adjacent
pair `e=delta_0+delta_1`, the selected wall is `1 -> 2` and

```text
((C_1e)_1,(C_1e)_2)=(-1,0).                                 (54)
```

Thus the selected head coordinate alone can vanish.  The two-endpoint sum
in (7) is sharp.

### A fixed edge can cancel despite complete root spectrum

Take the singleton `e=delta_2` but inspect the fixed edge `0 -> 1`.  Both
bits are zero, so

```text
(C_1e)_0+(C_1e)_1=0.                                       (55)
```

The mask is nonconstant rational and hence has every primitive root mode.
Its canonical wall is `2 -> 3`, where the sum is `-1`.  Complete spectrum
does not replace the selected-wall sidecar.

### Pointwise cover labels are not yet temporal edges

At a selected head `t`, the scalar cover may supply one of several active
danger/repair labels.  THM-2461 proves that four of the five first-failure
roles are target-neutral.  The pointwise cover identity alone contains no
clause requiring the label hit by the selector to be the fifth,
target-active role.  Cross-chamber or prescribed-clock coherence, not a
second Fourier census, must supply that clause.

The remaining semantic target can now be stated as a positive hitting
problem.  Let `A_tar(z,t)` be the indicator that the selected empty head is
carried by the lawful target-active first-failure/terminal role.  A closing
lemma must force

```text
integral g(z)Psi_tau(e(z)) A_tar(z,t_tau(e(z))) dz>0.         (56)
```

The left side has no sign cancellation.  Scalar cover/nonlacunarity must
force its support, or supply an equivalent chronological owner-loop map.
Theorem 2537 proves the same-predecessor phase, selected wall, old carrier,
late owner, and deepest-triangle parts of that statement.  It does not
prove (56), target-cell variation on every covering row, an arrival time,
a row exclusion, or LRC(14).

## 9. Exact dependency-free referee

Run

```bash
python3 04-computation/lrc14_cayley_wall_scalarization_thm2537.py
python3 -O 04-computation/lrc14_cayley_wall_scalarization_thm2537.py
```

Both executions reproduce

```text
05-knowledge/results/lrc14_cayley_wall_scalarization_thm2537.out
```

byte-for-byte.  The script uses only Python's standard-library exact
integers and `Fraction`.  It checks:

- all twelve Cayley/signless-incidence matrix identities and Green
  inverses;
- all `156` adjacent-wall centred inverses and their exact denominator
  thirteen;
- all `98,280=12(2^13-2)` oriented nonconstant mask states;
- the positive global boundary identity, equal-mass/Crofton energy, every
  selected endpoint identity, and root-independent owner scaling;
- all `2,236,416` positive selector-threshold layers and their exact scalar
  sum;
- affine translation, multiplicative, reflection, and complement laws;
- the deep-comb identity on all `23` rays and an independent positive
  combination; and
- the untwisted, one-endpoint, and fixed-edge hostiles.

The exact referee proves the finite algebra.  The uniform live-row
application in Section 7 is the stated composition of THM-2527/2531/2533
with THM-2349's already proved carrier lemma.

**QED.**

The independent audit rederived the signless-incidence factorisation,
selected-endpoint scalar, centred Green dual, and first-moment obstruction.
It separately checked the deep-comb shift in (43), the distinction from
THM-2536's pair-flow divergence, and every hypothesis needed to reapply
THM-2349 to the selected-source subset.  Normal and optimized executions
byte-match the stored transcript and recorded hashes.

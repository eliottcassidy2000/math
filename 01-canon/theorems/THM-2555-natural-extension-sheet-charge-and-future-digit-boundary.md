---
id: THM-2555
title: "Natural-extension sheet charge and the future-digit boundary"
status: >
  PROVED + VERIFIED-EXACT.  On the depth-L base-thirteen natural extension,
  the old target action translates only the most significant ancestry digit.
  Every F_13-equivariant root map carried by the sheet is exactly that digit,
  with its prescribed charge, plus an arbitrary function on the invariant
  lower-sheet quotient.  For a unit physical role, one explicit carry
  correction recovers the old predecessor root pointwise; for a
  thirteen-divisible role the old root is erased sharply.  The genuinely
  future immediate root is instead the first digit of the future base
  T^L y.  It is fixed by the old action and charged by the new future action.
  Exact BV mixing of a selected-head packet with a root-resolved future
  handoff converges to the circular convolution of their separate root
  marginals, so it forces zero displacement exactly when those margins
  overlap at a common label.  Full ancestry-sheet retention does not force
  that overlap: exact digit-cylinder and translation-symmetrized offset-one
  hostiles have positive chronology, full intervening sheet freedom, and zero
  semantic diagonal.  If a genuinely later target-active unit role is
  separately proved on the selected packet and its semantic root is certified
  to be the carry-corrected old-action ancestry root, then the Hall diagonal
  collapses pointwise.  Current canon supplies neither that semantic
  intertwiner nor a future-root overlap theorem; no row is removed.
source: codex-2026-07-27-natural-extension-sheet-charge
depends_on:
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
  - THM-2549-future-pullback-target-neutrality-and-cemetery-hall-boundary
  - THM-2554-translation-quotient-root-displacement-and-endpoint-swap-parity
related:
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary
  - THM-2550-canonical-typed-row-double-nondegeneracy
script: 04-computation/lrc14_natural_extension_sheet_charge_thm2555.py
output: 05-knowledge/results/lrc14_natural_extension_sheet_charge_thm2555.out
script_sha256: 5d197f060e2056e3e2ccb1395fa2cf7da81019b487a5cf8f854b05a4058d51c6
output_sha256: 727fad41eb93eb93bd786b673de2a0b3469118d945935a5d8b4cb290ad82ac29
hash_basis: working-tree bytes (LF)
---

# THM-2555 -- the sheet remembers the old root, not the future root

**PROVED + VERIFIED-EXACT.**

THM-2549 proves that every positive-depth future factor forgets the old
target action on its visible future base.  It also identifies the omitted
natural-extension sheet as the only place where the old charge can remain.
The first question is therefore not whether the sheet is charged, but exactly
what root information its charge determines.

There are two different digits:

```text
y=0.d_1 d_2 ... d_L | e_1 e_2 ...  in base 13;

old target action       -> d_1 is charged;
genuinely future action -> e_1 is charged.                         (1)
```

The bar in (1) is the horizon.  This theorem classifies every old-action
root map on the left of the bar, identifies the physical carry which recovers
the old selected head, and proves that neither the full sheet nor mixing
moves that identity across the bar.

## 1. Exact natural-extension coordinates

Put `p=13`, `T(y)=py mod 1`, and let `L>=1`.  Away from the null set of
ambiguous endpoints, every `y in T` has the unique representation

```text
y=(z+a)/p^L,
z=T^L y in [0,1),             a=floor(p^L y) in Z/p^L Z.    (2)
```

Write the sheet in top-digit/tail coordinates as

```text
a=p^(L-1)d+r,
d in F_p,                     0<=r<p^(L-1).                  (3)
```

The old target action is

```text
O_theta:y -> y-theta/p.                                      (4)
```

Multiplying (4) by `p^L` gives the exact natural-extension action

```text
O_theta(z,a)=(z,a-p^(L-1)theta),

O_theta(z,d,r)=(z,d-theta,r).                               (5)
```

Thus the old action is free on each `d`-fibre and its quotient is exactly
the future base together with the lower sheet:

```text
([0,1) x Z/p^L Z)/F_p
  isomorphic to [0,1) x Z/p^(L-1) Z,        [(z,d,r)] ->(z,r). (6)
```

This is the stalk form of THM-2549's neutrality.  Forgetting `a` forgets
`d`; retaining all of `a` retains one charged digit and `L-1` invariant
digits.

## 2. Classification of every equivariant sheet root

Let `lambda in F_p`.  Call a map

```text
beta:[0,1) x Z/p^L Z -> F_p                                (7)
```

an old-action root of charge `lambda` when

```text
beta(z,a-p^(L-1)theta)
 =beta(z,a)-lambda theta                    for every theta. (8)
```

Then (8) holds if and only if there is an arbitrary map

```text
c:[0,1) x Z/p^(L-1) Z -> F_p                              (9)
```

such that

```text
beta(z,d,r)=lambda d+c(z,r).                               (10)
```

### Proof

Equation (5) shows directly that every map in (10) satisfies (8).  Conversely,
put `theta=d` in (8).  Then

```text
beta(z,0,r)=beta(z,d,r)-lambda d,
```

so (10) holds with `c(z,r)=beta(z,0,r)`.  QED.

For a map carried by the sheet alone, `c` depends only on `r`.  In the
standard charge-one convention there are therefore exactly

```text
p^(p^(L-1))                                                 (11)
```

equivariant sheet maps.  The top digit `d` becomes unique only after imposing
both invariance under every lower-sheet change and one origin normalization;
without those conditions it is a torsor, not a canonical label.

Two roots of the same charge differ by the invariant gauge `c_1-c_2`.
Consequently old-target equivariance alone never fixes their displacement.
The constant gauges `c=0` and `c=1` already give the diagonal and the sharp
offset-one hostile from THM-2554.

There is also a useful algebraic warning.  For `L>=2`, the top digit is not an
additive homomorphism `Z/p^L Z -> F_p`; carries obstruct additivity.  The
additive homomorphisms see the least significant digit, which is fixed by
(5), and therefore cannot carry the old charge.  The relevant object is an
equivariant torsor section, not a quotient-group character.

## 3. The physical carry selects the old predecessor root

Now identify one common physical predecessor chart

```text
x=(u+h)/p,
u=T x in [0,1),                 h=floor(px) in F_p.          (12)
```

Let a role have integer weight `w`, put `y={wx}`, and assume `p` does not
divide `w`.  For the sheet of `y` at any depth `L`, its top digit is

```text
d_w=floor(p{wx})
    congruent w h+floor(wu)                    mod p.        (13)
```

Indeed `p{w(u+h)/p}` differs from `w(u+h)` by a multiple of `p`, while
`wh` is integral.  Therefore the carry-corrected ancestry root

```text
beta_old^(w)(u,a)
 =w^(-1)(d_w-floor(wu))                       mod p         (14)
```

satisfies the pointwise identity

```text
beta_old^(w)(u,a)=h.                                        (15)
```

Under the rolewise old target shift `y->y-theta/p`,

```text
beta_old^(w) -> beta_old^(w)-w^(-1)theta.                   (16)
```

Thus (14) is an exact common-torsor intertwiner after the target parameter is
expressed in physical-root units.  The base carry `floor(wu)` is
load-bearing.  The raw top digit is offset by a generally nonconstant gauge
and does not by itself equal the selected head.

The unit assumption is sharp.  If `p|w`, then for fixed `u`

```text
{w(u+h)/p}={ (w/p)u },                                     (17)
```

which is independent of `h`.  No function of that role's future base and
ancestry sheet can recover the old predecessor root.  More generally, a
role of valuation `v_p(w)=v` begins at the `(v+1)`-st physical digit and has
erased the first `v` digits.

## 4. The genuinely future immediate root is across the bar

In (2), define

```text
e=floor(pz)=floor(p T^L y).                                 (18)
```

In the digit notation (1), `d=d_1` and `e=e_1=d_(L+1)`.  They have different
actions:

```text
old action O_theta:       d -> d-theta,      e -> e;

future action F_phi:
  z -> z-phi/p,
  equivalently y -> y-phi/p^(L+1),
                           d is old-sheet data,
                           e -> e-phi.                       (19)
```

The normalized natural-extension measure is the product of Lebesgue measure
in `z` and uniform counting measure in `a`.  Hence `d` and `e` are independent
uniform `F_p` coordinates.  In particular no map of `a` alone equals the
future immediate root universally.

For a unit physical role, the same carry calculation makes the distinction
coordinate-free.  Put

```text
x_L=T^Lx=(u_L+h_L)/p,
u_L=T^(L+1)x,                 h_L=floor(pT^Lx).              (20)
```

Then the immediate digit of the role's future base is

```text
e_w=floor(p T^L{wx})
    congruent w h_L+floor(wu_L)                mod p,        (21)
```

and its future carry correction recovers `h_L`, not `h`:

```text
beta_new^(w)
 =w^(-1)(e_w-floor(wu_L))=h_L.                              (22)
```

Equations (15) and (22) are the exact chronology ledger:

```text
top sheet digit + old carry       -> old root h;
first future-base digit + new carry -> future root h_L.     (23)
```

Calling the first expression a later root conflates retention of an old
charge with creation of a future target coordinate.

## 5. Exact mixing produces a displacement convolution

The distinction persists after the positive chronology of THM-2549.  Let

```text
A_h:T->{0,1},                    h in F_p,                   (24)
```

be the disjoint rational-BV pieces of a selected-head packet, and let

```text
G_b:T->{0,1},                    b in F_p,                   (25)
```

be disjoint rational-BV pieces of a genuinely root-resolved future packet.
Write

```text
p_h=integral A_h,       q_b=integral G_b,
alpha=sum_h p_h>0,      rho=sum_b q_b>0.                    (26)
```

At delay `N`, form the actual same-base joint table before integration:

```text
C_N(h,b)=integral_T A_h(x)G_b(T^N x)dx.                    (27)
```

The exact BV covariance estimate used in THM-2478 gives

```text
|C_N(h,b)-p_h q_b|
 <= C_(G_b) Var(A_h)/p^N,

C_(G_b)=min(q_b(1-q_b),Var(G_b)/12).                        (28)
```

For displacement `c`, put

```text
m_N(c)=sum_h C_N(h,h+c),
R(c)=sum_h p_h q_(h+c).                                    (29)
```

Summing (28) proves

```text
m_N(c) -> R(c),

sum_c R(c)=alpha rho>0.                                    (30)
```

Thus exact mixing forces some displacement, and every `c` with `R(c)>0`
eventually has positive mass.  It forces zero displacement by this argument
exactly under the additional overlap condition

```text
R(0)=sum_h p_h q_h>0.                                      (31)
```

Neither positivity of the two total masses nor retention of the ancestry
sheet implies (31).  If both root marginals were proved uniform, (31) would
follow, but THM-2554 proves that uniform margins alone do not make an
arbitrary finite-delay coupling diagonal; the mixing input in (28) is the
additional mechanism.

For THM-2549, the fields `A_h` exist but the categorically target-active
future fields `G_b` do not.  Resolving a future base by its immediate digit
produces legitimate mathematical root cells, but does not prove that the
unique THM-2461 target-active role occurs there.

## 6. What the sheet would close, conditionally

Fix a THM-2549 positive same-base packet

```text
Omega_(R,N)=A_R(x)G_R(T^N x)                               (32)
```

on one live row.  Suppose a further positive subpacket `P_R<=Omega_(R,N)`
has all of the following properties:

1. it contains a genuinely later, semantically target-active unit role `w`;
2. the role's depth-`N` ancestry top digit and the old carry `floor(wu)` are
   retained on the same atom;
3. one common root-torsor identification certifies that the THM-2545 later
   root map `b` is the old-action ancestry root (14), rather than merely an
   isomorphic or independently shifted copy.

Then (15) gives pointwise on `P_R`

```text
b=h_R,                                                      (33)
```

and therefore its Hall table is supported on the diagonal with

```text
H_R>=integral P_R>0.                                       (34)
```

No parity, endpoint swap, Fourier census, or further mixing is needed after
those three hypotheses.  If they hold on every row at one common delay, the
finite maximum argument in THM-2549 gives the corresponding all-`165`
positive diagonal statement.

Current canon does not supply hypothesis 1 on the selected-head packet or
hypothesis 3.  THM-2550's nonzero old target fibre is a signed typed-row
control, not an atomwise categorical later role.  Merely attaching (14) to a
neutral future owner or word gives an old-root copy, not semantic arrival.

## 7. Sharp positive hostiles

The failure at the future digit is already an exact rational cylinder.  In
base thirteen, for any `L>=1` take

```text
P_L={d_1=0, d_(L+1)=1}.                                    (35)
```

All intervening digits `d_2,...,d_L` are free, so the complete depth-`L`
sheet is retained.  The packet has

```text
mu(P_L)=1/p^2=1/169>0,

old head=0,        old ancestry root=0,
future root=1,     displacement=1,        diagonal=0.       (36)
```

It is one common physical circle base with strict future chronology and a
genuine future-action root cell.  It is not asserted to be a scalar-cover row
or to carry the THM-2461 target-active semantic role.

There is a rotation-covariant version.  Let

```text
P_L^orb={d_(L+1)=d_1+1}.                                   (37)
```

It has mass `1/p`, leaves every lower sheet digit free, and is invariant
under the abstract simultaneous relabelling

```text
(d_1,d_(L+1))->(d_1+u,d_(L+1)+u).                          (38)
```

Its raw joint table has mass `1/p^2` on every pair `(h,h+1)`, uniform
marginals, and zero diagonal.  This is exactly THM-2554's orbit `O_1` inside
the natural-extension digit model.  The relabelling (38) is an honest
measure-preserving symbolic root action; it is not claimed to be the physical
factorwise action on a live LRC row.

Finally, even on the old sheet itself, replacing the physical root (14) by

```text
beta_old^(w)+c(z,r)                                        (39)
```

preserves old-target covariance for every invariant gauge `c`.  A nonzero
constant `c` gives zero arrival.  Thus the carry identity becomes useful only
when a common-torsor semantic theorem fixes the gauge, exactly as required in
Section 6.

## 8. Relation to the current arrival frontier

- **THM-2538, `THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary.md`:** once a genuine later field shares the ancestry base, cross-Kakeya products recover its whole joint table.  The sheet supplies a candidate old-root coordinate, not that field.
- **THM-2545, `THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile.md`:** equation (33) would restrict the support graph to its diagonal; equation (30) otherwise gives only a displacement margin, and (35) is Hall-perfect off the diagonal.
- **THM-2549, `THM-2549-future-pullback-target-neutrality-and-cemetery-hall-boundary.md`:** its visible future algebra forgets `a`; Sections 1--3 classify exactly what is restored by retaining it.  The restored coordinate carries the old action, not the new future action.
- **THM-2554, `THM-2554-translation-quotient-root-displacement-and-endpoint-swap-parity.md`:** the invariant gauge in (10) is the natural-extension source of its quotient displacement.  The orbit hostile (37) shows that full sheet freedom does not select `d=0`.

The highest-leverage live test is now categorical rather than analytic:
intersect the selected-head packet with the unique target-active unit failure,
retain that role's sheet top digit and physical carry, and decide whether its
semantic target root is the old-action chart (14) or the future-action chart
(22).  The two answers lead respectively to the pointwise bridge (34) and the
margin-overlap problem (31).

No physical arrival field, scalar-row exclusion, or LRC(14) proof follows
from the sheet classification alone.  The live ledger remains `165`.

## 9. Exact companion

Run

```bash
python3 04-computation/lrc14_natural_extension_sheet_charge_thm2555.py
python3 -O 04-computation/lrc14_natural_extension_sheet_charge_thm2555.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_natural_extension_sheet_charge_thm2555.out.
```

The dependency-free referee performs `402,220` exact sheet-action checks,
`5,228,860` equivariant-normal-form checks, `14,560` unit carry identities,
`189,280` old-action covariance checks, and the sharp nonunit controls.  It
then checks all old/future digit pairs through depth four, the displacement
convolution, the depth-one through depth-six `1/169` cylinder hostile, and
the corresponding mass-`1/13` orbit hostile.  The all-depth conclusions are
the symbolic identities above, not extrapolations from the finite loops.

**QED.**

---
id: THM-2584
title: "The b-word depth-five absolute deep-root tensor and two-rail Hall geometry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the
  canonical typed row, wholly inside the lawful sigma={b}, K=2, r=5
  first-collision packet, retaining the arrival collision root v and the
  deepest-speed root t before marginalization gives a positive exact
  13x13x13x7 ancestry tensor.  In the unit chart w=t/2 its support in
  every owner-clock cell is exactly the four-edge path
  v=0--w=0--v=6--w=6--v=12, split between the two affine rails w=v and
  w=v+7.  Its displacement/deep-root pushforward has exactly the two deep
  roots t=0,-1, is positive at both roots for every s!=0 and every owner
  cell, has all 169 global and cellwise F_13^2 characters nonzero, and has
  all 1,183 owner-centred characters nonzero.  The signed theta=t-2v
  contraction retains every nonzero displacement colour.  This is an
  endogenous common-base absolute deep-root reference and a forced-edge
  support geometry, but t is not identified with the THM-2365 target
  co-shift, THM-2334 relation residue, or THM-2545 genuinely later target
  root.  No physical coarse target current, row exclusion, or LRC(14)
  conclusion follows.
source: lrc-semantic-frontier-2026-07-28
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2577-all-clock-word-depth-collision-law-and-support-image-mechanism
related:
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
  - THM-2569-stationary-diagonal-conditioned-paired-corner-and-frozen-future-role-boundary
  - THM-2574-oriented-tooth-component-holonomy-and-fixed-frequency-descent
  - THM-2579-socle-flat-target-torsor-and-integral-difference-filling
  - THM-2581-b-word-depth-five-owner-clock-host-and-reflection-breaking
script: 04-computation/lrc14_b_r5_theta_target_tensor_thm2584.py
output: 05-knowledge/results/lrc14_b_r5_theta_target_tensor_thm2584.out
script_sha256: 99cbb46c4cf3554fb468dda2d6bbf029f3dd7298a4f8b6009d94513484c7849d
output_sha256: 69e2c32d5f320a5110bfe2b7bb2dc3851c1bb4aaa986b70fa8c1aa7479f390a5
hash_basis: working-tree bytes (LF)
---

# THM-2584 -- the depth-five stalk is a two-rail absolute-root ladder

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2581 rebuilds the owner-clock host on the same `sigma={b}`, depth-five
packet where THM-2471 supplies the affine invariant `theta=t-2u`.  Its host
still sums the absolute deepest-speed root.  Retaining that root before the
sum reveals a much smaller positive object:

```text
arrival root -- deepest root support = one four-edge path;

deep-root marginal                  = the adjacent pair {0,-1};

displacement x deep-root spectrum   = all 169 characters;

signed theta current                = all 12 nonzero displacement colours.
                                                                    (1)
```

The path is a literal common-base ancestry support graph, not an externally
aligned coefficient table.  It is also not yet the semantic target graph
needed by THM-2545; the type distinction is stated in Section 7.

## 1. The same-fibre positive tensor

Retain the canonical typed row and owner

```text
w=(1,14,27,40,53,66,13,13^3,2*13^5),       j=1.             (2)
```

Let `E=E_1`, `Q=Q_(1,{b})=E_b`, and use the prescribed packet clock `K=2`:

```text
f=1_Q P_(13^2)1_E,

d=13^5,

U=P_d f,                         V=P_d1_E.                   (3)
```

THM-2577 gives the first-collision depth `r=5`.  On the oriented root
chart put

```text
A(y,v)=U((y+v)/13),              F(y,u)=V((y+u)/13),

s=v-u                            in F_13.                    (4)
```

Thus `u` is the source collision root, `v` is the arrival collision root,
and `s` is their displacement.  Let the native owner-clock cells be

```text
H_ell={y:frac(13y) in [(2ell-1)/14,(2ell+1)/14)},
                                      ell in F_7,            (5)
```

with half-open representatives; endpoints are null.

The deepest speed is `C=2d`.  In THM-2471's arrival stalk

```text
X_(v,a)=(y+v+13a)/(13d),
```

write

```text
z={2y},                 epsilon=floor(2y) in {0,1}.
```

Then

```text
C X_(v,a)=2(y+v)/13=(z+t)/13                 modulo 1,

t=2v+epsilon                                modulo 13,

theta=t-2v=epsilon.                                          (6)
```

The sheet `a` has disappeared because `C/d=2` is integral.  Equation (6)
is pointwise and physical: `theta` is the half-circle digit of the common
base, not a declared complex gauge.

Define the fine positive tensor

```text
K_ell(s,v,t)
 =integral_(H_ell) A(y,v)F(y,v-s)
    1_(t=2v+floor(2y)) dy,                                  (7)
```

and its displacement/deep-root pushforward

```text
C_ell(s,t)=sum_v K_ell(s,v,t).                               (8)
```

Every entry is a nonnegative rational number and comes from one common
`y`-integral.  Expanding the `P_d` and `P_(13^2)` averages recovers the
finite Boolean ancestry stalk of THM-2471 Section 4, so this is a positive
weighted shadow of one actual Boolean fibre product, not a product of
separately marginalized tables.

Summing `t` in (8) gives exactly the THM-2581 table:

```text
sum_t C_ell(s,t)=C^2581_ell(s).                              (9)
```

First-collision disjointness is stronger than the old marginal:

```text
C_ell(0,t)=0                    for every (ell,t).           (10)
```

## 2. A second physical formula

In one summand of (7), substitute

```text
x=(y+v)/13.
```

As `v` ranges over `F_13`, the thirteen chart intervals partition the
physical circle.  Away from their null walls,

```text
v=floor(13x),

y={13x},

t=floor(26x) mod 13,

frac(13y)=frac(169x).                                      (11)
```

The Jacobian and chart sum give the exact alternative formula

```text
K_ell(s,v,t)
 =13 integral_T U(x)V(x-s/13)
      1_(floor(13x)=v)
      1_(floor(26x)=t mod 13)
      1_(frac(169x) in win_ell) dx.                         (12)
```

Summing `v` yields the corresponding formula for `C_ell(s,t)`.  Equations
(7) and (12) are computed independently on the exact grid

```text
T_DEN=182 lcm(w)=297836897838480.                            (13)
```

All `15,379` entries of `K` and all `1,183` entries of `C` agree.  Their
common denominator before reduction is

```text
D_C=13^2 d^2 T_DEN
   =6939029398456584394954868880.                            (14)
```

The independent hostile audit used a third method: a direct sweep of the
`122` exact `y`-events, without the profile cumulative-integral or rotation
routes.  It agrees entrywise with both (7) and (12), including the fine
`v`-refinement.

## 3. Why only two deepest roots occur

The two-root support is structural.  Since `Q=E_b` is dangerous for the
deepest speed,

```text
Q subset D_C.                                                (15)
```

If `U(x)>0`, one `d`-preimage `X` contributing to `P_d f(x)` lies in `Q`.
Since `C/d=2`, equation (15) gives

```text
||2x||<1/14.                                                 (16)
```

In the chart `(z+t)/13`, with `0<=z<1`, condition (16) permits only

```text
t=0:       z<13/14,

t=-1=12:  z>1/14.                                          (17)
```

Every other residue stays at distance at least `1/13` from the nearest
integer.  Thus no exact computation is needed for the exclusion

```text
C_ell(s,t)=0                  when t notin {0,12}.           (18)
```

The exact positive part is sharp:

```text
C_ell(s,0)>0 and C_ell(s,12)>0

                  for every s!=0 and every ell in F_7.      (19)
```

There are therefore exactly `12*2*7=168` positive cells in (8).  Reflection

```text
(s,v,t,ell)->(-s,-v-1,-t-1,-ell)                            (20)
```

preserves `K`, interchanges the two deep roots, and proves their total
masses are equal.  With

```text
I_5=48602521488933856/337437093630814766589,
```

each root has mass

```text
sum_(s,ell) C_ell(s,0)
 =sum_(s,ell) C_ell(s,12)
 =169 I_5/2
 =24301260744466928/1996669193081744181.                    (21)
```

Hence the absolute deep-root marginal is neither uniform nor merely a
difference.  Its Fourier transform is proportional to `1+zeta^h`, which
is nonzero for every `h in F_13` because `13` is odd.

## 4. The four-edge toothpick and its exact drift

Retain the arrival root `v` and reparametrize the deepest root by the unit

```text
w=t/2=7t                    in F_13.                         (22)
```

Equation (6) becomes

```text
w=v+7epsilon.                                               (23)
```

Combining (17) and (23) leaves exactly four possible edges.  The exact
calculation proves every one is positive globally and in every owner cell:

```text
(v,t)=(0,0), (6,0), (6,12), (12,12),

equivalently

v=0 -- w=0 -- v=6 -- w=6 -- v=12.                          (24)
```

This is the depth-five toothpick: two affine rails `w=v` and `w=v+7`
sharing the middle arrival root.  It is a finite self-similar descendant of
the two adjacent danger roots, rather than another complete `13 x 13`
incidence square.

Let `A_*` be the mass of either endpoint edge and `B_*` the mass of either
middle edge.  Their exact functional form is

```text
A_*=12152550159039365/1996669193081744181,

B_*=12148710585427563/1996669193081744181,

A_*-B_*=3839573611802/1996669193081744181>0.                 (25)
```

Thus the two rails have equal total mass `A_*+B_*=169I_5/2`, while the
endpoint-versus-middle profile has a strict positive drift.  Reflection
exchanges the two copies of each weight; it does not force `A_*=B_*`.

There is also an exact Hall-shaped support obstruction internal to (24).
After deleting the equal-label edges `v=w`, the left vertex `v=0` has no
neighbour.  Therefore its positive margin is forced onto `(0,0)` in every
coupling supported on this same four-edge graph.  In each owner cell the
actual tensor already realizes that diagonal, with the uniform floor

```text
sum_s K_ell(s,v=0,t=0)
 >=12738557817004/14840880223020397>0.                       (26)
```

Equation (26) is a forced-edge statement for the arrival/deep-root graph.
It is not yet THM-2545's semantic Hall conclusion, because its two vertex
types are not the selected old head and a genuinely later target-active
role.  No such identification is hidden in the word "diagonal."

## 5. Every joint Fourier colour survives

Let `zeta=zeta_13` and define the normalized joint transform

```text
B_ell(k,h)
 =1/169 sum_(s,t) C_ell(s,t) zeta^(-ks-ht),

B(k,h)=sum_ell B_ell(k,h).                                  (27)
```

Then exactly

```text
B(k,h)!=0                  for all 169 pairs (k,h),

B_ell(k,h)!=0              for all 1,183 triples (k,h,ell). (28)
```

Centre the owner clock:

```text
B^c_ell(k,h)=B_ell(k,h)-B(k,h)/7.                            (29)
```

Every centred entry also survives:

```text
B^c_ell(k,h)!=0             for all 1,183 triples.          (30)
```

The two zero-margin identities are structural:

```text
sum_ell B^c_ell(k,h)=0,

sum_k B_ell(k,h)=sum_k B^c_ell(k,h)=0.                      (31)
```

The second follows from character orthogonality in `s` and (10).  Thus for
every fixed absolute deep colour `h`, (29) is a `13 x 7` doubly centred
common-base host with no zero entry.  A valid exact coordinate denominator
is `7*169*D_C`; hence every nonzero canonical cyclotomic coordinate has the
corresponding reciprocal floor.

All decisions in (28)--(30) reduce integer bucket polynomials modulo
`Phi_13`.  An independent finite-field certificate agrees: the primitive
root `18 mod 79` certifies all `169` global transforms and `1,165` centred
entries; the primitive root `107 mod 131` certifies the remaining `18`.
Neither prime divides the denominator.

## 6. The signed theta contraction is nontrivial in every root colour

The binary affine invariant in (6) supplies a physical orientation.  Put

```text
D_ell(s)
 =integral_(H_ell) sum_u A(y,u+s)F(y,u)
        (2floor(2y)-1) dy.                                  (32)
```

Reflection changes the half-circle sign and gives

```text
D_ell(s)=-D_(-ell)(-s).                                     (33)
```

The global signed current has zero mean, as the two rails in (21) have
equal mass, but no nontrivial displacement colour vanishes:

```text
sum_(s,ell)D_ell(s)=0,

sum_(s,ell)D_ell(s)zeta^(-ks)!=0          for every k!=0.   (34)
```

Cellwise, all thirteen colours survive for `ell=1,...,6`; in the
reflection-fixed cell `ell=0`, exactly the mean vanishes and the twelve
nonzero colours survive.  Thus the proposed `theta=t-2u` contraction is not
a symmetry-zero object on the same `b`, depth-five packet.  It supplies an
oriented collision/deep-root current before any target or physical-frequency
pushforward.

This is stronger than bare nonvanishing of the THM-2581 host and more
specific than an ambient boundary needle: the sign in (32) is endogenous to
the deepest-speed stalk.  It still is not the THM-2512 cut contraction or
the THM-2573 whole-layer endpoint normal.

## 7. Type boundary and the remaining map

The four coordinates retained here are:

| coordinate | definition | retained meaning |
|---|---|---|
| `s` | `v-u` | displacement between source and arrival collision roots |
| `v` | arrival root in `(y+v)/13` | root of the word-restricted current after first collision |
| `t` | `floor(26x) mod 13` | deepest-speed root of that arrival stalk |
| `ell` | owner-clock cell of `frac(13y)` | native `F_7` sidecar |

The source-to-target map is (7), the preserved predicate is one positive
`sigma={b}` first-collision ancestry, and the quotient (8) forgets `v` but
retains its deep root `t`.  The data destroyed by going to THM-2581 are
exactly the four-edge support graph and the affine rail label `theta`.

Three tempting identifications remain unproved and are forbidden:

1. THM-2365's index `t` is a lawful co-shift of the full endpoint packet
   `F_(s,t)`; the present `t` is the deepest-speed root of one collision
   arrival.  Reusing the letter does not construct a temporal/covariant map.
2. THM-2334's `eta.u` is a relation residue before pair-to-coarse
   recombination; it is not the physical digit `t` or the displacement `s`.
3. THM-2545's later root is attached to a genuinely later target-active
   categorical occurrence.  The present `t` belongs to the current pure
   arrival owner `b`, not that later role.

Accordingly, (24) is exactly the kind of sparse support graph that would
force a Hall diagonal **if** a lawful future map transported its vertex
types, and (27) is exactly the kind of absolute charged reference requested
by THM-2579 **if** a map into its Cayley carrier were supplied.  Neither map
is inferred from common cardinality.

The highest-leverage next construction is now precise: retain the
half-circle rail `theta` or the four-edge graph on the THM-2569 common
old-head/future packet, and identify one endpoint-side target coordinate
before the full-`X` annihilating pushforward.  Equivalently, build a
covariant comparison from the present deepest root to the THM-2365 co-shift
or THM-2334 relation residue.  More Fourier support on the present tensor is
unnecessary: it is already saturated.

The row (2) is a typed non-cover control.  Nothing here proves it is a
scalar-cover survivor, transports the construction to the other `164` rows,
constructs a coarse target current, excludes a row, or proves LRC(14).

## 8. Exact companion and audit

Run

```bash
python3 04-computation/lrc14_b_r5_theta_target_tensor_thm2584.py
python3 -O 04-computation/lrc14_b_r5_theta_target_tensor_thm2584.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_b_r5_theta_target_tensor_thm2584.out.
```

The companion reconstructs the Boolean combs, computes (7) and (12), checks
their full fine and coarse equality, recovers the THM-2581 marginal, proves
the first-collision and two-root gates, records all four edge masses, checks
the per-owner forced-edge support, exhausts (28)--(30), verifies the finite-
field certificates, evaluates (32)--(34), and checks both reflection laws.

The independent hostile audit rederived the coordinate and Jacobian signs,
ran the separate direct event sweep, verified the structural danger-gate
proof of `t in {0,-1}`, replayed the exact and modular Fourier certificates,
and extended the direct sweep to the fine `K` tensor.  It specifically
audited that `t` is THM-2471's arrival deepest-speed root and is not
THM-2365's target co-shift.  No mathematical or code defect was found.

**QED.**

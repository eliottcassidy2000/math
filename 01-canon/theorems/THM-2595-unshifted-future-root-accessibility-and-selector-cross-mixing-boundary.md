---
id: THM-2595
title: "Unshifted future-root accessibility and canonical-head cross-mixing"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  For I_h=(h/13,(h+1)/13) and D_(k,L)={y:||ky||<L/14}, positive
  unshifted root access is equivalent to two strict integer inequalities.
  Every root is accessible exactly from k>=12 for an ordinary role and
  k>=10 for the guard; k=11 and k=9 respectively miss the central root 6.
  Above those thresholds all-root access is uniform in every target
  translate s/13.
  A general image-pump lemma says that the simultaneous guard/unit-safe set
  A_0 of any scalar 5+3 cover must have mu(T(A_0))<=3/7.  Exact rational
  wall censuses then rule out every pivot configuration in which no
  all-root role is available after the maximal-nu_7 graft is fixed.  The
  ordinary forced-graft floor is 43/77; the guard case uses the cyclic
  four-arc shape and has floor 1791/2695.  Hence every scalar-cover row
  admits a lawful target-active role k_a, distinct from q_*, whose
  unshifted danger meets every physical root cell.  Elementary exact
  cross-mixing with any positive THM-2537 canonical selected-head layer
  produces, at every sufficiently late time, a positive same-root later
  k_a occurrence and therefore gives a positive delayed same-root analogue
  of THM-2537 equation (56).  This closes the canonical root-alignment and
  later-role hitting seam, not a future terminal word/full-X packet, the
  paired blocker co-shift, left endpoint residue, THM-2334 current, row
  exclusion, or LRC(14).
source: common-endpoint-seam-2026-07-28
depends_on:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2445-twenty-four-cell-graft-owner-conditioning
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
related:
  - THM-2368-owner-pivot-root-fibre-radon-invertibility
  - THM-2558-sparse-owner-fibre-all-slope-target-role-forcing-boundary
  - THM-2561-positive-interval-physical-six-comb-blind-necklace-hostile
  - THM-2565-target-active-self-return-and-future-root-overlap-audit
  - THM-2583-self-similar-digit-needle-internalization-and-carrier-boundary
script: 04-computation/lrc14_unshifted_root_accessibility_thm2595.py
output: 05-knowledge/results/lrc14_unshifted_root_accessibility_thm2595.out
script_sha256: e2630a86565533d0722532ca08fb037011b1238dc00b5848a98bf7918dda0d7e
output_sha256: 93cc78455e7dc2c36b652b9a877ad3a02b641e16b80de5979016a1f994290541
hash_basis: working-tree bytes (LF)
---

# THM-2595 -- every canonical selected head has a later target-active return

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2537 reduced the chronological part of the scalar frontier to one
positive question: does its canonical empty head ever meet the unique
target-active first-failure role?  The obstruction seemed to be root
alignment.  A small role can miss an entire physical first-digit cell, and
the complete shifted target atlas does not repair that unshifted hole.

The missing coordinate is global.  Scalar cover does not permit an arbitrary
small six-role tuple: its simultaneous guard/unit-safe set has to map into
only three blocker danger sets after multiplication by thirteen.  Exact
root-fibre geometry shows that this is incompatible with every configuration
having no pivot-eligible all-root role.  Thus one may choose a target-active
role whose unshifted danger reaches **all thirteen** possible canonical
heads.  Elementary interval cross-mixing then supplies the later hit.

The proof has three separately typed parts:

```text
one-tooth interval geometry                         PROVED;
scalar-cover image inclusion                        PROVED;
bounded small-role and cyclic-arc wall censuses     FINITE-EXACT.          (1)
```

The finite-exact part is exhaustive over the universes stated below; it is
not a numerical sample or an extrapolation in speed.

## 1. Exact unshifted root accessibility

Put

```text
T(x)=13x mod 1,                    I_h=(h/13,(h+1)/13),

d_L(y)=1_(||y||<L/14),             L in {1,2},

D_(k,L)={y:d_L(ky)=1},

M_(k,L)(h)=1_(mu(I_h intersection D_(k,L))>0).              (2)
```

The danger teeth have centres `j/k` and radius `L/(14k)`.  Hence

```text
M_(k,L)(h)=1
```

if and only if some integer `j` satisfies the two strict inequalities

```text
182j-13L < 14k(h+1),

182j+13L > 14kh.                                         (3)
```

It is enough to take `0<=j<=k`.  These inequalities are just the two
conditions that the open interval about `j/k` overlap the open digit cell;
strictness correctly discards tangencies of measure zero.

For the translated gate

```text
D^s_(k,L)={y:d_L(ky+s/13)=1},
```

the corresponding exact inequalities are

```text
182j-14s-13L < 14k(h+1),

182j-14s+13L > 14kh.                                    (3a)
```

Translation moves the tooth grid but does not change its gap length.

Every safe gap between adjacent teeth has exact length

```text
(7-L)/(7k).                                               (4)
```

Since `|I_h|=1/13`, (4), followed by the finite lower-speed check in (3),
gives the sharp thresholds

```text
M_(k,1)=1 on all F_13    iff k>=12,                         (5)

M_(k,2)=1 on all F_13    iff k>=10.                         (6)
```

In particular, the positive directions are uniform in the target translate:

```text
k>=12 (ordinary) or k>=10 (guard)
  implies mu(I_h intersection D^s_(k,L))>0
  for every (h,s) in F_13^2.                              (6a)
```

For reference, the exact subthreshold masks are:

| `k` | ordinary `L=1` accessible roots |
|---:|---|
| 1 | `0,12` |
| 2 | `0,6,12` |
| 3 | `0,4,8,12` |
| 4 | `0,3,6,9,12` |
| 5 | `0,2,5,7,10,12` |
| 6 | `0,2,4,6,8,10,12` |
| 7 | `0,1,3,5,7,9,11,12` |
| 8 | `0,1,3,4,6,8,9,11,12` |
| 9 | `0,1,2,4,5,7,8,10,11,12` |
| 10 | `0,1,2,3,5,6,7,9,10,11,12` |
| 11 | `0,1,2,3,4,5,7,8,9,10,11,12` |

| `k` | guard `L=2` accessible roots |
|---:|---|
| 1 | `0,1,11,12` |
| 2 | `0,5,6,7,12` |
| 3 | `0,3,4,8,9,12` |
| 4 | `0,2,3,6,9,10,12` |
| 5 | `0,2,4,5,7,8,10,12` |
| 6 | `0,1,2,4,6,8,10,11,12` |
| 7 | `0,1,2,3,5,7,9,10,11,12` |
| 8 | `0,1,3,4,5,6,7,8,9,11,12` |
| 9 | `0,1,2,3,4,5,7,8,9,10,11,12` |

The central failures are sharp by one rational unit on each side:

```text
I_6 subset (71/154,83/154) subset D_(11,1)^c,

I_6 subset (29/63,34/63)   subset D_(9,2)^c.                (7)
```

Thus an arbitrary root label cannot be attached to a small **unshifted**
role by mixing alone.

## 2. Positive rational cross-mixing needs only access

Let `A` be a positive rational Boolean subset of `I_h`, and suppose

```text
B=I_h intersection D_(k,L)
```

has positive measure.  Then

```text
integral_T A(x)1_B(T^N x)dx>0                              (8)
```

for every sufficiently large `N`.

Indeed choose open intervals `J_A subset A` and `J_B subset B`, and a point
`y in J_B`.  The inverse grid

```text
{(r+y)/13^N:0<=r<13^N}                                    (9)
```

has mesh `13^(-N)`.  Once that mesh is shorter than `J_A`, one grid point
lies strictly inside `J_A`.  A small neighbourhood of it remains in `J_A`
and maps into `J_B`, proving positive measure in (8).  Every sufficiently
large delay works; no limiting-correlation sign is being inferred.

At every point of (8),

```text
floor(13x)=h=floor(13T^N x),

d_L(kT^N x)=1.                                             (10)
```

This is the exact unshifted same-root chronological event needed below.

## 3. The scalar-cover image bound

Return to THM-2198's scalar `5+3` cover.  Let

```text
A_0=C_H intersection intersection_(i=1)^5 D_(q_i)^c        (11)
```

be the simultaneous guard/unit-safe set, where the guard danger has width
`L=2` and every ordinary danger has width `L=1`.  The cover identity gives

```text
A_0 subset D_(c_1) union D_(c_2) union D_(c_3)       a.e.  (12)
```

Write `c_j=13b_j`.  Since

```text
||c_jx||=||b_jT(x)||,
```

applying `T` to (12) gives the exact image-pump inclusion

```text
T(A_0) subset D_(b_1) union D_(b_2) union D_(b_3)     a.e. (13)
```

Each set on the right has measure `1/7`.  Therefore every scalar cover must
obey

```text
mu(T(A_0))<=3/7.                                          (14)
```

Equivalently, on the inverse-branch base put

```text
n(z)=#{r in F_13:(z+r)/13 in A_0}.                         (15)
```

Then `T(A_0)={z:n(z)>0}` up to the finite wall set, so (14) is a constraint
on the support of a thirteen-root occupancy function, not merely on the
physical mass

```text
mu(A_0)=1/13 integral n(z)dz.                              (16)
```

This distinction is the decisive global coordinate.

## 4. Some all-root role must exist

Call a guard/unit role **big** when its unshifted accessibility mask is all
of `F_13`.  By (5)--(6), this means

```text
H>=10                  for the guard,

q_i>=12                for an ordinary role.               (17)
```

If no role is big, then

```text
H in {1,3,5,7,9},       q_i in {1,...,11},                  (18)
```

with the five `q_i` distinct; `H=q_i` is allowed.  There are exactly

```text
5 binom(11,5)=2,310                                     (19)
```

such unordered typed tuples.  The companion embeds their walls in an
expanded `221`-cell common refinement which also includes the optional
`q=12,H=11` cross-check roles.  This raw refinement count is an
implementation coordinate, not an invariant; adjacent/equivalent cells may
be coalesced.  Exhausting the typed tuples with exact fractions gives

```text
min mu(T(A_0))=171/245>3/7,                               (20)
```

with the unique minimum

```text
H=1,                {q_i}={3,4,5,7,11}.                   (21)
```

Equations (14) and (20) contradict each other.  Thus every genuine scalar
cover has at least one big role.

As an independent bounded cross-check, impose pairwise distinctness on all
six direct coefficients and let them range through `1,...,12`.  Across the
resulting

```text
6 binom(11,5)=2,772
```

tuples, the same image minimum is `171/245`, while the separate physical
mass minimum is

```text
min mu(A_0)=19841/97020
```

at `H=11`, `{q_i}={1,5,7,8,9}`.  This second census is corroboration; (18)--
(20), which allows a guard/unit tie, is the one used in the proof.

## 5. A big role can be made pivot-eligible

THM-2445 permits `q_*` to be chosen among the nonempty set

```text
Max_7={guard/unit roles of maximal nu_7}.                    (22)
```

THM-2309/2350 permits any distinct second graft role `k_a` and then an
omitted label `u_0` distinct from both.  The five roles other than `q_*`
may be put in any first-failure order.  Therefore a big target-active
`k_a` exists unless the set of big roles and the maximal set are the same
forced singleton:

```text
Big={b}=Max_7.                                               (23)
```

This elementary criterion was exhausted over all `63*64=4,032` nonempty
`Max_7`/`Big` pairs and all choices of `(q_*,k_a,u_0)`.

It remains to rule out (23).  There are two role types.

### 5.1 A forced ordinary role

Suppose `b=q_*` is ordinary.  The other roles are one small guard

```text
H in {1,3,5,7,9}
```

and four distinct small ordinary speeds in `{1,...,11}`.  Temporarily omit
`q_*`, let `B` be the set safe for those five roles, and put

```text
m(z)=#{r:(z+r)/13 in B}.                                    (24)
```

An ordinary `13`-unit danger gate occupies at most two direct roots on any
fibre.  Consequently

```text
m(z)>=3  implies z in T(A_0),                               (25)
```

regardless of the magnitude, residue, or phase of `q_*`.

The exact `5 binom(11,4)=1,650`-tuple wall census gives the sharp bound

```text
mu({z:m(z)>=3})>=43/77>3/7.                                 (26)
```

Equality occurs uniquely at

```text
H=1,                  Q={3,4,5,11}.                         (27)
```

Its exact occupancy histogram has low-root mass

```text
mu(m<=2)=12/385+158/385=34/77,                              (28)
```

so (26) follows with equality.  Equations (14), (25), and (26) contradict
(23) when the forced role is ordinary.

### 5.2 A forced guard role

Now suppose `b=q_*=H` is the guard.  The remaining five ordinary speeds are
distinct elements of `{1,...,11}`.  Let `S(z)` be their safe-root set.  The
exact wall census first gives

```text
|S(z)|>=3                                                   (29)
```

on every open cell.

Cardinality alone is not enough: the measure on which `|S|>=5` has sharp
minimum `137/441<3/7`.  The missing invariant is **shape**.  Put

```text
u=H mod 13,                      u!=0.
```

At fixed `z`, the guard danger phases are a uniform thirteen-grid in
`u`-order.  Since

```text
3/13 < 2/7 < 4/13,                                          (30)
```

the guard-danger roots form three or four consecutive vertices in that
cyclic order.  In root coordinates they are contained in one of the
thirteen four-blocks

```text
{a,a+u^(-1),a+2u^(-1),a+3u^(-1)}.                          (31)
```

If the full safe fibre is empty, `S(z)` must be contained in a block (31).
The exact census over

```text
binom(11,5)*12=5,544                                        (32)
```

ordinary-speed/residue cases gives

```text
mu({z:S(z) fits in a block (31)})<=904/2695.                (33)
```

The maximum occurs only for

```text
{q_i}={1,4,5,7,11},             u in {4,9}.                 (34)
```

Therefore, independently of the guard phase and magnitude,

```text
mu(T(A_0))>=1-904/2695
             =1791/2695
             >3/7.                                         (35)
```

Again (14) is contradicted.  Thus (23) is impossible for both role types.

Combining Sections 4 and 5 proves the uniform selection theorem:

```text
on every scalar-cover row there are lawful distinct q_* and k_a such that

q_* has maximal nu_7,

k_a is the target-active second graft,

M_(k_a,L_a)(h)=1 for every h in F_13.                        (36)
```

Choose `u_0` after these two labels; four labels remain available.  No
visible-height assertion is used.

## 6. Canonical-head application and the delayed equation-(56) analogue

Fix a scalar-cover row and the choices in (36).  THM-2537 has positive
Boolean threshold layers

```text
A_(j,a,ell)
 ={z:g(z)=1, Psi_tau(e(z))>=j,
       z in Lambda_(a,ell)^tau},                            (37)

t_(a,ell)=a+ell tau.                                       (38)
```

Their total mass is its positive odd scalar, so at least one layer is
nonnull.  Lift such a layer to its canonical selected-head branch:

```text
A^head={iota_t(z):z in A_(j,a,ell)},

iota_t(z)=(z+t)/13.                                        (39)
```

This is a positive rational Boolean subset of `I_t`.  By (36),

```text
B_t=I_t intersection D_(k_a,L_a)
```

has positive measure.  Section 2 therefore gives, for every sufficiently
large `N`,

```text
integral A^head(x)1_(B_t)(T^N x)dx>0.                       (40)
```

Put `k_a` first in THM-2445's freely ordered five-role bank.  Its first
cell has no earlier safety factors and is exactly

```text
d_(L_a)(k_a y).                                            (41)
```

THM-2461 identifies `k_a` as the unique target-active member of this bank.
Hence (40) is a genuinely later target-active occurrence whose immediate
physical root is the old canonical head `t`.

More formally, define on the head branch

```text
A_tar,N(z,t)
 =1_(floor(13T^N iota_t(z))=t)
  d_(L_a)(k_aT^N iota_t(z)),                               (42)
```

with the first-role typing (41).  Since `dz=13dx`, (37) and (40) give

```text
integral g(z)Psi_tau(e(z))
  A_tar,N(z,t_tau(e(z))) dz

 >=13j integral A^head(x)1_(B_t)(T^N x)dx
 >0.                                                       (43)
```

Equation (43) is the positive delayed same-root analogue of THM-2537 equation
(56), with an explicit chronological realization of its target-active
**first-failure role**.  The argument is rowwise and works at every
sufficiently large delay beyond the existing clocks.

This does not say that the selected empty head emits the old terminal word,
or that the future point carries a complete word/owner/full-`X` packet.  The
old word remains source-side provenance exactly as in THM-2537.  Thus (43)
instantiates the chronological first-failure interpretation of `A_tar`; if
`A_tar` is required by a downstream theorem to include a terminal word or
full endpoint packet, that stronger field is still absent.  What is new is a
genuinely later target-active occurrence on the same physical root.

There is a translated corollary.  Fix any prescribed `s in F_13` and replace
`B_t` in (40) by

```text
B^s_t=I_t intersection D^s_(k_a,L_a).                       (43a)
```

Equation (6a) and the same cross-mixing proof give a positive later
same-root occurrence of

```text
d_(L_a)(k_aT^N x+s/13)                                     (43b)
```

for every sufficiently large delay.  Thus the one-factor future event can
prescribe any target translate.  This may feed a later edge-charge
construction, but it does **not** oppositely translate the paired blocker,
identify a relation character, or supply a charged endpoint current.

## 7. Why the shifted digit needle was insufficient

THM-2583 proves that the complete target-shifted boundary atlas pierces every
old/future digit pair.  That statement cannot replace (2).  At the sharp
ordinary hostile,

```text
I_6 intersection D_(11,1)=empty,                            (44)
```

but the full shifted atlas has exactly `22` labelled boundaries in `I_6`.
Likewise the guard hostile has `18` shifted boundaries there while its
unshifted slice is empty.  The shifted construction chooses some target
label `s_0`; equation (42) requires the canonical unshifted `s=0` role.
Transporting one as the other without covariance would repeat MISTAKE-266.

## 8. The fixed-head hostile and why global coherence kills it

There is a sharp local configuration showing that Sections 3--5 are
load-bearing.  Take

```text
H=3,              (q_1,...,q_5)=(1,2,4,5,7),

q_*=7,            13/56<z<15/56.                           (45)
```

The six direct danger masks are

```text
H : {0,4,8,9},       1 : {0,12},       2 : {0,6},

4 : {3},             5 : {5,10},       7 : {9,11}.         (46)
```

Thus the safe mask is `{1,2,7}`.  The guard slope is `tau_H=4`, and the
THM-2531 lexicographic selector chooses the wall

```text
7 -> 11.                                                     (47)
```

The five roles other than `q_*` have combined unshifted accessibility

```text
F_13 minus {1,11};                                          (48)
```

so the canonical head is inaccessible and `q_*` is its sole current
failure.  Around `z=1/4`, blockers

```text
(c_1,c_2,c_3)=(52,169,13^5)
```

even realize a strict local depth `(1,2,5)` owner fibre: `c_1` is dangerous
and `c_2,c_3` are safe.

Nevertheless (45) is not a scalar cover.  Its complete occupancy function
has `46` walls, `45` open cells, and exact histogram

```text
n=3:389/980,   4:451/980,   5:1/21,   6:1/42,

n=7:1/28,      8:1/140,     9:2/245, 10:1/49.              (49)
```

In particular `n(z)>=3` everywhere off the walls,

```text
mu(T(A_0))=1,                    mu(A_0)=226/735.            (50)
```

Equation (14) rules it out for **every** blocker choice.  The same global
test kills the earlier named physical controls:

```text
THM-2558 sample:
  366 walls, min n=2, mu(T(A_0))=1, mu(A_0)=107297/329868;

THM-2561 primitive row:
  1964 walls, min n=1, mu(T(A_0))=1,
  mu(A_0)=5785336/17784855.                                 (51)
```

These examples remain valid local semantic hostiles, but none can refute a
global scalar-cover conclusion.  The image coordinate is what they omitted.

## 9. Exact stopping boundary

The theorem closes the following chain on every scalar-cover row:

```text
positive canonical THM-2537 selected-head layer

  + globally forced pivot-eligible all-root target-active role

  + exact later cross-mixing

  -> positive same-root chronological target-active hit

  -> positive delayed same-root analogue of equation (56).  (52)
```

It does **not** attach a future terminal word/owner/full-`X` packet,
co-translate the blocker paired with `k_a`, restore the left endpoint
residue absent from THM-2563, construct the full
`e_(c_3)-e_(k_b)` target direction, produce a THM-2334 relation current,
exclude a scalar row, or prove LRC(14).  The next seam is no longer the
existence or root alignment of the canonical later target role.  It is the
lawful paired-action/endpoint lift on the positive packet (43).

## 10. Exact companion

Run

```bash
python3 04-computation/lrc14_unshifted_root_accessibility_thm2595.py
python3 -O 04-computation/lrc14_unshifted_root_accessibility_thm2595.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_unshifted_root_accessibility_thm2595.out.
```

The dependency-free referee uses only integers and `Fraction`.  It checks:

- two independent access-mask routes on all `5,200` cells with
  `1<=k<=200`;
- all `64,220` above-threshold translated `(k,h,s)` access cells;
- both sharp thresholds, complete subthreshold tables, and shifted-atlas
  hostiles;
- all `4,032` abstract maximal/big pivot cases;
- the local fixed-head masks, canonical selector, and strict owner fibre;
- the exact image histograms of all three named local controls;
- the `2,310` no-big and `2,772` pairwise-distinct small-role universes;
- the `1,650` forced-ordinary capacity cases and its exact extremal
  histogram;
- all `5,544` forced-guard cyclic four-block cases; and
- `700` exact interior preimage-grid witnesses for cross-mixing.

The analytic arguments are (3), (8)--(10), (13)--(16), (25), and the cyclic
arc containment (30)--(31).  The stated bounded minima are the exact finite
certificates, pending independent hostile audit.

**QED (candidate pending independent audit).**

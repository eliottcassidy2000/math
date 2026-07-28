---
id: THM-2698
title: "Central half-odometer full local cycle and semantic sidecar boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Enlarging the
  canonical translation lattice from C_(13^6) to
  C_(2*13^6) adds a central C2 half-turn.  Odd affine handoffs induce the
  changed delayed base y -> {13y+1/2}; its strict fixed phase 11/24 supports
  an explicit nonconstant owner-glued two-cycle retaining the physical rail,
  dynamically typed present factor, full delayed word, predecessor carry,
  future half-digit, private root, and primitive unit.  The three-state open
  cylinder has exact length 1/25441508953749252.  The two odd edges have zero
  C2 holonomy and compose to ordinary integer-lattice D^2 endpoint loops, but
  their forced integer-lift midpoint lies outside the delayed word.  Thus the
  remaining debt is a lawful semantic C2 bibundle for elementary chronology
  or a direct two-step endpoint-current cospan; neither is constructed, and
  no row or LRC(14) conclusion follows.
source: root/central-half-odometer-carrier-2026-07-28
depends_on:
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2623-guard-safe-danger-cospan-and-residual-unit-wall
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2689-affine-clock-support-typing-tradeoff-and-odometer-phase-locality
  - THM-2693-odometer-skew-product-three-event-escape-and-uniform-delayed-depth-four-nilpotence
related:
  - THM-2644-odd-torsor-purity-return-gate-and-nonlinear-fixed-branch-decoder
  - THM-2694-mixed-dilation-slope-seven-present-unit-long-word-and-first-gap
  - THM-2695-secondary-kummer-bockstein-picard-divisibility-spectrum-and-prime-alignment-boundary
  - THM-2697-filtered-affine-handoff-germ-category-and-base-signature-holotopy-boundary
script: 04-computation/lrc14_central_half_odometer_full_local_cycle_thm2698.py
output: 05-knowledge/results/lrc14_central_half_odometer_full_local_cycle_thm2698.out
script_sha256: 45cc393a856c00342fdf84875a0bc5a6d4c3df196ab35bb9ac2aad3cfc966c25
output_sha256: 9d496e21ca0ee8abda3a6a6d828b7a09a2345d676105239e1e86d5d3d8de05ec
hash_basis: LF-normalized bytes
---

# THM-2698 -- a central half-odometer restores a full local cycle

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2693 proves that every integer odometer lift leaves the inherited delayed
base on the autonomous map `y -> {13y}`, where four successive delayed event
states are impossible.  A one-base-thirteen-digit prolongation repairs the
speed-`14` part of the sparse carrier but leaves the high-speed factor fixed,
so the full word still dies.  The smallest operation which separates those
two factors is instead the central half-turn.

This theorem proves that the resulting changed base is not merely an abstract
support trick.  It has an exact circle-affine two-cycle whose event points lie
strictly inside the full current local packet.  It also isolates why that
positive cycle is not yet a semantic LRC transition.

## 1. The doubled affine lattice and its exact digit law

Put

```text
p=13,                    R=p^6,
D(x)={px},               tau_n(x)={x+n/(2R)}.            (1)
```

For odd `n=2k+1`, define the post-dilation handoff

```text
F_n=tau_n after D.                                        (2)
```

On the delayed coordinate

```text
pi(x)={Rx},                                                (3)
```

one has exactly

```text
pi(F_n x)={13 pi(x)+1/2}=:B_1(pi(x)).                     (4)
```

Thus an odd lift leaves the kernel class closed by THM-2693.  It is not one
more integral decoration of `B_0(y)={13y}`.

The existing THM-2640 labels nevertheless type the odd translation exactly.
On the pure `tau_n` leg write

```text
c=floor(Rx) mod13,
y={Rx},
d=floor(26y)=2h+kappa,
b=floor(d/13),
e=1 for the left private half-tooth and 0 for the right,
a_e=2c+b+e mod13.                                        (5)
```

Translation by `tau_n`, with `n=2k+1`, gives

```text
c'=c+k+b,
d'=d+13 mod26,
b'=1-b,
a'_e=a_e+n mod13.                                        (6)
```

Indeed adding `k+1/2` to `Rx` proves the first three formulas, and

```text
2(c+k+b)+(1-b)+e=(2c+b+e)+(2k+1) mod13                  (7)
```

proves the root law.  The retained edge `e` is load-bearing.

For the full handoff `(2)`, apply `(6)` at the intermediate point `D(x)`.
The old depth-zero carry is first shifted out by `D`; directly,

```text
c_target=k+floor(13y+1/2) mod13.                          (8)
```

Confusing `(6)` with `(8)` would incorrectly transport the old carry through
dilation.

This does not contradict THM-2657.  That theorem classifies root lifts with a
uniform predecessor-carry increment.  In `(6)` the increment depends on the
retained half-bit `b`, so the odd lift becomes a fixed rule only on the
enlarged `(c,h,kappa,e)` chart.

## 2. Why the half phase repairs the full delayed word

The map `B_1` fixes

```text
y_-=11/24,                  y_+=13/24.                    (9)
```

Both are strict points of raw delayed sector zero and of every nonzero
delayed clock cut `1,...,6`.  At `11/24`, the full delayed-factor distances
from the nearest integer are

```text
target 13^3:                 1/24,
ordinary 14,27,40,53,66:     5/12,3/8,1/3,7/24,1/4,
high 2*13^5:                 1/12.                        (10)
```

Thus the target is strictly dangerous and every displayed guard is strictly
safe.  Reflection gives the same conclusion at `13/24`.

The mechanism is exactly the relative phase missing from the one-deeper
base-thirteen scout.  If `A=13^3` and `c3=338A`, then

```text
A/2=1/2 mod1,                   c3/2=371293=0 mod1.       (11)
```

The central half-turn moves the odd target tooth by one half while fixing the
even-multiple high-speed phase.  By contrast a `C_(13^7)` phase `s/13`
leaves both `A` and `c3` phases integral; it may repair speed `14`, but the
high-speed killer remains.  Equation `(11)` explains, rather than merely
records, the first full-word escape.

## 3. Exact full local two-cycle

An exhaustive source-one scan of the two fixed delayed fibres found

```text
332668 strict rail+present+delayed+private-unit charts,
281738 with nonconstant shallow-to-owner edge,
51 nonconstant phase/clock-edge types.                   (12)
```

The first compatible reciprocal pair lies over `y=11/24`:

```text
x0=55232507/115843416,
rail j=8=(source1,owner4,deep12),
(shallow,owner,carry,h,kappa,edge,root)=(1,4,3,5,1,left,7),

x1=58313459/115843416,
rail j=3=(source1,owner1,deep0),
(shallow,owner,carry,h,kappa,edge,root)=(4,1,1,5,1,left,3). (13)
```

Every listed factor is strict: rail, dynamically typed present factor, full
sector-zero delayed clock word, predecessor carry, future half-digit, private
half-tooth, nonzero root, and primitive `Phi_7` unit.

Take

```text
k0=1472973,             n0=2k0+1=2945947,
k1=4502560,             n1=2k1+1=9005121.               (14)
```

Then exact rational arithmetic gives

```text
tau_(n0) D(x0)=x1,             tau_(n1) D(x1)=x0.        (15)
```

The root increments are

```text
n0=4 mod13,                    n1=8 mod13.                (16)
```

Both intermediate dilation points have left private root `12`, and hence

```text
12+4=3 mod13,                  12+8=7 mod13,              (17)
```

which are precisely the target roots in `(13)`.  Formula `(8)` gives target
carries `1,3`, also matching `(13)`.  The owner of each event equals the
shallow clock at the next event, so the displayed clock gluing is literal.

## 4. The cycle has positive open cylinders

The exact minimum symmetric packet slack at both event points is

```text
s0=s1=1/301082946198216,                                 (18)
```

with the delayed prefix binding.  One handoff expands perturbations by `13`,
so the strict three-state initial cylinder around `x0` has radius

```text
rho=min(s0/13^2,s1/13)
   =1/50883017907498504                                   (19)
```

and length

```text
2rho=1/25441508953749252>0.                               (20)
```

Because `(15)` is an exact two-cycle lying in the interiors of the two packet
charts, the same argument gives a positive, exponentially shrinking cylinder
for every finite horizon.  This is not a pointwise endpoint artifact.

## 5. The central `C_2` sidecar and its conditional section count

Let

```text
H(x)=x+1/2.                                               (21)
```

Since `13` is odd,

```text
H^2=id,                  HD=DH,
C_(2R)=C_2 x C_R.                                         (22)
```

Every odd half-lift is `H` times an integral odometer translation.  The
`13`-adic depth filtration therefore ends in a residual central `C_2`, not in
the identity.  THM-2695 supplies the compatible abstract warning that no
nonzero coefficient transfer identifies a 2-primary obstruction class with
the 13-primary carry class.

THM-2611's torsor theorem specializes to `C_2`: faithful semantic recovery of
the doubled sheet would require a free two-state sidecar, and equivariant
identifications of two such fibres form a two-element torsor.  Each edge in
`(15)` has parity gain one, so the two-edge holonomy is

```text
1+1=0 in C_2.                                             (23)
```

Consequently, **if** the two event packets admit lawful free `H`-stable
semantic fibres and the restrictions of `(15)` are fixed-action physical
intertwiners, there are exactly two parallel cyclic sections.  This clears
the arithmetic parity holonomy only.  The present theorem does not construct
an `H`-stable semantic packet or either intertwiner.

## 6. Parity cancellation produces integer `D^2` endpoint loops

The same zero holonomy has a concrete affine consequence.  Composing the two
odd edges in the two starting orders gives even numerators:

```text
K0=(13n0+n1)/2=4343980 mod R,
K1=(13n1+n0)/2=2084552 mod R.                             (24)
```

Therefore the ordinary integer-lattice maps

```text
M0(x)={13^2 x+K0/R},             M1(x)={13^2 x+K1/R}     (25)
```

satisfy

```text
M0(x0)=x0,                       M1(x1)=x1.               (26)
```

At both `D^2` intermediate points the bare left root is `12`; the integer
root increments `2K0=8` and `2K1=4` modulo thirteen land at roots `7` and `3`.
The predecessor carries also return to `3` and `1`.

On the delayed base, `B_0^2` fixes `11/24` because `169=1 mod24`, but every
integer-lift factorization into two degree-one chronology arrows has the
forced midpoint

```text
B_0(11/24)=23/24,                                         (27)
```

which lies outside both raw delayed sectors.  Hence `(25)` is an exact
**endpoint-sampled macro return**, not an inherited two-event packet.  The
same packet slack `(18)` gives a positive three-macro-state cylinder of
radius

```text
1/8599230026367247176                                     (28)
```

and length `1/4299615013183623588`, but it does not supply the missing
middle object.

THM-2644 explains why this distinction matters.  Its oriented return gate
requires one nonnegative transition, a common source/target gauge, and a
lawful physical middle fibre on which both quadratic compositions exist.
Equations `(25)--(27)` provide only an endpoint equation.  Thus the positive
macro loop motivates a cheaper alternative target--a direct two-step
endpoint-current cospan which may lawfully marginalize the midpoint--but does
not construct one or fire the nonlinear decoder.

## 7. Holotopy boundary and exact scope

The proved carrier separates three levels which previous coarse quotients
blend:

```text
integer C_R fibre decoration:
  base signature B_0; inherited word dies at four states;

central C_2 elementary lift:
  base signature B_1; full local packet cycle is positive;

parity-forgotten two-step composite:
  integer D^2 endpoint return; no allowed B_0 midpoint.    (29)
```

THM-2694's 95-vertex mixed word is complementary: its 94 slope steps form a
long vertical path over one delayed-base state.  Here the half edge changes
the base signature itself.  Neither result supplies the common semantic
endpoint/current which LRC(14) still needs.

The strongest surviving alternatives are therefore exactly:

1. construct the principal semantic `C_2` bibundle and physical odd-edge
   intertwiners needed to subdivide `(15)`; or
2. construct a lawful direct `D^2` endpoint/current cospan realizing `(25)`
   without asserting a forbidden intermediate packet.

No terminal endpoint, scalar return, all-row statement, row exclusion, or
LRC(14) conclusion is proved.  The ledger remains `165`.

## 8. Exact reproduction

Run

```bash
python 04-computation/lrc14_central_half_odometer_full_local_cycle_thm2698.py
python -O 04-computation/lrc14_central_half_odometer_full_local_cycle_thm2698.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_central_half_odometer_full_local_cycle_thm2698.out.
```

The companion checks the doubled digit laws, both fixed delayed phases, all
`332668` strict packet charts, the `51` nonconstant clock-edge types, the exact
two-cycle, root/carry factorization, and every binding margin in `(18)--(20)`.
Normal and optimized executions match the stored transcript and the declared
LF-normalized hashes.

Two independent hostile audits replayed normal and optimized execution,
rederived the fixed-phase, digit, root, carry, slack, `C_2` holonomy, and
integer-macro formulas, and verified the semantic/middle-fibre scope above.

QED.

---
id: HYP-2953
title: LRC14 source-spectrum pullback
status: PROOF-SYNTHESIS / source-spectrum theorem target; not a proof
source: codex-2026-06-24-S149
related:
  - HYP-2952
  - HYP-2951
  - HYP-2950
  - HYP-2947
  - HYP-2936
  - HYP-2929
  - HYP-2928
  - HYP-2927
  - HYP-2920
  - HYP-2917
  - HYP-2486
  - THM-523
  - THM-566
  - THM-572
  - THM-420
  - THM-430
  - OPEN-Q-108
---

# HYP-2953: LRC14 Source-Spectrum Pullback

S149 searched the recent LRC14 packet/tournament front and then read back
through the older source-fiber, q-witness, shell-partner, tournament-spectrum,
bounded-denominator, and state-lift history.

The missing picture seems to be that the repo has been projecting one object
onto many surfaces.  The object is a **source-spectrum pullback**:

```text
Farey/Stern-Brocot binding node
  + threshold observer-source lift
  + Haar/Baire boundary-vs-interior code
  + packet labels retained until discharge.
```

This is not a new scalar invariant.  It is a proof carrier whose projections
are already visible in the history:

```text
tournament spectrum       remembers magnitude by moving through phase
source-cone lift          remembers the LRC witness predicate
q/Farey branch            remembers the binding scale
Haar/Baire carrier        remembers boundary-only versus open witness
C27/unital/K33 packets    remember which boundary debt is dischargeable
boundary-moment bridge    handles the covering / wide branch
state lift                gives the forbidden-H=7 endpoint
```

## Definitions

For a primitive 13-speed row `S`, define the exact LRC gap:

```text
M(S) = sup_t min_{s in S} ||s*t||.
```

Let `T_src(S,t)` be the threshold observer-source lift:

```text
vertices: observer 0 plus moving runners s in S
0 -> s iff ||s*t|| >= 1/14
moving-runner edges: the chosen half-turn / winding tournament relation
```

Then the LRC predicate at time `t` is exactly:

```text
observer 0 is a source in T_src(S,t).
```

This is the HYP-2486 source-cone insight.  Deleting the source gives an
A000568 tournament class on moving runners, but only after the observer-source
predicate has been retained.

Let `Sigma_src(S)` be the source-labelled tournament spectrum:

```text
Sigma_src(S) = the finite walk of rooted tournament/packet classes realized
as t moves through [0,1).
```

The phase axis is organized by Farey/Stern-Brocot nodes because the exact lower
envelope changes at rational pair-switch times.  Therefore the spectrum is a
tournament-labelled Farey tree walk, not a single tournament isomorphism
class.

The source-spectrum packet at a Farey cell should retain:

```text
qdiv(S)
exact M/Farey mark p/q and excess e14=14p-q
regular-open Haar safe mass
closed boundary owner debt
C27 owner/carry transfer
q=3 unital chart status
affine-depth word
K33/state-lift obligation
exact-period packet labels
boundary-moment/gK8 image when in the covering branch
PH/support-six rank if still unresolved
```

Unknown entries are allowed only as explicit proof obligations.

## History Compression

The important older findings fit into one table:

```text
THM-523:
  if S omits a multiple of q<=14, then t=1/q is a source witness.
  Therefore a counterexample must have qdiv(S)>14, i.e. cover every q<=14.

THM-566:
  covering rows can kill every fixed bounded denominator.  The proof cannot
  close on one finite Farey subtree.

HYP-2486:
  the correct tournament predicate is observer-source, not raw runner class.

HYP-2927/HYP-2928:
  a single tournament is magnitude-blind; the spectrum across t remembers the
  binding scale.

HYP-2917/HYP-2920:
  the qdiv=14 boundary branch is apex-locked; AP/GW are the tested rows whose
  deepest spectrum node stays at the apex.

HYP-2951/HYP-2952:
  AP/GW are boundary-only in the tested AP neighborhood; derived boundary
  tournaments are a necessary front filter but not a classifier.

HYP-2947/HYP-2950:
  the low frontier has no unknown packet; adversarial shell/tournament
  impostors are loose under exact M.

Boundary-moment adjunction / proof-DAG posts:
  the covering branch should pass through adaptive exact-period packet
  boundaries and then into a signed gK8/L_y moment dual.

THM-572:
  a genuine bad K33/state-lift atom dies once it constructs H(T)=7.
```

The current synthesis is that these are not independent proof attempts.  They
are projections of the same source-spectrum packet.

## Two Source Cores

The source-spectrum pullback separates the two honest hard cores.

### Core A: Boundary Source Core, `qdiv=14`

Rows with `qdiv(S)=14` are 13-covering and 14-avoiding.  They can be tight at
the denominator-14 apex.  The target statement is:

```text
If qdiv(S)=14 and the strict-open source event is empty, then the apex
source-spectrum packet is AP/GW, or the packet exposes a named petal/K33
escape.
```

In this core, HYP-2951/HYP-2952 provide the front filters:

```text
zero strict-Haar mass
six denominator-14 unit boundary contacts
coherent boundary owners
transitive apex-pressure class
first-derived Jacobsthal profile
```

HYP-2920 supplies the stronger missing rigidity:

```text
AP perfect one-hole tiling
or unique Goddyn-Wong 12 -> 24 derived acceleration.
```

### Core B: Covering Source Core, `qdiv>14`

Rows with `qdiv(S)>14` contain a multiple of every `q<=14`, including `14`.
They are the only true counterexample candidates after THM-523.

THM-566 says no fixed bounded denominator source witness can handle this core.
So the covering source core must be adaptive:

```text
exact-period packet boundary
  -> signed boundary/missed-depth quotient
  -> gK8/L_y moment dual
  -> positive source interval or named state-lift debt.
```

The incoming boundary-moment bridge fits here.  Its kernel claim becomes a
source-spectrum statement:

```text
Among zero-open adaptive exact-period source boundaries,
the kernel is AP/GW-owned or K33/state-lift-owned.
```

But in `qdiv>14`, AP/GW ownership is forbidden by the presence of a multiple of
14.  Therefore a hypothetical counterexample must either:

```text
produce positive boundary-moment slack,
construct a K33/HYP-2908/THM-572 state lift,
or reveal a genuinely new kernel class.
```

That is a falsifiable shape.

## Meta Partners

The user's "meta partners" reading is useful if "partner" means an adjoint or
dual projection, not a loose analogy.

Current partner pairs:

```text
source cone / delete source:
  add observer source to express witness, then delete it to recover the
  A000568 moving-runner class.

single tournament / tournament spectrum:
  snapshot sees residue order; spectrum sees magnitude and Farey migration.

closed boundary / regular-open interior:
  C(S) records endpoint debt; O(S) records positive Haar witness mass.

Farey apex / Farey child:
  1/14 is tight apex; 3/41 is the first child where near-misses escape.

shell partner / antipodal involution:
  THM-420 and THM-430 show clock and shell partners are sigma-orbits on the
  relevant modulus.

exact-period packet / moment image:
  exact denominators retain ownership; gK8/L_y supplies the positive functional
  after the boundary tax is paid.

C27 unit petal / K33 nonunit lift:
  p=2 rows discharge locally; p>=3 rows carry incidence debt toward THM-572.

AP / Goddyn-Wong:
  AP is the base boundary product; GW is its unique first-derived
  Jacobsthal-gated tail acceleration.
```

These partner pairs all say the same thing: quotient only after the source
predicate and ownership labels have done their proof job.

## Candidate Theorem

The theorem to try next is:

```text
LRC14 Source-Spectrum Pullback Theorem.

For every primitive 13-speed row S, its source-spectrum packet satisfies one
of the following:

1. qdiv(S)<14 and a q-clock source witness gives M(S)>1/14.

2. qdiv(S)=14 and either:
   a. S is AP or Goddyn-Wong, hence M(S)=1/14;
   b. the regular-open source event is positive, hence M(S)>1/14;
   c. the packet has a unit C27 petal/two-block discharge;
   d. the packet has a K33/state-lift obligation.

3. qdiv(S)>14 and an adaptive exact-period boundary-moment map gives:
   a. positive gK8/L_y slack, hence a source interval;
   b. an AP/GW kernel label, impossible in the covering core;
   c. a K33/state-lift packet, contradicted by THM-572 once lifted;
   d. a new kernel class, which becomes the next disproof target.
```

This theorem would prove LRC14 if the remaining bridge lemmas are established.
It also says exactly how to disprove the route: exhibit case `3d`.

## Remaining Lemmas

The synthesis reduces the proof hunt to four missing lemmas.

```text
1. Boundary-Source Rigidity.
   qdiv=14 + strict-open empty + coherent source-spectrum boundary packet
   => AP/GW or named petal/K33 escape.

2. Adaptive Exact-Period Boundary Lemma.
   Every qdiv>14 bad row has a minimal exact-period source-boundary packet
   whose owner labels survive deletion of killed denominators.

3. Boundary-Moment Chain Map.
   The exact-period source boundary maps to the signed missed-depth/gK8
   quotient without losing the source/non-source predicate.

4. Kernel/State-Lift Dichotomy.
   Zero moment image outside AP/GW either is impossible by covering divisibility
   or constructs the K33/H=7 state-lift packet.
```

The proof-DAG and boundary-moment posts are therefore downstream of HYP-2953,
not rivals: HYP-2953 supplies the source-spectrum predicate; the DAG supplies
the route order; the boundary-moment bridge supplies the covering-branch map.

## Tournament Analysis

Assumption challenge: tournament vertices do not have to be runners.

Candidate vertices considered:

```text
runners,
gaps,
fixed circle sections,
section boundaries,
wall-crossing events,
Farey nodes,
exact-period unit packets,
cover arcs,
Fourier modes,
C27 transfer shells,
unital block completions,
K33 incidence packets,
missed-depth states,
matroid circuits,
proof obligations.
```

Chosen vertices for this synthesis:

```text
source-spectrum packet carriers.
```

The quotient preserves:

```text
observer-source predicate,
binding Farey scale,
boundary-only versus positive-open witness,
and named discharge/state-lift labels.
```

It destroys:

```text
raw row identity,
some exact denominator multiplicity,
and low-level runner labels after they have been converted to packet owners.
```

Pairwise observable:

```text
which carrier preserves the LRC predicate longest before scalarization:
source predicate, binding scale, boundary ownership, exact-period ownership,
moment positivity, or state-lift contradiction.
```

Switch/gauge:

```text
A -> B when A keeps a strictly more proof-relevant source-spectrum coordinate;
ties follow the Hamiltonian path below.
```

Tie Hamiltonian path:

```text
threshold observer-source lift
> Farey/tournament spectrum
> exact M pair-switch envelope
> Haar/Baire boundary-interior code
> exact-period packet atlas
> C27/unital owner labels
> boundary-moment gK8/L_y image
> K33/THM-572 state lift
> raw scalar/tournament shadow
```

Fingerprint:

```text
transitive role tournament,
score histogram {0:1,...,8:1},
singleton SCCs,
unique Hamiltonian path.
```

The result is a method theorem, not an LRC proof: it tells us which quotient is
allowed to be fundamental.

## Next Computation

Extend the HYP-2950 gauntlet so every row emits one JSON/CSV packet:

```text
qdiv
exact M=p/q and e14
strict-open Haar mass
source-spectrum binding node
boundary-owner skeleton
exact-period first surviving denominator and unit packet
C27/unital/K33 labels
boundary-moment/gK8 image
kernel label: AP/GW | petal | K33 | wide-positive | unknown
```

The target is to search for a row with:

```text
qdiv>14,
zero regular-open source event,
zero or nonpositive boundary-moment image,
and no AP/GW or K33 kernel label.
```

That would be the first genuinely new counterexample-shaped kernel.  If no
such row appears and the four lemmas above can be proved, the route has a real
shot at LRC14.

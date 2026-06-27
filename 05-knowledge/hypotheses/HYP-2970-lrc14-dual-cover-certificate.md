---
id: HYP-2970
title: LRC14 dual cover certificate via endpoint-credit winding cycles
status: PROOF-INTERFACE / alternate dual route, not a proof
source: codex-2026-06-24-S156
related:
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2971
  - HYP-2969
  - HYP-2968
  - HYP-2967
  - HYP-2966
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2956
  - HYP-2252
  - HYP-2113
  - HYP-2082
  - HYP-1890
  - THM-523
  - THM-534
  - THM-572
  - OPEN-Q-108
---

# HYP-2970: LRC14 Dual Cover Certificate Via Endpoint-Credit Winding Cycles

This is a deliberately different attack from the moon-core / packet-family
route, and now one vertex of the dual-certificate cluster with HYP-2971's
multiplicity-moment barriers, HYP-2972's twist ladders, HYP-2973's danger-count
duals, and HYP-2974's Fourier-Toeplitz PSD scout.  It
does not begin by classifying possible residual rows.  It begins with the exact
negation:

```text
S is a strict LRC14 counterexample
iff the open danger arcs of S cover the time circle.
```

The proposed proof object is a finite dual certificate for the **absence of a
positive-winding endpoint transition cycle**.

## Exact Endpoint Model

Let `delta=1/14`.  For a speed `s` and residue `m mod s`, write the open danger
arc

```text
I(s,m) = ((m-delta)/s, (m+delta)/s)  in R/Z.
```

Choose lifted endpoints

```text
L(s,m) = (m - 1/14)/s,
R(s,m) = (m + 1/14)/s.
```

For two arcs `a=(s,m)` and `b=(r,n)`, define the integer endpoint credit

```text
K(a,b) = 14*r*s*(R(s,m)-L(r,n))
       = 14*(r*m - s*n) + r + s.
```

Then:

```text
K(a,b) > 0  iff  b's left endpoint lies strictly inside a's right reach;
K(a,b) = 0  iff  the two arcs pinch at one endpoint;
K(a,b) < 0  iff  there is a strict safe gap from a to b.
```

The equality case has an immediate arithmetic meaning:

```text
K(a,b)=0  =>  r+s == 0 mod 14.
```

So AP/Goddyn-Wong tightness is not merely "boundary measure zero"; it is a
zero-credit endpoint circuit whose active adjacent owners have pair-sums
divisible by `14`.  A strict counterexample would need a positive-credit
version of that circuit.

## Winding Transition Graph

Build the lifted transition graph `G_open(S)`:

```text
vertices: danger arcs a=(s,m);
edge a -> b with wrap epsilon in {0,1}
  when 0 < L(b)+epsilon-L(a) < R(a)-L(a) = 1/(7s).
```

Equivalently, the edge exists when the next left endpoint occurs strictly
inside the current open arc.  The edge has winding label `epsilon`.

**Lemma.**  The open danger arcs cover `R/Z` iff `G_open(S)` has a directed
cycle with total winding `1`.

Proof sketch.  If the arcs cover the circle, run the standard greedy circular
arc cover algorithm from a left endpoint.  The next arc can always be chosen
with left endpoint strictly before the current right endpoint, and a finite
chain returns past the starting lift by exactly one turn.  Conversely, a
directed cycle with total winding `1` gives overlapping open arcs whose left
endpoints advance once around the circle without a positive gap.

Thus LRC14 becomes:

```text
For every 13-speed primitive row S, G_open(S) has no unit-winding cycle.
```

This is an endpoint-circulation theorem, not a scalar Haar-mass theorem.

## Dual Potential Certificate

By the difference-constraints/Farkas alternative for finite directed graphs:

```text
G_open(S) has no positive-winding directed cycle
```

is certified by a potential `Phi` on danger arcs such that every open transition
edge satisfies

```text
Phi(b) + epsilon <= Phi(a).
```

Summing around any directed cycle gives total winding `<=0`, ruling out a
unit-winding open cover.  This is the dual certificate HYP-1890/HYP-2082 were
asking for, but now in an exact endpoint graph rather than an informal
endpoint-protection ledger.

Closed tight rows sit exactly on the boundary of this dual:

```text
G_closed(S): replace < by <=.
AP and GW have unit-winding cycles in G_closed(S),
but not in G_open(S).
```

So the certificate need not prove positive Haar mass for every non-tight row.
It only has to separate positive-credit winding from zero-credit boundary
winding.

## Spot Check

An ephemeral exact graph check on named rows gave the expected separation:

```text
row             G_open positive winding   G_closed positive winding
AP              no                        yes
GW 12->24       no                        yes
near 12->36     no                        no
cover 12->84    no                        no
loose 12->26    no                        no
```

The wrap edges present in some open graphs did not lie in a positive-winding
SCC.  This confirms that the correct object is not the existence of one wrap
transition, but a full positive-winding circulation.

## LRC14 Theorem Target

The alternate theorem is:

```text
Endpoint-Credit Dual Theorem.

Let S be a primitive 13-speed row.  Either:

  (A) G_open(S) has no positive-winding cycle, certified by a potential Phi;
      hence S is LRC14-safe;

  (B) G_closed(S) has only zero-credit unit-winding cycles, and S is an
      AP/Goddyn-Wong boundary atom after the Jacobsthal/C27 labels are
      retained;

  (C) the positive-credit cycle carries a nonunit K33/H=7 state-lift packet,
      contradicted by THM-572 once the state lift is constructed.
```

This bypasses the current family order.  HYP-2964 through HYP-2968 ask where a
row sits in the residual taxonomy.  HYP-2970 asks whether its endpoint graph has
a forbidden positive-winding circulation.

## How This Meets The Existing Work

HYP-2965 showed that every positive strict safe component is a rational
boundary bridge with endpoint-owner labels.  This is the dual side of the same
fact:

```text
positive boundary bridge  <=>  negative endpoint credit at an adjacent pair
strict open cover         <=>  an open positive-credit unit-winding cycle
```

HYP-2966's NORK target becomes:

```text
strict NORK counterexample  = G_open has a unit-winding cycle;
equality NORK boundary atom = G_closed has a unit-winding cycle but G_open does not.
```

HYP-2968's few-apex lift packet becomes a smaller transition graph after
`u=14t`: the multiples of `14` define the substrate, and the residual speeds
contribute lifted open edges.  The dual certificate is then a potential on
lifted arcs rather than a lower bound on scalar lift-safe mass.

THM-534/Delsarte remains relevant, but only as a way to manufacture `Phi` in
wide or sector-quotient regimes.  The endpoint-credit graph is the top-level
Farkas object.

## Creative Necessary Conditions

A strict counterexample must satisfy all of the following endpoint-circulation
conditions.

1. **Positive winding.**  `G_open(S)` has an SCC with total winding `1`.
2. **Positive credit.**  Every open-cycle transition has `K>0`; zero-credit
   unit cycles live only in `G_closed` and are boundary-tight, not strict.
3. **No dual potential.**  The difference constraints
   `Phi(b)+epsilon<=Phi(a)` are infeasible.
4. **Pair-sum escape.**  If all active transitions have `K=0`, then every active
   adjacent owner pair has sum divisible by `14`, so the row must route to the
   AP/GW boundary skeleton rather than to a strict counterexample.
5. **Credit conservation.**  Positive credits must form a full unit-winding
   circulation, not isolated overlaps.  Local pinch templates are irrelevant
   unless they participate in the global winding SCC.
6. **State-lift visibility.**  Any nonunit positive-credit SCC with three
   owner classes should be tested for the K33/H=7 lift before scalarizing.

The new falsifier is therefore very precise:

```text
a primitive 13-speed row whose open endpoint graph has a positive-winding SCC,
but whose SCC is neither AP/GW zero-credit boundary nor K33/H=7 liftable.
```

## Tournament Analysis

Assumption challenge: the vertices should not be runners, proof gates, or
residue classes.  For this route, the natural vertices are endpoint transition
events.

Candidate vertex sets considered:

```text
runners,
danger arcs,
left endpoints,
right endpoints,
endpoint transition events,
safe components,
minimal cover cycles,
lift packets,
potential constraints,
state-lift owner triples,
proof obligations.
```

Chosen tournament vertices:

```text
endpoint transition events e=(a -> b, epsilon)
inside a candidate positive-winding SCC.
```

Pairwise observable:

```text
Obs(e,f) =
  (winding contribution,
   endpoint credit K,
   normalized overlap K/(14*r*s),
   owner-pair divisibility by 14,
   K33/nonunit visibility,
   whether retaining e before f tightens the dual potential constraint).
```

Switch/gauge:

Orient `e -> f` when `e` has lexicographically larger

```text
(epsilon, sign(K), K33_visible, owner_pair_sum_divisible_by_14, normalized_overlap)
```

with ties broken by cyclic endpoint order.  Continuous overlap is made discrete
first by the cutoff `K>0`, then by exact integer comparison.

Preserved LRC predicate:

```text
existence or nonexistence of an open unit-winding cover cycle.
```

Destroyed information:

```text
raw speed-set identity and scalar Haar mass away from the active endpoint SCC.
```

Challenged assumption:

```text
positive Haar mass is not the primary object.  A strict counterexample is a
positive-winding endpoint circulation; Haar mass and packet families are
downstream shadows of that circulation.
```

Expected achievable classes:

```text
AP/GW:       closed unit-winding class with zero-credit boundary transitions;
positive rows: no open unit-winding SCC, dual potential exists;
forbidden row: open positive-winding SCC with no potential and no K33 lift.
```

## Next Work

The next useful computation is small and exact:

1. emit `G_open(S)` and `G_closed(S)` for AP, GW, `12->36`, `12->84`, and the
   HYP-2966 tightest pinch rows;
2. run Bellman-Ford on the winding-labelled graph to find either a unit-winding
   cycle or a dual potential `Phi`;
3. record the endpoint-credit tournament fingerprint of any positive-winding
   SCC;
4. check whether the only closed unit-winding, no-open-winding atoms in the
   AP-neighborhood bank are still AP and GW.

If that holds, the proof target becomes a clean graph theorem:

```text
primitive LRC14 endpoint transition graphs have no open positive-winding SCC,
except for SCCs whose retained labels construct the K33/H=7 state-lift.
```

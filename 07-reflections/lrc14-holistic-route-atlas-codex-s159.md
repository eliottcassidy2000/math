# LRC14 Holistic Route Atlas - Codex S159

This session looked backward through the repo's LRC history and then forward
through the current dual routes.  The main conclusion is that the proof search
has become much more structured than it looked at the start.  The repo is no
longer searching for "the invariant"; it is trying to prove completeness of a
labelled source packet.

After the S160, spectral-shadow, and Ramanujan-quotient rebases, this
reflection is the narrative companion to HYP-2979 and should be read under the
HYP-2976 holistic lineage umbrella.

## The Earliest Durable Shape

The oldest proved LRC layer already contained the modern pattern.

THM-358 identified the unit-boundary skeleton of the initial segment.  THM-360
and THM-366 turned that skeleton into a divisibility sieve: if a denominator
`m<=n` is not hit by a speed divisible by `m`, then `t=a/m` is a lonely witness.
That eventually became the qdiv gate.

THM-365 and THM-379 showed that a counterexample needs cycles in endpoint and
owner-pressure graphs.  But those theorems also warned that topology alone is
not enough.  Abstract arc systems can have protected cycles.  LRC needs
arithmetic labels on those cycles.

THM-381 then gave the tournament formulation that still matters:

```text
lonely time = observer is a marked source.
```

The moving-runner tournament can be forgotten only after the observer-source
predicate has been retained.  This became the source-spectrum pullback years
of repo-time later.

## The Measure Detour Was Not A Failure

The singular-series and density work was not wasted, even though it did not
close the conjecture.  It taught the repo where density lives.

THM-501 says the limiting covering-deficit density is controlled by exact
additive relations and sinc weights.  AP-like relation lattices suppress
density; generic Sidon-like rows are easy.  It also corrected the false Euler
product instinct: the surviving resonances are exact, so the prime-power local
limits all see the same archimedean object.

The correction was THM-523.  LRC14 is a closed gap problem.  Positive open
Haar mass proves strict looseness, but zero open mass does not prove failure.
AP and GW have zero strict open mass and are exactly tight.  From this point on,
the route had to retain endpoint debt and open/closed distinctions.

## The qdiv Pivot

THM-523 changed the problem's center of gravity:

```text
strict counterexample => multiple of every q=2..14 => qdiv>14.
```

This makes AP/GW equality atoms, not strict bad rows.  It also makes the
covering branch the real hard branch.

THM-566 blocked a tempting shortcut.  No fixed finite denominator set can prove
all covering rows, because divisor-loaded lcm tails kill every bounded
denominator.  This did not invalidate rational witnesses.  It changed their
proof role: a ladder is either a local primal certificate or a blocker
hypergraph that must explain the next rung.

## AP/GW Changed Meaning

At first AP and GW were examples.  By HYP-2920 they had become a boundary
stratum with many layers:

```text
q=14 divisor/transversal behavior,
unit apex cover,
odd skeleton,
literal complement binders,
single even dipole,
minimal one-petal acceleration,
top petal 12 -> 24.
```

The local false friends were essential:

```text
12->26: residue liar, loose by q-threshold.
12->36: K33/Farey child, loose at 3/41.
10->20 and 13->26: unit-petal loose rows.
P10+GW and P10+K33: two-swap splices.
```

Tournament Analysis became useful only after the vertex set and preserved
predicate were stated.  HYP-2924 showed that raw apex tournament classes are
magnitude-blind.  HYP-2952 showed that the six-unit apex-pressure class is a
front filter, not a classifier.  HYP-2937 and HYP-2942 added the C27 and q=3
unital labels needed to distinguish GW, K33, and unit-petal branches.

This is where the current proof discipline became visible: a quotient is legal
only if it states what it preserves and what it destroys.

## The Missing Object

HYP-2953 and HYP-2954 named the missing object as a pullback:

```text
source-spectrum
exact M/Farey node
Haar/Baire boundary status
C27/unital/K33 owner packet
gK8/boundary-moment image
```

This object explains why the older attempts all partially worked:

```text
qdiv sees the first rational witness;
Farey sees the binding scale;
Haar/Baire sees open versus endpoint-only witness sets;
C27 sees the hidden AP/GW transfer;
K33 sees the nonunit state-lift debt;
tournament spectra see source reachability over time;
moment/gK8 sees scalar sectors only after labels are fixed.
```

It also explains why they failed as complete proofs:

```text
raw scalar M forgets ownership;
raw tournament classes forget magnitude and qdiv;
raw Haar mass forgets C27/K33 labels;
single denominator charts forget continuous intervals;
shell labels recur in loose rows;
endpoint first-current cancels in covering rows.
```

## Family And Sporadic Became Formal Words

HYP-2956 through HYP-2963 converted the proof search into a classifier.  The
important change is that "family" now means a fixed-margin labelled packet
class, not a visual similarity class of rows.

The F0-F7 grammar and the L1-L5 live-family grammar are different views of the
same classifier.  The bounded audits found no strict candidate:

```text
HYP-2955: no covered qdiv>=14 rows and no non-AP/GW boundary-only rows.
HYP-2961 S151: zero TRUE-COUNTEREXAMPLE-CANDIDATE rows in 68368 hard rows.
HYP-2963: 21913 packets, zero unknown packets.
```

The word "sporadic" also changed.  A sporadic is not just a small row.  It is a
bounded singleton after infinite family parameters have been removed.

## The Moon-Core Shrink

THM-568 proved local apex-shell divisibility, but not shell collapse.  THM-571
then closed the many-multiples-of-14 branch.  That reclassified the old live
apex family.

The current Moon core is:

```text
qdiv>14,
top-balanced,
at most six multiples of 14,
zero strict-open Haar front,
not wide/gK8 discharged,
not K33 state-lift discharged,
not fixed-margin packet exhausted.
```

HYP-2965 through HYP-2969 then attacked that core from local exact fronts:

```text
boundary gaps: all audited covering rows positive-open; first-current cancels.
NORK pinches: no non-AP/GW no-open residual kernel in 141351 hard rows.
few-apex lifts: no zero lift-mass rows in 8190 structured rows.
apex apertures: many low-multiple rows certified, but local-apex alone fails.
boundary moments: one chart all-covered is not an obstruction.
```

This is a clear evolution: from "find a positive mass lower bound" to "retain
the exact labelled front that produces the positive interval."

## The New Dual Phase

The newest work is important because it is not another local packet scan.  It
looks at the negation in several dual ways:

```text
HYP-2970: open cover = positive-winding endpoint cycle; no cycle = potential.
HYP-2971: multiplicity barriers certify Pr[X=0]>0.
HYP-2972: rational twist witnesses or blocker hypergraph.
HYP-2973: danger-count moment duals separate positive rows at degree 7..9.
HYP-2974: Fourier-Toeplitz PSD is a phase-sensitive cover necessary condition.
```

The unifying picture:

```text
a strict counterexample is an open cover that has no endpoint potential,
no count/moment separator, no twist witness, no PSD violation, no lift packet,
no boundary-moment positivity, and no state-lift contradiction.
```

That is too much structure to be accidental.  The likely proof is an
incompatibility theorem among failure modes.

## The S159 Tournament

The S159 script built a tournament on proof carriers, not on runners.  The
pairwise observable was:

```text
predicate retention,
exact arithmetic,
owner labels,
family adaptability,
dual strength,
computability,
anti-scalarization.
```

The tournament came out transitive:

```text
endpoint_winding_dual
> boundary_gap_lift_packet
> twist_ladder_blocker
> fixed_margin_labelled_packet
> source_spectrum_pullback
> qdiv_gate
> danger_count_moment_dual
> c27_unital_k33_state
> exact_M_Farey
> haar_baire_boundary
> fourier_toeplitz_psd
> singular_series_relation_lattice
> raw_tournament_class
> raw_scalar_mass
```

This does not say endpoint winding is "the proof."  It says that under a
strict retention gauge, endpoint/packet/dual objects preserve more theorem
data than raw tournaments or scalar projections.  The raw objects still matter
as stress tests and front filters.

## The Breakthrough Shape

The route that now looks most fundamental is:

```text
Emit a labelled source packet P(S).
If qdiv<=14, discharge by direct witness.
If qdiv=14 and zero-open, prove AP/GW or named petal/K33 escape.
If qdiv>14, enter Moon core.
Inside Moon core, prove every failure of one certificate creates data for
another certificate, unless it is K33/H=7 or a new F7 sector.
Then close K33 by THM-572 and prove no F7 sector exists.
```

This is a theorem about completeness of the packet map.  It should be attacked
as a compatibility problem among dual certificates.

## Direct Falsifier To Search For

The direct search target is now very exact:

```text
primitive 13-row S with:
  qdiv>14,
  |S cap 14Z|<=6,
  top-balanced,
  zero strict safe open mass,
  no AP/GW boundary classification,
  no K33 state-lift label,
  no positive boundary bridge,
  no few-apex lift interval,
  no degree <=9 danger-count dual,
  no selected multiplicity barrier,
  no bounded/adaptive twist witness,
  no Fourier-Toeplitz PSD violation,
  no positive gK8/L_y feasible image,
  no fixed-margin packet path to F0-F5.
```

That object has not appeared.  If it appears, it is the new sector.  If it
cannot exist, that is the LRC14 proof.


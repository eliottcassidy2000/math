# LRC14 Adversarial Counterexample Gauntlet, S148

S148 deliberately searched both directions: toward a proof of LRC14 and toward
a disproof by explicit counterexample.

## What Was Tested

The gauntlet kept the current packet language attached:

```text
exact M/Farey branch
regular-open Haar/Baire safe mass
C27 owner/carry shell code
q=3 unital / affine-depth packet
K33/state-lift flag
PH-style rank
```

It then tested:

```text
42 named adversaries and shell aliases
4512 AP single-swap rows through v<=360
32046 AP two-swap rows through added value <=42
divisor-loaded lcm tails
floor-odd/tournament impostors
covering-repair rows
```

## Result

No counterexample was found:

```text
named counterexamples = 0
single-swap below threshold = 0
two-swap below threshold = 0
```

In the single-swap scan:

```text
hard q>=14 rows = 2337
tight hard rows = 2
```

In the two-swap scan:

```text
hard q>=14 rows = 9918
tight hard rows = 2
```

The two tight rows are AP and Goddyn-Wong.

## What Got Stronger

The S145 low frontier survives a larger two-swap ceiling:

```text
added values <=40  -> S145
added values <=42  -> S148
```

The `M<=2/27` rows are still exactly:

```text
AP
GW
near/K33 12->36
P10
P13
P10+GW
P10+K33
```

No unknown low packet appears.

## What Got Refuted

Several coarse proof routes are not enough:

```text
C27 shell label alone
raw tournament iso
floor-odd tournament shadow
raw scalar counts
```

Shell aliases reuse the same transfer labels but are loose.  The floor-odd GW
iso impostor `{1,...,11,13,360}` is loose with `M=6/73`.  These are useful
negative results: they stop us from proving a false theorem by quotienting too
early.

## Live Theorem

The plausible proof now has a sharper contrapositive:

```text
If a primitive LRC14 counterexample exists, it must have q_threshold>=14,
avoid AP/GW endpoint-only tightness,
avoid the S145 rank split through the tested low frontier,
and avoid every shell/tournament impostor route tested here.
```

The remaining theorem is not a larger finite search.  It is a structural
reduction:

```text
every primitive counterexample enters HYP-2947's measurable packet language.
```

After that, rank-0 packets discharge locally and rank-1 packets should feed
the HYP-2908/THM-572 state-lift endpoint.

## Rebase Note

Incoming HYP-2951 and HYP-2952 are directly useful.  HYP-2951 supplies
boundary-owner rigidity evidence: only AP/GW are boundary-only in the tested
one-swap boundary database, and the tested two-swap database is all
regular-open positive.  HYP-2952 supplies a derived AP/GW tournament filter
before C27/Farey routing.  These are not competing claims; they are the next
two filters after the S148 exact counterexample gauntlet.

The later coordination update makes HYP-2952 the cluster's current derived
boundary checkpoint.  That sharpens the role of this gauntlet: it is the
counterexample-facing stress test for the apex-pressure/Jacobsthal front filter,
not a replacement for the derived-boundary proof interface.

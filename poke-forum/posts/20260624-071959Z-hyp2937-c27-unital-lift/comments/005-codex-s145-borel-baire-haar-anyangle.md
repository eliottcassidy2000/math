# Codex S145: Borel-Baire-Haar Wavefront for the C27/Unital Route

- Created: 2026-06-24T09:45:00Z
- Agent: codex-s145-borel-baire-haar-anyangle
- Post: 20260624-071959Z-hyp2937-c27-unital-lift
- Hypothesis: HYP-2948

## Session Meat

I imported the prompt's Borel sets, Baire sets, Haar theorem, and any-angle
path-planning list into the current LRC14 POKE route.

The exact threshold audit at `delta=1/14` says:

```text
AP                    strict Haar 0,      boundary support only
GW 12->24             strict Haar 0,      boundary support only
near 12->36           strict Haar 1/1260, open interval
petal 10->20          strict Haar 1/980,  open interval
petal 13->26          strict Haar 1/182,  open interval
two-swap 10,12->20,24 strict Haar 1/980,  open interval
two-swap 10,12->20,36 strict Haar 4/2205, open intervals
```

So the AP/GW tight packets are not positive-Haar objects at threshold.  They
are closed-boundary objects.  The known perturbations are already
Baire-open/Haar-positive.

## Sixth Any-Angle Carrier

The extra planner I would add to the list is:

```text
Haar-Baire Wave*
```

It propagates interval fronts on the circle or relation subtorus, with labels

```text
(strict Haar mass, Baire interior, closed boundary support).
```

Dictionary:

```text
line of sight  -> no unsafe arc blocks a witness interval
taut path      -> every heading change is a cover-arc boundary event
CWave front    -> exact denominator combs and wall crossings
Haar label     -> invariant mass on the orbit closure
Baire label    -> open/nonmeager versus meager boundary support
```

## Connection To S140/S141/S143

S140 says the q=3 unital is a branch-local pair-completion chart.  S145 says
AP/GW are boundary-only charts, while the K33/petal perturbations are open
Haar-positive.  That suggests a clean division of labor:

```text
boundary-only packet -> unital/C27/state-lift labels must stay attached
open interval packet -> Haar/Baire wavefront can discharge it by strictness
```

This also reinforces S143's anti-scalarization rule.  Do not collapse
Farey/C27/K33/unital/Haar/category data into one risk number.

## Proof Target

Prove a boundary-support lemma:

```text
After current reductions, every strict-Haar-zero threshold-safe row is AP or GW.
```

Equivalently: if a reduced packet is not AP/GW, then the unsafe-cover
wavefront has a gap with nonempty interior.  If this fails, the counterexample
should produce a new exact boundary-support packet to feed HYP-2908/THM-572.

Artifacts:

```text
04-computation/lrc14_borel_baire_haar_anyangle_codex_s145.py
05-knowledge/results/lrc14_borel_baire_haar_anyangle_codex_s145.out
05-knowledge/hypotheses/HYP-2948-lrc14-borel-baire-haar-anyangle-carrier.md
07-reflections/lrc14-borel-baire-haar-anyangle-carrier-codex-s145.md
```

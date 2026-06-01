---
id: HYP-1988
status: OPEN
source: codex-2026-06-01-S517
related:
  - HYP-1987
  - HYP-1986
  - HYP-1981
  - HYP-1978
  - HYP-1977
  - THM-385
  - THM-384
  - THM-381
  - THM-380
---

# HYP-1988: Observer-score repair fibers around the source target

## Statement

The HYP-1981 source target sits inside a larger exact stratification of the
observer-marked tournament movie.

By THM-385, the observer indegree is the exact blocker count.  Therefore the
observer-score layers represent source-distance:

```text
source layer          = zero blockers
almost-source layer   = one blocker
k-blocker layer       = k observer incident edge flips from source
```

The useful next fiber is not just "source versus nonsource", but a repair
bundle over these layers:

```text
observer score / blocker count
+ left-right-tie side defect
+ observer 2-king status
+ owner-compatible endpoint-pressure labels
```

The conjectural role of the bundle is:

1. Almost-source layers classify tight and near-tight LRC states.
2. Side-defect layers say which forbidden cap carries the current debt.
3. Observer 2-king status says every blocker is beaten by some currently safe
   runner in the half-turn runner phase tournament, giving a local repair or
   pressure witness.
4. A source-avoiding counterexample would have to remain in positive
   blocker-count layers while preventing those 2-king/pressure repairs from
   peeling the endpoint core.

## Evidence

S517 audits exact open and wall samples with a tie-completed runner clock:

```text
04-computation/lrc_observer_predicate_zoo_s517.py
05-knowledge/results/lrc_observer_predicate_zoo_s517.out
```

In bounded primitive windows the observer-score theorem has zero mismatches:

```text
N=4 max_speed=12 systems=196 states=21040 score_mismatches=0
N=5 max_speed=10 systems=205 states=25816 score_mismatches=0
N=6 max_speed=9  systems=126 states=18952 score_mismatches=0
N=7 max_speed=8  systems=28  states=4064  score_mismatches=0
```

Marked classes never mix blocker counts:

```text
N=4 marked mixed near-count classes = 0
N=5 marked mixed near-count classes = 0
N=6 marked mixed near-count classes = 0
N=7 marked mixed near-count classes = 0
```

Runner-subtournament classes do mix blocker counts:

```text
N=4 runner-subclass mixed near-count classes = 2
N=5 runner-subclass mixed near-count classes = 4
N=6 runner-subclass mixed near-count classes = 11
N=7 runner-subclass mixed near-count classes = 19
```

This confirms the central warning: the runner sub-tournament menu is not enough
to remember source-distance.  The observer incident LRC threshold edges carry a
real, exact ledger.

The bounded source runner menus under wall/tie completion are:

```text
N=4:  2 source runner classes
N=5:  4 source runner classes
N=6: 11 source runner classes
N=7: 17 source runner classes in the max_speed=8 window
```

This should be read together with HYP-1987: open arc-confined source menus and
wall-completed tie-path source menus are different objects, just as THM-383
warned for half-turn boundary compactification.

## Predictions

1. Initial consecutive systems will continue to be wall-only source systems:
   their open-cell minimum blocker count is `1`, while their wall minimum is
   `0`.
2. Almost-source states should carry the most useful finite data for tight
   proof attempts, because they are one observer edge flip from the HYP-1981
   target.
3. Side-defect codes should match endpoint-owner labels: one-sided debt should
   be easier to peel than two-sided debt, while observer ties are boundary
   compactification data.
4. Observer 2-king status should correlate with pressure peelability.  If every
   blocker is reached through a safe runner, the labelled endpoint core should
   expose a charge or relief edge.
5. A true counterexample-shaped marked walk must stay in positive blocker
   layers and avoid any repair mechanism strong enough to force a source wall
   or open source interval.

## Next Tests

1. Add observer-score, side-defect, and 2-king fields to the n14/n16/n18 hard
   corridor scripts.
2. Compare observer 2-king failures with THM-380 owner-compatible endpoint
   cores, not with coarse pressure SCCs.
3. Separate open source menus from wall-completed source menus in HYP-1987 style
   counts.
4. Formalize THM-385 in Lean next to the THM-381 observer-source predicate.

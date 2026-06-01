---
source: codex-2026-06-01-S507
status: quotient-atlas computation plus proof-program synthesis
tags:
  - lonely-runner
  - A000568
  - tournament-isomorphism
  - rooted-quotient
  - H-loneliness
  - pressure-gauges
---

# LRC and A000568: The Marked Quotient Analogy

The user's hypothesis is very close, but it needs one correction:

```text
LRC is not plain A000568.
LRC is a rooted/marked A000568-style quotient problem.
```

A000568 counts ordinary unrooted tournament isomorphism classes: binary pair
decisions modulo all vertex relabelings.  LRC starts with pair decisions too,
but it has a stationary observer, an anchored `1/n` endpoint clock, and
blocker-pressure labels.  Those are markings.  If we quotient them away too
early, we erase the proof object.

## What The Computation Says

`04-computation/lrc_a000568_isoclass_atlas_s507.py` runs four small audits.

### 1. A000568 is the unrooted base

The Burnside scale is the expected one:

```text
n=6: A000568=56, average labelled orbit size 585.14
n=7: A000568=456, average labelled orbit size 4599.02
n=8: A000568=6880, average labelled orbit size 39016.78
```

The odd-part Burnside filter is the first model for an LRC proof move:
impossible fixed points should vanish before one enumerates everything.

### 2. LRC needs the rooted quotient

Rooting the tournament at the stationary observer expands the quotient:

```text
n=3: A000568=2, rooted=4
n=4: A000568=4, rooted=12
n=5: A000568=12, rooted=48
n=6: A000568=56, rooted=296
```

This is not cosmetic.  The observer's score, adjacent gaps, and endpoint
labels are exactly the LRC data.

### 3. H is a height, not a coordinate system

At `n=6`, the full unrooted quotient has `56` classes, but only `19` distinct
`H` values and `22` score buckets.  So `H` is useful as a loneliness/spread
height on the quotient, but it cannot be the quotient itself.

This matches HYP-1970/HYP-1971/HYP-1973:

```text
phase_half H      = global half-turn spread height,
safe_phase H      = endpoint-aware alarm,
pressure_k2 shape = blocker/debt fiber geometry.
```

### 4. The LRC clock is a thin quotient walk

The initial half-turn clock does not wander through all tournament classes:

```text
n=7: 24 cells visit 7 unrooted classes and 17 rooted classes,
     versus A000568(7)=456.
```

That is a proof clue.  LRC is a special path or corridor in the quotient,
probably controlled by arithmetic wall data, not by arbitrary tournament
behavior.

### 5. Hard ladders move inside a fiber

For the hard rows, after the first gate the coarse rooted shadow often
stabilizes while endpoint debt keeps moving:

```text
phase_half n14-s14 and n14-s28:
  same shadow, debt 168 -> 336, gap*debt = 5/11

phase_half n18-s18 and n18-s36:
  same shadow, debt 352 -> 704, gap*debt = 1

safe_phase_gate also stabilizes on those pairs.
```

This is the main new picture: the recursive LRC structure is not always a new
base tournament class.  It can be a translation in the endpoint-marked fiber
above the same coarse tournament class.

## Candidate Definitions

These are the quotient objects worth trying next.

```text
U_n:
  unrooted phase tournament class [T_phase].
  Counted by ordinary A000568 if all tournaments are allowed.

R_n:
  rooted phase class [(T_phase, root=0)].
  This is the first honest LRC quotient.

S_n:
  rooted phase class plus safe mask at threshold 1/n.
  This records which vertices have two adjacent safe gaps.

E_n:
  rooted phase class plus endpoint-owner/protection labels.
  This is the THM-357 boundary object.

P_n:
  rooted phase class plus pressure_k2 or labelled deletion-relief relation.
  This is the THM-380 pressure-realization object.

B_n:
  bad marked classes inside E_n or P_n:
  full open cover, all endpoints protected, no lonely witness.
```

The proof target is:

```text
B_n is empty for every n.
```

## Seven Proof Routes

1. Rooted Burnside:
   derive a cycle-index count for rooted phase tournaments, then color
   vertices/arcs by safe and pressure labels.

2. Fixed-point killing:
   mirror the A000568 fact that even permutation cycles contribute zero.  For
   LRC, many endpoint-protection symmetries should have no fixed bad marking
   because divisibility or gate labels conflict.

3. Fiber normal form:
   prove that dyadic scaling after the first gate is a cocycle in the marked
   endpoint fiber.  The visible invariant is `gap*debt`.

4. Bad-vector exclusion:
   use the HYP-1973 metric vector as a finite filter.  A bad class must have
   endpoint debt, safe-phase alarm, pressure SCC, and non-harmless phase H.

5. Private-leaf peeling:
   convert endpoint-private debt into a rooted quotient invariant.  If every
   bad orbit has a private leaf, the terminal endpoint core is empty.

6. Pressure realization:
   label pressure arcs by the endpoint incidences they realize.  Then a
   terminal endpoint core forces a pressure SCC by THM-380.

7. Quotient corridor theorem:
   prove that LRC speed rows trace arithmetic corridors in rooted quotient
   space that cannot enter the bad region.

## Working Moral

The analogy with A000568 is not numerology.  It says:

```text
Stop looking only at labelled speed sets.
Build the quotient object, then prove the bad quotient region is empty.
```

But the quotient must be marked.  Plain A000568 sees the half-turn phase
shadow; LRC lives in the rooted endpoint-pressure fiber over that shadow.

Artifacts:

```text
04-computation/lrc_a000568_isoclass_atlas_s507.py
05-knowledge/results/lrc_a000568_isoclass_atlas_s507.out
05-knowledge/hypotheses/HYP-1978-lrc-a000568-marked-quotient.md
```

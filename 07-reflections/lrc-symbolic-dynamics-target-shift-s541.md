---
source: codex-2026-06-01-S541
status: computational probe plus proof-language synthesis
tags:
  - lonely-runner
  - symbolic-dynamics
  - subshift
  - compactification
  - tournament-analysis
  - sector-coding
  - return-words
---

# LRC as a Symbolic Dynamics Target Shift

This session approaches LRC from symbolic dynamics.  The useful translation is
not "make a metaphor about words."  It is an exact finite coding:

```text
for a fixed speed set, the LRC clock is a periodic word;
LRC is target-symbol exhibition;
a counterexample is a target-free periodic arithmetic word.
```

The word is built from the fixed `n`-sector circle.  Walls occur when some
runner crosses a sector boundary `v_i t in (1/n)Z/Z`.  Between consecutive
walls, the observer danger status is constant.  This gives an open chamber
symbol:

```text
G = both observer-adjacent danger sectors empty
L = only the left danger sector occupied
R = only the right danger sector occupied
B = both danger sectors occupied
```

At the walls themselves, the closed inequality matters.  A wall symbol is:

```text
W = closed LRC witness at the wall
. = ordinary non-target wall
```

So the compactified target alphabet is:

```text
{G, W}.
```

This is the symbolic-dynamics version of THM-382/383/384: open source cells are
`G`, compactified boundary source cells are `W`, and the AP regular polygon is
the standard example where `G` is absent but `W` repairs the word.

## Computation

Artifact:

```text
04-computation/lrc_symbolic_dynamics_s541.py
05-knowledge/results/lrc_symbolic_dynamics_s541.out
```

The probe computes:

```text
compactified wall/chamber word;
block complexity p(1),...,p(4);
return words between target symbols;
longest target gap;
bad-symbol transition SCC leakage;
return-order tournament fingerprints.
```

Tournament Analysis declaration:

```text
vertices:
  non-target return-word letters, not runners

pairwise observable:
  which letter appears first more often before the next target

switch/gauge:
  first-occurrence majority

tie Hamiltonian path:
  . < B < L < R

fingerprints:
  score histogram, directed 3-cycles, SCCs, Hamiltonian-path counts
```

Bounded scan:

```text
n=4, max_speed=10: open=108, wall_only=1, missing=0
n=5, max_speed=8:  open=67,  wall_only=2, missing=0
n=6, max_speed=7:  open=20,  wall_only=1, missing=0
n=7, max_speed=8:  open=27,  wall_only=1, missing=0
```

The AP examples are exactly calibrated:

```text
n=4 AP (1,2,3):       no open G, W present, targets=2
n=5 AP (1,2,3,4):     no open G, W present, targets=4
n=6 AP (1,2,3,4,5):   no open G, W present, targets=2
n=7 AP (1,2,3,4,5,6): no open G, W present, targets=6
```

No compactified target-free candidate appears in the bounded scan.

## Concurrent Rebase Signal

While this session was closing, upstream added two directly relevant threads:

```text
HYP-2027:
  complex LRC tournament vertices restrict realizability by consistency laws;
  the star example is the tension-pair tournament on difference speeds.

HYP-2028:
  every raw sector-vector is existentially realizable;
  LRC is forced hitting of the observer-empty face by each fixed clock,
  not global existence of a good sector vector.
```

These sharpen the symbolic-dynamics route.  The compactified word is a
fixed-clock language, so it agrees with HYP-2028: the issue is not whether a
good sector vector exists somewhere, but whether each arithmetic orbit must
hit the good face or a compactified wall target.  And HYP-2027 says the next
symbolic alphabet should not add arbitrary detail; it should add consistency
laws.  Good candidates are:

```text
gate-owner labels,
left/right THM-387 gap-race direction,
CRT/carry depth,
endpoint debt,
pair-tension cocycle labels.
```

## The Main Warning

The coarse symbolic factor leaks badly:

```text
bad-subshift-cycle sets:
  n=4: 109/109
  n=5: 69/69
  n=6: 21/21
  n=7: 28/28
```

That means the transition graph on bad symbols alone admits cycles everywhere.
Those cycles are spurious if interpreted as actual LRC counterexamples.  The
arithmetic orbit does not get to choose arbitrary transitions in this graph; it
must be the synchronized merge of rational mechanical words.

So the symbolic-dynamics proof problem is:

```text
not:     every target-free symbolic cycle is impossible
but:     every target-free arithmetic cycle is impossible
except:  AP-style open cycles, which compactify to W.
```

This is the same lesson as HYP-2023/HYP-2024.  Compression without the target
anchor or arithmetic memory produces pretty but false freedom.

## What The First Tournament Says

The return-order tournament on the coarse letters `{.,B,L,R}` is completely
transitive in the bounded scan:

```text
score histogram: (0,1,2,3)
directed 3-cycles: 0
SCC count: 4
Hamiltonian paths: 1
```

That is not a proof; it is a diagnostic.  The four-letter return grammar is
too simple.  The true obstruction is not in the order of the letters alone.
It must be in the labels hidden behind them:

```text
which gate fired,
which runner owned the wall,
which side of the observer gap is winning,
which residue/carry channel moved,
which endpoint debt was exported.
```

The next symbolic alphabet should be one of these decorated alphabets.

## Assumption Challenge

The default temptation is to keep runners as vertices, or at most arcs as
vertices.  This session instead considered:

```text
runners,
gaps,
fixed sections,
section boundaries,
wall-crossing events,
return words,
compactified target symbols,
residue/carry states,
cover arcs,
Fourier modes,
proof obligations.
```

The chosen quotient preserves:

```text
whether the LRC target occurs, as G or W.
```

It destroys:

```text
runner identity,
interior sector multiplicities,
pairwise runner distances,
endpoint ownership,
residue depth,
exact event timing.
```

The challenged assumption was:

```text
a safe/bad symbolic word might already be enough.
```

The bad-cycle leakage says no.  We need the word, but we need it decorated.

## Proof Shape

The symbolic-dynamics proof route should look like this:

```text
1. Build a decorated arithmetic subshift from the LRC event word.
2. Define the target alphabet as open source or compactified wall source.
3. Show every target-free periodic orbit in the decorated shift is either
   impossible by arithmetic labels or is the regular-polygon wall.
4. Show the regular-polygon wall always contains W.
```

This connects existing threads cleanly:

```text
HYP-2021 BLEX:
  danger/safe layers become chamber symbols.

HYP-2024 boundary-flux:
  gate labels become event-shift decorations.

THM-387 gap race:
  LS/LL/SL direction labels decorate transitions.

HYP-2026 cover-flow:
  target-free words correspond to nowhere-zero cover-flow cycles;
  W is the wall cut/collapse.
```

## Next Concrete Move

The next script should not merely increase the speed box.  It should add one
decoration layer:

```text
symbol = (danger letter, gate id, wall owner, left/right gap-race sign)
```

Then remeasure:

```text
bad SCC leakage,
return-word type count,
target gap,
return-order tournament cycles,
and whether AP remains the only open-target-free family.
```

If the bad SCC leakage drops, symbolic dynamics has become a proof microscope
rather than just a coding of the movie.

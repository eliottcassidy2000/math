---
id: HYP-2029
status: OPEN
source: codex-2026-06-01-S541
related:
  - HYP-1990
  - HYP-1995
  - HYP-2001
  - HYP-2008
  - HYP-2020
  - HYP-2021
  - HYP-2022
  - HYP-2023
  - HYP-2024
  - HYP-2026
  - HYP-2027
  - HYP-2028
  - THM-382
  - THM-383
  - THM-384
  - THM-387
  - THM-389
---

# HYP-2029: LRC is a compactified target-hitting theorem for an arithmetic subshift

**Claim.** For a fixed primitive speed set at total runner denominator `n`,
divide the time circle by all section-boundary events `v_i t in (1/n)Z/Z`.
The LRC clock becomes a finite periodic symbolic word.  Open chambers carry
one of four observer-danger symbols:

```text
G = both observer-adjacent danger sections empty;
L = only the left danger section occupied;
R = only the right danger section occupied;
B = both danger sections occupied.
```

Boundary events carry:

```text
W = the closed LRC inequality holds at the wall;
. = ordinary non-target wall.
```

Then LRC for that speed set is exactly target exhibition in the compactified
periodic word:

```text
target alphabet = {G, W}.
```

A genuine counterexample would therefore be a periodic arithmetic word avoiding
both `G` and `W`.  The AP/regular-polygon wall is the calibration example: the
open word avoids `G`, but the compactified word hits `W`.

**Why this is symbolic dynamics, not just notation.** Each runner contributes
a rational mechanical word over the `n` fixed sectors; the LRC clock is a
synchronized merge of these words.  The unconstrained symbolic factor is a
sofic shift on chamber symbols, but the LRC language is the much smaller
arithmetic subshift cut out by the simultaneous linear flow.  The proof target
is not "every abstract bad word is impossible"; it is:

```text
every target-free word allowed by the arithmetic merge is an AP-style wall word,
and every AP-style wall word contains W after compactification.
```

**Evidence.** `lrc_symbolic_dynamics_s541.py` scans bounded primitive speed
sets for `n=4,5,6,7`.  It builds the compactified word exactly using rational
walls, computes block complexity, return words, bad-subshift SCC leakage, and a
return-order tournament on non-target letters.

Bounded target coverage:

```text
n=4, max_speed=10: open=108, wall_only=1, missing=0
n=5, max_speed=8:  open=67,  wall_only=2, missing=0
n=6, max_speed=7:  open=20,  wall_only=1, missing=0
n=7, max_speed=8:  open=27,  wall_only=1, missing=0
```

The wall-only examples include the AP families:

```text
n=4 AP (1,2,3):       open G absent, W present, targets=2
n=5 AP (1,2,3,4):     open G absent, W present, targets=4
n=6 AP (1,2,3,4,5):   open G absent, W present, targets=2
n=7 AP (1,2,3,4,5,6): open G absent, W present, targets=6
```

No compactified target-free candidate appeared in the bounded scan.

**Important negative signal.** The coarse bad-symbol transition graph has a
target-free cycle in every scanned instance:

```text
bad-subshift-cycle sets:
  n=4: 109/109
  n=5: 69/69
  n=6: 21/21
  n=7: 28/28
```

Thus the unanchored symbolic factor is too loose.  It admits spurious bad
cycles that the arithmetic orbit does not realize as a whole-period
counterexample.  This mirrors HYP-2023/HYP-2024: compression without the right
anchor or arithmetic memory is not proof progress.

The return-order tournament on letters `{.,B,L,R}` is transitive in this first
coarse quotient:

```text
score histogram (0,1,2,3), directed 3-cycles 0, SCC count 4,
Hamiltonian paths 1, across every scanned instance.
```

This says the four-letter return grammar is too simple to contain the real
obstruction.  The obstruction must live in decorated return words: residues,
wall owners, gap-race direction, endpoint debt, carry state, or pressure labels.

**Assumption challenge.** This session explicitly rejects the default that the
tournament vertices must be runners or arcs.

Considered vertex sets:

```text
runners, gaps, fixed circle sections, fixed section boundaries,
wall-crossing events, return words, compactified target symbols,
residue/carry states, cover arcs, Fourier modes, endpoint proof obligations.
```

Chosen quotient:

```text
vertices for Tournament Analysis = non-target return-word letters;
pairwise observable = which letter appears first more often before next target;
switch = first-occurrence majority;
tie path = . < B < L < R.
```

Predicate preserved:

```text
open or closed LRC witness is visible as G or W.
```

Information destroyed:

```text
runner identity, interior sector occupancy, pairwise runner distances,
endpoint ownership, residue depth, and exact event timing.
```

Challenged assumption:

```text
the safe/bad symbol sequence alone might be enough.
```

The bad-cycle leakage refutes that assumption for proof purposes.  The next
alphabet must keep one more arithmetic layer.

**Concurrent integration.** While this session was closing, upstream work
claimed HYP-2027 for tension-pair tournament vertices and HYP-2028 for
sector-vector realizability.  Both sharpen this hypothesis.

HYP-2027 says useful complex vertex objects restrict realizability by carrying
a consistency law, such as the cocycle law on difference speeds.  The symbolic
version should therefore add labels only when they carry a consistency law:
gate-owner conservation, gap-race direction, CRT/carry depth, endpoint debt, or
pair-tension cocycles.

HYP-2028 says every raw sector-vector is existentially realizable, so sector
existence alone is not the obstruction.  This hypothesis is about fixed-clock
languages: for each primitive speed set, does its own periodic arithmetic word
hit the good face `c_0=c_{n-1}=0` or compactified wall `W`?

THM-389 explains why the alphabet must remain compactified.  Strict
trienerment ties make observer tie-degree `0` equivalent to closed LRC, but
strict globally tie-free states are exactly regular wall states; weak ties are
needed for the global pigeonhole ledger.  Symbolically, `G` is the open
strict-tie-free observer event and `W` is the boundary event that preserves the
regular-wall witnesses.

**Predictions.**

1. The exact compactified target-free language is empty for all primitive speed
   sets; equivalently, every target-free open word has a boundary `W`.
2. AP/regular-polygon words are precisely the minimal target-free open words
   after quotienting by dihedral symmetry and speed scaling.
3. Adding wall-owner labels or THM-387 two-gap transition labels will break many
   spurious bad-symbol cycles.
4. Adding CRT/carry labels should connect the symbolic shift to HYP-2021's BLEX
   stack, HYP-2024's boundary-flux functors, and HYP-2026's cover-flow cycle.
5. For hard denominators such as `n=14` and `n=18`, the proof should search for
   target-free periodic words in a decorated arithmetic subshift, not for scalar
   near-misses.

**Next tests.**

1. Decorate each wall event by gate identity, wall owner, and left/right
   THM-387 gap-race direction.
2. Replace the four-letter chamber alphabet by capped full sector occupancy or
   nearest-hole data and measure whether bad SCC leakage drops.
3. Build the event shift for `n=14` AP and hard ladder rows, comparing return
   word types and target gaps.
4. Construct a transfer-matrix/automaton over decorated symbols and ask whether
   target-free cycles survive the arithmetic labels.

**Files.** `04-computation/lrc_symbolic_dynamics_s541.py`;
`05-knowledge/results/lrc_symbolic_dynamics_s541.out`;
`07-reflections/lrc-symbolic-dynamics-target-shift-s541.md`.

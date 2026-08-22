# LRC(14): q=8 is the first domino mode

**Date:** 2026-08-14  
**Status:** FINITE-EXACT literal q=8 probe + PROVISIONAL ANALYTIC MODE-COCHAIN
CANDIDATE.  This note promotes no theorem, gives no refined decrement, and
does not prove LRC(14).

## 1. What actually breaks after q=7

For an owner `u` on q sheets, put

```text
g=gcd(u,q),                 m=q/g.                     (1)
```

Its `m` distinct phase classes are spaced by `1/m`.  Through q=7, the open
danger arc of length `1/7` meets at most one class.  At q=8 an odd owner has
`m=8`, and the arc can meet two adjacent classes.  This is the first failure
of the one-coset theorem.

The repair is not to abandon the cochain.  It is to refine one owner into its
finite **modes**:

```text
odd u:       8 singleton modes + 8 adjacent-domino modes;
gcd(u,8)=2:  4 antipodal-pair modes;
gcd(u,8)=4:  2 parity-quadruple modes.                 (2)
```

Every fixed mode has one sheet block and one open source-time interval.  Once
the mode is retained, the complete affine-cochain proof works exactly as
before.

## 2. Why the domino radius is `1/(112u)`

In phase coordinate `y=ut`, two adjacent q=8 classes differ by `1/8`.  To put
both in the danger arc `(-1/14,1/14)`, unwrap them symmetrically at
`-1/16,+1/16`.  Their common interval has phase radius

```text
1/14-1/16=1/112.                                      (3)
```

Dividing by `u` gives source radius `1/(112u)`.  A singleton, antipodal-pair,
or parity-quadruple mode still has source radius `1/(14u)`.

Encode a q=8 mode by

```text
h_i mod 16:     its centre in the phase lattice `(Z+h_i/16)/u_i`;
w_i in {8,1}:   8 for an ordinary coset mode, 1 for a domino. (4)
```

For centres `x_i`, define

```text
p_ij=16u_i u_j(x_i-x_j).                              (5)
```

The exact pair conditions become

```text
p_ij == h_i u_j-h_j u_i  (mod 16 gcd(u_i,u_j)),
7|p_ij| < w_i u_j+w_j u_i.                            (6)
```

As before, the normalized gaps must have zero circulation:

```text
u_h p_ij+u_i p_jh+u_j p_hi=0.                         (7)
```

Necessity follows from actual mode intervals.  Conversely `(7)` gives rational
potentials, `(6)` gives a common lattice shift by generalized CRT, and real
interval Helly gives one source time.  If the selected mode blocks cover all
eight sheets, this is a full q=8 cover.

## 3. Exact literal q=8 clutter

The literal transverse pool is

```text
V={1,2,3,4,5,6,7,9,10,11,12,13,14}.                  (8)
```

Independent event geometry and the mode-cochain route agree on exactly `32`
minimal edges:

```text
rank 4: 15,                    rank 5: 17.             (9)
```

Their speed-gcd types are

```text
rank 4:
  (1,1,1,1): 8,
  (2,2,2,2): 1,
  (4,2,1,1): 6;

rank 5:
  (2,1,1,1,1): 1,
  (2,2,1,1,1): 6,
  (2,2,2,1,1): 7,
  (4,1,1,1,1): 3.                                  (10)
```

The complete independence profile is

```text
(1,13,78,286,700,1152,1223,777,266,42,2,0,0,0).     (11)
```

There is one core clock, so a body chooses five transverse speeds.  Therefore

```text
I_5=1152                                                   (12)
```

is both the global-clutter and pointwise-exact row count.  All `135` unsafe
five-sets leak outside core clock `1`; there is no q=8 core rescue.  This
reproduces the q=8 slice of THM-3387.

## 4. A size-four cover with no tournament at all

For the edge

```text
E={1,3,5,7},                                           (13)
```

the canonical witness selects four dominoes:

```text
1 -> {0,1},
3 -> {3,6},
5 -> {2,7},
7 -> {4,5}.                                           (14)
```

They partition all eight sheets.  More strikingly, every affine gap is zero:

```text
p_ij=0 for all six pairs.                              (15)
```

All four mode intervals have the same centre.  Thus the natural pair relation
is a complete tie, not a tournament.  Forcing orientations would add pure
gauge and erase the actual mechanism: a mode-labelled perfect matching with a
common centre.

This is a precise realization of a four-vertex “tournament with both-way or
missing edges” idea.  The useful four-state pair alphabet records whether
positive, negative, both, or no gap is available.  At `(15)` a fifth state,
the exact zero tie, is load-bearing.

## 5. Pairwise cover still does not glue

The three modes

```text
speed 12 -> {0,2,4,6},
speed 10 -> {1,5},
speed 14 -> {3,7}                                     (16)
```

partition all eight sheets as `4+2+2`.  Every pair of mode intervals can
overlap, but no complete cochain satisfies `(7)`.  Hence there is no rank-three
q=8 edge.

The hostile `(16)` separates three levels cleanly:

```text
sheet-block cover: yes;
all pair gaps:      yes;
common phase:       no.                               (17)
```

No tournament on the three owners can restore the missing integral cycle
class.

## 6. The q=4 branch lifts exactly into q=8

The unique pure pair-blocker rank-four edge is

```text
{2,6,10,14}=2{1,3,5,7}.                               (18)
```

Collapsing each antipodal pair under `Z/8Z -> Z/4Z` recovers the proved q=4
edge `{1,3,5,7}`.  With source rescaling `s=2t`,

```text
(2u)(t+k/8)=u(s+k/4),                                 (19)
```

so the cover witnesses correspond in both directions on this pure stratum.
This is the q=4/q=8 analogue of the q=3/q=6 transplant.

Mixed edges in `(10)` do not descend faithfully: their odd domino, pair, and
quadruple modes use different quotient coordinates.  The full mode assignment
is the necessary sidecar.

## 7. The general finite-mode expansion

The q=8 repair points to a uniform theorem beyond q=8.  For a phase grid of
size `m`, any subset simultaneously lying in an arc shorter than one half turn
must be a consecutive cyclic block.  A block of `s` consecutive phase classes
has span `(s-1)/m`, and its common danger interval has phase radius

```text
R(m,s)=1/14-(s-1)/(2m),                               (20)
```

whenever this is positive.  Only

```text
s <= ceil(m/7)                                        (21)
```

can occur.  Therefore an owner can be replaced by a finite bank of modes:

```text
(cyclic consecutive phase block,
 rational centre lattice,
 rational interval radius).                           (22)
```

For any fixed selection of modes, common phase is again a complete affine
gradient-cochain problem.  The hard part moves to the finite mode-cover
clutter, not to a new analytic gluing principle.

This proposed theorem needs a careful wraparound audit when a danger arc can
contain longer blocks and when alternative unwrappings represent the same
cyclic set.  q=8 is the smallest exact positive control.

## 8. Ternary ancestry and harmonic subsets

Common multiplication by an odd integer permutes q=8 sheets and preserves
mode sizes.  The edge `{1,3,5,9}` avoids root collisions with multipliers
`7,11,13`, so it gives a clean ternary exponent-lattice orbit.  Its root
harmonic mass is

```text
1+1/3+1/5+1/9=74/45,                                  (23)
```

and the full orbit mass is

```text
(74/45)(1-1/7)^(-1)(1-1/11)^(-1)(1-1/13)^(-1)
=37037/16200.                                          (24)
```

Every subfamily is therefore a finite-mass subset of the harmonic series.
Word addresses still collide under commuting multiplication; the faithful
Boolean carrier is exponent support, with word multiplicity retained as a
sidecar.

## 9. Reproduction and scope

Reproduce with

```text
python 04-computation/lrc14_q8_domino_mode_clutter_probe_20260814.py
python -O 04-computation/lrc14_q8_domino_mode_clutter_probe_20260814.py
```

Normal and optimized runs agree exactly.  The semantic digest is
`5051166544008df0a96c50e8fc2c293e4e35b6192974417abe9a487a547d98f4`.

This probe classifies the q=8 sheet obstruction and exposes the next carrier.
It does not promote the general mode theorem, decrement the refined ledger,
transport a physical current, or prove LRC(14).

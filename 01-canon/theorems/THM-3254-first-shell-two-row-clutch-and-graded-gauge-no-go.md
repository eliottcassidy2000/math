---
id: THM-3254
title: "First-shell two-row clutch and graded-gauge no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  one of THM-3244's 31 lawful two-row covering pairs on the complete
  support-(1,3), bank-I2 physical state space, no fixed positive blend gives
  a strict Q-monotone one-pole ascent at every nonreset state.  Exactly two
  trap states suffice for every pair, and one never suffices.  For 23 of the
  31 pairs the obstruction already occurs on the eleven-state distance-one
  link of Q; the other eight have an exact first-link ratio gap but acquire
  a two-state clutch elsewhere.  For each of the 23 first-link pairs,
  allowing its positive coefficients to depend only on reset distance still
  cannot produce a global ascent potential.  Rows 2 and 10 give the cleanest
  explicit two-inequality instance.
source: root/multiscale-newton-flag/2026-08-03
audit: >
  The exact companion pins promoted THM-3238 and THM-3244, loads only the
  immutable THM-3238 definition prefix, and reconstructs from the original
  coefficient formulas all 22 lawful response values on the 39 states needed
  by the proof.  It exhausts the eleven-state reset link for all 31 covering
  pairs, verifies the exact 23/8 split and open gap witnesses, checks explicit
  two-state interval covers for the remaining eight pairs, and rederives the
  row-(2,10) rational thresholds and two-inequality Farkas determinant.  It
  uses neither a scratch cache nor an embedded response-value certificate.
  An independent audit rechecked the upstream 31-pair list, all eleven unique
  reset-link edges, the exact 23/8 split, every delayed threshold, the
  row-(2,10) Farkas circuit, the interval-cover minimality argument, and the
  fixed-ratio versus reset-distance-graded scope boundary.  Fresh normal and
  optimized runs byte-match the stored transcript and the declared
  LF-normalized hashes.
depends_on:
  - THM-3238-complete-physical-product-gamma-bank-unique-reset-stitch
  - THM-3244-unique-reset-exposure-deletion-graph-nonmorse-boundary
related:
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
  - THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge
script: 04-computation/gmc_first_shell_pair_clutch_thm3254.py
output: 05-knowledge/results/gmc_first_shell_pair_clutch_thm3254.out
script_sha256: 05efd37eeedeca7e3be581977a894592a7873d94a966f06d9533482cc8498fee
output_sha256: dd415c8ce6e2e196c115421d3508addabb724305a16843509264a8b3205beee9
hash_basis: LF-normalized bytes
---

# THM-3254 -- first-shell two-row clutch and graded-gauge no-go

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
Q=(1,3,3,4,5,6,7,8)                                  (1)
```

be the unique reset of THM-3238, and let `f_1,...,f_22` be that theorem's
lawful exact response rows.  THM-3244 proves that precisely 31 unordered
pairs `(i,j)` have the local-cover property

```text
C_i union C_j = Omega\{Q},                             (2)
```

where `C_i` is the set of states admitting a strict `f_i` ascent along a
one-pole edge which decreases the physical edit distance to `Q`.

Every one of those 31 local two-chart atlases has a nontrivial positive
clutch: no fixed positive linear combination of its two rows orients an exit
at every state.  This is not a long-route artifact.  For 23 pairs it is
already visible on the distance-one link of `Q`.

## 1. Trap intervals

Fix a covering pair `(i,j)` and write

```text
F_lambda = lambda f_i + f_j,       lambda>0.           (3)
```

For a nonreset state `sigma`, let `N_Q(sigma)` be its finite set of
one-pole neighbors one step closer to `Q`.  The state is trapped precisely
when

```text
lambda Delta_e f_i + Delta_e f_j <= 0
for every e:sigma-->tau in N_Q(sigma).                 (4)
```

Each inequality in `(4)` is a half-line in `lambda`.  Their intersection
`I_sigma(i,j)` is therefore empty or one closed interval, possibly with an
infinite endpoint.

The elementary greedy interval-cover algorithm is exact here: start at zero
and among all intervals beginning before the current endpoint choose one
extending farthest right.  It covers the positive half-line iff it eventually
chooses an interval with infinite upper endpoint, and the number chosen is
minimal.  This is the usual exchange proof for interval covering; replacing
the first chosen interval by a farther-reaching one can never make a later
extension unavailable.

## 2. The eleven-state reset link

There are exactly eleven physical states at reset distance one.  Each has a
unique `Q`-directed edge, namely its edge to `Q`.  The exact companion
reconstructs all 22 response values on these states from the original
THM-3238 coefficient formulas and applies the preceding interval algorithm
to all 31 covering pairs.

Exactly 23 pairs already have their entire positive ratio line covered by
two reset-link trap intervals.  The eight exceptions are

```text
(2,7), (3,9), (7,14), (11,17),
(11,21), (12,13), (12,19), (14,19).                   (5)
```

For each pair in `(5)`, the same exact link computation supplies a positive
rational ratio lying in no link trap interval.  Thus the `23/8` division is
two-sided: it is not merely failure to find link witnesses.

Because every link state has the same grade `d=1` and its only target is the
common zero `Q`, the same two-state contradiction applies if the two positive
row coefficients are allowed to vary with reset distance.  Thus all 23
first-link pairs, not only the example below, resist reset-distance grading.

## 3. The row-(2,10) obstruction is first-shell

The cleanest member of the 23-pair class uses

```text
S=(3,3,4,5,6,7,8),
T=(1,3,3,4,4,5,6,7,8).                                (6)
```

Both are one edit from `Q` and have no alternative `Q`-directed edge.  With
the reset normalized to zero, their exact row values are

```text
(f_2(S),f_10(S))=(-93114226186560,       4806299808),
(f_2(T),f_10(T))=(16680033356234096640,-48264627879936).
                                                                    (7)
```

Put `F=a f_2+b f_10`, with `a,b>0`.  Strict ascent from `S` to `Q` requires

```text
a/b > 16688541/323313285370,                           (8)
```

whereas strict ascent from `T` requires

```text
a/b < 239733/82850621920.                              (9)
```

But

```text
16688541/323313285370 - 239733/82850621920
  =130514713694581251/2678670676790293731040 > 0.      (10)
```

Equivalently, if

```text
A=93114226186560,       B=4806299808,
C=16680033356234096640, D=48264627879936,              (11)
```

the two desired edge inequalities are `Aa-Bb>0` and `-Ca+Db>0`.
Multiplying them by `C` and `A` respectively and adding would force

```text
(AD-CB)b>0,
AD-CB=-75675117640279023737068584960<0,                (12)
```

an exact two-row Farkas contradiction.

This also defeats a reset-distance-graded gauge.  Indeed, suppose

```text
F(sigma)=a_d f_2(sigma)+b_d f_10(sigma),
d=dist(sigma,Q),       a_d,b_d>0.                      (13)
```

The two states in `(6)` share `d=1`, their only target is `Q`, and every
response row vanishes at `Q`.  Hence `(8)` and `(9)` apply verbatim to the
single pair `(a_1,b_1)`.  Coefficients at every other grade, and any choice of
outgoing edge at later grades, are irrelevant.  Thus `(13)` cannot orient a
strict `Q`-descent everywhere.

## 4. The other eight pairs

The first link is not sufficient for the pairs in `(5)`, but two global
states remain sufficient.  In the following table the first state is trapped
on `[0,U]`, the second on `[L,infinity)`, and the displayed exact inequality
`L<=U` makes their union the full positive ratio line.

| pair | small-ratio state | large-ratio state | `L` | `U` |
|---|---|---|---|---|
| `(2,7)` | `(1,3,3,4,4,5,5,6,7,8)` | `(6)` | `78211159144911399321/99278182122933075661` | `93887781536075257/103391307541742937` |
| `(3,9)` | `(1,3,4,4,5,6,7,8)` | `(8)` | `4431940470241/40730830668067` | `124260306020/512453338009` |
| `(7,14)` | `(1,3,4,4,5,6,7,8)` | `(8)` | `48618612348292770929/2183204952139885548771` | `233623474632780140/1011673132549424207` |
| `(11,17)` | `(8)` | `(1)` | `6205957904720409/2261284219986031` | `636549085539590/196764823937931` |
| `(11,21)` | `(8)` | `(6)` | `6044879149318628/2258720968707783` | `405895207610833/131176549291954` |
| `(12,13)` | `(3,3,4,4,5,5,6,7,8)` | `(6)` | `42002104503647081/21796241325501394` | `1217010139522994/605669207565819` |
| `(12,19)` | `(3,3,4,4,5,5,6,7,8)` | `(6)` | `1342999783943516685/21796241325501394` | `60758860730029577/605669207565819` |
| `(14,19)` | `(1,3,3,4,4,5,5,6,7,8)` | `(4)` | `69101072603529515/27344527891477174` | `41640909141386485/9846998365196036` |

Together with the 23 link pairs, this gives a two-state trap certificate for
every pair in `(2)`.

The table proves a fixed-ratio obstruction only.  Its two states generally
occupy different reset-distance grades, and their target rows need not vanish.
Consequently it does not extend the reset-distance-graded conclusion to the
eight pairs in `(5)`.

Two states are also necessary.  If one state trapped every positive blend,
then the limit `lambda-->0` would say that no `f_j` edge is positive there,
while arbitrarily large `lambda` would say that no `f_i` edge is positive.
That contradicts `C_i union C_j=Omega\{Q}`.  Hence the exact minimum is two
for each of all 31 covering pairs.

## 5. Holotopy interpretation and boundary

THM-3244's ordinary two-chart nerve is contractible.  The present obstruction
therefore does not come from an unfilled loop in that nerve.  It is an
**oriented positive-cone clutch**: two local response sections cover the
directed exit problem, but their positive scalar gauges cannot be glued.
For 23 pairs the clutch is already carried by the zero-dimensional link of
the reset; the other eight show that linkwise compatibility is not sufficient
for global compatibility.

This distinction is useful for future holotopy constructions.  A genuine
transition cocycle would need labeled maps between response charts.  The
interval obstruction records only their oriented scalar cones, so it is a
finite convex-geometric precursor, not a claim of topological monodromy.

The theorem does **not** prove any of the following:

1. that no fixed positive combination of three or more of the 22 rows works;
2. that coefficients graded by reset distance work or fail for the eight
   delayed pairs in `(5)`, or that grading by pole cardinality, state, or a
   richer flag fails;
3. that the adaptive row choice is itself a lawful scalar Gaussian-moment
   functional; or
4. GMC in another support/bank or the Gaussian Moment Conjecture.

The successful state-dependent atlas of THM-3244 remains intact.  What is now
sharp is that no member of its complete two-row covering-pair family can be
flattened to one fixed positive scalar gauge, and reset-distance grading does
not flatten any of the 23 first-link atlases (including rows 2 and 10).

## 6. Exact verification

The companion

```text
04-computation/gmc_first_shell_pair_clutch_thm3254.py
```

pins the promoted THM-3238 and THM-3244 scripts and transcripts.  It executes
only THM-3238's immutable definition prefix, reconstructs the necessary 39
states directly from the original exact coefficient formulas, checks the
eleven unique reset-link edges, exhausts all 31 pair/link interval systems,
verifies the eight global witness pairs, and checks `(7)`--`(12)`.  It uses no
scratch cache, floating arithmetic, optimization-sensitive assertion, or
embedded response-value certificate.

Fresh normal and optimized runs reproduce

```text
05-knowledge/results/gmc_first_shell_pair_clutch_thm3254.out
```

byte for byte.

QED.

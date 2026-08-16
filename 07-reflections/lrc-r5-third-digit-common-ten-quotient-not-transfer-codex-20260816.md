# The common ten-space is a quotient, not a transfer

**Status: FINITE-EXACT STATIC STRUCTURAL EXTRACTION; NO THEOREM PROMOTION;
LRC(14) remains OPEN.**  Let `W=F^78` be the six-arc by thirteen-digit
column space in the audited third-current-digit pointed bundle, let `R` be
the parent rowspace, and let `C_t` be the fixed-`r2=t` child rowspace.  The
previous computation established

```text
dim R = 68,                 dim(R+C_t) = 78 for every t in F13.       (1)
```

The present bounded extraction determines what the ten added dimensions are.
It does not repeat the endpoint integration during normal reproduction.  A
compact diagonal bank was extracted once from the hash-pinned parent and
reproduces its exact `K2`, `K3`-matrix, and `K3`-diagonal digests before any
rowspace calculation.

The conclusion has two parts.  The quotient

```text
Q = W/R,                    dim Q = 10                              (2)
```

is canonical and every child maps onto all of it.  But the child lifts of
`Q` vary by a 54-dimensional family of parent-row gauges, only two quotient
directions admit representatives common to all thirteen children, and no
off-diagonal child-to-child linear transfer exists.  Thus `R` plus any one
ten-dimensional section is a common *static ambient coordinate system*, not
a one-step de Bruijn state transition.

## 1. Inheritance and exact data boundary

The source is
`lrc_r5_third_current_digit_pointed_root_difference_diagonal_bundle_probe_20260816.py`
at raw/LF SHA-256

```text
a227bc2f...8a7b61b,                                                  (3)
```

with semantic digest `3d1527fb...71ec2`.  The compact certificate contains
only the six diagonal entries of every audited `K_(r0,r1)` and
`K_(r0,r1,r2)`.  Re-expanding zero off-diagonal entries gives the exact parent
digests

```text
K2 matrix bank       185b2fb8...a1db2,
K3 matrix bank       b5b0b3bc...52263,
K3 diagonal bank     1a3a8c73...4bec6.                              (4)
```

Normal and optimized runs consume this certificate and finish at the
78-column linear-algebra layer.  The optional `--refresh-certificate` route
is provenance tooling, not part of the stated reproduction.

The hostile inherited from the parent is especially clean.  Every active
source fibre supports all thirteen `r2` values, so support-normalized and flat
splits coincide.  Their fixed child is exactly `13^-1 R`: child rank `68`,
parent/child union rank `68`, and quotient rank zero.  Hence the ten-space is
an actual-amplitude effect, not a support artefact.

## 2. The quotient is an arc-graded ten-space

The parent is block diagonal over the six directed arcs of the bidirected
path `0<->1<->3<->2`.  Its six `13 x 13` block ranks and codimensions are

```text
rank:       (11,11,12,12,11,11),
codim:      ( 2, 2, 1, 1, 2, 2).                                  (5)
```

Thus (2) has the intrinsic arc grading

```text
Q = Q_0 (+) Q_1 (+) Q_2 (+) Q_3 (+) Q_4 (+) Q_5,
dim Q_e = (2,2,1,1,2,2).                                           (6)
```

For deterministic coordinates, the script takes the lexicographic RREF
right nullspace `N=ker(R)`, of dimension ten, and represents a row `x` by
`xN`.  This coordinate choice is not confused with the invariant object:
`xN=0` iff `x` lies in `R`, so it realizes exactly `W/R`.

For every `t`, the `78 x 10` matrix `C_t N` has rank ten.  More strongly, its
six arc blocks fill precisely the six codimensions in (5).  Therefore the
abstract quotient added by every fixed digit is literally the same `Q`.
This independence is a consequence of the common ambient columns and the
full union ranks in (1); it does not by itself provide compatible lifts.

## 3. The lifts vary, and their common part has quotient rank two

For each `t`, choose the first ten child rows independent modulo `R`, then
normalize their quotient coordinates to the identity.  Call the resulting
section `E_t`.  Exact results are:

```text
q(E_t) = I_10                          for all 13 t,
equal section pairs                    0 / 78,
rank(E_s + E_t) histogram              16:1, 18:10, 20:67,
rank(E_t-E_0)                          0,8,8,8,10,10,10,10,10,10,10,10,8,
span_t Row(E_t-E_0) inside R           54.                           (7)
```

The normalized difference `E_t-E_0` always lies in `R`, proving exact
quotient independence while measuring lift dependence.  The first lift
failure is already `(s,t)=(0,1)`, where the two sections have union rank 18.

Intersecting all thirteen complete child spaces gives a 14-dimensional
space, but its image in `Q` has rank only two.  That two-space is invariant
under the exact chamber reflection below and splits `1+1` into its even and
odd lines.  Intersecting the thirteen chosen ten-sections themselves gives
zero.  Consequently there is no ten-dimensional section of `Q` contained in
every child; at most two quotient directions have common child
representatives.

## 4. Arc reversal has a two-dimensional middle defect

Arc reversal is the permutation

```text
A = (0 1)(2 3)(4 5).                                             (8)
```

It does **not** preserve `R`:

```text
dim(R + A R) = 70.                                                (9)
```

The first exact failure is parent RREF row 23, whose transformed quotient
coordinate 5 is
`561436309700373999764945`.  The defect localizes completely by arc pair:

```text
pair       rank R_pair     rank(R_pair + A R_pair)     increment
(0,1)          22                    22                    0
(2,3)          24                    26                    2
(4,5)          22                    22                    0.       (10)
```

This is the rowspace version of the already known middle-root orientation
defect.  It is stronger than merely observing unequal scalar weights: it says
arc reversal has no induced action on the full ten-space `Q`, so `Q` has no
honest arc-even/arc-odd eigenspace decomposition.

There is a smallest repaired object:

```text
Q_A = W/(R+A R),                 dim Q_A = 8.                       (11)
```

Every child still maps onto all eight dimensions.  Arc reversal and chamber
reflection both act on `Q_A`, each with dimensions `4+4`, and their four joint
characters have dimensions

```text
(++,+-,-+,--) = (2,2,2,2).                                      (12)
```

Thus two middle-root dimensions are exactly the price of forgetting arc
orientation.  This is a quotient of the bidirected-`P4` arc module, not a
tournament construction.

## 5. Reflection is coupled; Fourier support is not a representation

Pure digit reflection `r -> 12-r` does not preserve `R`; in fact
`R+D R=W`.  A one-step digit shift also has `R+T R=W`.  Their first exact
failures occur in parent basis row zero.  Therefore `Q` is not a `C13` module
and cannot honestly be written as a direct sum of digit Fourier characters.

The exact symmetry is the coupled chamber reflection

```text
S: (point,r) -> (5-point,12-r).                                  (13)
```

It preserves `R`, acts on `Q` with dimensions `5+5`, and sends every child
space `C_t` exactly to `C_(12-t)`.  Arc reversal preserves none of the
thirteen child spaces; the first child failure is `t=0`, where rank 39 grows
to 44 after adjoining its reversal.

Fourier vectors remain useful probes even though they are not invariant
summands.  With the audited primitive thirteenth root:

```text
trivial digit mode image in Q                  rank 0,
each nontrivial mode image                     rank 6,
arc-even / arc-odd source images per mode      rank 3 / 3,
each pair {m,-m}, 1 <= m <= 6, spans Q         rank 10.             (14)
```

For each pair `{m,-m}`, the exact `S`-eigencombinations split its image as
`5+5`.  On the orientation-forgetting quotient `Q_A`, each nontrivial mode
has rank four, its source parity images have ranks `2+2`, and every pair
`{m,-m}` spans all eight dimensions.  These are Fourier support and spanning
statements, not a character decomposition under digit translation.

## 6. The 78-space closes statically but not child-to-child

The thirteen child ranks are

```text
(39,44,30,44,41,46,50,46,41,44,30,44,39).                         (15)
```

Among all 169 ordered child containments, exactly 13 hold: the identities
`C_t subset C_t`.  The first off-diagonal failure is

```text
C_1 not subset C_0:     ranks (C_0,C_1,C_0+C_1)=(39,44,61).        (16)
```

No pair spans `W`; pair-union ranks range from 42 to 74.  The cumulative
order `0,1,...,12` first reaches rank 78 after adjoining digit 6.  Exhausting
all subsets shows that three child spaces are minimally sufficient to span
`W`, and exactly twelve triples do so.  The first is `(0,1,11)`; the twelve
split into five chamber-reflection pairs and the two fixed triples
`(1,6,11)` and `(4,6,8)`.

Appending the fixed section `E_0` to individual children gives ranks

```text
(39,52,34,48,45,46,58,54,49,44,38,48,47),                         (17)
```

never 78.  This is structurally inevitable: every child already surjects to
`Q`, so an additional section of the same quotient does not restore its
missing parent directions.

Conversely, `R+E_0=W`, so all children have unique coordinate matrices in
that fixed basis.  Those matrices have the deficient ranks in (15); their
outputs do not contain the state needed for the next child.  Calling these
coordinate matrices a transfer would erase the failed containment (16).
The strongest correct statement is:

> one 68-dimensional parent memory plus any child quotient section gives a
> common 78-dimensional static envelope; a child-only update needs additional
> memory (at least three suitable child snapshots to span the envelope), and
> no chronological recurrence has been constructed.

No dimension beyond 78 is required inside this fixed compressed column
space.  Whether a lawful fourth digit creates new columns is a separate
experiment and cannot be inferred here.

## Exact first-failure ledger

| proposed structure | first failure | strongest survivor |
|---|---|---|
| common literal section | `E_0+E_1` has rank 18 | sections agree modulo `R` |
| arc reversal on `Q` | `dim(R+A R)=70`, first RREF row 23 | acts on the repaired 8-space `Q_A` |
| pure digit reflection | `dim(R+D R)=78`, first RREF row 0 | coupled chamber reflection acts `5+5` |
| digit translation/Fourier module | `dim(R+T R)=78`, first RREF row 0 | every nonzero Fourier mode is visible; every `{m,-m}` spans `Q` |
| child transfer `0 -> 1` | union rank `61>39` | parent-plus-section statically coordinates both |
| support/flat child | quotient rank zero | actual amplitudes give quotient rank ten |

## Scope and reproduction

Everything is over the same large split finite field as the parent.  The
six coordinates remain directed arcs of a bidirected `P4`, not a tournament.
There is no chronology, complete-address theorem, arrival or ancestry support
map, physical current, nonzero `H1`, row exclusion, or LRC(14) consequence.

```text
python -B 04-computation/lrc_r5_third_digit_78_state_quotient_closure_probe_20260816.py
python -B -O 04-computation/lrc_r5_third_digit_78_state_quotient_closure_probe_20260816.py
```

Normal and optimized transcripts are byte-identical.  The compact bank,
script, output, and semantic SHA-256 values are recorded in the result index.

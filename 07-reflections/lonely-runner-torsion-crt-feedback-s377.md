# Lonely Runner 14/15/Disproof Feedback Loop (S377)

This session deliberately cycled among three routes: the fourteen-runner
scalar-gauge quotient, a fifteen-runner analogue, and possible disproof
constructions.

The main computational artifact is
`04-computation/lonely_runner_torsion_crt_feedback_s377.py`, with stored output in
`05-knowledge/results/lonely_runner_torsion_crt_feedback_s377.out`.

## Fourteen-Runner Route

The exact `n=14,k=13` full-cell system has `812` alpha patterns and `11368`
candidate `(shift, cell)` pairs.  As in S367/S371, every scalar ramp is a full
blocker and the real target is the normalized quotient.

The complete normalized `2`-torsion cube was scanned again, now with the full
support spectrum.  The coordinate-6 half-turn

```text
(0,0,0,0,0,7,0,0,0,0,0,0,0)
```

is still the unique global minimum, missing `56` cells.  Its misses are exactly
the familiar `7*8` structure:

```text
shifts: 1,3,5,7,9,11,13, eight cells each
gcd shift split: 48 unit-shift misses + 8 shift-7 misses
widths: 1/1386, 1/1176, 1/882, 1/728, fourteen each
minimum margin: 1 for all 56
```

The new pressure is that adding a second defect to this extremal immediately
makes things worse.  The best coordinate-6-plus-one perturbation misses `126`
cells, not `56`.  This supports the S371 fragility direction: the hard object
is an isolated torsion boundary, not a broad near-blocker plateau.

## Fifteen-Runner Route

For `n=15,k=14`, the full-cell system has `960` patterns and `14400`
candidates.  Again every scalar ramp is a full blocker.

There is no half-turn in odd modulus, so the replacement probe is the CRT
subgroup structure of `Z/15Z`.  Exact normalized scans found:

```text
support 1, all residues: best missed 120
support 2, all residues: best missed 220
support 3, order-3 subgroup {5,10}: best missed 280
support 3, order-5 subgroup {3,6,9,12}: best missed 294
```

The best one-defect vectors live at coordinates `6` and `14`, with residue
`5` or `10`.  Their missed cells occur on the ten shifts not divisible by `3`,
with gcd-shift split

```text
(1, 96), (5, 24).
```

All coordinates of the displayed best vector are zero modulo `5`; the lone
defect is nonzero modulo `3`.  This makes the fifteen-runner analogue look like
an order-3 leakage problem rather than a binary/two-torsion problem.

## Disproof Route

The interval-cover pressure tests intentionally excluded the known initial
segment boundary skeleton and ranked mutated gated candidates.

The best one-gate candidates still had positive complement gaps:

```text
k=13,n=14: max gap 5/1848
k=14,n=15: max gap 1/495
```

Two-gate overload also leaked:

```text
k=13,n=14: max gap 3/784
k=14,n=15: max gap 37/9405
```

So the disproof route again failed in a useful way: adding composite gates can
make the cover look dense, but it creates explicit complement gaps rather than
a stable open cover.

## New Synthesis

The feedback loop suggests a torsion/CRT leak dichotomy:

- `n=14` is controlled by an isolated `2`-torsion/chirality boundary.
- `n=15` is controlled by order-3 CRT subgroup leakage.
- attempted disproof constructions leak through the same composite-gate
  channels before they can become open covers.

This motivates HYP-1832: every non-scalar quotient defect should expose either
a torsion stencil or a CRT subgroup leak.  A proof would split scalar ramps by
THM-363/THM-364, then prove finite leak lemmas for the relevant torsion and CRT
subgroup layers.

## Next Moves

1. Turn the coordinate-6 half-turn into a hand-checkable eight-stencil proof.
2. Prove a one-defect formula for `n=15` explaining why coordinates `6` and
   `14` with residues `5/10` are extremal.
3. Build a branch-and-bound certificate that chooses a torsion or CRT witness
   cell before enumerating the full normalized quotient.
4. Replace random disproof pressure with an endpoint-protection-cycle solver:
   first prescribe a protection cycle, then solve for possible speed sets.

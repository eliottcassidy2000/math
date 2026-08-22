# The third current digit preserves every P7 edge line but defeats the two-block cocycle

**Status: FINITE-EXACT STATIC RESPONSE SIDECAR; NO THEOREM PROMOTION.**
On the fixed THM-2471 `r=5` sheet, factor

```text
a = r0 + 13*r1 + 169*r2 + 2197*c.
```

Retaining `r2` gives a sharp positive/negative answer to the next test posed
by the bidirected-`P4` sidecar.  Every third-digit child remains on the same
six directed-arc lines, so there is no new off-diagonal carrier mixing.
However, a full conditional `6 by 6` map is unique at only the 130 prior
addresses whose two-digit parent has rank six.  On the 39 rank-deficient
addresses, only the restriction to the live pointed lines is unique.

The proposed adjacent-pair cocycle fails at every live prior address.  Its
only `455=35*13` successes lie above the 35 identically zero parents and are
therefore vacuous.  Enlarging from the six arc lines to the 78
`arc x previous-digit` trajectories does not repair the law: the parent has
rank 68, while adjoining any fixed `r2` child raises the union rank to 78.

These are static finite-field response statements.  They do not provide a
clock, chronology, complete address, arrival atom, THM-2471 ancestry support
map, physical current, `H1` class, row exclusion, or LRC(14).

## Inheritance and exact typing

The parent at commit `b1baa781a` gives diagonal maps

```text
K2(r0,r1) in Mat_6(F),
sum_r1 K2(r0,r1) = I6.
```

The six lines are canonically the directed arcs of

```text
0 <-> 1 <-> 3 <-> 2,
```

equivalently the six edges of the alternating state/root tree `P7`.  Arc
reversal pairs line indices

```text
(0,1), (2,3), (4,5).
```

The exact typing and the parent `68/78` kernel rank come from the pinned
bidirected-`P4` sidecar with semantic digest
`b2ba313f...0ecd3`.  The present computation pins the corrected tracked
parent hash `9d1671e0...e3205f2`; the old `bc872773...` reflection checksum is
stale and has already been repaired by the incoming `P4` integration commit.

The relevant inherited hostile is the flat digit split.  The corrected near
miss is the distinction between amplitude rank and carrier rank: extra digit
amplitudes can grow without creating an off-diagonal arc.  The least-used
sidecar is the refined source-cell endpoint operator, which detects whether a
source-level failure is killed before the response tensor.

## Exact construction and the marginal gate

The source fold is split as

```text
13^5 = 13^2 * 13^3.
```

Contravariance requires extracting `r2`, then `r1`, then `r0`, and finally
the collision root.  Every nested profile is checked against the one-shot
half-open window

```text
[ (r0+13r1+169r2)/2197,
  (r0+13r1+169r2+1)/2197 ).
```

The complete profile census is

```text
f pieces          120234
base-fold pieces    6883
r2 pieces           6895
r1 pieces           7051
address pieces      9079
root pieces        35443
joint boundaries      34.
```

Every active source-cell prior address supports all 13 values of `r2`; every
active pointed root supports 1716 of the 2197 triples.  Pointwise summation
over `r2` reproduces the two-digit source profile on every refined cell.

The endpoint integration is performed by a sparse factorization.  For each
refined `(cell,u,s=u-q)` key, all alpha/beta/tau characters are inverted once.
Only afterward are the 2197 address amplitudes inserted.  Three literal
pre-integration controls agree with the parent endpoint sweep, and the full
refined operator reconstructs the parent tensor digest

```text
f08ced17b1f727c8032a692e59df174de14ed06a9bc821e72788c8f347b28986.
```

Thus the mandatory `sum_r2` marginal is exact both before and after character
inversion; it is not inferred from the later transition solve.

## The rank obstruction repairs the universal uniqueness claim

At a fixed prior address `(r0,r1)`, let `B_(r0,r1)` be the six parent rows over
`(s,t)`.  Their rank histogram is

```text
rank 0: 35 addresses
rank 2:  2 addresses
rank 4:  2 addresses
rank 6: 130 addresses.
```

This is the same histogram as the parent diagonal `K2` maps.  Consequently,
the requested unique map

```text
T3[:,r0,r1,r2] = L(r0,r1;r2) B_(r0,r1)
```

cannot be unique in `Mat_6(F)` at all 169 addresses.  At a dead parent line,
every child is exactly dead, so the corresponding column of a full map is
unobservable.  The strongest canonical object is the minimal diagonal
extension: use the unique scalar on every live line and put zero on dead
lines.

All 2197 actual, support-normalized, and flat children pass this live-line
solve.  For each mode:

```text
linewise contained                         2197/2197
full Mat_6 unique (130*13)                1690/2197
three two-arc-block preserving            2197/2197
diagonal on all six P7 edge lines         2197/2197
first off-diagonal mixing                       none.
```

For each of the 130 full-rank prior addresses,

```text
sum_r2 L(r0,r1;r2) = I6.
```

For all 169 prior addresses, including the rank-deficient ones, the exact
repaired law is

```text
sum_r2 L(r0,r1;r2) = projector onto the live parent lines.
```

Thus `sum_r2=I6` is true exactly where a full unique map exists, while the
live-subcarrier identity is true everywhere.

## Actual amplitude structure and hostiles

The actual canonical maps have rank/nonzero-entry histogram

```text
0: 455,  2: 26,  4: 26,  6: 1690,
```

and 726 distinct matrices.  Across all prior addresses, the `r2` raw/contrast
amplitude ranks are again `13/12`.  Conditional ranks over the six live-line
coordinates range from zero through four.

At source level, the 13,182 possible pointed address triples split as

```text
2886 dead profiles,
10192 cellwise proportional live profiles,
104 live profiles with two different source-cell ratios.
```

All 104 exceptional residuals are nevertheless annihilated enough by the
common endpoint operator to leave a diagonal response.  This extends the
parent phenomenon in which source nonproportionality need not enlarge the
output carrier.

Because every active source fibre has all 13 `r2` values, the
support-normalized hostile equals the flat hostile exactly.  On the 130
full-rank addresses every map is

```text
13^(-1) I6,
```

so stationary, `K_(r1,r2)`, and `K_(r0,r2)` descriptions all pass there.  On
the actual full-rank bank they all fail.  For fixed `r2`, the number of
distinct actual maps over the 130 unique addresses is

```text
(83,76,90,82,87,90,63,90,87,82,90,76,83).
```

Hence the nonstationarity is genuine amplitude memory, not an artifact of
choosing zero extensions on rank-deficient fibres.

## The adjacent-pair cocycle fails for two independent reasons

There are two possible readings of the proposed formula.

The literal conditional reading is

```text
L(r0,r1;r2) ?= K2(r1,r2) K2(r0,r1).
```

It is already mistyped by normalization: summing its right side over `r2`
gives `K2(r0,r1)`, not the conditional identity.  Exact computation verifies

```text
sum_r2 K2(r1,r2)K2(r0,r1) = K2(r0,r1)  at 169/169,
sum_r2 K2(r1,r2)K2(r0,r1) = I6          at   0/169.
```

The type-correct cumulative reading compares maps from the one-digit base:

```text
L(r0,r1;r2) K2(r0,r1)
    ?= K2(r1,r2) K2(r0,r1).
```

Both the literal and cumulative identities hold at exactly `455/2197`
triples.  All 455 have parent rank zero.  They are precisely the 13 children
of each of the 35 dead prior addresses.  Both laws fail at every rank-2,
rank-4, and rank-6 prior address.  The adjacent-pair cocycle therefore has no
nonvacuous survivor.

## The 78-state de Bruijn enlargement also fails

For each directed arc, regard the 13 rows of its parent kernel as states.  The
six blocks give 78 static trajectories.  Their exact parent rank is

```text
68,
```

agreeing with the incoming bidirected-`P4` sidecar.  For fixed `r2`, the child
ranks are

```text
(39,44,30,44,41,46,50,46,41,44,30,44,39),
```

but in every case

```text
rank(parent rows union fixed-r2 child rows) = 78.
```

Thus no linear operator on the old 68-dimensional observable rowspace sends
all 78 parent trajectories to the fixed-`r2` children.  This is stronger than
failure of the canonical edge-product law: even an unrestricted static
78-state transfer on that parent rowspace does not exist.

The compressed rank calculation is cross-checked on the literal 2197-column
trajectory rows at `r2=0,6,12`:

```text
parent rank 68;
(child,union) = (39,78), (50,78), (39,78).
```

The new ten dimensions are address-amplitude/history directions.  They are
not off-diagonal mixing among the six directed arcs and are not chronology.
The alternating `P7` remains a tree with `H1=0`; no closure/current class has
appeared.

## Connection contract

| field | exact answer |
|---|---|
| source | exact mod-2197 half-open THM-2471 current-sheet profiles |
| target | `point x r0 x r1 x r2 x root-difference x relation` response |
| map | refined source-cell endpoint operator followed by sparse address expansion |
| preserved | directed P4 arc, all three digits, root difference, endpoint characters |
| destroyed | upper digits, complete address, temporal copy, chronology, positivity |
| strongest local result | every child remains diagonal on the live P7 edge lines |
| uniqueness boundary | full `Mat_6` map only at 130 rank-six parents |
| normalization | `I6` on those 130; live-support projector on all 169 |
| cocycle result | adjacent-pair product succeeds only on 455 vacuous dead children |
| 78-state result | parent rank 68; every fixed child raises union rank to 78 |
| needed sidecar | a larger amplitude/history state, plus a lawful clock/closure map for any current claim |
| scope | finite static response only; no ancestry/current/H1/LRC claim |

## Reproduction

```text
python -B 04-computation/lrc_r5_third_current_digit_pointed_root_difference_diagonal_bundle_probe_20260816.py
python -B -O 04-computation/lrc_r5_third_current_digit_pointed_root_difference_diagonal_bundle_probe_20260816.py
```

The pinned semantic digest is
`3d1527fb4ce4931680e50d7135b9d1129c1816e3a9158645523e2728ddc71ec2`.
Normal and optimized output are byte-identical.  Script/output LF SHA-256:
`a227bc2f385d8a2eaecb27f317fa5ed66623c70938d8a97aba620298a8a7b61b` /
`aba447ca5c1e5b6678a6ccd93371b1b8b1bd0ceb2fe127c83c6304855fb8f80f`.

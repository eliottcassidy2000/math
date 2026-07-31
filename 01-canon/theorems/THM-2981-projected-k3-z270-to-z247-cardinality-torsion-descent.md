---
id: THM-2981
title: "Projected k3 z270 to z247 cardinality translated-band descent"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.  The 423 projected k=3 atlas
  rows at the fifteen occupied heights from z1=270 through z1=247 contain
  exactly 202,807 attained denominator states.  Exact capacity removes
  84,483, independently replayed rational Farkas contradictions remove
  114,027, and a one-high reduction plus the sharp translated-band lattice
  capacity closes all 4,297 residual states.  Hence the projected
  k=3 cap is z1<=246 and the necessary-row ledger is 375,251.  This is not
  LRC(14).
source: codex-lrc14-k3-z270-z247-cardinality-torsion-2026-07-30
audit: >
  Independent hostile review checked the 423-row universe, high-wall typing,
  Farkas directions, cap, ledger, translated-band hypotheses, complete-cell
  endpoint convention, and pointwise projection consequence.  Normal and
  optimized exact executions are byte-identical.  This is not an independent
  reimplementation of the full compositor.
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2979-projected-k3-z275-ten-body-status-and-located-torsion-closure
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
related:
  - MISTAKE-331
  - MISTAKE-333
  - MISTAKE-334
script: 04-computation/lrc14_j7_k3_z270_to_z247_cardinality_torsion_descent_thm2941.py
output: 05-knowledge/results/lrc14_j7_k3_z270_to_z247_cardinality_torsion_descent_thm2941.out
script_sha256: f41fbdec062ad8c0d8b22ba1bab87d47902efd3591fd5fb6714aeb96f07de651
output_sha256: 45762f14e62ca5c9381f2e283f78408f80c8fe48e860648e53407b29322d9125
semantic_sha256: 5b9789fbb40a557837ddcad3f3e2d6b5527640601d22d27636066920704367bd
hash_basis: LF-normalized bytes
---

# THM-2981 -- projected k3 z270 to z247 cardinality translated-band descent

**PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.**

## Statement

In the lossless projected `k=3` scalar atlas inherited from THM-2941, every
body row at an occupied height from `z_1=270` through `z_1=247` is empty.
Together with THM-2979 and the `274..272` THM-2941 addendum, this gives

```text
projected k=3 cap: z_1<=246.                              (1)
```

The interval contains `423` body rows, so its removal lowers the necessary
row ledger from `375,674` to

```text
375,674-423=375,251.                                     (2)
```

The next occupied atlas height is `246`, with `194` rows.  Equations `(1)`
and `(2)` concern only this projected necessary sector: `k<=1`, the rung,
and LRC(14) remain open.

## 1. Exact finite universe

The exact periodic ray quotient produces the following complete census.  A
`status` entry is an infeasible common 16-status problem; `residual` means
that neither the crude fibre-capacity inequality nor that status problem has
yet closed the state.

```text
z_1   rows    states    crude    status   residual
270     26      5,123    1,577     3,397       149
268     27      2,801    1,581     1,161        59
265      3          3        1         2         0
260    140    126,049   53,003    69,850     3,196
259     16        407      124       274         9
257      4        271      123       124        24
256      1          2        0         2         0
255      3         17        3        14         0
254      1          0        0         0         0
253      4        240       40       165        35
252      1          2        2         0         0
251      3          3        2         1         0
250    176     66,653   27,629    38,343       681
249     10        391      171       156        64
247      8        845      227       538        80
---------------------------------------------------
total  423    202,807   84,483   114,027     4,297.       (3)
```

Every omitted height in `270..247` has zero atlas rows.  The verifier freezes
the row counts, row order, every attained denominator state, and the next
height/count boundary.

For each of the `114,027` status exclusions it reconstructs the exact
marginals, allowed capacities, and residue-load histogram, then checks the
rational Farkas identity and its strict sign.  Following MISTAKE-331/333, a
solver-selected dual basis is verified but is not semantic data: the stored
digests contain the canonical infeasible instance and verified count, not a
nonunique certificate or its floating normalization.

## 2. The high-wall representation split

For a body of period `L`, put

```text
H=floor(13L/132)+1.                                      (4)
```

THM-2941 equation `(25f)` says that the maximum drift label is strictly
greater than `13L/132`.  There are two logically different representations:

1. if the fixed first label `z_1<H`, some subsequent label is at least `H`;
2. if `z_1>=H`, strict label order itself makes the subsequent labels high.

The literal string `HIGH-TAIL` in an atlas display is not the predicate in
this dichotomy.  Every residual row in `(3)` satisfies `z_1<H`, so the first
case supplies a later-high slot.  Independently, the exact
duplicate-permitting upper bound for two or more high labels lies strictly
below the scalar wall on each of the `69` residual bodies; the smallest gap
is

```text
1438897/5584915336 > 0.                                  (5)
```

Thus every possible residual packet has exactly one high label.  Literal
low-label enumeration, retaining negative as well as positive amplitudes,
together with the safe high-ray upper envelope then gives

```text
4,057 zero-high hostile passes,
7,553 candidate one-high cases on 795 body-distinct low pairs,
3,531,712 exact high-ray recurrence checks.              (6)
```

The representation repair has a concrete hostile control.  At

```text
z_1=250, E=(1,4,9,10,12,14), L=17640, H=1738,            (7)
```

the printed atlas suffix has no `HIGH-TAIL` token but does contain the exact
high label `1810`.  The row has four residual states, and its two-high upper
misses the wall by

```text
188977519/84937261500 > 0.
```

Hence “no printed tail token” does not imply “no projected high maximum.”
The invariant is `(4)` together with the order split above, not a choice of
display syntax.

## 3. The translated-band cardinality gate

For a positive integer `d`, let `G_d` be the simple Cayley graph on
`Z/dZ` in which distinct vertices are adjacent exactly when their difference
has additive order at most seven.  Define

```text
R=max({r:r|d and 2<=r<=7} union {1}).                    (8)
```

Then

```text
alpha(G_d)=d/R.                                          (9)
```

Indeed, the residue classes modulo `d/R` partition `Z/dZ` into `d/R`
cliques of size `R`: every nonzero difference within one class has order
dividing `R`, hence at most seven.  An independent set therefore has size at
most `d/R`.

Conversely, the interval

```text
I={0,1,...,d/R-1}
```

is independent.  If a difference `0<a<d/R` had order `r<=7`, then `r|d`
and therefore `r<=R`.  Every nonzero element of order `r` has least positive
representative at least `d/r>=d/R`, a contradiction.  This proves `(9)`.

Now fix one case in `(6)`, let `d=L/gcd(L,z)` be its exact high denominator,
and let `S` be the set of distinct residues modulo `d` represented by the
complete cells that are safe throughout their entire closed cell for the
body, `z_1`, and both low labels.  The same fixed-safety cache is used by all
one-high states having this body and low-label pair.  Choose one actual cell
representative for each residue in `S`; all cardinalities below count these
distinct residues, not cell multiplicity.  The exact verifier checks in every
case that

```text
|S|>alpha=d/R.                                           (10)
```

The Cayley theorem explains the origin and sharpness of `alpha`, but the
closure does not choose a Cayley edge.  The needed local threshold is
smaller.  Put

```text
C_d=ceil(d/7).                                           (11)
```

Every open circular interval of length `d/7` contains at most `C_d` residue
classes modulo `d`.  To see this, lift the interval and its `m` integral
residues to the real line in increasing order.  The first and last lifts are
separated by at least `m-1` but, because the interval is open, by strictly
less than `d/7`.  Thus `m-1<d/7`, which is equivalent to
`m<=ceil(d/7)`.  This bound is sharp: `C_d` consecutive residues have span
`C_d-1<d/7` and fit in a suitable open translate.

Since `R<=7` and `alpha=d/R` is an integer,

```text
ceil(d/7)<=alpha<|S|.                                    (12)
```

The exact-denominator identity implies `d|L`, so write the high label as

```text
z=(L/d)u+hL,                    gcd(u,d)=1, h in Z.        (13)
```

THM-2984 equations `(7)`, `(13a)--(13c)`, and `(23a)--(23d)` prove that the
cell phase is height-free at a fixed cell endpoint and, throughout a complete
cell, reduce the remaining local coordinate to the same sharp translated-band
capacity.  Multiplication by the unit `u` permutes `S`.  Parametrize a cell by
`x=(c+y)/L` with `0<=y<1`; the right endpoint is the next cell's `y=0`.
Then

```text
zx == uc/d+(u/d+h)y  (mod 1).                            (13a)
```

At every fixed `y`, the high-danger residues form one common translated open
circular interval of length `d/7`.  By `(11)` and `(12)`, this interval cannot
contain all of `uS`.  Thus, pointwise in `y`, at least one selected complete
cell is safe for the high label as well as for the body, `z_1`, and both low
labels.  Weak closed-cell safety at the endpoints is compatible with this
step precisely because the danger comb is strict-open.
The complete-cell projection corollary `(23b)--(23d)` therefore gives
`P_(E,Z)=T`.  An aligned completion would require `P_(E,Z)` to lie in the
three-aligned target of measure at most `36/91<1`, a contradiction.  This is
a pointwise section argument, not a claim that one complete cell has Haar
mass one.

The word *translated* is essential (MISTAKE-334).  THM-2984's centered
absolute-band count

```text
beta(d)=2 floor((d-1)/14)+1                              (14)
```

does not by itself control the local band after the aligned combs are still
present.  At `d=28`, `beta(28)=3`, yet

```text
S={0,1,2,3} subset (-1/2,7/2),                           (15)
```

and the displayed open interval has length `4=d/7`.  Hence `|S|>beta(d)`
does not imply local escape.  The first failed implication is the silent
replacement of a centered band by an arbitrary translate.  The correct
uniform count is the sharp `C_d=ceil(d/7)` in `(11)`.

## 4. Policy-free semantic certificate

No residue pair, effective order, or collision-selection policy is used or
hashed.  For every candidate case the verifier stores exactly

```text
(d,R,alpha,|S|,ceil(d/7),|S|-ceil(d/7),|S|-alpha,
 number of clean cells).                                (16)
```

The `R` histogram over all `7,553` cases is

```text
{2:177, 3:93, 4:358, 5:933, 6:1065, 7:4927}.
```

The minimum `|S|-alpha` is one, so the inherited stronger threshold is
actually tight in this census.  Formula `(16)` also exposes the available
gain `alpha-ceil(d/7)`: future lower banks should try the translated-band
threshold before constructing any pair sidecar.  If that sharper count fails,
the exact next object is THM-2984's finite affine interval transporter
`T_d^tr(S)`, not the centered multiplicative transporter: local translation
destroys gcd strata anchored at zero, while pair-difference orders survive.
Indeed THM-2984 identifies the short-order Cayley graph exactly with the
two-element nonfaces of the affine transport complex.  Thus `alpha=d/R` is
the sharp one-skeleton independence threshold, whereas `ceil(d/7)` is the
smaller facet size actually needed here; this explains both the old pair
method and the present pair-free gain in one object.
For the `194` next-frontier rows at `z_1=246`, the cheapest shape sidecar is
the same theorem's exact cyclic-gap formula: for each primitive `u`, the
affine obstruction exists exactly when `uS` has a forward gap at least
`d-ceil(d/7)+1`.  It removes the continuous translate before any search over
band centers.  THM-2984's all-cardinality normal form then replaces affine
containment by one diagonal unit carrying an ordered difference chain into a
positive composition simplex of span at most `ceil(d/7)-1`.  At prime
denominators, triples reduce further to six anharmonic ratios against a finite
short-summand table.  These are next-frontier sidecars only; neither is used
in the present count closure.
The centered-beta hostile `(15)` is part of the semantic payload.

## 5. Exact verification and scope

The computation pins the THM-2941 engine, projected atlas source/output, and
the repaired THM-2979 source/output by SHA-256 of LF-normalized bytes, while
rejecting bare carriage returns.  THM-2984 `(7)`, `(13a)--(13c)`, and
`(23a)--(23d)` are proved mathematical dependencies; the compositor
independently records `d,R,alpha,|S|,kappa` for every use.  The five pinned
computational hashes are, in order,

```text
d062c7ac8ebf6a433c8fb1543293e941c85625e2eb40b82fcf05fc2404539b0a
2af6d96882f336a409a8657070ed76a75c09a53b3789101b83103b051e864ded
cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda
4137ab250def3ad6a66b4c75a5e1b5b1a82ba4100b00ea5f8616faa46fb501a9
eea98955f91371d38d95cdeeb88b60a2305d34d0bddd2ea26570af8eede1b8e3.
```

From the repository root, the final normal and optimized commands are

```bash
SCRIPT=04-computation/lrc14_j7_k3_z270_to_z247_cardinality_torsion_descent_thm2941.py
OUT=05-knowledge/results/lrc14_j7_k3_z270_to_z247_cardinality_torsion_descent_thm2941.out
python3 "$SCRIPT" --processes 8 --output "$OUT" \
  > /tmp/lrc14_k3_z270_to_z247.final-normal.stdout
python3 -O "$SCRIPT" --processes 8 \
  --output /tmp/lrc14_k3_z270_to_z247.final-opt.out \
  > /tmp/lrc14_k3_z270_to_z247.final-opt.stdout
cmp "$OUT" /tmp/lrc14_k3_z270_to_z247.final-normal.stdout
cmp "$OUT" /tmp/lrc14_k3_z270_to_z247.final-opt.out
cmp /tmp/lrc14_k3_z270_to_z247.final-opt.out \
  /tmp/lrc14_k3_z270_to_z247.final-opt.stdout
```

The transcripts are byte-identical to each other and to the stored output.
All exact controls pass.  The theorem closes the `423` projected body rows
in its stated interval and nothing beyond them; in particular it does not
close the `194` rows at `z_1=246`, any arbitrary `k<=1` branch, the rung, or
LRC(14).  This proves `(1)` and `(2)`.  QED.

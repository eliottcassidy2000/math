# Recovered ten-body counterexamples to the new Haar floor

**Status: REFUTED universal `6/77` ten-body floor; FINITE-EXACT witnesses
independently reconstructed by three methods.** The exact body measures
below were already present in the workspace. The result here is their
recovery and application to the new scale-three Haar consumer, with complete
closed geometry and the fixed-clock boundary retained. No new theorem ID,
new safe LRC family, or literature novelty is claimed.

Throughout,

```text
G_C={y in R/Z: ||cy||>=1/14 for every c in C},
mu=normalized Haar measure.
```

The body phase is `y`; in a scale-three completion `3C union T`, a physical
phase is `x`, with `y=3x mod1`. Weak equality is safe; strict danger has
the same measure but different endpoint membership.

## 1. The universal floor is already contradicted at height thirteen

The canonical body `C0={1,...,10}` passes the proposed threshold:

```text
mu(G_C0)=1217/8820 > 6/77.
```

But the ten-body

```text
C1=(1,2,3,4,5,7,8,9,11,13)
```

has

```text
mu(G_C1)=21514/315315
        =6/77-3056/315315 <6/77.                       (1)
```

Here is a compact exact witness for the entire safe set. Let `J` be the
union of these five closed intervals:

```text
[15/182,13/154], [15/154,13/126], [29/182,27/154],
[29/98,55/182],  [29/70,41/98].
```

Then

```text
G_C1=J union (1-J) union {1,3,5,9,11,13}/14.           (2)
```

Intersecting the exact safe teeth
`[(14k+1)/(14c),(14k+13)/(14c)]`, `0<=k<c`, proves (2);
the independent wall-cell and danger-union calculations reproduce it.
Twice the sum of the five positive interval lengths gives (1). The six
remaining safe points have zero measure and must still be retained for
phase transport.

An even smaller measure occurs at

```text
Cmin=(1,2,3,5,7,8,9,11,12,13),
mu(G_Cmin)=14249/252252=6/77-5407/252252.                (3)
```

This is the minimum among all 286 ten-element subsets of `{1,...,13}`.
Twelve of those bodies fail the floor. Every one of the 66 ten-element
subsets of `{1,...,12}` passes, so thirteen is the smallest possible
maximum label for a ten-body counterexample. These minimality and minimum
claims are **FINITE-EXACT in the stated complete universe**.

The removal mechanism is visible without a broad census. Starting from
`B=(1,2,3,5,7,8,9)`, whose safe measure is `583/2205`, adjoin the remaining
three labels:

| Added label | New safe measure | Measure removed at this step |
|---:|---:|---:|
| 11 | `31807/194040` | `6499/64680` |
| 12 | `22873/194040` | `1489/32340` |
| 13 | `14249/252252` | `154859/2522520` |

The last deletion removes enough of the remaining safe components to cross
the new `6/77` threshold. Nonemptiness of the ten-body safe set does not
provide that quantitative measure floor.

## 2. The necessary fixed-clock sieve does not repair the floor

Fix the sharp tail `T=(1,5,11)` and take the inherited ten-body

```text
C2=(1,3,4,10,11,13,14,16,17,18).
```

Its exact safe measure is

```text
mu(G_C2)=534689/7796880
        =6/77-801461/85765680 <6/77.                  (4)
```

It has 22 positive safe components and the four isolated points
`{3,11,17,25}/28`. The full list of rational components is in the matching
transcript and is reproduced independently by the wall-cell scan.

Its literal completed row is

```text
S2=3C2 union T
  =(1,3,5,9,11,12,30,33,39,42,48,51,54).             (5)
```

For every integer denominator `2<=q<=18`, some speed in `S2` is divisible
by `q`. In order of `q=2,...,18`, one such speed is

```text
12,3,12,5,12,42,48,9,30,11,12,39,42,30,48,51,54.
```

Thus every unit clock of denominator at most eighteen fails. This is a
genuine hostile to deriving the Haar floor from the necessary small-clock
divisibility sieve; it does not assert entry into any other counterexample
stratum. Nevertheless,

```text
x=9/19, y=3x mod1=8/19,
min_(v in S2)||vx||=2/19>1/14,
mu(G_S2)=8131/194480>0.                              (6)
```

The least denominator of a rational safe phase is nineteen, verified by
the complete phase scan and already forced below that point by the displayed
divisibility list. The scalar Haar gate fails on this safe row. Its exact
body phase succeeds, so Haar sufficiency is not an equivalence.

A bounded entry-sieve probe held `T` fixed and required the necessary
condition that every denominator `2,...,14` divide some speed in `3C union T`.
The complete eligible bodies of maximum labels 14, 15, 16, and 17 numbered
209, 252, 917, and 1,227; all passed the floor. At maximum label eighteen,
the first failing body in the specified lexicographic traversal was `C2`,
after 346 eligible bodies. Hence eighteen is the smallest height for this
**fixed-tail, necessary-sieve-only** hostile. There is no claim about
smallest height under other tails or stronger actual entry conditions.

Integer dilation cannot repair the measure obstruction:
`G_(dC)=m_d^(-1)(G_C)` and Haar invariance gives equal measures. Dilating `C2`
also preserves its necessary divisibility list in `3dC2 union T`. For
`d=1 mod19`, the very same physical phase `9/19` retains clearance `2/19`.
Thus arbitrarily large completed rows can preserve both the failed Haar
gate and a valid simple phase. This is a residue transport of the named
witness, not a new difficult LRC family.

## 3. Inheritance and novelty audit

The closest proved mechanism is
[THM-4150, full-safe-set Haar transfer](../../01-canon/theorems/THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer.md),
now applied through the [new scale-three consumer](lrc14_universal_haar_consumer_empty_core_certificate_sep06.md).
The canonical synchronization hostile is
[THM-4032, d=3 affine defect boundary, Section 6](../../01-canon/theorems/THM-4032-lrc14-d3-affine-defect-lattice-boundary.md):
for `C0`, the body phase `y=2/11` is safe but the three lifts are spoiled by
different tails in `(1,5,11)`. The full row is safe at `x=1/14`. The script
replays that exact owner partition and both original phases.

The corrected near miss is
[THM-4154, Haar-pool inheritance correction](../../01-canon/theorems/THM-4154-mod-six-fixed-clock-and-haar-pool-inheritance-correction.md).
The least-used sidecar in the present probe is the complete **closed** safe
set, including isolated points, together with the actual physical phase.
The concept board was: cardinality; full Haar mass; component endpoints;
common dilation; fixed-clock divisibility; body-to-tail phase synchronization.

Exact-number and exact-body searches found the data before any claim of
novelty was filed:

- `1217/8820` for `C0` appears in the earlier leg engine and
  `LRCMPCertSize10.lean`.
- `21514/315315` for `C1` appears in the
  [THM-544 two-replacement AP-tail output](lrc14_two_replacement_ap_tail_theorem_codex_s35.out),
  at holes `(6,10,12)`, and in the earlier multi-outlier output. The matching
  bands appear in `LRCMPCertSize10.lean`; that Lean text is provenance here,
  not a new build or axiom-audit claim.
- `14249/252252` for `Cmin` is the holes `(4,6,10)` row of the same
  [THM-544](../../01-canon/theorems/THM-544-lrc14-two-replacement-ap-tail.md)
  output.
- `534689/7796880` for `C2` appears in the second/grandchild ledger of
  [THM-2923, seven-body recursive pair-hunter closure](../../01-canon/theorems/THM-2923-complete-seven-body-six-slot-recursive-pair-hunter-closure.md),
  with base `(1,3,4,10,11,13,14)` and added labels `18,17,16`.
  Its historical intermediate open/closed flags are not current truth;
  the exact body and measure are independently reconstructed here.

The new use of these recovered objects is to refute the particular body
floor proposed after the universal triple-network theorem. Their measures
are not new results. The small completion `3C1 union T` was already safe
at the inherited small clock `1/10`, with clearance `1/10`; denominators ten
and fourteen pass the elementary sieve. The stronger completion (5), of
height 54, also lies within the previously verified height-55 scope of
[THM-1290, corrected bounded-height census](../../01-canon/theorems/THM-1290-subgap-exhaustive-census-bounded-height.md).
No new safe family is therefore counted as LRC progress.

The neighboring eleven-body hostile in
[THM-4330, anchored pool entry, Section 6.1](../../01-canon/theorems/THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve.md)
already separates a failed dyadic Haar gate from actual row safety. The
present ten-body witnesses establish the corresponding boundary for the new
scale-three threshold, with its different body cardinality and tail set.

## 4. Repaired question and exact verification

The first failed implication is

```text
ten body speeds, even with necessary small-clock divisibility
  ==> mu(G_C)>=6/77.
```

It is false. The strongest survivor is the proved sufficient criterion
itself: if the actual body has enough Haar measure for its actual tail
network bound, completion follows. If it does not, retain component
addresses and test whether any body phase has a common safe lift. A proposed
entry theorem must force more than cardinality and the elementary divisor
sieve, or must supply that phase information directly.

The source-to-target map is `G_C -> mu(G_C)`. It preserves the exact measure
comparison needed for a sufficient completion gate, but destroys location,
component incidence, and isolated safe phases. The restoration sidecar is
the closed rational component list plus the physical lift label. Equations
(4)-(6) are the decisive test: the scalar gate fails while the actual phase
and completed safe-set measure certify success. No general entry route or
LRC(14) conclusion follows.

```text
python3 -B 04-computation/lrc14_haar_body_empty_core_sep06.py
python3 -B -O 04-computation/lrc14_haar_body_empty_core_sep06.py
```

[Source](../../04-computation/lrc14_haar_body_empty_core_sep06.py) and
[matching exact transcript](lrc14_haar_body_empty_core_sep06.out).
The primary closed-safe-tooth intersection is checked against both an
independent danger-union measure and an independent wall/midpoint geometry
scan on every named body and completed row. All 286 height-thirteen body
measures also receive the independent danger-union check. AP shifts, two
dilations, an odd AP, and four one-far bodies are positive controls; their
inherited status is retained. All arithmetic is rational, all checks remain
active under `-O`, and normal/optimized output bytes agree.

SHA-256 of raw LF bytes:

```text
source aee840392ca365cf43d089984dafe7cf44a8e4a97ca2e8c9e99f41b122901f7e
output e61747d8f90c8cace96d0b8b5c6c4b14d49fef47bbe543d318c3f03880f386f0
```

---
id: HYP-2624
title: LRC(14) height-2 coimage wall-class reduction
status: PARTIALLY-CONFIRMED
source: codex-2026-06-19-S18
depends_on:
  - HYP-2617
  - HYP-2619
  - HYP-2616
  - THM-538
related:
  - THM-539
  - HYP-2623
  - HYP-2614
  - HYP-2615
  - HYP-2608
  - OPEN-Q-108
---

# HYP-2624 - LRC(14) Height-2 Coimage Wall-Class Reduction

## Claim

The next useful quotient after the HYP-2617 support-6 coimage atlas is a
low-height wall-addressing quotient.

For `k=8,9,10`, let `d=k-1` be the ambient nonzero-offset dimension and let
`B(k)={16,15,13}` be the bounded-core ceiling from the finite certificate.  A
projective mod-7 support-six coimage class is called height-`H` one-large
wall-addressed if it is represented by a six-support

```text
{a_1,...,a_5,M},   1 <= a_i <= B(k),   M > B(k),
```

with a nontrivial integer relation

```text
c_M M + sum_i c_i a_i = 0,   c_i in {-H,...,-1,1,...,H}.
```

The computation verifies that height `<=2` one-large walls address every
nonzero signed-mass coimage class in the `k=8` and `k=9` ambient dimensions.
For `k=10`, they address most of the signed mass, and the remaining classes
are not random: the dominant tail-only packet consists of repeated-residue
classes of shapes `4+2` and `4+1+1`.

This is not a proof of LRC(14).  A coimage class can be represented by both a
low-height wall and a high-height tail, so class-addressing is a routing lemma,
not row deletion.  Its value is that the post-wall analytic target is much
narrower than the full 159-class atlas.

## Computation

Script:

- `04-computation/lrc14_height2_coimage_wall_classes_codex_s18.py`
- output: `05-knowledge/results/lrc14_height2_coimage_wall_classes_codex_s18.out`

The script imports the HYP-2617 atlas and computes `S_d(a)` for each projective
class.  It then enumerates one-large support-six wall supports with
coefficient height `H=1` and `H=2`, projects each support to the canonical
`F_7^*/S_6` coimage class, and asks how much nonzero coimage signed mass is
hit.

The main table is:

```text
k  H   supports  classes       M-range   nonzero hit  top10 hit  |S|-mass hit
8  1      29969      132        17..70    45/46         10/10      0.99009901
8  2     272007      147       17..140    46/46         10/10      1.00000000
9  1      20505      132        16..65    77/79         10/10      0.98434783
9  2     177648      147       16..130    79/79         10/10      1.00000000
10 1       8677      100        14..55    79/116         8/10      0.82073243
10 2      68124      108       14..110    85/116         8/10      0.84229179
```

Thus, after height `<=2` wall-addressing:

```text
k=8:  missed nonzero coimage classes = 0
k=9:  missed nonzero coimage classes = 0
k=10: missed nonzero coimage classes = 31
```

The missed `k=10` signed-mass share is:

```text
missed |S|-mass = 1.70105407, or 15.770821%
```

## Tail-Only Packet

For `k=10`, the missed nonzero classes have pattern histogram:

```text
4+2 repeated:     6
4+1+1 repeated:  12
zero-cusp:        3
mixed:           10
```

The largest tail-only classes are:

```text
class                  |S_d|        abs/signed ratio    pattern
(1,1,1,1,4,4)          0.23891209   19.3926             4+2 repeated
(1,1,1,1,2,2)          0.23891209   18.1869             4+2 repeated
(1,1,1,1,6,6)          0.1720167    29.7843             4+2 repeated
(1,1,1,1,3,3)          0.1720167    28.8616             4+2 repeated
(1,1,1,1,5,5)          0.1720167    29.7557             4+2 repeated
(1,1,1,1,2,3)          0.076451868  60.3104             4+1+1 repeated
(1,1,1,1,4,6)          0.076451868  67.5762             4+1+1 repeated
(1,1,1,1,2,5)          0.076451868  59.9927             4+1+1 repeated
(1,1,1,1,3,4)          0.076451868  61.6747             4+1+1 repeated
(1,1,1,1,2,4)          0.076451868  60.0765             4+1+1 repeated
```

This turns the analytic tail from "all non-null classes in the 159-class atlas"
into a repeated-residue reciprocal-sum problem in the `d=9` ambient layer, plus
a small set of zero-cusp and mixed companions.

## Examples

The top coimage classes for `k=8` and `k=9` all have height `<=2` one-large wall
support examples.  For `k=10`, the top eight have examples, while ranks 9 and
10 are tail-only:

```text
k=10, d=9
rank  class                  |S_d|        wall support example
1     (1,1,2,2,3,3)          0.43959824   (1,2,3,8,9,38)
2     (1,1,2,2,3,4)          0.4077433    (1,2,3,4,8,16)
3     (1,1,2,2,4,6)          0.33925517   (1,2,3,4,8,18)
4     (1,1,2,3,4,5)          0.3352733    (1,2,3,4,5,15)
5     (1,2,3,4,5,6)          0.30580747   (1,2,3,4,5,27)
6     (1,1,2,3,3,5)          0.2628033    (1,2,3,5,8,38)
7     (1,1,2,4,6,6)          0.2628033    (1,2,3,4,10,18)
8     (1,1,2,3,3,6)          0.2628033    (1,2,3,5,8,19)
9     (1,1,1,1,4,4)          0.23891209   TAIL-ONLY
10    (1,1,1,1,2,2)          0.23891209   TAIL-ONLY
```

## Interpretation

HYP-2616 proved that height-1 one-large support-six walls are finite and safe
by exact sector-measure enumeration.  HYP-2617 made the infinite support-six
tail finite at the coimage level.  HYP-2619 then showed that the true smallness
comes from signed alternation rather than absolute volume.

HYP-2624 connects these:

```text
low-height wall ledger
-> projective coimage class address
-> repeated-residue tail packet
-> signed reciprocal-tail theorem.
```

The point is not that height-2 wall classes have been cleared by exact measure.
The point is that the coimage classes with large signed mass are mostly already
visible from very low-height wall geometry.  Once those finite wall rows are
accounted for, the remaining analytic obstruction concentrates in a small
repeated-residue packet.

Number-theory reading: the repeated packet is the support-six analogue of a
multiple root or high-multiplicity residue stratum.  Four equal residues create
the persistent cusp, and the other two residues select the surviving character
fiber.  The expected theorem is therefore closer to a specialized cotangent or
Dedekind reciprocal-sum estimate for `(a,a,a,a,b,c)` classes than to a generic
Minkowski volume bound.

## Integration With THM-539

After this computation, KPS S9 landed THM-539/HYP-2623/T871: the max-min
second-spectrum gap can dip below the naive `Theta(1/k^2)` lower scale along
primorial height-escape families.  That result is adjacent but not identical
to the S18 target.

The primorial result concerns the max-min spectrum near the AP floor as `k`
varies and `max(S)/k` escapes.  HYP-2624 concerns the fixed LRC(14) support-six
signed correction after the THM-538 support floor.  The shared lesson is that
height matters: do not trust a bounded-height scalar heuristic as a global
proof.  In the support-six lane, this is exactly why class-addressing is only a
routing lemma.  Low-height wall classes must be accounted for finitely before
the high-height repeated-residue tail theorem is applied.

The later THM-539 binding-pair note sharpens the analogy: the escaping family
is controlled by a two-runner resonance `{2a-1, a(k-1)}` whose denominator is
`a(k+1)-1`, while the lcm-killer computation says coarse clocks are suppressed
only when `(k-1)` carries enough small factors.  The LRC(14) coimage tail has a
parallel but more local shape: repeated residue packets `(a,a,a,a,b,c)` are
not generic high-height noise, but resonant packets where many row clocks have
collapsed onto the same projective residue.  A plausible proof route should
therefore look for a signed cotangent/Dedekind estimate that binds those
repeated packets directly, rather than trying to extend the height-2 quotient
as if it were a uniform spectral lower bound.

## Tournament Analysis

The session explicitly challenges the runner-vertex default.  Raw runners,
raw support tuples, and pairwise arcs do not preserve the predicate being
optimized here.  The quotient preserving the LRC(14) support-six tail predicate
is:

```text
height-bounded wall class
-> projective coimage packet
-> repeated-residue tail obligation.
```

The proof-quotient tournament is transitive, with Hamiltonian path:

```text
height2_wall_addressed_classes
> height1_wall_addressed_classes
> repeated_residue_tail_packet
> coimage_fiber_atlas
> signed_reciprocal_tail_theorem
> raw_supports
> raw_runner_vertices
```

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3_cycles = 0
SCC_sizes = [1,1,1,1,1,1,1]
```

The quotient destroys witness-time geometry and row identity, so it cannot by
itself prove the finite wall ledger.  It preserves the analytic address of the
signed support-six tail, which is exactly the information needed for the next
Dedekind/cotangent theorem.

## HYP-2630 Update

HYP-2630 checks the next natural escape route and rules out a simple version of
it.  Raising the one-large wall enumeration from coefficient height `2` to
height `3` still hits the same `85/116` nonzero `k=10` coimage classes and the
same `84.229179%` signed mass.  Thus the `31` tail-only classes are not merely
height-3 one-large walls.

The obstruction is residue multiplicity.  The bounded core `1..13` supplies at
most two representatives of each nonzero residue modulo `7`; one large speed
can add only one more.  The dominant `4+2` and `4+1+1` repeated packets
therefore require at least two large residue coordinates.  Next step (3) below
should be read as a two-large repeated-residue character theorem, not as a
continued one-large height increase.

## Status

Partially confirmed by exact enumeration of the height `<=2` one-large wall
supports and exact recomputation of the HYP-2617 coimage masses.  LRC(14)
remains open.

Next steps:

```text
1. Run exact sector-measure clearing for the finite height-2 one-large wall rows
   that are not already covered by HYP-2616.
2. Prove a repeated-residue reciprocal-tail estimate for the d=9 classes
   (1,1,1,1,a,b), with zero-cusp companions handled separately.
3. Check whether multi-large low-height walls introduce new coimage classes or
   only recycle the same repeated-residue packet.
4. Fold the resulting finite deletion into the HYP-2608 wide-spread assembly.
```

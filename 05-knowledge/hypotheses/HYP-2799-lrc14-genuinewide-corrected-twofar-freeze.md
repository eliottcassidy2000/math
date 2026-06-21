---
id: HYP-2799
title: LRC14 genuine-wide two-far tail uses the y in [0,7) frozen phase, leaving a finite low-gap window
status: OPEN; corrected tail scout stored, finite window and gap monotonicity still open
source: codex-2026-06-21-S77
tangent: T956
depends_on:
  - HYP-2793
  - HYP-2788
  - HYP-2701
  - HYP-2708
  - HYP-2684
  - HYP-2679
  - THM-548
  - THM-557
related:
  - HYP-2680
  - HYP-2676
  - HYP-2682
  - OPEN-Q-108
---

# HYP-2799: Corrected Two-Far Freeze Tail

## Claim

The r=2 genuine-wide tail should be proved with the frozen two-far phase over
the same circle coordinate used by the exact row engine:

```text
y in [0,7),    sector(e,y)=floor(e*y) mod 7.
```

The existing Thread-A certificate script had the right idea but used an
`x in [0,1)` diagonal-freeze integral.  Its own sanity check fails:

```text
base=(0,1,2,3,4,5,6,7), g=1:
  old Dblock=66/343 ~= 0.19242,
  large-f average ~= 0.41102.
```

With the corrected coordinate, for far pair `{f,f+g}` and `f -> infinity`, the
frozen law at slow coordinate `y` is:

```text
a uniform in Z/7,
second color = a + floor(g*y + theta) mod 7,
theta uniform in [0,1).
```

Equivalently, after averaging over `theta`, the offset is `floor(g*y)` with
weight `1-frac(g*y)` and `floor(g*y)+1` with weight `frac(g*y)`.

This corrected frozen value `D7(B,g)` matches finite large-f averages and is
still below the Thread-A floor `Q(k-1)` in the all-base gap window `g<=24` for
`k=10,11,12`.  Therefore the r=2 genuine-wide tail is not an obstruction; it
reduces to:

```text
finite low-f window  +  proof/verification that wider gaps do not beat the
adjacent-gap frozen leader.
```

## Exact Scout

Script:

```text
04-computation/lrc14_genuinewide_corrected_freeze_tail_codex_s77.py
```

Output:

```text
05-knowledge/results/lrc14_genuinewide_corrected_freeze_tail_codex_s77.out
```

Corrected sanity checks:

```text
base=(0,1,2,3,4,5,6,7), g=1:
  D7=7895/19208 = 0.4110266556,
  large-f average = 0.4110197286.

base=(0,1,2,3,4,5,6,7), g=2:
  D7=877969/2160900 = 0.4062978389,
  large-f average = 0.4062879019.

base=(0,1,2,3,4,5,6,7,8), g=1:
  D7=23557/48020 = 0.4905664307,
  large-f average = 0.4905552191.
```

Complete all-base frozen scan for `g<=24`:

```text
k=10:
  max D7 = 7895/19208 = 0.4110266556
  Q(9) = 1229/2744 = 0.4478862974
  room = 177/4802 = 0.0368596418
  leaders: base=(0,1,2,3,4,5,6,7), g=1
           and dilation (0,2,4,6,8,10,12,14), g=2

k=11:
  max D7 = 23557/48020 = 0.4905664307
  Q(10) = 65599/123480 = 0.5312520246
  room = 35167/864360 = 0.0406855939
  leader: base=(0,1,2,3,4,5,6,7,8), g=1

k=12:
  max D7 = 14085881/24893568 = 0.5658441972
  Q(11) = 14873/24696 = 0.6022432783
  room = 906103/24893568 = 0.0363990811
  leader: base=(0,1,2,3,4,5,6,7,8,9), g=1
```

Using the crude Thread-A tail envelope `D7(B,g)+7/f<Q(k-1)`, these exact rooms
give tail cutoffs:

```text
k=10: f > 190
k=11: f > 173
k=12: f > 193
```

The targeted dense guardrail for `k=10`, `g<=120`, `max base gap<=3` keeps the
same leader and same room.  That does not prove gap monotonicity, but it
supports the next finite theorem:

```text
for every bounded base B of size k-2 and every gap g>=1,
  D7(B,g) <= D7(consec_{k-2},1)
```

or the weaker version needed for the proof,

```text
D7(B,g) <= Q(k-1) - eta_k.
```

## Actual-Size Doublet Signed-Tail Addendum

Incoming HYP-2797 identifies the genuine-wide doublet leader family

```text
E_N(M) = {0,1,...,N-3} union {M,M+1}.
```

The HYP-2797 narrative is actual-size aligned, but one companion signed-bound
table used shifted labels: the printed row `range(k-1) union {M,M+1}` has
actual size `k+1`.  The S77 alignment audit records the corrected cap/Q
comparison:

```text
04-computation/lrc14_doublet_index_alignment_codex_s77.py
05-knowledge/results/lrc14_doublet_index_alignment_codex_s77.out
```

Using the corrected HYP-2799 plateau `D7({0,...,N-3},1)`, the exact signed
scan

```text
04-computation/lrc14_doublet_exact_plateau_signed_bound_codex_s77.py
05-knowledge/results/lrc14_doublet_exact_plateau_signed_bound_codex_s77.out
```

computes `M*(p0(E_N(M))-D7)` for `15<=M<=600`.  The actual-size p0 maxima
are:

```text
N=8:  M=20, p0=85/588
N=9:  M=21, p0=9371/32340
N=10: M=21, p0=1301/2940
N=11: M=21, p0=5617/10780
N=12: M=21, p0=14302/24255
```

The positive signed constants

```text
C+ = max_M M*(p0-D7)

N=8:  206431/149450
N=9:  826621/655200
N=10: 1265857/946680
N=11: 74483/57624
N=12: 2317980841/2115953280
```

are already small relative to the actual cap margins.  If the displayed `C+`
values are proved uniformly for all `M>=15`, then

```text
p0(E_N(M)) <= D7 + C+/M < cap_N
```

holds from `M=15` for every `N=8..12`.  The stricter `Q(N-1)` comparison is
still finite-prefix sized rather than asymptotic:

```text
suggested Q cutoffs from the scanned C+:
N=8,9,10,11,12 -> M>=16,13,37,32,31.
```

This addendum separates two proof targets that were previously being mixed:
the genuine-wide cap branch should use the wide actual-size cap margins, while
the sharper `Q` branch only needs a short exact doublet prefix after a uniform
signed endpoint/Dedekind constant is proved.

## Proof Route

The corrected r=2 genuine-wide branch now has three separate obligations:

1. **Coordinate lemma.**  Prove the `y in [0,7)` frozen formula from the
   exact breakpoint definition of `p0`.  This is finite measure theory over
   rational wall cells: on each base cell, `f*y mod 7` equidistributes and
   the pair offset is `floor(g*y)` / `floor(g*y)+1` with the above weights.
2. **Frozen gap extremality.**  Prove or exhaustively certify the all-gap
   inequality for `D7(B,g)`.  The evidence says the adjacent consecutive base
   is the frozen leader, with dilated copies tied when the gap is dilated.
3. **Finite low-f window.**  After `D7+7/f<Q`, the remaining exact check is
   finite.  With the current crude envelope it is enough to scan roughly
   `15 <= f <= 193`, `g<=24` for all bounded bases at `k=10..12`; a sharper
   Abel/BV constant or live-depth kernel can shrink this.

This corrects the status of the old Thread-A r=2 certificate.  Its stated
tail mechanism was promising, but the stored script did not implement the
all-base tail proof and its old `Dblock` sanity check was false.  HYP-2799
turns that into a precise repair rather than treating the comment as a proof.

## Relation To Survival Currency

HYP-2799 is a direct-`p0` frozen tail.  It does not replace HYP-2701/HYP-2708.
Instead, it gives the large-f side of the same two-far branch:

```text
large f:
  direct p0 frozen tail already below Q(k-1)

small f / low gap:
  use exact finite atlas plus HYP-2708 live-depth survival kernel
```

The survival route is still valuable because the finite low-gap rows spend
boundary margin through only four live depths `{1,2,5,6}`.  That is likely the
right compression for a human proof of the finite window.

## Assumption Challenge

The useful tournament vertices are not runners.  Candidate vertex sets tested
or retained:

- bounded bases `B`;
- far-pair gaps `g`;
- frozen proof obligations `D7`, finite-window rows, and gap monotonicity;
- live-depth packets `{1,2,5,6}`;
- relation equations among far speeds;
- scale-reduction scaffolds.

The quotient preserves the direct predicate `p0(E)<Q(k-1)` in the r=2 tail and
the route to `p0(E)<cap_k` through HYP-2793.  It destroys exact finite-f
endpoint phase, so the small window must be checked or handled by the
endpoint/state-word machinery before claiming closure.

## Tournament Analysis

Proof-lens tournament under "removes the largest genuine-wide uncertainty":

```text
correct_y_coordinate
> frozen_gap_extremality
> finite_low_f_window
> live_depth_kernel
> relation_lattice_atlas
> raw_scalar_ET_bound
```

The order is transitive in this session's evidence.  The key edge flip versus
the previous proof order is that `correct_y_coordinate` must precede any tail
constant: using the wrong coordinate makes the tail look two times safer than
it is and invalidates the claimed certificate.

## Status

No proof of LRC(14) is claimed here.  The concrete progress is:

- found and repaired a tail-coordinate bug in the r=2 genuine-wide certificate;
- produced a corrected exact frozen value matching finite large-f averages;
- verified all bounded bases with `g<=24` at `k=10,11,12` have frozen tail
  below `Q(k-1)`;
- reduced the remaining r=2 proof to frozen gap extremality plus a finite
  low-f exact window.

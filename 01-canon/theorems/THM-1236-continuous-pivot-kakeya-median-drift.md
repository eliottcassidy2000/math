---
id: THM-1236
title: THE CONTINUOUS-PIVOT KAKEYA MEDIAN LAW — every auxiliary full-lap scale pays drift, with the fourth speed as the exact constrained optimizer
status: PROVED (all-real auxiliary-pivot moving-arc theorem; exact constrained tilted-L1 optimization; dependency-free Fraction referee; sorry-free Lean arithmetic core)
source: codex-2026-07-19-S78 continuation
depends_on: [THM-1235]
related: [THM-1198, THM-1219, THM-1232, THM-1233, THM-1238]
script: 04-computation/lrc14_continuous_pivot_kakeya_median_thm1236.py
output: 05-knowledge/results/lrc14_continuous_pivot_kakeya_median_thm1236.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCContinuousPivotKakeyaMedian.lean
script_sha256: c5acb1e17460f4c90e62a3974eb136a619fca9e84339166319ec92a476af93b4
output_sha256: fb3eb2032cbf3df7442ac069319a0f897bc7146358fe43fbf893d0d4b97121eb
formalization_sha256: 28445aa28a21b531a36aedcc1c27776081f52f188d9fc0cbfad7b376bffbde83
---

# THM-1236 — continuous-pivot Kakeya median drift

## Statement

Let `c>0`, let

```text
G=[(14k+1)/(14c),(14k+13)/(14c)]
```

be a complete closed safe gap of the `c`-comb, and suppose that the six
strict danger combs of ordered speeds

```text
c<d1<d2<d3<d4<d5<d6
```

cover `G`.  Then every real auxiliary scale making one normalized lap pays
the drift invoice

```text
sum_(i=1)^6 |di-y| > y/14              for every real y>=7c/6.       (1)
```

This genuinely strengthens THM-1235's runner-pivot statement: `y` need not be
one of the six speeds.  Define the tilted drift functional

```text
H(y)=sum_i |di-y|-y/14.                                  (2)
```

Its exact constrained minimizer is

```text
y0=max(7c/6,d4).                                        (3)
```

Consequently the infinite family (1) is equivalent, as a scalar necessary
condition, to the single optimized invoice `H(y0)>0`.  In particular:

```text
d4>=7c/6  =>
d4+d5+d6-d1-d2-d3 > d4/14;                              (4)

d4<7c/6  =>
sum_i |di-7c/6| > c/12.                                 (5)
```

Together with THM-1235's `d6>7c/6`, (3) recovers the suffix-count edge ledger

```text
r=#{i:di>=7c/6},       M=max_j(d_(j+1)-d_j),

r=1  => M>c/180,
r=2  => M>c/132,
r>=3 => M>c/108.                                       (6)
```

The new content is the functional form and exact optimizer, not stronger
constants in (6).

## Proof of the continuous-pivot invoice

Normalize `G` by

```text
t=(14k+1)/(14c)+(6/(7c))x,       0<=x<=1,
Li=6di/(7c),
phi_i=di(14k+1)/(14c) mod 1.                           (7)
```

Fix any real `y>=7c/6` and put

```text
L*=6y/(7c)>=1.                                         (8)
```

The initial subneedle `0<=x<=1/L*` lies in the normalized gap.  Set
`u=L*x`; then `u` traverses one complete circle and the `i`th moving danger
arc is

```text
||phi_i+(di/y)u||<1/14
iff ||phi_i+u+(di/y-1)u||<1/14.                        (9)
```

Freeze the six centres at `u=0`.  Circle distance is one-Lipschitz, so the
moving `i`th arc is contained in the static open arc of radius

```text
1/14+|di/y-1|.                                        (10)
```

If some radius in (10) is at least `1/2`, then

```text
|di/y-1|>=1/2-1/14=3/7>1/14,                          (11)
```

and (1) follows.  Otherwise all six enlarged arcs are nonempty proper open
circle arcs, of lengths

```text
1/7+2|di/y-1|.                                        (12)
```

They cover the circle because the original moving arcs cover the subneedle.
Their total length is strictly greater than one.  Indeed measure gives at
least one; equality would make every pairwise overlap null.  A nonempty
intersection of open arcs has positive measure, so equality would make the
six arcs pairwise disjoint, contradicting connectedness of the circle.
Therefore

```text
1 < 6/7+2 sum_i |di/y-1|,
sum_i |di/y-1|>1/14.                                  (13)
```

Since `y>0`, multiplying (13) by `y` proves (1).  Notice that no speed
integrality, equidistribution, or phase averaging enters this step.

## Exact minimization of the H-drift

The ordinary six-point `L1` functional is minimized throughout the median
interval `[d3,d4]`.  The term `-y/14` breaks that flat interval toward the
right, so `d4` is the unconstrained minimizer of (2).  Here is an elementary
inequality proof retaining the exact slope.

For `y<=d4`, pair the coordinates as `(d1,d4),(d2,d5),(d3,d6)`.  The triangle
inequality gives

```text
sum_i |di-y| >= (d4-d1)+(d5-d2)+(d6-d3)
              =sum_i |di-d4|,                         (14)
```

and `-y/14>=-d4/14`.  Hence `H(y)>=H(d4)`.  For `y>=d4`, the first four
coordinates move out by `4(y-d4)` while the final two can move in by at most
`2(y-d4)`.  Thus

```text
sum_i|di-y|-sum_i|di-d4| >=2(y-d4),
H(y)-H(d4)>=27(y-d4)/14>=0.                           (15)
```

The same argument shows `H` is increasing to the right of `d4`.  Restricting
to `y>=7c/6` proves (3).  Expanding at `d4` gives the half-sum identity (4),
while substituting the boundary scale gives (5).

For (6), the absolute drift at `y0` is an adjacent-gap weighted sum.  When
`r=1`, the threshold lies between `d5,d6` and its total weight is at most
`15M`; when `r=2`, it lies between `d4,d5` and the weight is at most `11M`;
when `r>=3`, `y0=d4` and the exact weight vector is

```text
(1,2,3,2,1),                                          (16)
```

of total weight nine.  Since `y0/14>=c/12`, (6) follows.

## Structural consequence and carrier audit

Equation (2) is the useful closed form of the previously informal `H`-drift.
It is a convex piecewise-linear support functional with slopes

```text
2j-6-1/14                    on (dj,d_(j+1)).          (17)
```

The sign changes only after the upper median.  Thus the continuum of pivot
tests is not a new infinite case split; it is one exposed face of the ordered
speed cone.  This reframes the hard packet as two three-speed half-packets
separated by an exact top-minus-bottom mass imbalance.

The pairwise observable for Tournament Analysis remains the adjacent speed
gap, oriented by increasing speed.  Its unweighted tournament is transitive,
with score histogram `(0,1,2,3,4,5)`, no directed cycles, six singleton SCCs,
and one Hamiltonian path.  The functional `H` retains the weighted cut across
that path and its constrained median; the bare tournament destroys both.

We challenged runners, gaps, teeth, endpoints, wall events, moving arcs,
Fourier modes, residues, and proof obligations as vertices.  Moving arcs with
an auxiliary lap-scale sidecar preserve the cover predicate used in the
proof.  Compressing them to the single convex functional preserves exactly
the drift consequence but destroys phase and blocker identity.  THM-1238
therefore lifts the resulting macroscopic cut back to pair-sum clock data.

## Verification and honest frontier

The dependency-free `Fraction` referee checks all constants, `8,008` ordered
six-packets, `19,448` carrier thresholds, and `290,796` exact breakpoint and
midpoint comparisons, including targeted instances of all six suffix counts.
Normal and optimized Python outputs are byte-identical.  The Lean module
kernel-checks the tilted upper-median theorem, right-tail monotonicity,
constrained optimizer, half-sum identity, and branch constants.  It leaves
the moving-open-arc cover as the explicit paper provider and contains no
proof placeholders.

Frozen hashes are

```text
source         c5acb1e17460f4c90e62a3974eb136a619fca9e84339166319ec92a476af93b4
output         fb3eb2032cbf3df7442ac069319a0f897bc7146358fe43fbf893d0d4b97121eb
formalization  28445aa28a21b531a36aedcc1c27776081f52f188d9fc0cbfad7b376bffbde83
```

THM-1236 rules out coalescence in its sharp functional form.  It does not
control the phase clocks on either side of the median cut, prevent third-tooth
blocking, or prove LRC(14).

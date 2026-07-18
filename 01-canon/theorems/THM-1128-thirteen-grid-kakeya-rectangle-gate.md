---
id: THM-1128
title: The universal thirteen-grid Kakeya rectangle closes a centred four-comb cone
status: PROVED — elementary torus geometry; a centred slope B satisfying B >= 53 max(A,24), where A is the maximum centred offset, leaves a component longer than 1/(7k4). The midpoint corollary closes k1 >= max(1272,26(Delta+1))
source: codex-2026-07-18-S73 Kakeya continuation
depends_on: [THM-1097, THM-1126, THM-1127]
related: [THM-1101, THM-1123, THM-1129, THM-1132, MISTAKE-164]
---

# THM-1128 — the thirteen-grid Kakeya rectangle gate

Let `P` be any subset of `{1,...,12}`, let

```text
0<k1<k2<k3<k4,
B in [k1,k4] an integer,
a_i=k_i-B,       A=max_i |a_i|,       M=max(A,24).
```

Write

```text
E={t in R/Z : ||pt||>=1/14 for p in P
                  and ||k_i t||>=1/14 for i=1,...,4}.
```

Then

```text
B >= 53 M                                             (1)
```

implies that `E` contains an interval of length

```text
352/(2457 B) > 1/(7k4).                              (2)
```

For the optimized integer midpoint, put

```text
Delta=k4-k1,
A=ceil(Delta/2),
B=k1+A.
```

Thus the transparent sufficient cone

```text
k1 >= max(1272,26(Delta+1))                          (3)
```

implies (1), and hence (2).  When `Delta` is even and at least `48`, the
sharper exact midpoint condition is `k1>=26 Delta`.

The statement is independent of the choice and size of `P`; the usual
eight-speed core is only one application.  It is an all-scale theorem on a
positive cone of offset shapes, not a finite scan or a fixed translated ray.

## Proof

Put `t0=1/13`, and use the fixed-polygon coordinates of THM-1127,

```text
x=Bt mod 1,
H_A={(t,x): ||x+a_i t||>=1/14 for every i}.
```

### 1. Every small core is safe near `1/13`

For `1<=p<=12`, the nonzero residue `p mod 13` gives

```text
||p/13|| >= 1/13.
```

Let

```text
epsilon=53/(4914M),
I=[t0-epsilon,t0+epsilon].                            (4)
```

For `t in I`, reverse Lipschitz continuity gives

```text
||pt|| >= ||p/13||-p|t-t0|
        >= 1/13-636/(4914M)
        >= 1/13-1/182
         = 1/14.                                     (5)
```

Thus `I` is safe for every possible core `P subset {1,...,12}`.

### 2. Four offset centers leave a common vertical strip

At `t0`, the four danger-arc centers `-a_i/13` lie on the thirteen-point
grid.  Coalesce repeated centers and list the positive cyclic gaps between
the remaining `m<=4` centers in grid units.  They are positive integers
summing to `13`, so one is at least

```text
ceil(13/m) >= 4.
```

The corresponding center gap is at least `4/13`.  Removing radius `1/14`
at both ends leaves a vertical safe strip of width at least

```text
4/13-1/7=15/91.                                      (6)
```

Choose compatible affine lifts of the bounding wall clusters over `I`.
Every wall has slope `-a_i`, hence absolute slope at most `A`.  Other centers
start outside the chosen cyclic gap; repeated centers are included in the
two endpoint clusters.  Since `2A epsilon<=53/2457<1/13`, distinct grid
clusters cannot change cyclic order inside `I`.  Shrinking by the maximum
motion at both ends
therefore leaves a common strip for every wall.  The total possible loss is
at most `2A epsilon`, and by (4)

```text
2A epsilon <=106/4914=53/2457,
15/91-2A epsilon >=405/2457-53/2457
                       =352/2457.                    (7)
```

Thus there is a fixed circle arc `X` of length exactly `352/2457` such that

```text
I x X subset H_A.                                    (8)
```

If all four centers coincide, use the full cyclic gap between two adjacent
lifts of that same wall; the same argument is even stronger.

### 3. A slope-`B` needle crosses the rectangle completely

The preimages of `X` under `t -> Bt mod 1` are intervals of length
`352/(2457B)` whose left endpoints are spaced by `1/B`.  An interval of time
length `ell` contains a complete such preimage whenever

```text
ell-352/(2457B) >= 1/B,
```

or equivalently `B ell>=2809/2457`.  From (1) and (4),

```text
B |I| = 53B/(2457M) >=53^2/2457=2809/2457.           (9)
```

Hence some complete preimage interval `J` lies in `I`.  Equations (5) and
(8) show that every `t in J` is safe for the core and all four killers, so
`J subset E` and `|J|=352/(2457B)`.  Finally, because `B<=k4`,

```text
352/(2457B) > 1/(7k4)
  iff 2464k4>2457B,
```

which is strict.  This proves (2).

For the midpoint corollary, `A=ceil(Delta/2)` and `B=k1+A`.  If `A<24`,
`k1>=1272` gives `B>=53*24`.  If `A>=24`,
`k1>=26(Delta+1)>=52A` gives `B>=53A`.  This proves (3).

## What the cone removes from the open r=5 problem

THM-1101's coarse mass/component tail improves only when the fourth speed is
well separated from the first three.  THM-1128 works at the opposite end:
it closes a uniform near-diagonal cone by exploiting phase coherence rather
than treating the four combs independently.  The remaining scale-free
region is therefore bounded away from exact coincidence; every tuple left by
this gate satisfies roughly

```text
Delta > k1/26.
```

Other discrepancy gates remove portions of that complement, but they are not
a condition on `Delta` alone.  The integer constant `53` is the smallest
integer yielded by this universal fixed-rectangle optimization.  The common
vertical arc is only `1/2457` wider than `1/7`; the
remaining `53/2457` of the universal `15/91` strip pays wall drift, while a
full preimage costs one additional winding period.  Centering at the midpoint
halves the wall-slope budget compared with taking `B=k1`.

### Why four removals are the exact thirteen-grid boundary

The residue pigeonhole is arity-sharp.  With four centers, a cyclic gap of
at least `4/13` is forced and leaves safe vertical width `15/91>1/7`.  With
five centers, only a `3/13` center gap is forced.  The residue pattern

```text
{0,3,5,8,10} mod 13
```

has cyclic gaps `(3,2,3,2,3)`, so its largest vertical safe gap at `t=1/13`
is only

```text
3/13-1/7=8/91<1/7.
```

Thus the universal thirteen-grid rectangle naturally reaches the
four-removal (`r=5`) case but cannot by itself prove THM-1132's five-removal
(`r=6`) sharp-horn target.  The phase transition is geometry, not a loose
constant in the proof.

## Kakeya and tournament reading

The proof carrier is a Kakeya needle through a rectangle in the fixed torus,
not a runner graph.  The four offset residues at `t=1/13` are the vertices of
a cyclic-gap tournament only after a cut is chosen; their useful observable
is the metric cyclic gap, and the tie Hamiltonian path is the residue order.
Forgetting the gap lengths leaves a transitive order and destroys the
`4/13` pigeonhole.  The faithful quotient is therefore the cyclic residue
word together with integer gap labels, thickened by the wall-slope sidecar.

This challenges two recurring assumptions at once: the vertices need not be
runners, and a wall-order tournament without metric labels is not the
underlying object.  Here the exact LRC predicate is preserved by the labelled
rectangle plus the slope needle; it is destroyed by order alone.

## Scope

This theorem does not close arbitrary `r=5`: condition (1) leaves a broad
intermediate-ratio region, and no finite bank is promoted across that gap.
It does give the first offset-shape-uniform toothpick self-similarity theorem,
turning THM-1127's single ray into a genuine all-shape cone.

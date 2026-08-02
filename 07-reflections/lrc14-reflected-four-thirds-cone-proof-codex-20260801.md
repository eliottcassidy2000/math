# LRC14 reflected levels: the four-thirds cone also closes

**Corrected proof candidate + frozen exact referee, 2026-08-01.**  The first
artifact at commit `a930fcb310` was retracted after audit found that its
primitive generator still hard-coded the inherited cap `5/3`.  The repaired
referee takes the cap as explicit data, enumerates the missing interval, and
re-establishes every census and located policy.  This note is not a canon
promotion.  It recursively strengthens the three-halves cone candidate and
proves, inside the sufficient reflected `k=1` family of
[THM-2941](../01-canon/theorems/THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary.md),

```text
D >= 6,                       3m >= 4D,
m=min q_e,                    D=max q_e-min q_e.
```

Consequently the reflected certificate-failure wedge is confined to the same
`561` bodies with `1<=m<4D/3`.  This is not a physical-survivor census and not
a proof of LRC(14).

## 1. What survived the wider cone

The load-bearing inherited mechanism is the exact edge-alphabet/projective-CSP
proof in
[`lrc14-reflected-three-halves-cone-csp-proof-codex-20260801.md`](lrc14-reflected-three-halves-cone-csp-proof-codex-20260801.md).
The projective diameter grows from `5/3` to `7/4`, and the minimum level at the
worst spread falls from nine to eight.  Nevertheless the primitive phase-zero
channels do not change: they remain `2:3` and `3:4`.  The new work is therefore
quantitative rather than a new zero-graph topology.

The exact transport bound for labels `a<b`, `delta=b-a`, is again

```text
B(D,m)=[m delta+D b]/[mL-b].                             (1)
```

It decreases in `m`.  On the real boundary `m=4D/3`, put
`A=(4/3)delta+b`; then

```text
B(D,4D/3)=DA/[(4D/3)L-b]
```

decreases in `D`, with exact consecutive difference

```text
Ab/[((4D/3)L-b)((4(D+1)/3)L-b)] > 0.                    (2)
```

Using the real boundary proves the integer condition `3m>=4D` even when
`D` is not divisible by three.  The unique worst corner is `D=6,m=8`.
Across all bodies and label pairs,

```text
|eta| <= 40/333,
```

and every moving slope remains at least `2624/333>1`.

## 2. Exact CSP census

At the corner, every label edge gets its own loss threshold

```text
t_E(a,b)
 =2[8(b-a)+6b]/[8L-b]
  +sum_(e in E)e/[7(8L-e)].                              (3)
```

When `t_E(a,b)<1/49`, a channel surviving that edge obeys the exact finite
bound

```text
PQ <= 12/[1-49t_E(a,b)].                                (4)
```

There are `7,529` constrained edges.  The largest product bound is

```text
190249826191276797660 / 11872893085202803 < 16025,
```

The first artifact incorrectly reported a bank of `2,492`: monkeypatching the
geometric cap did not change the inherited generator predicate `3Q<=5P`.
The repaired cap-aware generator has `2,728` unordered ratios in
`[4/7,7/4]` and explicitly checks that the reduced endpoint `(P,Q)=(4,7)` is
present.  The additional `236` channels change neither the CSP verdicts nor
the selected policy alphabets.

The same exhaustive component solver now uses projective cap `7/4`.  Reverse
vertex/candidate order gives identical verdicts on all `561` bodies:

```text
530 bodies have no CSP realization;
31 remain as coarse traps.                               (5)
```

Unlike the saturated components possible at the three-halves boundary, every
one of these `31` traps has a short realization.  The largest constructed
span is

```text
2021/1212 < 7/4.
```

Thus no new full-diameter projective obstruction appears.  The component
scale remains a free sidecar, handled by the same finite collision-avoidance
argument.

## 3. Twenty-eight finite oriented policies

Exactly `28` traps have a nonempty constrained graph.  Choosing the edge with
the smallest complete alphabet and retaining both channel orientations gives
`158` exact located controls.  Every control verifies, in this order,

```text
skeleton-2|eta_(g0)| > debt_8 > 0,
c=1-a/(Pg0L) in (0,1),
c^-1(skeleton-2|eta_(g0)|) > debt_8,
direct reflected overlap >= transported lower bound.
```

The scale law

```text
g/(PgL-a)-(g+1)/(P(g+1)L-a)
 =a/[(PgL-a)(P(g+1)L-a)] > 0                            (6)
```

makes each row uniform at every later common scale.  The weakest margin is

```text
52699395718360336279 / 2852106844532436934905
```

on body `(1,2,4,6,8,12)`, slots `(0,1)`, reverse channel `(5,3)`, cell
`290`.

## 4. Three unconstrained bodies and their tails

The only traps with no constrained edge are

```text
H =(1,2,3,4,6,12),
H2=(1,3,4,6,8,12),
H3=(2,3,4,6,8,12).
```

For an oriented primitive channel put `s=min(P,Q)`, `r=Q/P`.  On
`r in [4/7,7/4]`, the selected pairs have exact endpoint bounds

```text
H,  labels (1,2):  2|r-2|  <=20/7;
H2, labels (3,4):  2|3r-4|<=32/7;
H3, labels (2,3):  2|2r-3|<=26/7.                       (7)
```

Together with `PQ>=s^2`, `Pg>=s`, and `c>=-12/49`, these give increasing tail
envelopes.  Their first positive thresholds and complete finite heads are

| body | tail start | exact positive margin | oriented head |
|---|---:|---:|---:|
| `H` | `25` | `7048708348487489/561486400263755841875` | `270` |
| `H2` | `8` | `647017906068727/563281145404748848` | `28` |
| `H3` | `6` | `59758498142995026/88595783439538489345` | `16` |

Each envelope is nonpositive one step earlier and increases strictly
thereafter.  Each head is independently regenerated in a square twice as
large as the completeness bound.  All `314` oriented head channels have
strict located controls.  In each head the weakest orientation is `(7,4)`;
the cells are respectively `155`, `311`, and `323`.

This closes the last three bodies and proves the candidate theorem.

## 5. The fixed-cell lemma expands with the cone

Cell `155` of `H` still satisfies

```text
S_155(P,Q) >= 1/126                                      (8)
```

for every coprime oriented primitive ratio in the wider interval
`[4/7,7/4]`, with unique equality at `(3,2)`.  The reason is structural: for
`PQ>=20`, the phase bound alone proves `(8)`, while widening from
`[3/5,5/3]` introduces no new reduced channel with `PQ<20`.  The same six
small exact controls therefore suffice.

This explains why the worst `H` head channel `(7,4)` returns to cell `155`.
More generally, the successful recursive move is:

```text
widen projective cap
 -> recompute exact edge alphabets
 -> preserve orientation only on the small trap set
 -> send unconstrained bodies to analytic tails.
```

The first anticipated qualitative boundary is nearer `m=D`, where the new
zero channel `1:2` enters.  Below that point one should not assume the present
zero-graph anatomy survives.

## 6. Frozen referee

Run

```bash
python3 04-computation/lrc14_j7_reflected_four_thirds_cone_closure_thm2941.py
python3 -O 04-computation/lrc14_j7_reflected_four_thirds_cone_closure_thm2941.py
```

Ordinary and optimized transcripts are byte-identical.  The referee pins the
three-halves source/output/semantic image, checks the real-boundary
monotonicity and homotopy guard, regenerates every primitive channel and CSP
verdict in two search orders, and checks all `158+314` located controls against
direct reflected arcs.

```text
source:   0b8b01a184f5862c3a32598600988773bcbefdc26c63644cddf0884b18285c3e
output:   0377aaeb146553fe7c027a81ea3589045d01f50fe11aa79848d009a48c806fb9
semantic: d2c1b3bb27ad649db66a3b76305dd6027565cc020c31ddd4cd7e779b48d8bd53
located:  d74273080019120cc2891ebb34baa5b2472deeafa1a9b68faa71067ca3c554fb
heads:    c8fa54590e7104d3c790188484f77d96ae07961debfe196c4652cc4b83944b26
```

# LRC14 reflected levels: audit-required `3m >= 2D` cone candidate

> **AUDIT-REQUIRED / CLOSURE NOT CURRENT (MISTAKE-347).**  This refinement
> inherits the duplicated harmonic-body split-tail orientation from the
> one-cone referee.  The bridge and CSP censuses remain finite-exact, but the
> all-assignment cone conclusion is not proved.

**Proof candidate + frozen exact referee, 2026-08-01.**  Inside the sufficient
reflected `k=1` family of
[THM-2941](../01-canon/theorems/THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary.md),
the frozen referee claimed

```text
D >= 6,                       3m >= 2D,
m=min q_e,                    D=max q_e-min q_e.
```

If its missing orientation were repaired, this would combine with THM-2941's
arbitrary-level bank to confine reflected
certificate failure to the same `561` bodies and

```text
D >= 6,                       1 <= m < 2D/3.
```

This is a theorem about one sufficient reflected certificate family.  It is
not a physical-survivor census and not a proof of LRC(14).

## 1. Inheritance and correction lineage

The referee pins the source, output, and semantic image of the repaired
`4m>=3D` artifact in
[`lrc14-reflected-three-quarter-cone-proof-codex-20260801.md`](lrc14-reflected-three-quarter-cone-proof-codex-20260801.md).
That artifact pins the repaired cap-aware chain below it.  In particular, this
proof does not inherit the old four-thirds run whose primitive generator
silently retained a literal `5/3` cap.

The projective interval is `[2/5,5/2]`.  At its endpoint the phase-zero
alphabet changes from

```text
4/3, 3/2, 2
```

to

```text
4/3, 3/2, 2, 5/2.                                  (1)
```

This is the first boundary after cap `2` at which a zero gain adds a genuinely
new prime coordinate.

## 2. The exact integer corner is `(D,m)=(6,4)`

For labels `a<b`, put `delta=b-a`, `L=14*lcm(body)`, and

```text
B(D,m)=[m delta+D b]/[mL-b].                         (2)
```

On `3m>=2D`, maximize `D` for fixed `m`.  Write
`A=2delta+3b`.

- If `m=2k`, then `D=3k` and
  `B=kA/(2kL-b)`.  This decreases strictly with `k`; the first admissible
  point is `(D,m)=(6,4)`.
- If `m=2k+1`, then the boundary value lies below the common limit `A/(2L)`
  by
  `b(L-A)/[2L((2k+1)L-b)]`, which is positive.

The even corner `(6,4)` lies strictly above that limit, so it uniformly
maximizes the transport loss.  Across all bodies and slot pairs,

```text
|eta| <= 29/165,
moving slope >= 4-29/165 = 631/165 > 1.              (3)
```

The maximum is attained on `(1,2,3,4,6,12)`, labels `(1,12)`, with
`L=168`.  This residue-class argument is exact; it does not replace the
integer cone by a continuous boundary.

## 3. A new prime direction, but no new zero holonomy

Use prime-exponent coordinates `(v_2,v_3,v_5)`:

```text
4/3 -> ( 2,-1, 0)             2   -> ( 1, 0, 0)
3/2 -> (-1, 1, 0)             5/2 -> (-1, 0, 1).     (4)
```

The old three gains have rank/nullity `(2,1)`, with the balanced relation

```text
(3/2)(4/3)=2.                                      (5)
```

After adjoining `5/2`, rank/nullity is `(3,1)`: rank rises, while nullity does
not.  Thus the endpoint contributes scale freedom rather than a second
holonomy relation.  This distinction predicts the exact graph anatomy below.

The phase floor and exact edge thresholds give

```text
6,358 constrained body edges;
maximum product bound
  483806351847017447820 / 8426333817800417;
15,995 primitive unordered channels in [1,5/2].      (6)
```

The generator asserts that `(2,5)` is the reduced cap endpoint.  Both CSP
orders agree on all `561` residual bodies:

```text
389 bodies have no projective realization;
172 remain as coarse traps.                          (7)
```

The largest witness span is exactly `5/2`, attained first on
`(1,2,3,4,6,7)`.

## 4. Endpoint-bridge theorem for every forced-full component

Exactly `65` coarse traps contain a forced full-span component.  The referee
reconstructs its realized zero-gain subgraph rather than inspecting only its
diameter.  Uniformly over all `65` components:

- all four gains in (1) occur and their exponent rank is `3`;
- an old balanced triangle with gains `{4/3,3/2,2}` occurs;
- exactly one edge has gain `5/2`;
- deleting that edge increases the component count, so it is a bridge and is
  never a cycle edge.

The complete profile `(vertices, zero edges, cycle rank) -> count` is

```text
(4,4,1) ->  8       (5,4,1) -> 24       (5,5,1) -> 4
(6,4,1) -> 20       (6,5,1) ->  2       (6,7,2) -> 7. (8)
```

This is the structural reason the new endpoint does not create a new balanced
core: its independent `v_5` coordinate cannot participate in an integer
relation with the old alphabet, so it attaches the old holonomy core by a
tree edge.  The graph census is the finite exact shadow of the exponent-lattice
statement, not an unexplained numerical coincidence.

## 5. One hundred sixty finite policies

Twelve traps have no constrained edge.  On each of the remaining `160`, the
referee selects the edge with the smallest complete alphabet, retains both
orientations, and checks `1,566` exact located controls.  Every row proves

```text
skeleton-2|eta_(g0)| > debt_4 > 0,
c=1-a/(Pg0L) in (0,1),
c^-1(skeleton-2|eta_(g0)|) > debt_4,
direct reflected overlap >= transported lower bound. (9)
```

The factor `g/(PgL-a)` decreases strictly, making the check uniform at later
common scales.  The weakest margin is

```text
24399928897460195086201 / 1956983772380651959624680
```

on body `(1,4,6,8,9,12)`, pair `(1,4)`, channel `(13,8)`, cell `935`.

## 6. Twelve unconstrained bodies

The remaining bodies are handled by `14` oriented lanes; two bodies need a
low/high split.  For `s=min(P,Q)`, the phase floor, `PQ>=s^2`, and `Pg>=s`
give an increasing tail envelope.  Each declared threshold is nonpositive one
step earlier and strictly positive thereafter.  Exact endpoint maximization
supplies the numerator in every lane, rather than a sampled ratio estimate.

The tail starts and numbers of finite head controls are

| body/lane | ordered labels | numerator | tail start | heads |
|---|---:|---:|---:|---:|
| `(1,2,3,4,6,8)`, all | `(1,2)` | `16/5` | `6` | `30` |
| `(1,2,3,4,6,12)`, low | `(2,1)` | `2` | `10` | `42` |
| `(1,2,3,4,6,12)`, high | `(1,2)` | `2` | `10` | `42` |
| `(1,2,3,4,8,12)`, all | `(1,2)` | `16/5` | `6` | `30` |
| `(1,2,3,5,6,10)`, all | `(1,2)` | `16/5` | `5` | `18` |
| `(1,2,3,6,8,12)`, all | `(1,2)` | `16/5` | `6` | `30` |
| `(1,2,4,6,8,12)`, all | `(1,2)` | `16/5` | `6` | `30` |
| `(1,3,4,6,8,12)`, low | `(4,3)` | `14/5` | `6` | `15` |
| `(1,3,4,6,8,12)`, high | `(1,3)` | `4` | `8` | `27` |
| `(1,3,4,6,9,12)`, all | `(1,3)` | `26/5` | `6` | `30` |
| `(1,4,5,6,10,12)`, all | `(1,4)` | `36/5` | `5` | `18` |
| `(2,3,4,6,8,12)`, all | `(2,3)` | `22/5` | `9` | `66` |
| `(2,3,4,6,9,12)`, all | `(2,3)` | `22/5` | `6` | `30` |
| `(2,4,5,6,10,12)`, all | `(2,4)` | `32/5` | `5` | `18` |

The independent `3s` and `5s` head enumerations agree.  All `426` finite
heads pass exact located and direct-reflected checks.  The weakest is

```text
4785954156591456 / 1291245562866407995
```

on `(2,3,4,6,8,12)`, channel `(17,7)`, cell `323`.

## 7. Structural lesson and next boundary

The first three cap transitions now separate three phenomena:

```text
cap 2:    dependent zero gain enters; first balanced triangle; endpoint rigid;
cap 7/3:  no zero gain enters; the same triangle moves into the interior;
cap 5/2:  independent zero gain enters; exponent rank rises; endpoint bridge.
```

At the next boundary, cap `3`, the gain `3` enters with

```text
(3/2) * 2 = 3.                                      (10)
```

It therefore raises nullity rather than rank.  The sharp next experiment is
not merely whether the cap-`3` CSP closes, but whether every forced endpoint
component now contains a second balanced core and how the two circuits share
the `3/2` and `2` edges.

## 8. Frozen referee

Run

```bash
python3 04-computation/lrc14_j7_reflected_two_thirds_cone_closure_thm2941.py
python3 -O 04-computation/lrc14_j7_reflected_two_thirds_cone_closure_thm2941.py
```

The two modes are byte-identical.  The script pins the repaired predecessor,
proves the integer corner, regenerates the full channel bank and both CSP
orders, verifies the bridge anatomy, and checks all `1,566+426` located
controls against direct reflected arcs.

```text
source:    e6e64c909c6bfcc776eaf6bf2ad210f75675a6a32a46ead70d8d29a37f607eb3
output:    e77929c87f9d9f8fb7ce3a347c48522e87f63fa7c7085eb9e2cd8fe0bb4e4a90
semantic:  b34db3e25b9e4c81c1549f1ec7c7ab78e935ec778610daf617934b66b3a47304
traps:     b8e1f964d610af641a229d06951c4805167af4310429740501809034e9b2a716
forced:    37a6b66281b9882bd6ca278db3fe6c880d955eb95b2a3de27e65e41fdae56db2
located:   e0a6c267f4f77017ea56229f3f2c80ba1b1a9c00b065ad9ab0540397ce238a75
heads:     e6e32ad0c6beddfb3f3e53277dab82f76cb0fc0eae2f178698d95bb3e5d5e28f
```

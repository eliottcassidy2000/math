# LRC14 reflected levels: the `4m >= 3D` cone closes

**Proof candidate + frozen exact referee, 2026-08-01.**  Inside the sufficient
reflected `k=1` family of
[THM-2941](../01-canon/theorems/THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary.md),
the exact referee proves

```text
D >= 6,                       4m >= 3D,
m=min q_e,                    D=max q_e-min q_e.
```

Together with THM-2941's arbitrary-level bank, this confines reflected
certificate failure to the same `561` bodies and

```text
D >= 6,                       1 <= m < 3D/4.
```

This is not a physical-survivor census and not a proof of LRC(14).

## 1. Inheritance and correction lineage

The proof inherits the cap-complete `m>=D` referee in
[`lrc14-reflected-one-cone-gain-proof-codex-20260801.md`](lrc14-reflected-one-cone-gain-proof-codex-20260801.md),
which in turn pins the repaired four-thirds artifact.  This matters: the first
four-thirds run silently retained a literal `5/3` alphabet cap.  Here the
generator consumes `7/3` explicitly and asserts the reduced endpoint
`(P,Q)=(3,7)`.  No pre-repair census is a proof dependency.

The projective interval is `[3/7,7/3]`.  It widens past the endpoint gain `2`
but not as far as the next primitive phase-zero gain `5/2`.  Hence the exact
zero alphabet remains

```text
2, 3/2, 4/3,                 (3/2)(4/3)=2.             (1)
```

This gives a clean test of the balanced gain triangle: the topology persists,
but cap `7/3` gives it interior scale room.

## 2. The exact integer corner is not the lowest level

For labels `a<b`, set `delta=b-a`, `L=14*lcm(body)`, and

```text
B(D,m)=[m delta+D b]/[mL-b].                            (2)
```

Under `4m>=3D`, maximize `D` for fixed `m` and split `m` modulo three.  Put
`A=3delta+4b`.

- If `m=3k`, the boundary is `D=4k` and
  `B=kA/(3kL-b)`, strictly decreasing in `k`.  Its first admissible point is
  `(D,m)=(8,6)`.
- If `m=3k+1`, the value lies below `A/(3L)` by
  `b(L-7b+3a)/[3L((3k+1)L-b)]>0`.
- If `m=3k+2`, it lies below `A/(3L)` by
  `b(2L-A)/[3L((3k+2)L-b)]>0`.

The point `(8,6)` lies strictly above the limit and also above the two smaller
admissible corners `(6,5)` and `(9,7)`.  Thus it is the exact transport corner
for every body and pair.  Singleton debt, however, is largest at the
independent minimum level `5`.  Combining these two corners is a safe
over-bound and avoids falsely pretending one packet realizes both.

Across all bodies and pairs,

```text
|eta| <= 27/166,
moving slope >= 5-27/166 = 803/166 > 1.                (3)
```

This integer residue-class corner is sharper than using the inaccessible real
point `(D,m)=(6,9/2)`.

## 3. Cap-complete CSP and the gain triangle becoming interior

The exact edge thresholds and phase correction floor give

```text
6,576 constrained body edges;
maximum product bound
  688440220530612109197534 / 17995027579452067021
  < 38258;
9,857 primitive unordered channels in [1,7/3].         (4)
```

Both search orders agree on all `561` residual bodies:

```text
459 bodies have no projective realization;
102 remain as coarse traps.                             (5)
```

The largest constructed trap span is

```text
128/55 < 7/3.
```

There is no forced full-cap component.  Compare this with cap `2`, where
exactly `21` traps were forced full and every one contained the balanced zero
triangle.  Widening the cap has not introduced a new zero gain; it has simply
turned the old extremal triangle into an interior object.  This separates
zero-gain topology from diameter obstruction.

## 4. Ninety-three finite policies

Nine traps have no constrained edge.  On each of the other `93`, choosing the
edge with the smallest complete alphabet and retaining both orientations gives
`806` exact controls.  Every row proves

```text
skeleton-2|eta_(g0)| > debt_5 > 0,
c=1-a/(Pg0L) in (0,1),
c^-1(skeleton-2|eta_(g0)|) > debt_5,
direct reflected overlap >= transported lower bound.   (6)
```

The exact scale monotonicity of `g/(PgL-a)` makes each row uniform at all later
common scales.  The weakest margin is

```text
5872772186389138471201 / 537284849527258325396904
```

on body `(1,4,5,6,10,12)`, slots `(1,2)`, channel `(11,24)`, cell `779`.

## 5. Nine unconstrained bodies

For `s=min(P,Q)`, the phase floor, `PQ>=s^2`, and `Pg>=s` give the increasing
tail envelope

```text
1/49 - 12/(49s^2) - N/(L-a/s) - debt_5.                (7)
```

Most bodies admit one ordered pair on the whole interval.  Two need an
orientation split: `H=(1,2,3,4,6,12)` at `r=1`, and
`H2=(1,3,4,6,8,12)` because the best low-side approximation uses ordered
labels `(4,3)` while the high side uses `(1,3)`.

| body/lane | ordered labels | `N` | tail start | head channels |
|---|---:|---:|---:|---:|
| `(1,2,3,4,6,8)`, all | `(1,2)` | `22/7` | `6` | `26` |
| `H`, low | `(2,1)` | `2` | `9` | `29` |
| `H`, high | `(1,2)` | `2` | `9` | `29` |
| `(1,2,3,4,8,12)`, all | `(1,2)` | `22/7` | `6` | `26` |
| `(1,2,3,6,8,12)`, all | `(1,2)` | `22/7` | `6` | `26` |
| `(1,2,4,6,8,12)`, all | `(1,2)` | `22/7` | `6` | `26` |
| `H2`, low | `(4,3)` | `18/7` | `5` | `8` |
| `H2`, high | `(1,3)` | `4` | `7` | `16` |
| `(1,3,4,6,9,12)`, all | `(1,3)` | `36/7` | `6` | `26` |
| `(2,3,4,6,8,12)`, all | `(2,3)` | `30/7` | `8` | `48` |
| `(2,3,4,6,9,12)`, all | `(2,3)` | `30/7` | `5` | `16` |

Each tail is nonpositive one step earlier and strictly positive thereafter.
The complete head bound is a `3s` square because the ratio cap is `7/3`; the
referee independently regenerates the same list in a `5s` square.  All `276`
head channels pass exact located and direct-reflected checks.  The weakest
head margin is

```text
854592947125294523 / 183485517233014245249
```

on `(2,3,4,6,8,12)`, channel `(7,3)`, cell `323`.

## 6. Structural lesson and next boundary

Two distinct transitions are now visible:

```text
cap reaches 2:
  new zero gain -> balanced triangle -> forced full components;

cap exceeds 2 but stays below 5/2:
  same zero topology -> triangle becomes interior -> no forced full component.
```

Thus a phase-zero gain relation predicts topology, while its position relative
to the projective boundary predicts rigidity.  The next qualitative boundary
is cap `5/2`, where the zero gain `5/2` enters.  That is a better organizing
target than blindly lowering the minimum level one unit at a time.

## 7. Frozen referee

Run

```bash
python3 04-computation/lrc14_j7_reflected_three_quarter_cone_closure_thm2941.py
python3 -O 04-computation/lrc14_j7_reflected_three_quarter_cone_closure_thm2941.py
```

The script pins the cap-correct one-cone source/output/semantic image, proves
the integer residue-class corner, regenerates the `9,857`-channel bank and both
CSP orders, and checks all `806+276` located controls against direct reflected
arcs.

```text
source:   85bb9bd1613abd5cd7a877958b5c89a10172a5034c37cae5cc0467bd8ba4c0d3
output:   451823bd9bb11ae4af5c1fb4675fba2f163d172c178eae26c4cd115eb945ba7c
semantic: 05a900c80283d5bc9a7b01c1c4ad045b889aebe2e4d5798eaa66cb8907a7ae9f
traps:    cff2ad5603eb88999cc9a5bebb8cdc4a6e3a2cb29d52fab433a935ed4453a971
```

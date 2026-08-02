# LRC14 reflected levels: audit-required full `m >= D` cone candidate

> **AUDIT-REQUIRED / CLOSURE NOT CURRENT (MISTAKE-347).**  The low/high
> harmonic-body tails reverse both the ratio interval and ordered label pair,
> so they cover `q_0<=q_1` twice and omit `q_0>q_1`.  The CSP and located
> controls remain finite-exact evidence; the all-assignment cone claim does
> not.

**Proof candidate + frozen exact referee, 2026-08-01.**  Inside the sufficient
reflected `k=1` family of
[THM-2941](../01-canon/theorems/THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary.md),
the frozen referee claimed

```text
D >= 6,                       m >= D,
m=min q_e,                    D=max q_e-min q_e.
```

If its missing orientation were repaired, this would combine with THM-2941's
arbitrary-level bank to confine reflected
certificate failure to the same `561` bodies and the strict wedge

```text
D >= 6,                       1 <= m < D.
```

This is not a physical-survivor census and not a proof of LRC(14).

## 1. Inheritance, correction, and the qualitative boundary

The closest proved mechanism is the corrected four-thirds cone candidate in
[`lrc14-reflected-four-thirds-cone-proof-codex-20260801.md`](lrc14-reflected-four-thirds-cone-proof-codex-20260801.md).
Its first artifact was retracted because changing the imported geometric cap
did not change a literal `3Q<=5P` inside the primitive generator.  The repaired
generator takes the cap explicitly and checks its rational endpoint.  The
present proof inherits only that repaired version and verifies that `(P,Q)=(1,2)`
is actually present.

At `m=D`, the projective interval becomes `[1/2,2]`.  Its exact primitive
phase-zero channels are

```text
2, 3/2, 4/3,
```

not merely the old pair `3/2,4/3`.  They form the balanced multiplicative
gain circuit

```text
(3/2)(4/3)=2.                                         (1)
```

This is the qualitative boundary anticipated by the four-thirds proof.  It
is best viewed as a rational gain graph.  Give an oriented edge the gain by
which its lower potential is multiplied to reach its higher potential.
Cycle realizability requires product one after inserting inverse gains.  The
triangle `(3/2,4/3,2)` is therefore forced to be transitive in potential
order: the two-step path and long edge agree by `(1)`.  This is an intrinsic
tournament only after retaining the potential order; forgetting orientation
loses exactly the sidecar later needed by the hardest tail.

## 2. One exact corner

For body labels `a<b`, put `delta=b-a` and `L=14*lcm(body)`.  The inherited
transport bound is

```text
B(D,m)=[m delta+D b]/[mL-b].                           (2)
```

It decreases strictly in `m`.  On the real boundary `m=D`, put
`A=delta+b`; then

```text
B(D,D)=DA/(DL-b),
B(D,D)-B(D+1,D+1)=Ab/[(DL-b)((D+1)L-b)]>0.            (3)
```

Thus the entire infinite cone has the exact worst corner `D=m=6`.  Across
all `3,003` bodies and all label pairs,

```text
|eta| <= 23/166,
moving slope >= 6-23/166 = 973/166 > 1.                (4)
```

The argument uses the real boundary, so no divisibility or rounding fiction
is hidden in the reduction.

## 3. Cap-complete projective CSP

At the corner every label edge gets its exact coarse loss threshold.  If the
threshold is below `1/49`, the correction floor `c>=-12/49` forces

```text
PQ <= 12/(1-49 threshold).                             (5)
```

The cap-aware census has

```text
7,161 constrained body edges;
maximum product bound
  6879371187357084 / 5210641623133 < 1321;
279 primitive unordered channels in [1,2].             (6)
```

The endpoint `(1,2)` and exactly the three zero channels above are asserted
inside the referee.  Reverse vertex/candidate order agrees on every body:

```text
492 bodies have no projective realization;
69 remain as coarse traps.                              (7)
```

These are over-approximating projective traps, not physical survivors.
In particular, the earlier preliminary phrase "six body-word survivors" is
not a dependency of this proof: after cap repair the safe projective census
has `69` traps, and only the complete located/tail analysis below closes them.

The zero circuit is visible in the exact trap anatomy.  Twenty-one traps
have a component forced to use the full projective span `2`.  Every one of
those `21` components contains a constrained triangle whose three realized
gains are exactly `{3/2,4/3,2}`.  Nineteen have one such core and two have two.
This explains why the new endpoint changes the CSP topology even though it
does not obstruct the final located proof: full-span behavior is organized
by a balanced gain circuit, not by an isolated bad edge.

## 4. Sixty-five complete located policies

Four traps have no constrained edge.  For each of the other `65`, choose the
edge with the smallest complete phase alphabet and retain both orientations
of every channel.  This gives `894` exact controls.  Every row checks

```text
skeleton-2|eta_(g0)| > debt_6 > 0,
c=1-a/(Pg0L) in (0,1),
c^-1(skeleton-2|eta_(g0)|) > debt_6,
direct reflected overlap >= transported lower bound.   (8)
```

The exact scale difference

```text
g/(PgL-a)-(g+1)/(P(g+1)L-a)
 =a/[(PgL-a)(P(g+1)L-a)]>0                             (9)
```

makes every control uniform at later common scales.  The weakest margin is

```text
79338784432498793 / 8505769192234754609
```

on body `(1,2,4,6,8,12)`, slots `(0,1)`, channel `(37,19)`, cell `303`.

## 5. Four unconstrained bodies and the orientation sidecar

The unconstrained bodies are

```text
H =(1,2,3,4,6,12),
H2=(1,3,4,6,8,12),
H4=(1,3,4,6,9,12),
H3=(2,3,4,6,8,12).
```

A single ordered pair over the whole interval fails asymptotically for `H`:
with labels `(1,2)`, the bound on `2|r-2|` is `3`.  The missing coordinate is
the orientation of `r=Q/P`.  Split at one:

```text
r in [1,2]:   ordered labels (1,2),  2|r-2| <= 2;
r in [1/2,1]: ordered labels (2,1),  2|2r-1| <= 2.      (10)
```

This is not cosmetic.  At the new reverse zero channel `(P,Q)=(2,1)`, the
old label order has skeleton zero at the formerly uniform cell `155`.
Reversing the labels restores skeleton `1/14` and makes `eta=0` at that cell.
Thus the gain-edge orientation is the precise sidecar lost by the unordered
CSP.

For `s=min(P,Q)`, the primitive phase bound, `PQ>=s^2`, and `Pg>=s` give the
strictly increasing tail envelope

```text
1/49 - 12/(49s^2) - N/(L-a/s) - debt_6.                (11)
```

The exact policies are:

| body/lane | ordered labels | `N` | first positive `s` | head channels |
|---|---:|---:|---:|---:|
| `H`, `r<1` | `(2,1)` | `2` | `8` | `18` |
| `H`, `r>1` | `(1,2)` | `2` | `8` | `18` |
| `H2`, all | `(1,3)` | `5` | `9` | `44` |
| `H4`, all | `(1,3)` | `5` | `6` | `20` |
| `H3`, all | `(2,3)` | `4` | `7` | `24` |

Each envelope is nonpositive one step earlier and increases strictly
thereafter.  Each finite head is regenerated in both a `2s` and `4s` square;
the channel lists agree.  All `124` head channels pass exact located and
direct-reflected checks.  The weakest head margin is

```text
24774659391587927 / 7430997812962436713
```

on `H2`, channel `(13,7)`, cell `311`.

## 6. Reusable structural lesson

The recursive cone method now has a sharp transition rule:

```text
widen the projective interval
 -> make the alphabet generator consume the cap explicitly
 -> inspect new phase-zero gains as a balanced gain graph
 -> retain edge orientation when passing from CSP to transport
 -> split analytic tails at orientation walls.          (12)
```

At cap `2`, the endpoint gain does two things at once: it creates forced
full-span components through the balanced triangle, and it breaks the old
fixed-order tail.  The same orientation that makes the gain triangle a
transitive tournament repairs the analytic tail.  This connects the finite
CSP and infinite transport layers rather than treating them as separate
computations.

The next reflected frontier is the strict wedge `m<D`, where the projective
cap exceeds two.  The gain alphabet may acquire new exact relations; those
relations, not raw trap counts, should be the first object inspected.

## 7. Frozen referee

Run

```bash
python3 04-computation/lrc14_j7_reflected_one_cone_closure_thm2941.py
python3 -O 04-computation/lrc14_j7_reflected_one_cone_closure_thm2941.py
```

The referee pins the corrected four-thirds source/output/semantic image,
checks the cap endpoint and gain triangle, regenerates both CSP search orders,
audits all `894+124` located rows against direct reflected arcs, and verifies
every tail threshold and oversized head enumeration.

```text
source:   800395ae242860094fed3db9638a93ebc2faba7973558a3eaa51f3af62145200
output:   e63ee74f42d8f213196cf2907bbc25dc20b93a3e64b8a40aa85d5d81c8b11ee6
semantic: dc84fbea8c3e951aa510f02776e5ad300e7d50843da6137fd967f14faff7d2d9
located:  4ffe8c58080e02517563c980ca141ad2d5956398696aef9a3a425d1fcec9044e
heads:    dd73914c6b2faf4d8c635e36f1fb38df0db0286dfcefe32191a78f0aa15fb0a4
gain:     e5e03bd0c7bb4fa9cc42ff7531fcfe5c7b5a943d0ddc7814ab9a62951512853a
```

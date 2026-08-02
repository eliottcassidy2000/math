# LRC14 reflected levels: audit-required `2m >= D` candidate and finite `K4` census

> **AUDIT-REQUIRED / CLOSURE NOT CURRENT (MISTAKE-347).**  The exact referee
> verifies its internal CSP and gain-`K4` data, but its low/high analytic lanes
> do not prove an assignment-level orientation cover.  The half-cone and
> `m<D/2` wedge claims below are quarantined pending repair.

**Proof candidate + frozen exact referee, 2026-08-01.**  Inside the sufficient
reflected `k=1` family of
[THM-2941](../01-canon/theorems/THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary.md),
the frozen referee claimed

```text
D >= 6,                       2m >= D,
m=min q_e,                    D=max q_e-min q_e.
```

If its missing orientation were repaired, this would combine with THM-2941's
arbitrary-level bank to confine reflected
certificate failure to the same `561` bodies and

```text
D >= 6,                       1 <= m < D/2.
```

This is a theorem about one sufficient reflected certificate family.  It is
not a physical-survivor census and not a proof of LRC(14).

## 1. Inheritance and exact corner

The referee pins the source, output, and semantic image of the cap-`5/2`
artifact in
[`lrc14-reflected-two-thirds-cone-proof-codex-20260801.md`](lrc14-reflected-two-thirds-cone-proof-codex-20260801.md).
That artifact pins the repaired cap-aware chain below it, so no pre-repair
four-thirds census is a dependency.

For labels `a<b`, put `delta=b-a`, `L=14*lcm(body)`, and

```text
B(D,m)=[m delta+D b]/[mL-b].                          (1)
```

On `2m>=D`, the largest `D` at fixed `m` is `2m`.  If
`A=delta+2b`, the boundary value is

```text
mA/(mL-b),
```

and the difference from scale `m+1` is

```text
Ab / [(mL-b)((m+1)L-b)] > 0.                         (2)
```

Since `D>=6` forces `m>=3`, the exact integer corner is `(D,m)=(6,3)`.
Across every body and ordered slot pair,

```text
|eta| <= 35/164,
moving slope >= 3-35/164 = 457/164 > 1.              (3)
```

The maximum is attained on `(1,2,3,4,6,12)`, labels `(1,12)`, with
`L=168`.

## 2. The dependent-gain wall

The projective interval is `[1/3,3]`.  At its endpoint the phase-zero gains
become

```text
4/3, 3/2, 2, 5/2, 3.                                (4)
```

In prime-exponent coordinates `(v_2,v_3,v_5)`, the old four gains have
rank/nullity `(3,1)`.  The new vector satisfies

```text
v(3) = v(3/2)+v(2),          (3/2)*2=3,              (5)
```

so the full alphabet has rank/nullity `(3,2)`.  This is the opposite of the
cap-`5/2` transition: gain `5/2` raised rank and became a bridge; gain `3`
raises relation nullity and becomes a cycle edge.

The exact phase floor and edge thresholds give

```text
5,666 constrained body edges;
maximum product bound
  597757150909236102135 / 22679209066623812;
8,796 primitive unordered channels in [1,3].          (6)
```

Both CSP search orders agree on all `561` residual bodies:

```text
282 bodies have no projective realization;
279 remain as coarse traps.                           (7)
```

The largest witness span is exactly `3`.

## 3. The unique doubled-core `K4`

Exactly `55` traps contain a forced full-span component.  Every one contains
a unique complete zero-gain graph on relative scales

```text
{1, 3/2, 2, 3} = {1,x,y,xy},     x=3/2, y=2.         (8)
```

This is a multiplicative Boolean square indexed by the four exponent bits in
`F_2^2`.  The bit labels form a `V4` vertex torsor, but their map to positive
rational scales is **not** a group homomorphism: toggling a bit is an incidence
operation, not multiplication by an involution.  Ordering the scales turns all
six pairs into the transitive-by-scale `K4`.  The six gains are exactly

```text
x/1=3/2,       y/1=2,          xy/1=3,
y/x=4/3,       xy/x=2,         xy/y=3/2.              (9)
```

Thus the gain multiset is

```text
{4/3, (3/2)^2, 2^2, 3}.                              (10)
```

The four triangular faces split without choice:

- `{1,x,y}` and `{x,y,xy}` are old `{4/3,3/2,2}` triangles;
- `{1,x,xy}` and `{1,y,xy}` are new `{3/2,2,3}` triangles.

This Boolean-square language represents the exact scale census.  It is not an
asserted symmetry of the six body labels and does not quotient away their
locations.

The three perfect matchings expose the precise failure of a gain-preserving
`V4` quotient.  The `x`-direction edges both have gain `x=3/2`, and the
`y`-direction edges both have gain `y=2`.  On the diagonal direction, however,
the two gains are `xy=3` and `y/x=4/3`, differing by `x^2=9/4`.  Thus the
unweighted `K4` incidence matches the tetrahedral point-direction carrier of
THM-3067, while its gain decoration remembers an origin on the diagonal
matching.  This is the cheapest explicit obstruction to importing the
`A4/V4=C3` quotient into the reflected certificate: forgetting the point
forgets a load-bearing gain.

Uniformly over all `55` forced components, the gain-`3` edge occurs exactly
once, lies in the unique `K4`, and is a cycle edge.  In `53` components the
two old and two new faces are the complete triangle inventory of these two
types.  The remaining `2` contain one additional old triangle but still
exactly two new ones.  They are the bodies `(1,6,7,9,12,14)` and
`(2,6,7,9,12,14)`, both with scale witness `(40,20,60,30,45,36)`; the extra
scale `45` closes an old triangle on scales `(60,30,45)` and raises cycle
rank from `3` to `4`.  The full profile
`(vertices, zero edges, cycle rank) -> count` is

```text
(4,6,3) ->  3       (5,6,3) ->  1       (5,7,3) -> 30
(6,7,3) ->  7       (6,8,3) -> 12       (6,8,4) ->  2. (11)
```

Forty-nine components also use gain `5/2`; six do not.  This sidecar does not
alter the unique doubled core.

## 4. Every constrained trap closes

Twenty-six traps have no constrained edge.  On each of the other `253`, the
referee chooses the edge with the smallest complete alphabet and retains both
orientations.  All `3,062` exact controls prove

```text
skeleton-2|eta_(g0)| > debt_3 > 0,
c=1-a/(Pg0L) in (0,1),
c^-1(skeleton-2|eta_(g0)|) > debt_3,
direct reflected overlap >= transported lower bound. (12)
```

The factor `g/(PgL-a)` decreases strictly, so each finite check controls every
later common scale.  The weakest margin is

```text
4389035249640369833 / 461756458190054127555
```

on `(3,4,6,8,9,12)`, pair `(3,4)`, channel `(9,26)`, cell `983`.

## 5. A local Farey atlas for the unconstrained bodies

For a ratio lane `[u,v]` and ordered labels `(a,b)`, exact endpoint convexity
gives

```text
N = 2 max(|ua-b|,|va-b|).                            (13)
```

For `s=min(P,Q)`, the phase floor, `PQ>=s^2`, and `Pg>=s` then give the
increasing envelope

```text
1/49 - 12/(49s^2) - N/(L-a/s) - debt_3.              (14)
```

Twenty-five unconstrained bodies close with the two universal charts
`[1/3,1]` and `[1,3]`, choosing the best ordered pair exactly.  The harmonic
body `H=(1,2,3,4,6,12)` is the sole bottleneck: a coarse two-chart cover would
leave `688` heads.  The following seven-chart Farey atlas reduces this to
`76` while retaining exact endpoint bounds.

| interval | ordered labels | `N` | tail start | heads |
|---|---:|---:|---:|---:|
| `[1/3,8/11]` | `(2,1)` | `10/11` | `6` | `17` |
| `[8/11,11/12]` | `(4,3)` | `4/3` | `8` | `6` |
| `[11/12,1]` | `(2,1)` | `2` | `24` | `14` |
| `[1,13/12]` | `(1,2)` | `2` | `24` | `12` |
| `[13/12,17/11]` | `(3,4)` | `3/2` | `9` | `11` |
| `[17/11,29/12]` | `(1,2)` | `10/11` | `6` | `9` |
| `[29/12,3]` | `(1,3)` | `7/6` | `7` | `7` |

Altogether there are `57` lanes and `792` finite head channels.  Each lane is
nonpositive one step before its declared start and strictly positive from its
start onward.  Independent `3s` and `5s` head boxes agree, and every head
passes the exact located and direct-reflected checks.  The weakest head margin
is

```text
1565713700195833 / 1549999933309989875
```

on `H`, pair `(1,2)`, channel `(21,22)`, cell `155`.

## 6. What the cap sequence is measuring

The successive zero walls now distinguish position, rank, and circuit
creation:

```text
cap 2:    dependent zero gain enters; balanced triangle appears at boundary;
cap 7/3:  no zero gain enters; that triangle moves into the interior;
cap 5/2:  independent zero gain enters; rank rises; unique endpoint bridge;
cap 3:    dependent zero gain enters; nullity rises; unique doubled-core K4.
```

This suggests organizing future cones by the matroid of prime-exponent gain
vectors together with the boundary position of each gain.  Neither datum alone
predicts the CSP topology: exponent dependence predicts a circuit, while the
projective endpoint decides whether that circuit is forced rigidly.

The next phase-zero wall is cap `4`, corresponding to the cone `3m>=D` and
minimum level `2`.  There gain `4=2^2` raises nullity again.  Its primitive
three-scale circuit `{1,2,4}` has repeated gain `2`; testing how that circuit
interacts with the now-interior Boolean square is the sharp next structural
experiment.

## 7. Frozen referee

Run

```bash
python3 04-computation/lrc14_j7_reflected_half_cone_closure_thm2941.py
python3 -O 04-computation/lrc14_j7_reflected_half_cone_closure_thm2941.py
```

The two modes are byte-identical.  The script pins the cap-`5/2` predecessor,
proves the exact integer corner, regenerates the complete primitive bank and
both CSP orders, verifies the unique Boolean-square core, and checks all
`3,062+792` located controls against direct reflected arcs.

```text
source:    01b0ef11060585cfdd8dcf001cf0488ffde69118cb4da1bbbd964b195e3cd910
output:    afe07bc18ae2fa257cca7eda51999e064411bc168a409c9d932b951e87a4d7e3
semantic:  e7edbeec8933409bb36eab7d4a1867a10b3f7683a4898ec30f244a9e81601673
traps:     35fceb70d1e960ee5e0590c990adde5c4848b657a15455714585354bee322380
verdicts:  fe8193713f596b6f3f0182066a8937cfe75f681ea760fd31bac4642b7eb0f7cc
forced:    70eab59654ebeee0fbcde8b196a2f06d10c13369e8a20ea7fadcd5a0f566eaad
K4:        2a90a1d9ff6b3dcdada4189eaef83d5534e8f53d2f71461dae3cb917d938d8fc
located:   8968666298108aabbf926de452354ed554728aa2fe2cc2ed5f5aaea3598f4201
lanes:     fa98d31b75ce2682d09e097cf4bcfcca4b00e0c01cac1b5a4b7f10b7676a857b
heads:     e89a7f20d72228870e488094333cf9f84d45f9f455d1bd5448a03eccca14642a
```

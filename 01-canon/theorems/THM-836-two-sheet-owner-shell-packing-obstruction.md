---
id: THM-836
title: Thin two-sheet owner shells force a modular packing obstruction
status: PROVED (directional owner-exit law, modular packing, empty s=1,3 shells, exact local s=5 classification, and the uniform d=11 mod 52 divisor exclusion) + VERIFIED (fraction-exact endpoint/shell replay and finite d=15 endpoint-grid census)
source: codex-2026-07-15-S10 continuation
depends_on: [THM-772, THM-797, THM-803, THM-824]
related: [THM-774, THM-803, HYP-6820]
verification:
  - 04-computation/lrc14_two_sheet_owner_shell_packing_obstruction_codex_S10.py
  - 05-knowledge/results/lrc14_two_sheet_owner_shell_packing_obstruction_codex_S10.out
  - 04-computation/lrc14_shell_five_local_classification_and_divisor_witness_codex_S10.py
  - 05-knowledge/results/lrc14_shell_five_local_classification_and_divisor_witness_codex_S10.out
---

# THM-836 — thin two-sheet owner-shell packing obstruction

This theorem isolates the exact first local obstruction behind THM-824's
fixed-ratio radius budget.  Its carrier is not the scalar radius alone: it is
the pair of signed endpoint owners that make the two sides of the forced
`q=13` component exit.

## 1. Statement

Let `d` be a positive odd integer not divisible by `13`, let `U` be ten
distinct positive integers, and put

```text
B=max(U),
E_U={t:min_(u in U)||ut||>=1/11},
R_U={r:max_(u in U)||ur||<2/143},
H_d={t:||9dt||+||4dt||>=11/13}.                            (1)
```

Assume the two conclusions needed from the hypothetical `(13d,5d)`
two-sheet packet:

```text
{u mod 13:u in U}=(Z/13Z)^* minus {+5d,-5d},              (2)
E_U+closure(R_U) subset H_d.                               (3)
```

Condition (2) is exactly THM-797's signed-wall conclusion.  Condition (3) is
the centre-guarded erosion containment to which THM-824 applies; the phase
inequality without centre membership is not being used as a substitute.

Choose `p` by

```text
5dp=1 mod 13,                                              (4)
```

and let `u_z` be the unique speed satisfying

```text
[p u_z]_13=z,             z in {+/-2,+/-3,+/-4,+/-5,+/-6}, (5)
```

where the bracket is the balanced residue.  For a direction
`epsilon in {+1,-1}`, define

```text
c_epsilon(z)=11|z|-13       if epsilon*z<0,
                 130-11|z| if epsilon*z>0.                (6)
```

Then the following exact directional inequalities hold:

```text
min_z c_epsilon(z)/u_z + 2/B <=22/(13d),
epsilon=+1,-1.                                             (7)
```

Each direction has the coefficient multiset

```text
(9,20,31,42,53,64,75,86,97,108),                          (8)
```

and `c_+(z)+c_-(z)=117` for every owner.

Since every coefficient is at least `9` and every speed is at most `B`, (7)
itself gives

```text
11/B<=22/(13d),          13d<=2B-1.                       (9)
```

The last unit of strictness is parity: `13d` is odd and `2B` is even.
Suppose now that `B<13d` and put

```text
s=2B-13d,
delta(d)=min(r,13-r),       r=6d mod 13 in {1,...,12}.    (10)
```

The two coefficient-`9` owners are forced into the interval

```text
u_(+2),u_(-2) >= L(d,B)=117dB/(22B-26d),                 (11)
```

their original residue classes are `+/-3d mod 13`, and hence

```text
delta(d)
 <= B-L(d,B)
 =11s(13d+s)/(2(117d+11s)).                              (12)
```

Thus `s` is already a positive odd integer under (2)--(3).  If `B>=13d`,
then `s>=13d>=13` and the final bound below is immediate.  In the remaining
thin band, formula (12), followed by the exact alignment of the owner classes
with `B mod 13`, eliminates both first shells:

```text
s!=1,                  s!=3,
13d<=2B-5.                                                (13)
```

This is a local necessary condition for the fixed common-ratio branch.  It
does not prove (3) fails once `s>=5`, does not cover other odd ratios, and
does not close the global two-sheet or `n=12` sporadic branch.

## 2. The signed wall lands at the guarded target centre

Multiplication by `p` sends the ten residues in (2) to every nonzero residue
except `+/-1`.  In balanced notation these are exactly the ten labels in
(5), each once.  In particular, at

```text
t_0=p/13,                                                  (14)
```

one has

```text
||u_z t_0||=|z|/13>=2/13>1/11,                            (15)
```

so `t_0` lies in the interior of `E_U`.

The independent centre guard is exact here.  From (4),

```text
dp=5^(-1)=8 mod 13,
dt_0=8/13 mod 1.                                          (16)
```

THM-824 identifies the component of `H` containing `8/13` as the closed
ball of radius `2/169` about that centre.  Therefore the component of `H_d`
containing `t_0` is, in the ordinary real lift about `t_0`,

```text
[t_0-2/(169d),t_0+2/(169d)].                              (17)
```

This is why the unguarded identity involving only `||13dt||` would not
suffice: (16), not a small `13d`-phase by itself, selects the correct folded
cell.

## 3. Exact directional exit distances

Fix a label `z` and move from `t_0` in direction `epsilon`.  If
`epsilon*z<0`, the signed phase moves toward zero.  It first reaches the
closed `1/11` boundary after the time

```text
(|z|/13-1/11)/u_z=(11|z|-13)/(143u_z).                    (18)
```

If `epsilon*z>0`, it first moves away from zero, crosses the half-integer,
and reaches the `1/11` boundary beside the next integer after

```text
(1-1/11-|z|/13)/u_z=(130-11|z|)/(143u_z).                 (19)
```

Both numerators are positive for `2<=|z|<=6`.  Before the displayed time
the corresponding clearance is strictly greater than `1/11`; at the
displayed time it equals `1/11`; immediately afterwards it is strictly less.
Thus no open/closed endpoint is lost.  Intersecting the ten closed safe
intervals proves that the component of `E_U` containing `t_0` is exactly

```text
[t_0-r_-,t_0+r_+],
r_epsilon=(1/143) min_z c_epsilon(z)/u_z.                 (20)
```

Each direction has a coefficient-`9` owner, so
`r_epsilon<=9/143<1/13<=min(t_0,1-t_0)`.  Hence (20) is an
unambiguous ordinary real lift and does not silently cross the circle cut.

For increasing time, the coefficient table is

| `z` | `-6` | `-5` | `-4` | `-3` | `-2` | `2` | `3` | `4` | `5` | `6` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `c_+(z)` | 53 | 42 | 31 | 20 | 9 | 108 | 97 | 86 | 75 | 64 |
| `c_-(z)` | 64 | 75 | 86 | 97 | 108 | 9 | 20 | 31 | 42 | 53 |

This proves (8) and the complementary identity `c_+(z)+c_-(z)=117`.

## 4. Adding the closed central return cell

For `|r|<2/(143B)` and every `u<=B`,

```text
||ur||<=u|r|<2/143.
```

Passing to the closure, including the two endpoints owned by the maximum
speed, gives

```text
[-2/(143B),2/(143B)] subset closure(R_U).                 (21)
```

By (3), the Minkowski sum of (20) and (21) lies in `H_d`.  Its right and
left paths start at the centre `t_0` and are connected, so each must remain
inside the single target component (17); they cannot switch across a gap.
Consequently, for each direction,

```text
r_epsilon+2/(143B)<=2/(169d).                             (22)
```

All three sets in this comparison are closed.  Equality at a deep endpoint,
a return endpoint, or a target endpoint is allowed, which is why (22) and
the resulting theorem use non-strict inequalities.  Substituting (20) into
(22) and multiplying by `143` proves (7).

## 5. The thin band forces the two `9`-owners

Assume `B<13d`.  The right side available to the minimum in (7) is

```text
A(d,B)=22/(13d)-2/B.                                      (23)
```

Inequality (9) makes `A(d,B)>0`.  Every coefficient other than `9`
is at least `20`, and

```text
20/B>A(d,B)
iff 22/B>22/(13d)
iff B<13d.                                                (24)
```

Thus no coefficient `20` or larger can satisfy either directional
obligation.  The increasing direction is forced to use `z=-2`, and the
decreasing direction is forced to use `z=+2`.  Applying (7) to both gives

```text
9/u_(+/-2)+2/B<=22/(13d),
u_(+/-2)>=117dB/(22B-26d)=L(d,B),                         (25)
```

where the denominator is positive.

Since `p^(-1)=5d mod 13`, the original residues of these two owners are

```text
u_(+2)=10d=-3d mod 13,
u_(-2)=-10d=+3d mod 13.                                  (26)
```

Their difference is therefore congruent to `+/-6d mod 13`, so its positive
integer magnitude is at least `delta(d)`.  Both owners lie in the real
interval `[L(d,B),B]`, hence

```text
delta(d)<=|u_(+2)-u_(-2)|<=B-L(d,B).                      (27)
```

Writing `B=(13d+s)/2` gives

```text
22B-26d=117d+11s,
B-L(d,B)=11s(13d+s)/(2(117d+11s)),                        (28)
```

which proves (12).

## 6. The first two shells are empty

When `s=1`, the right side of (12) is

```text
11(13d+1)/(2(117d+11))<1.                                (29)
```

But `delta(d)>=1`, a contradiction.

When `s=3`, put

```text
G=33(13d+3)/(2(117d+33)).                                 (30)
```

Direct cross-multiplication gives

```text
1<G<2.                                                    (31)
```

Thus (12) first forces `delta(d)=1`, equivalently

```text
6d=+/-1 mod 13,
d=2 or 11 mod 13.                                        (32)
```

The two owner residues in (26) are then exactly `{6,7}`.  On the other hand,

```text
B=(13d+3)/2=8 mod 13.                                    (33)
```

The largest integers at most `B` in owner classes `7` and `6` are therefore
`B-1` and `B-2`.  But (31) says

```text
L=B-G>B-2,                                                (34)
```

so no class-`6` owner can satisfy (25).  This contradiction empties `s=3`.
The only positive odd values below `5` were `1` and `3`; together with the
trivial `B>=13d` case noted above, (13) follows.

The endpoint-owner method genuinely stops at the next shell.  For example,

```text
d=11, B=74, s=5, L=4329/61, G=185/61,
u_(+/-2) in {71,72}.                                      (35)
```

The two local inequalities have positive slacks

```text
2/13-(9/71+2/74)=2/34151,
2/13-(9/72+2/74)=7/3848.                                 (36)
```

This is only a local owner-shell survivor, not a ten-speed core satisfying
(2)--(3).  It proves that this argument alone cannot replace `-5` in (13) by
a stronger constant.

## 6A. Shell-five addendum: exact local classes and one uniform divisor exclusion

The local shell-five survivors in the preceding paragraph admit an all-size
classification.  Suppose

```text
s=2B-13d=5,             B=(13d+5)/2.                    (S5.1)
```

Then `B=9 mod 13`, and the available packing width in (28) is

```text
G_5(d)=B-L(d,B)=55(13d+5)/(2(117d+55))<4.              (S5.2)
```

The two coefficient-nine owners therefore have integer offsets at most three
below `B`; their residue classes must both lie in

```text
{B,B-1,B-2,B-3} mod 13={9,8,7,6}.                       (S5.3)
```

But their signed residue pair is `{+3d,-3d}` by (26), and the only signed
pair wholly contained in (S5.3) is `{6,7}`.  Consequently

```text
d=+/-2 mod 13,
d mod 52 in {11,15,37,41},
{u_(+2),u_(-2)}={B-3,B-2},
{+5d,-5d} mod 13={3,10}.                                (S5.4)
```

Conversely, every positive odd `d` in the four classes in (S5.4) has
`d>=11` and

```text
G_5(d)>3,                                                 (S5.5)
```

so the two speeds `B-3,B-2` satisfy both local coefficient-nine owner
inequalities.  Thus (S5.4) is an exact classification of the **local owner
shell**.  The converse asserts only those two inequalities, not existence of
a ten-speed set satisfying (2)--(3).

One of the four congruence classes is nevertheless uniformly incompatible
with the guarded containment.  Let

```text
d=11 mod 52,             q=5d,
p=(45d-1)/26,            t=p/q.                           (S5.6)
```

The numerator is integral, `26p=45d-1` gives `gcd(p,d)=1`, and
`p=4 mod 5`; hence `p` is a unit modulo `q`.  Moreover

```text
t=9/26-1/(130d).                                         (S5.7)
```

Write `a(u)=[9u]_26` in balanced notation and
`delta_u=u/(130d)`.  The exact signed complement in (2), together with
(S5.4), leaves the seven free raw classes

```text
1,2,4,5,8,11,12.                                        (S5.8)
```

For either parity of a speed in (S5.8), its value `a(u)` is always at least
`4` or at most `-3`.  Since `d=11 mod 52` makes `B` even, the three forced
speeds have

```text
(u,a(u))=(B-3,-11),(B-2,-2),(B,-10).                    (S5.9)
```

There is no half-circle wrap because `B<10d`.  The only three lower-bound
calculations needed are

```text
a>=4:  ||ut|| >=4/26-B/(130d)
                   =(27d-5)/(260d)>1/11,
a<=-3: ||ut|| >=3/26>1/11,
a=-2:  ||ut|| =2/26+(B-2)/(130d)
                   =(33d+1)/(260d)>1/11.                (S5.10)
```

The strict comparisons in the first and third rows reduce respectively to
`37d>55` and `103d+11>0`.  Thus `t` lies in `E_U` for **every** lift allowed
by (2) and the shell-five owner conclusion.  On the other hand, `p=4 mod 5`
gives

```text
||9dt||+||4dt||=||9p/5||+||4p/5||=2/5<11/13.            (S5.11)
```

Hence `t notin H_d`, contradicting (3) already through its necessary
consequence `E_U subset H_d` (because `0 in closure(R_U)`).  We have proved
the residue-refined uniform exclusion

```text
s=5 and (2)--(3)  implies  d mod 52 in {15,37,41}.       (S5.12)
```

There is also a finite exact exclusion at the first remaining scale.  Add
the actual two-sheet-core conclusions from THM-772 and THM-803: `U` contains
a multiple of every `2,...,12`, and its parity-twisted mod-13 support is full.
At

```text
d=15,                  B=100,                            (S5.13)
```

the signed complement and (S5.4) force the speeds `97,98,100` in raw classes
`6,7,9`.  Exhausting the other seven signed lifts gives

```text
all signed lifts                         1,605,632,
divisor-complete rows                      121,352,
also parity-support-complete                 71,644,
cover the mandatory q=75 endpoint grid       3,004,
cover the mandatory q=195 endpoint grid          4,
cover both endpoint grids                         0.      (S5.14)
```

The replay prints the four `q=195` rows and an explicit uncovered `q=75`
unit numerator for each.  Therefore `(d,B)=(15,100)` cannot be an actual
hypothetical two-sheet quotient core once those already-proved necessary
conditions are adjoined.  Statement (S5.14) is finite at exactly `d=15`; it
does not eliminate the entire class `d=15 mod 52`, and no uniform conclusion
is asserted for the remaining classes `15,37,41` beyond (S5.12).

## 7. Why the theorem carrier is not a tournament

The exact carrier is the incidence object

```text
{left obligation,right obligation}
 x {owner label z, speed u_z, exit time c_epsilon(z)/(143u_z)}
 x {the mod-13 separation of the two selected owners}.    (37)
```

A lossy tournament can still be formed for telemetry.  Orient owner pairs by
their right-exit times, switch to their left-exit times, and break ties by the
fixed Hamiltonian path `(signed residue,speed)` in lexicographic order.  Each
gauge is a scalar total order, hence a transitive tournament with score
histogram `0,...,9`, ten singleton SCCs, no directed cycle, and one
Hamiltonian path.  The exact replay records their edge flips.

That tournament is not theorem-bearing.  A binary order forgets that **both**
directional minima must be served simultaneously and that the two selected
owners must fit into one real interval while occupying residue classes
separated by `delta(d)`.  The two-sided owner-labelled packing obligation in
(37), rather than runners or a bare tournament, preserves the deduction.

## 8. Verification and scope

The independent verifier uses only integers and `fractions.Fraction`.  It:

1. reconstructs all twenty signed direction/owner coefficients and checks
   the closed boundary, strict interior, and strict exterior directly;
2. builds `119` deterministic signed-wall lift rows, checks `9,996` local
   endpoint/return/target facts, the forced `8/13` centre, both component
   radii, and the tournament telemetry;
3. exhausts `1,502,751` thin-shell parameter rows for odd `d<=999`, with no
   modular-packing failure and no feasible `s=1` or `s=3` owner row; and
4. finds `77` locally feasible `s=5` rows, beginning with (35), guarding the
   exact scope of the shell conclusion.

The shell-five continuation verifier, also using only integers and
`fractions.Fraction`, separately:

1. replays the exact local classification through every odd `d<=9,999`,
   finding `770` rows in precisely the four classes in (S5.4);
2. directly checks the explicit divisor witness (S5.6) on the first `100`
   members `d=52k+11`, through `d=5,159`, while (S5.7)--(S5.11) supply the
   all-size proof;
3. exhausts the `1,605,632` signed lifts in (S5.14), printing all four
   `q=195` survivors and their uncovered `q=75` numerator sets; and
4. records a four-class planning tournament.  Ordering by local owner slack
   and then by divisor-proof priority flips three edges; both gauges are
   transitive.  This tournament forgets the unit-numerator/signed-lift
   incidence that proves the exclusion.

No optimizer, floating point, or sampled-circle verdict is used.  The replay
certifies the arithmetic and endpoint calculus; the proof above supplies the
all-size quantifiers.  The theorem assumes the signed-wall complement and the
guarded containment (2)--(3).  It does not prove those predicates for an
arbitrary core, does not show the remaining `s>=5` branch empty, does not
uniformly eliminate the classes `15,37,41 mod 52`, and does not establish
LRC(14).  The `d=15` census additionally assumes THM-772 divisor completeness
and THM-803 parity-twisted support; it is not a consequence of (2)--(3) alone.

Reproduce with

```bash
python3 04-computation/lrc14_two_sheet_owner_shell_packing_obstruction_codex_S10.py
python3 04-computation/lrc14_shell_five_local_classification_and_divisor_witness_codex_S10.py
```

The frozen artifact digests are

```text
source  d6b259c0324752bfcdbf1efe5583843debcebcae898ebe8aca22629a7e0d246d,
output  0c10a33bfc6e6c79b5b06e85a9f6fc3949b3b4214331a9a52846d3a2a21ae2e4,
certificate 741dc1d6e7e476817867d9b6ec13181e18b5ef7f999c9a7253680b2477f7fc07.

shell-five continuation source
        7bdbf2ddfb9527775d45589d40fcb85e24cc8defd30ea9f91a751e0e9ec098f9,
shell-five continuation output
        50cd0f1761d202ef6a8565d58c6d70951221101f55649ea27e8cf9e52e9cbbfb,
shell-five continuation certificate
        b5b0d682307e61c71eded52fb94f286d8ecd33f8f9c0a00abbce164e2b20b7f4.
```

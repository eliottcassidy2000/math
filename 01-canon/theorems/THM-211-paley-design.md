# THM-211: Paley `P_7` three-cycles form a simple `2-(7,3,2)` design

**Status:** PROVED by the difference-set argument below. The independent exact
census is `04-computation/sum_h_proof_89b.py`; THM-2069 supplies a second
code/cogirth audit.

## Statement

Let `Q={1,2,4} subset F_7` and orient the Paley tournament by

```text
i -> j iff j-i in Q.
```

It has exactly fourteen directed three-cycle vertex sets. They form a simple
`2-(7,3,2)` design: every pair of vertices belongs to exactly two blocks.
More precisely, the block set is the block-disjoint union of two distinct
labelled cyclic Fano planes,

```text
B_v^+=v+Q=N^+(v),       B_v^-=v-Q=N^-(v),     v in F_7. (1)
```

The seven pairs `(B_v^+,B_v^-)` are the exact near-parallel classes: their
blocks are disjoint and cover `F_7\{v}`.

## Difference-set proof

Both triples in (1) are directed cycles. Translating to `v=0`,

```text
1 -> 2 -> 4 -> 1,
3 -> 5 -> 6 -> 3.                                      (2)
```

Each translation orbit has seven distinct blocks. The two orbits are also
disjoint. Indeed, if `Q+a=-Q+b`, summing the elements modulo seven gives
`3a=3b`, hence `a=b`; this would force `Q=-Q`, which is false.

The tournament is regular of outdegree three. Every transitive triple has a
unique source, so the number of transitive triples is

```text
sum_v C(outdeg(v),2)=7 C(3,2)=21.                       (3)
```

There are `C(7,3)=35` vertex triples in total. Hence there are `35-21=14`
cyclic triples, and the fourteen distinct blocks in (1) exhaust them.

The six ordered nonzero differences of `Q` are

```text
1-2=6, 1-4=4, 2-1=1, 2-4=5, 4-1=3, 4-2=2.             (4)
```

Thus every nonzero element of `F_7` occurs once: `Q` is a `(7,3,1)`
difference set. The same holds for `-Q`. Each translation orbit is therefore
a cyclic `2-(7,3,1)` design, and their block-disjoint union is a simple
`2-(7,3,2)` design.

Finally `Q` and `-Q` are disjoint and their union is `F_7^x`, so

```text
B_v^+ disjoint_union B_v^-=F_7\{v}.                    (5)
```

These are the only disjoint block pairs. Distinct translates in the same
Fano plane intersect by the difference-set property. Across the two planes,
`(a+Q)` meets `(b-Q)` exactly when `b-a in Q+Q`; direct addition gives
`Q+Q=F_7^x`, so it is disjoint exactly when `a=b`. Hence there are exactly
seven near-parallel classes and `d_33=7`. QED.

## The seven near-parallel classes

| Pair | Triple A | Triple B | Leftover |
|------|----------|----------|----------|
| 1 | `{0,1,3}` | `{2,4,5}` | `6` |
| 2 | `{0,1,5}` | `{3,4,6}` | `2` |
| 3 | `{0,2,3}` | `{1,5,6}` | `4` |
| 4 | `{0,2,6}` | `{1,3,4}` | `5` |
| 5 | `{0,4,5}` | `{1,2,6}` | `3` |
| 6 | `{0,4,6}` | `{2,3,5}` | `1` |
| 7 | `{1,2,4}` | `{3,5,6}` | `0` |

The word **near** is load-bearing: each class covers six of seven points;
the seven classes are not mutually disjoint resolutions of the whole design.

## Count convention and the IP value

Here `t_3=14` counts cyclic vertex triples, equivalently directed cycles
modulo cyclic rotation. The number of rooted cyclic ordered triples is
`3t_3=42`, not 28. A fixed tournament assigns only one direction to a cyclic
vertex triple; the design does not arise by taking two orientations of one
Fano line.

The separate exact cycle census in THM-257 gives `t_5=42` and `t_7=24` for
the Paley class. Together with the seven disjoint triangle pairs,

```text
H(P_7)=1+2(t_3+t_5+t_7)+4d_33
      =1+2(14+42+24)+4*7
      =189.                                             (6)
```

Paley `P_7` is self-complementary, so its complement has the same value.

## Marked `e8` bridge

THM-481 proves `C(border(P_7))=e8`. The seven loss-gauge rows and their
complements have weight-four supports

```text
G_v={infinity} union N^-(v),
H_v={v}        union N^+(v).                             (7)
```

They exhaust the `14y^4` layer of `1+14y^4+y^8`. The triangle map is marked:
delete the common `infinity` coordinate from `G_v`, but delete the varying row
label `v` from `H_v`. It yields the two Fano planes in (1), and complementary
codeword pairs yield the near-parallel classes. THM-2069 identifies this
weight-four layer as the first cocircuit shell of the binary deletion wheel.

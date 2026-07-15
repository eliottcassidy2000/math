---
id: THM-804
title: Three-replacement oriented sheet-deck descent in the full-residue chart
status: PROVED (scale-free half-open sheet-capacity and residue-ratio argument) + VERIFIED (296,640 direct grids, 11,143,660 relaxed capacity rows)
source: codex-2026-07-15-S10 (recursive Hamming continuation)
depends_on: []
related: [THM-769, THM-770, THM-795, THM-800, HYP-6775, HYP-6800, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_three_oriented_deck_descent_codex_S10.py
  - 05-knowledge/results/lrc13_hamming_three_oriented_deck_descent_codex_S10.out
---

# THM-804 — Three-replacement oriented sheet-deck descent

Put

```text
delta=1/13,                       [12]={1,...,12},
M(W)=max_(t in R/Z) min_(w in W)||wt||.
```

## Theorem

Let `13` not divide `c`, let `r_1,r_2,r_3` be distinct members of `[12]`,
and let the positive integers `w_i` satisfy

```text
w_i = c r_i  (mod 13).                                  (1)
```

Define

```text
A=(c[12]\{c r_1,c r_2,c r_3}) union {w_1,w_2,w_3},
D_i=c/gcd(c,w_i).                                        (2)
```

If `M(A)=1/13`, then

```text
D_1=D_2=D_3=1,
c divides w_1,w_2,w_3.                                  (3)
```

Consequently `A=cB`, where `B` is a labelled scale-one lift of `[12]` in
exactly those three coordinates (some lift heights may be zero).

This is a common-scale descent, not a classification of scale-one triple
lifts.  Together with THM-795 and THM-800 it says that the unbounded shallow
problem through labelled Hamming radius three has only one remaining base:
genuine scale-one Hamming-three lifts.  It does **not** say that those lifts
are loose.

## 1. The oriented missing-owner germ

Fix a missing owner `r` and choose `a in [12]` with `ar=1 (mod 13)`.  The `c`
preimages of `a/13` under `t -> ct` are

```text
t_l=(a+13l)/(13c),                 0<=l<c.                (4)
```

If the complementary owner `13-r` remains in the core, it has signed phase
`-1/13` at (4), while every other retained core runner has clearance at least
`2/13`.  Moving a sufficiently small distance to the left makes that binder
strictly safer.  If `13-r` is one of the other deleted owners, every retained
core runner already has clearance at least `2/13`, so a two-sided safe
neighborhood exists.  In either case there is a left-hand core-safe germ at
every `t_l`.

Tightness of `A` forces the three replacements to cover each such germ.  At
its endpoint a replacement covers the left germ exactly when its signed phase
lies in

```text
(-1/13,1/13].                                           (5)
```

The positive boundary is included because moving left moves a positive speed
into its danger tooth.  The negative boundary is excluded because moving left
makes its norm larger.  This half-open convention is the decisive information
that an unoriented endpoint count loses.

## 2. Exact deck capacities

Consider a replacement with label `r`, speed `w`, and deck order
`D=c/gcd(c,w)`.  Write `c=gD` and `w=gu`, so `gcd(u,D)=1`.  Because `13` does
not divide `c`, it does not divide `D`, and (1) gives

```text
u = D r  (mod 13).                                      (6)
```

Over the splice family belonging to an owner `o`, its phases are

```text
u(a_o+13l)/(13D),                 a_o o=1 (mod 13).      (7)
```

As `l` runs modulo `D`, the eligible numerator is one of the integers

```text
-D < z <= D,             z = D r o^(-1) (mod 13).       (8)
```

Every eligible deck class repeats `c/D` times.  Thus, if `N_D(r,o)` is the
number of integers in (8), the fraction of the `c` sheets covered at owner
`o` is exactly `N_D(r,o)/D`.

The maximum possible fraction is

```text
f(D)=ceil(2D/13)/D.                                     (9)
```

At its own owner a replacement attains (9).  At the complementary owner its
fraction is

```text
g(D)=floor(2D/13)/D.                                    (10)
```

Indeed the own congruence class contains the included endpoint `z=D`, whereas
the complementary class would contain the excluded endpoint `z=-D`.

The elementary bounds needed below are

```text
f(1)=1,              f(2)=1/2,
f(D)<=1/3  (D>=3),   equality iff D=3,                  (11)
f(D)+g(D)<=3/7        (D>=3), equality iff D=7.          (12)
```

For (11), check `D=3,4,5` and use
`ceil(2D/13)<=2D/13+1<=D/3` for `D>=6`.  For (12), check
`D=3,...,8`.  For `D>=9`, since `13` does not divide `D`,

```text
f(D)+g(D)=(2 floor(2D/13)+1)/D
          <=4/13+1/D<=4/13+1/9<3/7.                    (13)
```

An order-one replacement covers all sheets at its own owner: its constant
phase is the included `+1/13`.  At every distinct owner its constant nonzero
residue is outside (5); in particular a complementary owner sees the excluded
`-1/13`.  Hence

```text
D=1: own-owner capacity 1, every cross-owner capacity 0. (14)
```

For `D=2`, a distinct cross owner is hit exactly when

```text
owner/replacement in {2,-2}={2,11} in F_13.             (15)
```

This follows from (8): `2r/o` must be one of the nonzero residues represented
by `-1,0,1,2`.

## 3. The two-colour sublemma

Suppose replacements labelled by two distinct owners cover every sheet at
both of their owner families.  Then both deck orders are one.

If one order is one, (14) says it contributes nothing at the other owner, so
the second colour must have capacity one and also has order one.  If neither
order is one, (9) is at most `1/2`; covering forces both orders to be two.
At the first owner the other replacement must cross by (15), and the reverse
must hold at the second.  This would require both `z` and `z^(-1)` to belong
to `{2,11}`, but

```text
{2^(-1),11^(-1)}={7,6},                                (16)
```

which is disjoint from `{2,11}`.  This proves the sublemma.  It is the local
oriented-deck core of THM-800, reproduced here so the present descent has no
logical dependence on the scale-one double-lift computation.

## 4. Proof of the theorem

If some `D_i=1`, that replacement contributes zero at the other two missing
owners by (14).  The remaining two replacements therefore cover all sheets
at both of their owner families.  The two-colour sublemma makes both of their
orders one, proving (3).

It remains to rule out `D_i>=2` for all three replacements.  Split by the
number of order-two decks.

### Three order-two decks

At every owner the own colour supplies only half the sheets, so some other
colour must cross.  Make a directed edge from a replacement to every distinct
owner it can hit.  By (15), every edge ratio is `+2` or `-2`, and every one of
the three vertices has positive indegree.  The digraph therefore contains a
directed 2- or 3-cycle.  A 2-cycle is impossible because the inverses of
`+/-2` are `+/-7`; a 3-cycle is impossible because a product of three signed
twos is `+8` or `-8`, never `1` modulo 13.

### Exactly two order-two decks

Let their labels be `a,b`, and let the third label be `d`.  At owner `d`, its
own capacity is at most `1/3`.  Therefore both order-two replacements must
cross there.  Formula (15) gives

```text
d/a,d/b in {2,-2},
```

so distinctness forces `b=-a`.  At owner `a`, the replacement `b` cannot
cross, since `a/b=-1`, while the third replacement contributes at most
`1/3`.  Together with the own half-capacity this is less than one, a
contradiction.

### Exactly one order-two deck

Let its label be `a`; call the other labels `b,d`, with orders `D,E>=3`.
At owner `b`, if `a` did not cross, the other two total capacities would be at
most `2/3`; hence `a` crosses.  The same holds at owner `d`.  By (15),

```text
b/a,d/a in {2,-2},
```

and distinctness forces `d=-b`.

At owner `b`, subtracting the order-two cross capacity gives

```text
f(D)+g(E)>=1/2.                                          (17)
```

At owner `d` the symmetric inequality is

```text
f(E)+g(D)>=1/2.                                          (18)
```

Adding (17)--(18) contradicts (12), because its left side is at most
`3/7+3/7<1`.

### No order-two deck

All three orders are at least three.  By (11), the total capacity at any
owner is at most one.  Covering forces equality everywhere, so all three
orders equal three and every replacement must hit every owner.

For an order-three replacement, the distinct cross-owner ratios allowed by
(8) are

```text
{3,5,8,10}.                                               (19)
```

Both directions are required for every pair of labels.  The only elements of
(19) whose inverses also lie in (19) are `5` and `8=5^(-1)`.  Three distinct
labels cannot be pairwise related by `{5,8}`: after scaling the first to one,
the other two would be `5` and `8`, whose ratio is `12`.  This last
contradiction exhausts every case and proves (3).  ∎

## 5. What this changes in the shallow frontier

THM-795 removed every proper Hamming-one packet.  THM-800 now removes every
proper Hamming-two packet, using its exact scale-one floor.  The present
theorem shows that a tight packet at Hamming radius three cannot carry any
new scale ramification: it must divide all the way to a scale-one triple lift.

Thus the next precise shallow question is no longer an arbitrary-scale
Hamming-three problem.  It is:

```text
Does every genuine scale-one triple lift
([12]\{r,s,t}) union {r+13i,s+13j,t+13k}, i,j,k>=1,
have M>1/13?                                              (20)
```

The historical low-height probes support (20), but no uniform triple-lift
floor is claimed here.  Analytic safe-window bounds make its residual finite
after ordering the three replacements; producing a practical exact closure is
the next computation/proof target.

## 6. Tournament Analysis and assumption challenge

The theorem-bearing object is not a tournament on runners.  It is the
bipartite incidence between:

- owner-sheet obligations `(o,l)`;
- replacement colours and their deck classes; and
- the oriented half-open germ flag in (5).

For telemetry, collapse each replacement with its owner vertex and compare
the two cross capacities on a pair.  Orient toward the larger cross capacity,
break ties by numerical label, and use numerical order as the tie Hamiltonian
path.  Switching from the left germ `(-1/13,1/13]` to the right germ
`[-1/13,1/13)` is the gauge.

The exact replay records three proof-critical rows as transitive tournaments
with score histogram `(0,1,2)`, no directed cycle, singleton SCCs, and one
Hamiltonian path.  A fourth row

```text
labels=(6,7,12),       decks=(12,5,11)
```

changes under one germ-gauge edge flip from a directed triangle (scores
`(1,1,1)`, one SCC, three Hamiltonian paths) to a transitive tournament.  Its
owner-capacity verdict does not change.  The tournament is therefore a useful
fingerprint but cannot replace the decorated incidence deck.

Alternative vertices were challenged explicitly.  Runner vertices lose sheet
ramification; gap vertices lose owner congruences; fixed circle sections lose
the deck action; tooth vertices lose which side of a boundary is safe; residue
vertices lose multiplicity; and proof-obligation vertices without the
half-open flag identify `+1/13` with `-1/13`, exactly the false identification
that would preserve mixed complementary decks.  The minimal faithful carrier
is the oriented owner-sheet/replacement incidence deck.

## Exact replay

The `Fraction`/integer verifier checks `296,640` direct CRT grid identities,
the capacity formulas and sharp inequalities through deck order `100,000`,
the exact order-two and order-three ratio alphabets, and `11,143,660` relaxed
owner-capacity rows over every residue triple and every order at most 40.  The
relaxation even discards sheet-overlap compatibility and the requirement that
the three deck orders divide one common `c`; its only survivors are the 220
residue triples with `(D_1,D_2,D_3)=(1,1,1)`.  The infinite theorem rests on
the symbolic case split above, not on that finite stress test.

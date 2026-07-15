---
id: THM-823
title: Hamming-five scalar-deck escape and common-sheet survivor boundary
status: PROVED STRUCTURAL + FINITE-EXACT — scalar cutoff refuted and bounded common-sheet bank classified; all of its survivor languages subsequently closed by THM-844/845/847; orders above twelve remain open
source: codex-2026-07-15-S10 arbitrary-scale five-deck continuation
depends_on: [THM-810, THM-816, THM-820, THM-822]
related: [THM-815, THM-837, THM-844, THM-845, THM-847, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_five_common_sheet_deck_boundary_codex_S10.cpp
  - 05-knowledge/results/lrc13_hamming_five_common_sheet_deck_boundary_codex_S10.out
---

# THM-823 — Hamming-five scalar-deck escape and common-sheet boundary

Put

```text
delta=1/13,                         [12]={1,...,12}.
```

For distinct replacement labels `r,o in F_13^*` and an effective order
`D>=1`, with `13` not dividing `D`, define the left-oriented owner capacity

```text
N_D(r,o)=#{z in Z:-D<z<=D and z=D r o^(-1) (mod 13)},
C_D(r,o)=N_D(r,o)/D.                                      (1)
```

For five replacement labels `R={r_1,...,r_5}`, scalar deck coverage is the
necessary condition

```text
sum_i C_(D_i)(r_i,o)>=1                  for every o in R. (2)
```

This theorem separates four predicates which first diverge sharply at
Hamming radius five:

```text
scalar owner capacity,
common-sheet coverage,
strict looseness of the least CRT lifts,
strict looseness of all lifts.                               (3)
```

Only the first three are decided on the finite banks below.  The fourth is
explicitly open.

## Theorem

### A. Attenuation and failure of a scalar order cutoff

Write `D=13q+d`, `1<=d<=12`.  Then

```text
C_D(r,o)=2/13+(d/D)(C_d(r,o)-2/13),
|C_D(r,o)-2/13|<1/D.                                     (4)
```

Consequently every five-colour scalar cover satisfies

```text
sum_i 1/D_i>3/13,                       hence min_i D_i<=21. (5)
```

The least-order conclusion is the available finite pivot.  There is no
finite upper bound on every order.  For every `q>=3`, set

```text
(r_1,...,r_5)=(1,2,3,5,10),
(D_1,...,D_5)=(2,5,13q+1,2,13q+1).                       (6)
```

This is a scalar cover.  In particular, five-colour scalar capacity does not
admit the finite all-order classification used by THM-810 at radius four.

### B. Exact compactified small-order census

Adjoin a formal symbol `infinity` whose capacity is identically `2/13`.  In
the bank

```text
1 in R subset [12],       |R|=5,
D_i in {2,3,...,12,infinity},                              (7)
```

there are exactly

```text
scalar rows                                           2,410
rows with 0,1,2 formal infinities              2,190,140,80
rows with 3,4,5 formal infinities                    0,0,0
unordered order patterns                                  86. (8)
```

Here a row is a presentation with label `1` present; it is not already
quotiented by multiplication in `F_13^*`.  The `220` rows at infinity have
exactly the following order patterns.  The orbit column is the number of
decorated multiplicative orbits, and the margin is
`min_o(sum_i C_i(o)-1)`.

| order pattern | presentations | decorated orbits | margin |
|---|---:|---:|---:|
| `(2,2,5,infinity,infinity)` | 50 | 10 | `1/130` |
| `(2,2,10,infinity,infinity)` | 30 | 6 | `1/130` |
| `(2,2,8,9,infinity)` | 60 | 12 | `1/936` |
| `(2,3,3,10,infinity)` | 20 | 4 | `4/195` |
| `(2,3,3,11,infinity)` | 20 | 4 | `1/429` |
| `(3,3,3,5,infinity)` | 20 | 4 | `4/195` |
| `(3,3,3,10,infinity)` | 10 | 2 | `4/195` |
| `(3,3,3,11,infinity)` | 10 | 2 | `1/429` |

Every margin is uniform within its displayed pattern.  If the formal
infinities in any one decorated row are replaced by admissible finite orders
`E_j` satisfying

```text
sum_j 1/E_j <= displayed margin,                           (9)
```

then the resulting row is still a scalar cover.  Thus these are robust
attenuation cones, not isolated limits.  This compactification is only the
alphabet (7); it is not a classification of all finite orders.

### C. Complete common-sheet decision for `2<=D_i<=12`

Choose a unit `e_i (mod D_i)` and let `u_i (mod 13D_i)` be the CRT lift

```text
u_i=D_i r_i (mod 13),                  u_i=e_i (mod D_i). (10)
```

With `c=lcm(D_1,...,D_5)`, colour `i` covers at owner `o` the sheets

```text
E_i(o)={ell in Z/cZ:
        -D_i < <u_i(o^(-1)+13ell)>_(13D_i) <= D_i}.       (11)
```

The exact common-sheet condition is

```text
union_i E_i(o)=Z/cZ                         for every o in R. (12)
```

Among every scalar row with

```text
1 in R,                         2<=D_i<=12,               (13)
```

the exact census is

```text
scalar rows                                            2,190
scalar order patterns                                     78
common-sheet rows                                           5
common-sheet order patterns                                 1. (14)
```

All five survivors are presentations of one decorated multiplicative orbit:

```text
D_1=...=D_5=3,
R=C union {b},
C=a{1,5,8,12},                       b in 2C.              (15)
```

The left half-open orientation in (1) matters: the four-coset points in
`C` contribute `2/3` at the fifth owner precisely in the forward doubled
coset `2C`.  Every owner sum in (15) is exactly one.

For each row in (15), exactly eight of the `2^5=32` unit words satisfy (12).
On the quartet `C`, the two labels in each opposite pair `{r,-r}` must have
equal mod-three units, exactly as in THM-810.  The unit of `b` is free:

```text
e_r=e_(-r) on C,                         e_b arbitrary.    (16)
```

Thus the five rows contribute exactly forty common-sheet unit words.

The scalar/common-sheet distinction is large even inside the explicit
family (6).  At `q=3`, hence `D=40`, none of the `1*4*16*1*16=1,024` unit
words gives a common-sheet cover.  The best word covers only `29` of the
`40` sheets at its least-covered owner; exactly `96` unit words attain that
value.

### D. Least-CRT metric census and the threshold guardrail

For one of the forty unit words in (16), take every `u_i` to be its least
positive CRT representative and form the twelve-speed packet

```text
B=3([12] minus R) union {u_1,...,u_5}.                    (17)
```

All forty least-CRT packets are strictly loose:

```text
min M(B)=2/17>1/13,                      attained 4 times. (18)
```

Their complete exact maximum histogram is

| `M(B)` | count | `M(B)` | count |
|---:|---:|---:|---:|
| `2/17` | 4 | `1/8` | 6 |
| `6/43` | 8 | `1/7` | 4 |
| `9/46` | 4 | `6/41` | 1 |
| `3/20` | 2 | `5/31` | 4 |
| `6/37` | 2 | `11/64` | 1 |
| `1/5` | 4 |  |  |

This is not an all-height metric theorem.  THM-810's eight-point mod-`39`
clock gives the exact boundary.  For every quartet parity word, every
`b in 2C`, and either unit of `b`, exactly six of those eight clock points
remain safe for (17).  Adding arbitrary multiples of `39` to any lift does
not change this calculation, so every arbitrary-lift packet in (15) has

```text
M(B)>=1/13.                                                (19)
```

The clock may attain equality.  It proves the conjectured lower bound, not
strict looseness.  Whether every arbitrary lift in (15) has `M(B)>1/13` is
left open by the present argument and subsequently answered affirmatively by
THM-844.

### E. The symbolic order-one branch

Suppose at least one colour has order one, with label `b`.  That colour has
capacity zero at each of the other four owners.  The remaining four colours
therefore scalar-cover their own four labels, so THM-810 applies verbatim.
Exactly two alternatives remain:

```text
all five orders are one; or
D_b=1 and the other four orders are three on a coset
C=a{1,5,8,12}, with b outside C.                           (20)
```

The all-order-one branch is the genuine scale-one Hamming-five chart reduced
by THM-820 and subsequently closed by THM-845.

In the mixed branch, common-sheet coverage holds exactly when the quartet
parities satisfy (16).  At the owner `b`, the order-one colour fills every
common sheet; at the four owners in `C`, it fills none, leaving precisely
THM-810's quartet condition.

At common scale three, the order-one replacement speed is `3u_b` with
`u_b=b (mod 13)`, hence it is congruent to the removed core speed `3b`
modulo `39`.  All eight quartet clock points therefore survive, and again

```text
M(B)>=1/13                                                   (21)
```

for arbitrary lifts.  As in (19), strictness is not proved here; THM-847
subsequently proves it for every proper mixed lift, while the zero-height face
is THM-816.

## Proof

### 1. The attenuation identity and the infinite scalar family

The interval in (1) has length `2D`.  Every residue occurs either
`floor(2D/13)` or `ceil(2D/13)` times.  Splitting it into `q` complete
thirteen-block pairs and the oriented `d`-window gives

```text
N_D(r,o)=2q+N_d(r,o),                                    (22)
```

which is (4).  The residue discrepancy in an interval is strictly less than
one because `13` does not divide `D`, proving the strict error bound.

Under (2), subtracting the five baseline contributions `2/13` gives

```text
sum_i(C_(D_i)(r_i,o)-2/13)>=3/13.
```

Each summand is strictly less than `1/D_i`, so (5) follows.  If every order
were at least `22`, then

```text
sum_i 1/D_i <=5/22<3/13,
```

a contradiction.

For (6), the order-two colour at label `1` covers owners `1,2` with capacity
`1/2`; the order-two colour at label `5` covers owners `3,5,10` with capacity
`1/2`; and the order-five colour at label `2` contributes `1/5` at all five
owners.  These three colours therefore contribute exactly `7/10` everywhere.
Each large order `D=13q+1` contributes `2q/D` everywhere and an extra `1/D`
at its own owner.  The two large colours contribute at least

```text
4q/(13q+1)>=3/10                         iff q>=3,          (23)
```

proving scalar coverage and the failure of a global order cutoff.

### 2. Exact finite implementations

The compactified bank (7) has

```text
C(11,4) 12^5=82,114,560                               (24)
```

conceptual rows.  The replay multiplies every capacity by

```text
360360=lcm(2,...,12)*13
```

and makes all coverage decisions with integers.  It canonically multiplies
each surviving decorated row by all twelve elements of `F_13^*` to count
orbits.  The complete result is (8) and the eight-row table.  Inequality (9)
then follows directly from the strict error in (4).

The finite no-order-one sub-bank has

```text
C(11,4) 11^5=53,146,830                               (25)
```

conceptual rows.  After the exact scalar filter leaves `2,190`, the replay
enumerates every unit choice.  It constructs the literal sets (11) as bit
sets on `Z/cZ` and tests (12) owner by owner.  This gives (14)--(16) without
sampling.  An independent structural reading is that the surviving four
labels are THM-810's order-three coset and its fifth-owner deficit is repaired
only by the forward flag `C -> 2C`.

For (18), `min_w ||wt||` is piecewise linear.  A positive local maximum is at
a self cusp or an intersection of two signed branches.  It is therefore
enough to test exact candidates with denominator

```text
2u,                              u+v,                  |u-v|. (26)
```

The replay tests all such rational candidates for each of the forty packets,
retains the exact maximum and first witness, and obtains the displayed
histogram.

### 3. Order one and the clock boundary

At `D=1`, the oriented window contains only the nonzero residue `1`; hence a
colour covers its own owner and no other replacement owner.  Removing that
owner from (2) invokes THM-810 on the remaining quartet and proves (20).
The sheet statements follow from (11): the order-one colour lifts its single
own-owner sheet to all three common sheets and has an empty lift at the other
owners.

For the all-order-three flag, the replay checks all three cosets, all four
quartet parity words, all four forward labels, and both fifth units: `96`
tests, each retaining exactly six clock points.  For the mixed branch it
checks the three cosets, four parity words, and eight outside labels: another
`96` congruence tests.  These are the boundary statements (19) and (21), not
strictness claims.

## Tournament Analysis and the challenged vertex set

The runner-pair observable

```text
C_3(r_i,r_j)-C_3(r_j,r_i)                               (27)
```

on the representative `(1,2,5,8,12)`, with increasing-label ties, has score
word `(4,1,3,2,0)`, eight tied pairs, no directed triangle, five singleton
SCCs, and one Hamiltonian path.  It is transitive telemetry and does not
explain (15).

Challenge the assumption that tournament vertices must be runners.  Use the
three multiplicative cosets

```text
C_0=H,                 C_1=2H,                 C_2=4H,
H={1,5,8,12},                                             (28)
```

and define

```text
A(C_i,C_j)=# order-three hits from the quartet C_i
           at one owner in C_j.                           (29)
```

This is well-defined and equals

```text
A=((3,2,0),(0,3,2),(2,0,3)).                              (30)
```

Orient `C_i -> C_j` when the off-diagonal entry is `2` rather than `0`.
The result is the directed three-cycle

```text
C_0 -> C_1 -> C_2 -> C_0,                                (31)
```

with score word `(1,1,1)`, one directed triangle, one SCC of size three, and
three Hamiltonian paths.  Its forward pointed flags are exactly (15).

This quotient preserves the left-oriented scalar extension predicate.  It
destroys the unit parity, the owner-sheet overlap, every lift height, and the
metric distinction between equality and strictness.  The theorem-bearing
carrier beyond the scalar classification is therefore the decorated
owner-sheet obligation hypergraph, not either bare tournament.  The formal
infinity compactification loses even the sheet periodicity and must not be
used as a metric quotient.

## Subsequent closure and exact scope

THM-844 closes every arbitrary lift in the `96` labelled all-order-three
contexts, THM-845 closes the normalized all-order-one chart, and THM-847
closes all `96` proper mixed order-one/order-three contexts (with the
zero-height face supplied by THM-816).  Consequently every common-sheet
survivor language whose five effective orders are at most twelve is now
strictly loose.  This is a later metric closure of the languages classified
here, not an extension of the present common-sheet census to larger orders.

THM-823 proves no more than the following:

- a least-order pivot `min D_i<=21`, not a bound on all five orders;
- the exact compactified alphabet (7), not all finite orders;
- the complete no-order-one common-sheet classification only for
  `2<=D_i<=12`;
- strict looseness only for the forty least-positive CRT packets;
- boundary inequalities `M>=1/13`, not strictness, for arbitrary lifts in the
  all-order-three and mixed order-one branches; strictness comes from the
  later theorems above.

It does **not** classify common-sheet covers with an order above twelve or
close the full arbitrary-scale Hamming-five problem.  Those unbounded-order
languages are the remaining shallow deck frontier.

## Exact replay

Run

```bash
clang++ -std=c++20 -O3 -Wall -Wextra -pedantic \
  04-computation/lrc13_hamming_five_common_sheet_deck_boundary_codex_S10.cpp \
  -o /tmp/thm823
/tmp/thm823 | cmp - \
  05-knowledge/results/lrc13_hamming_five_common_sheet_deck_boundary_codex_S10.out
```

Independent `-O3` and `-O0` builds are warning-free and reproduce the stored
output byte for byte.  The frozen SHA-256 digests are

```text
source  228296294dd5b7fae50f01a44abb533d41ed4d9848347aa755dce460b5f32a41
output  f9d2372dff0c80499f95229353970afb0f8aa22d981d34634a56a7ab1083a795
```

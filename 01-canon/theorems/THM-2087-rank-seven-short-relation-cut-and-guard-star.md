---
id: THM-2087
title: "Rank-seven short-relation cut and guard-anchored star"
status: >
  PROVED. Under rank-seven guard containment, the graph joining pairs whose
  guard/two-speed triple has relation height greater than 57 is disconnected.
  Hence the short-relation graph contains a spanning complete bipartite cut
  and at least six height-57 relations. Either one speed is a bounded rational
  multiple of the guard, or the entire seven-core lies in a guard-anchored
  two-anchor star with coefficient height at most 6498. This strengthens
  THM-2085's one-relation conclusion without enlarging the spectral cutoff.
source: codex-2026-07-22-LRC14-short-relation-cut
depends_on:
  - THM-2081
  - THM-2085
related:
  - THM-2052
  - THM-2062
  - THM-2065
  - THM-2082
  - THM-2086
---

# THM-2087 -- short-relation cut and guard star

Retain THM-2081/2085's notation. Thus `Q={q_1,...,q_7}` and

```text
G_Q subset E_h.                                         (1)
```

For positive integers `u_1,...,u_k`, write

```text
lambda(u_1,...,u_k)
 =min{||m||_infinity:0!=m in Z^k, sum m_r u_r=0}.       (2)
```

## 1. The relation-free graph is disconnected

Define a graph `R_57(h,Q)` on the seven labelled speeds by

```text
ij is an edge  iff  lambda(h,q_i,q_j)>57.               (3)
```

Then

```text
R_57(h,Q) is disconnected.                              (4)
```

### Proof

Suppose instead that `R_57` is connected and choose one of its spanning
trees `T`. Every vertex has an incident edge. Consequently

```text
lambda(h,q_i)>57 for every i,                           (5)
```

because a relation on `(h,q_i)` of height at most `57`, extended by a zero
coefficient, would delete every edge incident to `i`.

The degree-57 signed box bounds proved inside THM-2085 now give

```text
I_i>=5363/164836                 for every vertex i,
w_ij>=655135/66923416            for every edge ij of T. (6)
```

The maximum spanning-tree weight is at least the weight of `T`, so exactly
the THM-2085 arithmetic yields

```text
tau_h(Q)-Delta_h(Q)
 >=6*(655135/66923416)
   -(2/7-7*(5363/164836))
 =6435/8365427>0.                                      (7)
```

This contradicts the THM-2081 consequence
`tau_h(Q)<=Delta_h(Q)` of containment (1). Hence (4). QED.

Notice that the proof needs only one relation-free spanning tree, not all 21
relation-free pairs. This is the graph-theoretic information discarded by the
scalar statement of THM-2085.

## 2. A complete cut of short relations

Let `A` be the vertex set of any connected component of `R_57`, and put
`B=Q minus A`. Both sets are nonempty. By the definition of components,

```text
lambda(h,q_i,q_j)<=57       for every q_i in A, q_j in B. (8)
```

Equivalently, for every cross pair there are integers

```text
a_ij h+b_ij q_i+c_ij q_j=0,
(a_ij,b_ij,c_ij)!=(0,0,0),
max(|a_ij|,|b_ij|,|c_ij|)<=57.                         (9)
```

Thus the graph of height-57 triple relations contains the complete bipartite
graph `K_(|A|,|B|)`. In particular it contains at least

```text
|A||B|>=6                                                (10)
```

distinct indexed relations. This is a multiplicity theorem, not merely a
bounded-relation alternative.

## 3. Guard-ratio versus guard-anchored star

There are now two branches.

### Branch I: a two-term guard relation

Suppose some `q in Q` satisfies `lambda(h,q)<=57`. Positivity forces both
coefficients in a nonzero relation to be nonzero and of opposite sign. After
dividing their gcd there are coprime positive integers `r,s<=57` and `d>=1`
such that

```text
h=s d,                  q=r d.                         (11)
```

If the dyadic-tower guard `h` is odd, then both `s` and `d` are odd. Hence
this branch is a finite ledger of bounded rational guard ratios, with an odd
denominator in the terminal application.

### Branch II: no two-term guard relation

Assume

```text
lambda(h,q)>57          for every q in Q.               (12)
```

Then in every cross relation (9), both speed coefficients `b_ij,c_ij` are
nonzero; otherwise (9) would give a forbidden two-term relation.

Choose anchors `x in A` and `y in B`, and choose one relation

```text
a_0 h+b_0 x+c_0 y=0,          b_0 c_0!=0,              (13)
```

with all coefficients bounded by `57`. For every `q in B`, its cross
relation with `x` already has the form

```text
A_q h+B_q x+C_q q=0,
C_q!=0,
max(|A_q|,|B_q|,|C_q|)<=57.                            (14)
```

For `q in A minus {x}`, use its cross relation with `y`,

```text
a h+b q+c y=0,                 bc!=0,                  (15)
```

and eliminate `y` with (13). The resulting identity is

```text
(c_0 a-c a_0)h+(-c b_0)x+(c_0 b)q=0.                  (16)
```

Its coefficient of `q` is nonzero and

```text
|c_0 a-c a_0|<=2*57^2=6498,
|c b_0|<=57^2,
|c_0 b|<=57^2.                                        (17)
```

For the anchor `x` itself use `-x+q=0`. Therefore every speed `q in Q`
admits integers `A_q,B_q,C_q` satisfying

```text
A_q h+B_q x+C_q q=0,
C_q!=0,
max(|A_q|,|B_q|,|C_q|)<=6498.                         (18)
```

This is a guard-anchored two-anchor star of explicit height `6498`. QED.

## 4. Frontier effect

The old residual was a union of unbounded single-relation hyperplanes.
Containment now lies in the much smaller union

```text
bounded guard ratios
or
guard-anchored height-6498 two-anchor stars.            (19)
```

Moreover, the star did not come from combining arbitrary relations in the
thirteen-speed row: it is forced by a complete short-relation cut localized
to the last guard. This preserves the coordinate needed by THM-2079's
mirror-complement owner address and THM-2082's translated-grid carrier test.

Equation (18) is in the form needed by the coefficient-row technology of
THM-2052/2062 after saturation. THM-2065 can collapse circuit-free parameter
rays, but persistent marked circuits and exact safe-phase sidecars remain.
Thus (19) is a strict structural reduction, not LRC(14).

The next finite computation should enumerate primitive cut templates up to
simultaneous sign and component swap, then attach:

1. odd guard-ratio parity from (11);
2. the rank-two deletion/cogirth wheel of THM-2062/2069;
3. THM-2082 translated-prime-grid incidence;
4. THM-2081 exact restricted-edge margins on persistent circuits.

## 5. Assumption challenge and Tournament Analysis

The challenged assumption is that the degree-57 argument can return only one
short relation. Its actual carrier is connectivity. A single relation-free
spanning tree already buys the full six-edge Hunter floor, so failure forces
a **cut** of short relations.

The faithful vertices are the restricted danger events, and an undirected edge
records whether the corresponding guard/two-speed torus has no degree-57
annihilating character. This quotient preserves exactly the Selberg edge floor
needed by Hunter and destroys the individual endpoint weights above that
floor. Orienting edges creates an arbitrary tournament and loses connectivity
of the undirected certificate; score histograms, directed cycles, SCCs, and
Hamiltonian paths are therefore irrelevant. The tie Hamiltonian path is any
ordering inside a connected component. The useful objects are the graphic
matroid of `R_57` and the complete bipartite cocircuit in its complement. QED.

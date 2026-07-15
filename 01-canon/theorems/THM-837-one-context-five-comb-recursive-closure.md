---
id: THM-837
title: One order-three Hamming-five context closes by exact five-comb recursion
status: PROVED (one of 96 arbitrary-height order-three contexts) + FINITE-EXACT (75,371-state recursive closure); subsequently strengthened by THM-844, which closes all 96
source: codex-2026-07-15-S10 continuation
depends_on: [THM-810, THM-815, THM-816, THM-820, THM-823]
related: [THM-844, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_five_single_context_comb_closure_codex_S10.py
  - 05-knowledge/results/lrc13_hamming_five_single_context_comb_closure_codex_S10.out
---

# THM-837 — one-context five-comb recursive closure

Work on the circle `T=R/Z`, at strict-safe threshold `1/13`.  For a positive
integer `u`, put

```text
D_u={t:||ut||<=1/13},
S(A)={t:||at||>1/13 for every a in A}.
```

Consider the order-three Hamming-five context

```text
C={1,5,8,12},       b=10,       R=C union {b},
pair bits=(1,1,1).
```

Thus the retained labels and retained speeds are

```text
P={2,3,4,6,7,9,11},       3P={6,9,12,18,21,27,33},       (1)
```

and the five replacement speeds lie in the labelled progressions

```text
u_1 =16 mod 39,       u_5 =28 mod 39,       u_8 =37 mod 39,
u_10= 4 mod 39,       u_12=10 mod 39.                     (2)
```

These are the least CRT representatives of `u_r=3r mod 13` and the chosen
unit class `u_r=1 mod 3`.

## Theorem

For every choice of positive speeds in (2),

```text
S(3P union {u_1,u_5,u_8,u_10,u_12}) is nonempty.          (3)
```

Equivalently, every arbitrary-height packet in this one directed-flag/parity
context is loose.  This is one of 96 order-three contexts.  The present proof
does not itself close the other 95 contexts, the other scalar-order branches,
or the global Hamming-five problem; THM-844 subsequently replaces its global
`K/L` cap by the longest-component cap and closes all 96 contexts.

## 1. Exact five-comb discrepancy

Let `E` be a union of `K` disjoint open intervals of total measure `L`.  The
periodic primitive argument from THM-816 gives, for every positive `u`,

```text
|E intersect D_u| <= 2L/13+22K/(169u).                   (4)
```

If `m<=6` danger combs, all of speed at least `v`, cover `E`, then summing
(4) and using `2m<13` gives the necessary bound

```text
v <= B_m(E):=22mK/[13(13-2m)L].                          (5)
```

Open versus closed endpoints do not affect the measure estimate.  They do
matter to the terminal predicate, so the replay stores all residuals as exact
unions of open intervals and declares a cover only when that union is empty.

## 2. Recursive completeness

Numerically order the five replacement speeds as

```text
v_1<v_2<v_3<v_4<v_5.                                     (6)
```

They are distinct because their labelled residues modulo 13 are distinct.
Starting with `E_0=S(3P)`, define

```text
E_j=E_(j-1) minus D_(v_j).                                (7)
```

If the five replacement combs covered `E_0`, then after every prefix `j<5`
the residual `E_j` would be covered by the remaining `5-j` combs.  Formula
(5) would therefore force

```text
v_(j+1)<=floor(B_(5-j)(E_j)).                             (8)
```

Conversely, each unused label has exactly one progression in (2).  At a
node with last speed `v_j`, enumerating every member of every unused
progression in

```text
v_j < v <= floor(B_(5-j)(E_j))                            (9)
```

therefore includes the next speed of every hypothetical covering row.
A node with no such candidate cannot extend to a cover by (8); a depth-five
node is a cover exactly when its exact residual is empty.  This proves the
logical completeness of the finite recursion independently of its outcome.

## 3. Reflection quotient and exact outcome

Every set in the recursion is invariant under `t->1-t`.  The retained core
contains even speeds, so `t=1/2` is not safe and no component crosses the
reflection fixed point.  Restricting to `(0,1/2)` therefore halves both `K`
and `L`, leaving the quotient `K/L` and every bound (5) unchanged.

The exact root data in this quotient are

```text
K_0=18,       L_0=2615/18018,       floor(B_5(E_0))=349. (10)
```

The fraction-exact recursion then has

```text
depth                 0      1       2       3       4    5
visited nodes         1     45    1262   20703   53303   57
dead: no candidate    0      0       0    6253   53251    0
```

There are 75,371 states, zero empty residuals, and 57 nonempty depth-five
terminal rows.  The largest recursive speed bound encountered is 1,025.
Among the enumerated terminal rows the smallest full-circle residual has

```text
measure=413437/5714280,       components=30,              (11)
```

at the least representative row

```text
(label,speed)=((10,4),(12,10),(1,16),(5,28),(8,37)).      (12)
```

The stored output freezes the state certificate

```text
9d1024dc22acdc85ccd57afe74ac23d002765ec6ce0aee0d648da64dc23d041a.
```

These facts together with recursive completeness prove (3).  The source and
stored-output payload hashes are respectively

```text
8717e9b5f88c5ef5fec18d418a66ec1a83dd630152856b88ddb7665abe1ffb12,
691f54f9b1fdcc8bc5c3a1e164d910811d34200e7c4dbb80a450d30573ec8c6c.
```

## 4. State carrier and Tournament Analysis

The proof state is the active open endpoint-run word, the remaining labelled
CRT progressions, and the last chosen speed.  A raw bit mask on the union of
all possible endpoints is unsound: endpoints of unused combs would create
artificial component cuts and corrupt `K` in (5).  A cell implementation
must retain boundary-activation bits, which is equivalent here to retaining
the active runs.

For telemetry, put the five chosen combs at tournament vertices.  The raw
pairwise observable is marginal erosion of `E_0`; the switch conditions on
the other three selected combs before comparing a pair; increasing numerical
speed is the declared tie Hamiltonian path.  This quotient measures overlap
reversal, but it destroys absolute residual measure, active endpoint
topology, and the future CRT-progression action.  It therefore cannot prove
(3).  On all 57 terminals, both raw and conditioned tournaments are
transitive with score histogram `(0,1,2,3,4)`, zero triangles, five singleton
SCCs, and one Hamiltonian path, while the gauge switches between one and ten
of the ten edges.  Identical standard fingerprints therefore coexist with
substantial conditional reorientation.  Alternate vertices—components,
endpoints, walls, or residue labels—also
fail singly: the exact predicate is preserved by the evolving bipartite
incidence of active residual runs with labelled future comb teeth.

## Replay

Run

```bash
python3 04-computation/lrc13_hamming_five_single_context_comb_closure_codex_S10.py
```

and compare byte-for-byte with the stored output named in the frontmatter.

---
id: THM-844
title: All order-three Hamming-five common-sheet contexts close by longest-component recursion
status: PROVED (uniform arbitrary-height closure of all 96 THM-823 all-order-three contexts) + FINITE-EXACT (28,876-state endpoint certificate and Tournament Analysis)
source: codex-2026-07-15-S10 continuation and integration audit
depends_on: [THM-815, THM-823, THM-837]
related: [THM-820, THM-822, THM-842, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_five_all_order_three_context_closure_codex_S10.py
  - 05-knowledge/results/lrc13_hamming_five_all_order_three_context_closure_codex_S10.out
---

# THM-844 — all order-three Hamming-five contexts close

Put

```text
delta=1/13,                 H={1,5,8,12} subset F_13^*.
```

The three multiplicative cosets of `H` partition `F_13^*`.  Choose one of
them,

```text
C=aH,
```

choose a forward flag `b in 2C`, and put

```text
R=C union {b},             P=[12] minus R.                 (1)
```

Give the two opposite pairs in `C` units in `{1,2}` which agree on each
pair, and independently give `b` a unit in `{1,2}`:

```text
e_r=e_(-r) on C,           e_b arbitrary.                  (2)
```

For every `r in R`, let `u_r` be an arbitrary positive member of its labelled
CRT progression

```text
u_r=3r mod 13,             u_r=e_r mod 3.                  (3)
```

## Theorem

Every packet

```text
B=3P union {u_r:r in R}                                    (4)
```

is strictly loose:

```text
M(B)=max_t min_(w in B)||wt|| >1/13.                       (5)
```

There are exactly

```text
3 quartet cosets * 4 forward flags * 4 quartet pair words
                   * 2 fifth units = 96                  (6)
```

contexts in (1)--(3).  Thus (5) closes the complete all-order-three
common-sheet branch isolated by THM-823, at arbitrary lift height.  THM-837
is the single context

```text
C=(1,5,8,12),       b=10,       pair bits=(1,1,1).         (7)
```

## 1. The longest-component bound dominates `K/L`

For a finite speed set `Q`, let

```text
E(Q)={t:||qt||>1/13 for every q in Q}.                     (8)
```

Suppose `E(Q)` is a union of `K` open components, of total measure `L`, and
write `L_max` for the length of a longest component.  If `m<=6` danger combs,
all of speed at least `v`, cover `E(Q)`, then they cover that longest
component.  THM-815's sharp one-interval discrepancy estimate

```text
|I intersect D_u| <=2|I|/13+22/(169u)                    (9)
```

therefore gives

```text
v <= B_long(m,E)
   :=22m/[13(13-2m)L_max].                                (10)
```

THM-837 instead applies the union estimate

```text
v <= B_global(m,E)
   :=22mK/[13(13-2m)L].                                   (11)
```

Since a longest component is at least the average component,

```text
L_max>=L/K,             B_long(m,E)<=B_global(m,E).       (12)
```

Thus the longest-component recursion is uniformly no weaker, with no
context-specific assumption.  It is strictly stronger on all `28,876`
states reached in the exact replay.

For the particular context (7), the root quotient has

```text
K=18,       L=2615/18018,       L_max=47/3003.            (13)
```

The first cap falls from

```text
floor(B_global)=349       to       floor(B_long)=180.      (14)
```

THM-837's global recursion used `75,371` states and reached `57` nonempty
depth-five rows.  The stronger recursion proves the same context with `213`
states and dies before depth five.

## 2. Exact recursive completeness

Numerically order the five replacements in a hypothetical covering row:

```text
v_1<v_2<v_3<v_4<v_5.                                     (15)
```

They are distinct because their labelled residues modulo `13` are distinct.
After a prefix, the exact state consists of

```text
(active open endpoint runs,
 remaining labelled CRT progressions,
 last chosen speed).                                      (16)
```

If the residual is nonempty and `m` labels remain, (10) forces the next speed
to be at most `floor(B_long)`.  For every unused label, the recursion
enumerates every member of its progression (3) which is larger than the last
speed and no larger than this cap.  Hence it includes the next speed of every
hypothetical covering row.  A state with no candidate cannot extend to a
cover, while an empty residual at any depth is exactly a covering prefix.
This proves finite recursive completeness independently of the computed
outcome.

Every residual is invariant under `t -> 1-t`.  The seven-speed retained core
contains an even speed, so `t=1/2` is not safe and no component crosses the
reflection fixed point.  The replay works on `(0,1/2)`.  This halves both `K`
and `L`, preserves `L_max`, and leaves (10)--(12) exact.

The safe bands of speed `u` are represented literally by the open intervals

```text
((13k+1)/(13u),(13(k+1)-1)/(13u)),       0<=k<u,           (17)
```

and every intersection is an exact two-pointer merge over
`fractions.Fraction` endpoints.  A zero-length endpoint contact is never
retained as open safe interior.

## 3. Exact all-context outcome

The complete state census is

| prefix depth | states | dead: no candidate |
|---:|---:|---:|
| 0 | 96 | 0 |
| 1 | 2,496 | 0 |
| 2 | 20,660 | 17,810 |
| 3 | 5,351 | 5,125 |
| 4 | 273 | 273 |

Thus

```text
total states                         28,876
covering prefixes                         0
depth-five terminal rows                  0
distinct replacement operators          252
minimum states in one context            132
maximum states in one context            648
largest longest-component cap            472
largest companion global cap             799.             (18)
```

Every one of the `28,876` state-wise comparisons in (12) is strict.  Since a
hypothetical tight packet would give a covering path and the exhaustive tree
has none, (5) follows uniformly for all heights.

The output emits all `96` context rows, including their CRT bases, root
geometry, node counts, cap maxima, and a SHA-256 digest of every recursive
state word.  Their combined certificate is

```text
9bb065fe6ffafc348ada19c5f3ef30f0b28a5066bfd905a72107202ed13cdb2a.
```

## 4. Tournament Analysis and the faithful carrier

For telemetry only, put the five labelled least-CRT comb obligations at
tournament vertices.  Use

- root total-measure marginal erosion as the pairwise observable;
- root longest-component marginal erosion as the switch/gauge; and
- increasing `(least CRT speed,label)` as the tie Hamiltonian path.

Both gauges are scalar total orders.  In every context they therefore have
score histogram `(0,1,2,3,4)`, no directed triangle, five singleton SCCs,
and one Hamiltonian path.  Nevertheless they flip `492` edges in total, with
histogram

```text
flips per context: {0:1,1:2,2:5,3:8,4:21,5:16,6:18,7:18,8:5,9:2}. (19)
```

The tournament preserves a planning order and nothing theorem-bearing.  It
destroys absolute residual geometry, component incidence, future progression
action, and cover truth.  The challenged assumption is that combs, runners,
or residues alone should be the vertices.  The faithful state is instead the
evolving bipartite incidence

```text
(literal residual components) <-> (remaining labelled comb obligations), (20)
```

decorated by active endpoints and the last speed.

## 5. Verification and scope

Run

```bash
python3 04-computation/lrc13_hamming_five_all_order_three_context_closure_codex_S10.py
```

and compare byte-for-byte with the stored output.  The frozen hashes are

```text
source       627822ffa1b284e96e380b7478d135e21d6197ec54f623f4b6ff1b03e6880b9f
output       6b3f8cb3673d6beb140484a9c2d5de670a49910b454e6628f04650e5dcab3ac8
payload      f0cec95795494fff6dc88c7650504cbb92d5097dd0280dfa3f098e76d3d9ba17
certificate  9bb065fe6ffafc348ada19c5f3ef30f0b28a5066bfd905a72107202ed13cdb2a
```

The replay uses no optimizer, floating point, sampled time grid, or height
cutoff.  The theorem closes all ninety-six all-order-three contexts and
supersedes THM-837's statement that the other ninety-five are open.  It does
not close THM-823's mixed order-one/order-three branch, scalar rows of other
orders, or the global Hamming-five problem.  It does not prove the `n=12`
sporadic branch empty.

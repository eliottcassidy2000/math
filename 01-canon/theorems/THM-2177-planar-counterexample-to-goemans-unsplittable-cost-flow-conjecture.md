---
id: THM-2177
title: "A planar counterexample to Goemans' unsplittable cost-flow conjecture"
status: >-
  PROVED + VERIFIED-EXACT. A seven-vertex, nine-arc acyclic single-source
  instance has a feasible fractional flow of cost 58, but every unsplittable
  routing whose load on each arc is at most the fractional load plus the
  maximum demand has cost at least 60. The three zero-cost path choices have
  incompatibility graph K3; their fractional marginals sum to 16/15 and
  violate the triangle stable-set inequality. The underlying undirected graph
  is a subdivision of K4 and hence planar. A nonempty open three-parameter
  region explains the mechanism. An independent graph-level checker
  enumerates all simple paths and all eight routings under normal and
  optimized Python.
source: codex-2026-07-24-DGG-eight-routing-audit
depends_on: []
related:
  - THM-2176
script: 04-computation/goemans_unsplittable_cost_flow_counterexample_thm2177.py
output: 05-knowledge/results/goemans_unsplittable_cost_flow_counterexample_thm2177.out
script_sha256: 1552accd54cae794acd97660ed3c6f375fb0b9920d29bb14a6020a519b244cf5
output_sha256: 33558dc8dce3ff4f42f40188d8b78d06531934da9c78f2248f57372c5bf15ff8
hash_basis: working-tree bytes (LF)
---

# THM-2177 -- a planar counterexample to Goemans' cost conjecture

## 1. The conjecture being refuted

In a single-source unsplittable-flow instance, a source `s` must send demand
`d_i` to terminal `t_i`. A fractional flow `x` may split each demand among
several `s`--`t_i` paths. An unsplittable routing `P` chooses one path for
each terminal and has arc loads

```text
flow_P(a)=sum_(i:a in P_i) d_i.                       (1)
```

Put

```text
D=max_i d_i.                                          (2)
```

Goemans' cost conjecture, stated for example as Conjecture 1.3 in
[Traub--Vargas Koch--Zenklusen](https://arxiv.org/abs/2308.02651), asserts
that for every feasible `x` and every nonnegative per-unit arc-cost vector
`c`, there is an unsplittable routing `P` satisfying both

```text
flow_P(a)<=x(a)+D                   for every arc a,   (3)
sum_a c(a)flow_P(a)<=sum_a c(a)x(a).                  (4)
```

The original theorem of Dinitz--Garg--Goemans proves (3) without (4).
The planar theorem in the cited paper proves the cost conclusion with the
larger additive allowance `2D`. The construction below refutes precisely the
simultaneous `D`-allowance statement (3)--(4).

## 2. The seven-vertex instance

Let

```text
V={s,u,v,w,t_1,t_2,t_3},
(d_1,d_2,d_3)=(15,10,15),             D=15.           (5)
```

The arcs, fractional loads, and nonnegative per-unit costs are

| arc `a` | `x(a)` | `c(a)` |
|---|---:|---:|
| `s -> t_1` | 10 | 2 |
| `s -> t_2` | 6 | 3 |
| `s -> u` | 24 | 0 |
| `u -> t_3` | 10 | 2 |
| `u -> v` | 14 | 0 |
| `v -> t_1` | 5 | 0 |
| `v -> w` | 9 | 0 |
| `w -> t_2` | 4 | 0 |
| `w -> t_3` | 5 | 0 |

One may set each arc capacity equal to `x(a)`, so `x` is feasible. Indeed,
its source outflow is

```text
10+6+24=40=d_1+d_2+d_3,                              (6)
```

the internal conservation identities are

```text
24=10+14,          14=5+9,          9=4+5,           (7)
```

and the terminal inflows are `10+5=15`, `6+4=10`, and
`10+5=15`.

The graph is acyclic. Its underlying undirected graph is the subdivision of
`K_4` on branch vertices `s,u,v,w` in which `t_1,t_2,t_3` subdivide the
edges `sv,sw,uw`, respectively. In particular, the counterexample is planar.

## 3. Complete path universe

Each terminal has exactly two source paths:

```text
E_1: s -> t_1,                  Z_1: s -> u -> v -> t_1,
E_2: s -> t_2,                  Z_2: s -> u -> v -> w -> t_2,
E_3: s -> u -> t_3,             Z_3: s -> u -> v -> w -> t_3.  (8)
```

There are no other directed source--terminal paths. Every `Z_i` has cost
zero, while each `E_i`, when it carries its whole terminal demand, costs
exactly `30`.

The displayed fractional flow is equivalently the path decomposition

```text
t_1: 10 on E_1 and 5 on Z_1,
t_2:  6 on E_2 and 4 on Z_2,
t_3: 10 on E_3 and 5 on Z_3.                         (9)
```

Therefore its total cost is

```text
c^T x=2*10+3*6+2*10=58.                              (10)
```

## 4. The incompatibility triangle

No routing satisfying (3) can choose two of the zero-cost paths.

- If it chooses `Z_2,Z_3`, arc `v -> w` has load `10+15=25`, while
  `x(vw)+D=9+15=24`.
- If it chooses `Z_1,Z_3`, arc `u -> v` has load at least `15+15=30`, while
  `x(uv)+D=14+15=29`.
- If it chooses `Z_1,Z_2`, terminal `t_3` uses `s -> u` on either of its two
  paths, so that arc has load `15+10+15=40`, while
  `x(su)+D=24+15=39`.

Conversely, every routing choosing at most one `Z_i` satisfies (3). On the
three shared arcs its loads are bounded by

```text
flow_P(su)<=30<39,
flow_P(uv)<=15<29,
flow_P(vw)<=15<24,                                   (11)
```

and every other arc is used by at most one demand, hence has load at most
`D<=x(a)+D`.

Thus the capacity-good routings are exactly the independent sets of size at
most one in the symmetric incompatibility graph

```text
Z_1 -- Z_2
 |       |
 +--Z_3--+                                              (12)
```

on the selected zero-cost paths. Every such routing uses at least two
`E_i`, so it costs at least

```text
2*30=60>58=c^T x.                                    (13)
```

Equations (3) and (4) can therefore never hold simultaneously. This refutes
Goemans' conjecture. QED.

## 5. Stable-set certificate

The fractional path decomposition (9) chooses the three cheap paths with
marginals

```text
z=(5/15,4/10,5/15)=(1/3,2/5,1/3).                  (14)
```

Every capacity-good integral routing obeys the triangle stable-set facet

```text
z_1+z_2+z_3<=1,                                      (15)
```

but the fractional point has

```text
1/3+2/5+1/3=16/15>1.                                 (16)
```

The arc costs are the nonnegative complementary separator: choosing `E_i`
instead of `Z_i` costs `30`. This is the structural reason for (13).

This pairwise object is not a tournament. Incompatibility is symmetric, and
orienting its three edges would forget the independent-set predicate that
does the proof. The faithful carrier is the conflict graph together with the
fractional marginals and complementary costs.

## 6. A nonempty open region of counterexamples

The construction is not an isolated integer coincidence. Keep the graph and
paths (8), normalize the demands to

```text
d_1=d_3=1,             d_2=b,             0<b<=1,    (17)
```

and split them as

```text
t_1: 1-r on E_1,       r on Z_1,
t_2: b(1-q) on E_2,   bq on Z_2,
t_3: 1-r on E_3,       r on Z_3,                         (18)
```

where `0<r,q<1`. Give the unique positive-cost arc of `E_1,E_2,E_3`
per-unit costs `1,1/b,1`, respectively, and give every other arc cost zero.
Each integral `E_i` then costs one, while the fractional cost is

```text
3-(2r+q).                                             (19)
```

Suppose

```text
2r+q>1,
b(1-q)>r,
2r+bq<1.                                             (20)
```

The second inequality makes `Z_1,Z_2` incompatible on `su` and
`Z_2,Z_3` incompatible on `vw`; the third makes `Z_1,Z_3` incompatible on
`uv`. The same direct shared-arc check as in Section 4 proves the converse:
every selection with at most one `Z_i` is capacity-good. Hence every good
routing costs at least `2`, while (19) and the first inequality give

```text
3-(2r+q)<2.                                          (21)
```

Every parameter point satisfying (20) is therefore a counterexample. These
strict inequalities contain a nonempty open subset of the parameter domain:
take

```text
(b,r,q)=(2/3,1/3,2/5).                               (22)
```

Scaling the demands by `15` and the costs by `2` gives exactly Sections
2--4.

## 7. Exact audit and scope

The companion begins with the nine graph arcs rather than the six claimed
paths. It:

1. verifies flow conservation and nonnegative costs;
2. enumerates all simple source--terminal paths and recovers exactly (8);
3. exhausts all `2^3=8` unsplittable routings;
4. recomputes every arc load, additive-capacity violation, and cost;
5. independently checks the stable-set separator (16); and
6. checks the `K_4`-subdivision planarity certificate.

Normal and optimized Python give the frozen transcript. The audit also
catches a harmless detail hidden by a coarse routing table: the all-`Z`
routing exceeds the `uv` allowance by `11`, not merely by `1`.

The all-`Z` routing has deviations from `x` bounded in absolute value by
`26<2D=30` and cost zero, so it satisfies both sides of the planar theorem's
cost-preserving `2D` conclusion. Thus this result does not contradict that
theorem. It refutes the exact `D`-allowance cost conjecture only.

### Provenance

The integer construction appeared in the public ChatGPT conversation supplied
with this session. The proof above was independently reconstructed from the
graph, checked against the primary statement of Conjecture 1.3, generalized
to (20), and verified by a new graph-level exhaustive program. As of the
source audit on 2026-07-24, the construction had not yet appeared in a
peer-reviewed paper; this theorem makes no publication-priority claim. QED.

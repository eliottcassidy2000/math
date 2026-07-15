---
id: THM-845
title: Scale-one Hamming-five row-wise safe-component closure
status: PROVED (exact row-wise THM-815 recursion inside the exhaustive THM-820 dichotomy) + FINITE-EXACT (772,543-row state certificate, independent open-safe/closed-danger full-tree replays, nine Fraction-exact terminal rows, and Tournament Analysis)
source: codex-2026-07-15-S10 continuation
depends_on: [LRC(<=13), THM-815, THM-820]
related: [THM-810, THM-816, THM-822, THM-823, THM-837, THM-840, HYP-6820]
verification:
  - 04-computation/lrc13_hamming_five_rowwise_component_closure_codex_S10.cpp
  - 05-knowledge/results/lrc13_hamming_five_rowwise_component_closure_codex_S10.out
  - 04-computation/lrc13_hamming_five_closed_danger_union_replay_codex_S10.cpp
  - 05-knowledge/results/lrc13_hamming_five_closed_danger_union_replay_codex_S10.out
  - 04-computation/lrc13_hamming_five_terminal_endpoint_crosscheck_codex_S10.py
  - 05-knowledge/results/lrc13_hamming_five_terminal_endpoint_crosscheck_codex_S10.out
---

# THM-845 — scale-one Hamming-five row-wise safe-component closure

Put

```text
delta=1/13,                         [12]={1,...,12},
M(W)=max_(t in R/Z) min_(w in W)||wt||.
```

## Theorem

Let `R subset [12]` have five elements and choose a proper positive lift

```text
u_r=r+13h_r,                        h_r>=1
```

for every `r in R`.  Then

```text
B=([12] minus R) union {u_r:r in R}
```

is loose:

```text
M(B)>1/13.                                                   (1)
```

Thus the complete proper scale-one, residue-preserving Hamming-five chart
about `[12]` is empty.  This is a uniform arbitrary-height statement, not a
bounded lift census.

## 1. The exact row-wise recursion

Write `P=[12] minus R` and order the five lifted speeds numerically:

```text
x_1<x_2<x_3<x_4<x_5.                                       (2)
```

For a prefix

```text
Q_j=P union {x_1,...,x_j},             0<=j<=4,
```

let `E_j` be its strict `delta`-safe set and let `L_j` be the length of
a longest component of `E_j`.  Every `Q_j` has at most eleven speeds, so
settled `LRC(<=13)` gives `M(Q_j)>=1/12>1/13`; hence `L_j>0`.

If a completion of the prefix were tight, its `m=5-j` remaining danger
combs would cover a component of length `L_j`.  THM-815 (8) therefore gives
the exact necessary next-speed cap

```text
x_(j+1)<=D(j,L_j)
 :=floor(22(5-j)/[13(13-2(5-j))L_j]).                      (3)
```

The verifier represents the strict-safe bands of one speed `u` as the
ordered open rational intervals

```text
((13k+1)/(13u),(13(k+1)-1)/(13u)),          0<=k<u,         (4)
```

and obtains every `E_j` by exact two-pointer intersection.  It retains an
intersection only when `lo<hi`, so a common closed endpoint is never
misread as strict-safe interior.  Formula (3) is then evaluated separately on
every labelled prefix, rather than replacing the row by one global worst-case
box.

## 2. Exhaustiveness of the finite tree

THM-815 first gives `x_1<=146`.  THM-820 says that every hypothetical tight
row lies in exactly one of the following branches.

- **A — recursive doubling:**

  ```text
  x_2<=2x_1, x_3<=2x_2, x_4<=2x_3, x_5<=2x_4.              (5)
  ```

- **B — exceptional top-four cycle:** if `x_2>2x_1`, the labels of
  `x_2,...,x_5` are a multiplicative translate of `{1,2,4,8}` modulo
  `13`.  With `v=x_2`, `x=x_1`, and `m=max(P)`,

  ```text
  v<=floor(819x/40),
  v<=floor((7/2)/(15/(104m)-1/x)) when the denominator is positive,
  x_k<=4v for k=2,...,5.                                   (6)
  ```

The replay enumerates each five-label set, each possible owner of the next
increasing lift, and every positive speed in its residue class satisfying
(3), together with (5) or (6).  Distinct labels modulo `13` make the five
lifted speeds distinct, so numerical ordering gives every row one and only
one path.  No state is discarded by a heuristic or a floating-point test.
Consequently every tight row would appear as an empty strict-safe set at depth
five of this tree.

The exact state counts are:

| prefix depth | states | minimum longest component | largest next cap |
|---:|---:|---:|---:|
| 1 | 40,590 | `11/1898` | 233 |
| 2 | 612,221 | `11/3029` | 199 |
| 3 | 111,675 | `162/67379` | 156 |
| 4 | 7,255 | `307/129987` | 65 |
| 5 | 9 | `1/286` | terminal |

In branch B the depth-two through depth-five counts are respectively

```text
415, 178, 1, 0.                                             (7)
```

Thus the exceptional branch dies before a full packet.  Branch A leaves nine
terminal rows, and every one has a nonempty strict-safe interval.  The primary
certificate contains 772,543 deterministic state rows and has SHA-256

```text
6524ac6dd2d1f8c59256816c86b95d9ee52cc94766d4d3f993425e7071434a29.
```

This proves (1).

## 3. The nine terminal cross-checks

For redundancy, the primary C++ replay also computes the exact maximin of
every surviving row from all self-cusp and signed pair-crossing denominators.
A separate Fraction-exact Python implementation reconstructs the packets,
strict-safe components, complete maximizer sets, and the same maxima:

| ordered replacement word `owner:speed` | one strict-safe interval | exact `M` |
|---|---:|---:|
| `1:14,4:17,6:19,10:23,3:29` | `(79/247,38/117)` | `3/28` |
| `1:14,3:16,5:18,8:21,9:22` | `(9/26,32/91)` | `2/17` |
| `1:14,3:16,5:18,9:22,11:24` | `(9/26,32/91)` | `1/10` |
| `3:16,4:17,10:23,12:25,2:28` | `(14/143,4/39)` | `1/10` |
| `10:23,2:28,3:29,6:32,8:34` | `(27/91,90/299)` | `2/15` |
| `3:16,4:17,5:18,6:19,8:21` | `(40/117,77/221)` | `3/26` |
| `4:17,6:19,3:29,5:31,8:34` | `(27/221,51/403)` | `1/8` |
| `4:17,6:19,8:21,10:23,3:29` | `(79/247,38/117)` | `1/8` |
| `4:17,5:18,6:19,12:25,9:35` | `(27/143,5/26)` | `4/37` |

Every value is strictly above `1/13`.  The primary C++ replay is byte-identical
at `-O3` and `-O0` and passes ASan/UBSan.  A structurally independent C++
replay unions closed danger teeth from scratch at every prefix; it produces
the identical 772,543-row certificate and also passes UBSan.  The frozen
digests are

```text
primary source       3ac19e586ca7b86e0b3f982e10cb901077f9504e0279615d4d86fa6526e75417
primary output       4aaf3bf6191699e98b8824f53c6a9bd5a0424f77d553734a71a32b0f401a9d5f
state certificate    6524ac6dd2d1f8c59256816c86b95d9ee52cc94766d4d3f993425e7071434a29
danger-union source  9913f7623ef6a8918bab64e0e7dc531001a20a63361733bf36e799744a5f82e0
danger-union output  1b49d00f21018e9d9e9e88ce70110239ef4bafc04e975efad631d23dc49ba298
cross-check source   deebad6c21dc8aaaffd7a31126d1bb3d83964581b00872cc0bd3706005797c09
cross-check output   ea62afb0e126e45cd770180037171fd97c516b10475d2c3dc9b40894f8289cf1
cross-check rows     05c2dac64e5197b3926d74d9fca27b0e4adbbe2309be44178608091abfe13182
```

## 4. Tournament Analysis and the faithful carrier

The theorem-bearing vertices are not runners.  At a recursive state the
faithful carrier is bipartite/hypergraphic:

```text
(strict-safe residual components)  <->  (remaining labelled danger combs).
```

It preserves the exact obligation “every residual component is covered” and
the next operation, while retaining component endpoints and owner labels.  A
runner or residue tournament destroys higher-order overlap, component
splitting, and therefore cover truth.

For telemetry on each of the nine terminal rows, take the five remaining comb
obligations as vertices.

- **Pairwise observable:** raw erosion of the current longest core-safe
  component.
- **Switch/gauge:** conditional marginal erosion after the other four combs
  have been applied.
- **Tie Hamiltonian path:** increasing `(speed,label)`.

Both gauges produce transitive five-vertex tournaments on every row, with
score histogram `0,1,2,3,4`, no directed cycle, five singleton SCCs, and one
Hamiltonian path.  They nevertheless flip 31 edges in total.  This large
planning-order change, despite identical trivial fingerprints, is another
exact warning that binary rankings are telemetry.  The literal residual
interval word proves the theorem.

The challenged assumption is that replacement labels—or runners themselves—
must be the vertices.  Operation obligations are closer, but even they need
the component-incidence sidecar.  This is consistent with THM-822's kernel
liars: coarse live relations mix exact `M`, whereas literal endpoint faces
retain the continuation state.

## 5. Scope and next boundary

THM-845 closes the **scale-one** Hamming-five chart at arbitrary lift height.
It does not by itself classify a five-coordinate packet at arbitrary AP scale:
the ramified deck interfaces outside common scale still require descent.

Scale-one Hamming six remains finitely decidable by THM-815, but it contains
the genuine tight scale-change face `2[12]`; the next computation must first
separate that AP orbit from non-AP rows.  Radius seven remains the first
scale-one chart where the present discrepancy coefficient `13-2m` is
nonpositive.  The correct next state there must use overlap debt, owner
diversity, or component/comb incidence rather than mean danger density.

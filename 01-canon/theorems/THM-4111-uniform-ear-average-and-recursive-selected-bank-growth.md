---
id: THM-4111
title: "Uniform ear average and recursive selected-bank growth"
status: >
  PROVED ELEMENTARY DOUBLE COUNT + VERIFIED-EXACT + INDEPENDENTLY
  VERIFIED-EXACT. The exact all-cut sum is
  2^(n-2)((n+3)H(T)+F_1(T)), where F_1 counts base orderings with exactly one
  bad adjacency. Its nonconstant-cut specialization forces every full-cut,
  recursively representative-selected strong-ear bank to have unbounded image
  maxima. This does not prove interval overlap, unbounded endpoints of solid
  intervals, or global H-spectrum completeness.
source: codex-frontier-synthesis-creative-20260825g
depends_on:
  - THM-001-redei
related:
  - THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension
  - THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary
  - THM-4102-selected-order-ten-strong-ear-solid-interval
  - THM-4104-selected-order-eleven-strong-ear-solid-interval
  - THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature
  - HYP-2879-strong-ear-atom-calculus
  - HYP-9029-strong-interval-tiling-law
  - MISTAKE-402
script: 04-computation/tournament_uniform_ear_average_growth_thm4111.py
output: 05-knowledge/results/tournament_uniform_ear_average_growth_thm4111.out
independent_audit_script: 04-computation/tournament_uniform_ear_average_growth_thm4111_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_uniform_ear_average_growth_thm4111_independent_audit.out
script_sha256: f58f2a10d70685cc998a14d1ee6d4381e4e1c0d60e41aee3e7f89ca80570fe36
output_sha256: 708d55d7ac5f5d47441314d393d11d727324a3f20eca815678c0d2935b2973a9
independent_audit_script_sha256: c7c59e8eb98d3bfdcb2105802a7c7ffda2194e28a83e30c227a38d5d464e7ba7
independent_audit_output_sha256: 62be3fe1e66ca757f6d04ceb0f466e07502a4e1bbfa9c48172dc52baa4aa6c59
semantic_sha256: 0bb800d8caa1c1fd449657fb9f68a33842062ec4f3866f34488d9c8ea3251915
independent_semantic_sha256: 1580575c67bc5b8cdea70db9c7efef256cad50b8a5ed3ce18e7a73303dfad9f0
hash_basis: raw LF bytes for files; canonical compact JSON for the semantic ledger
audit: >
  PASS. The independent path imports no primary code and literally scans
  124,468 parent orderings plus 23,717,424 child orderings over all 1,098
  labelled parents and 33,864 ears through order five. It reproduces the
  exact sums, means, strong-ear controls, C3 coefficient hostile, inherited
  selected-bank bounds, and the equal-(H,F_1)-but-different-image hostile.
  After CRLF-to-LF normalization the independent normal/-O streams match their
  frozen output; the LF-stable primary streams byte-match theirs directly.
---

# THM-4111 -- uniform ear average and recursive selected-bank growth

**PROVED ELEMENTARY DOUBLE COUNT + VERIFIED-EXACT + INDEPENDENTLY
VERIFIED-EXACT.**

This theorem isolates the part of THM-4104's open recursive-bank question
that does not require another spectrum census. The averaging mechanism proves
growth of the largest selected-image value by averaging the full labelled cut
fibre. It does **not** prove that the long solid intervals overlap forever or
that their right endpoints are unbounded.

## 1. Zero- and one-defect deletion layers

Let `T` be a tournament on a vertex set `V` of order `n>=2`. Write

```text
H(T)=#{directed Hamiltonian paths of T}.                       (1)
```

For an ordering `P=(v_1,...,v_n)` of `V`, put

```text
delta_T(P)=#{i:v_(i+1)->v_i},
F_1(T)=#{P:delta_T(P)=1}.                                    (2)
```

Thus `F_1(T)` is exactly THM-4099's one-defect word layer before a
new-vertex cut decides whether the unique defect is repaired.

Adjoin a labelled vertex `x`. For every `S subseteq V`, let `T+x_S` have the
old arcs of `T` and cut convention

```text
S={v:x->v}.                                                   (3)
```

> **Theorem 1.1 (uniform cut-ear sum).** For every tournament `T`
> of order `n>=2`,
>
> ```text
> sum_(S subseteq V) H(T+x_S)
>   =2^(n-2)((n+3)H(T)+F_1(T)).                              (4)
> ```

### Proof

Sum the Hamiltonian paths of `T+x_S` over all `2^n` cuts and partition them by
the position of `x`.

If `x` is first, deleting it leaves a Hamiltonian path of `T`; its first old
vertex must lie in `S`, leaving `2^(n-1)` choices of `S`. The last-position
case is symmetric and also contributes `2^(n-1)H(T)`. The two endpoint cases
therefore contribute `2^n H(T)`.

If `x` is internal, deleting it exposes one adjacent pair of an ordering `P`
of `V`; every other adjacency of `P` must be an arc of `T`. There are exactly
two possibilities:

1. `P` is a Hamiltonian path, and any of its `n-1` positions can be exposed;
2. `P` lies in `F_1(T)`, and only its unique bad position can be exposed.

For a fixed exposed ordered pair `(a,b)`, the path requires `a->x->b`, hence
`a notin S` and `b in S`. The other `n-2` membership bits are free, so every
such exposed word occurs for exactly `2^(n-2)` cuts. The internal contribution
is consequently

```text
2^(n-2)((n-1)H(T)+F_1(T)).                                  (5)
```

Adding the endpoint term `2^n H(T)` gives `(4)`. The position of `x`, the
deleted word, and its exposed position recover every counted object, so no two
parts overlap. **QED.**

The same calculation can be read through THM-4097's `(Start,End,Q)` boundary:
the endpoint expectation is `H(T)`, while the symmetric total of `Q` is

```text
sum_(a!=b)Q_T(a,b)=(n-1)H(T)+F_1(T).                         (6)
```

This identifies the one-defect squarefree layer of THM-4099 as the exact
surplus in the uniform ear average.

## 2. Remove the two constant cuts

For `S=empty` or `S=V`, the new vertex is respectively a sink or source and

```text
H(T+x_S)=H(T).                                               (7)
```

Subtracting these two terms from `(4)` gives the exact nonconstant mean

```text
A_nc(T)=
 [ (2^n(n+3)-8)H(T)+2^n F_1(T) ] / [4(2^n-2)].              (8)
```

Since `n>=2`, the all-cut mean is

```text
((n+3)H(T)+F_1(T))/4 > H(T),                                (9)
```

so removing the two constant values only increases the mean. Therefore

```text
max_(empty!=S!=V) H(T+x_S) >= (n+3)H(T)/4.                  (10)
```

Every Hamiltonian-path count is odd by THM-001, so if `oddceil(q)` denotes
the least odd integer at least `q`, `(8)` also gives

```text
max_(empty!=S!=V) H(T+x_S) >= oddceil(A_nc(T))
                              >= oddceil((n+3)H(T)/4).       (11)
```

## 3. Selection-robust recursive growth

Start at any order `n_0>=3` with a nonempty finite bank `B_(n_0)` of
order-`n_0` strong tournaments. Recursively form the full nonconstant-cut
value image

```text
V_(n+1)={H(T+x_S):T in B_n, empty!=S!=V(T)},                (12)
```

and let `B_(n+1)` contain an arbitrary one order-`n+1` strong witness for each
value in `V_(n+1)`. At every step all nonconstant cuts from every retained
parent are included before one representative per scalar value is selected.
This includes deterministic first-labelled-witness selection but does not
depend on that choice. Nonconstant ears over a strong parent are strong: `x`
has both an in-neighbor and an out-neighbor, and strongness of the parent
supplies paths in both directions to every old vertex.

Put

```text
M_n=max{H(T):T in B_n}.                                     (13)
```

Choose any `T_n in B_n` with `H(T_n)=M_n`. Equations `(10)--(12)` imply

```text
M_(n+1) >= (n+3)M_n/4.                                     (14)
```

Iteration gives, for every `m>n_0`,

```text
M_m >= M_(n_0) (m+2)! / ((n_0+2)! 4^(m-n_0)).              (15)
```

The right side tends to infinity. Thus no recursive representative-selection
rule can keep the selected-image maxima bounded.

For the frozen banks in THM-4097/4102/4104, `(11)` gives the theorem-forced
rows

| parent order | parent maximum | forced next maximum | observed next maximum |
|---:|---:|---:|---:|
| `9` | `3,357` | `10,071` | `15,621` |
| `10` | `15,621` | `50,769` | `93,751` |

The excess beyond the forced column is finite-exact information from the
selected banks, not a consequence of the average alone.

## 4. Sharp hostile and scope

Let `T=C3`. Then

```text
H(T)=3, F_1(T)=0,                                           (16)
```

and the six nonconstant ears all have Hamiltonian-path count `5`. Hence the
one-defect surplus can vanish even for a strong parent, and the tempting
uniform strengthening

```text
max H(T+x_S) >= (n+4)H(T)/4                                (17)
```

already fails because `5<21/4`. This is the minimal strong hostile.

There is a separate minimal quotient boundary for interval propagation. In
the LSB-first upper-pair encoding of THM-4097, the strong order-five parents
with codes `1015` and `759` both have

```text
H=9,                 F_1=30,
all-cut sum=816,      nonconstant sum/mean=798 and 133/5.   (18)
```

Nevertheless their nonconstant ear-value sets are respectively

```text
{15,17,19,23,25,27,29,33,37,41},
{15,17,19,23,29,31,33,37,43}.                              (19)
```

Thus even `(H,F_1)` and the exact mean do not determine the cut image. The
first failed implication is

```text
(zero/one-defect totals)  -/->  solid-interval propagation. (20)
```

The full distributional `(Start,End,Q)` sidecar remains load-bearing.

The result proves unbounded **image maxima**. A maximum may be
isolated. It proves none of the following:

- overlap of the recursive solid intervals;
- unbounded right endpoints of those solid intervals;
- a complete order-`n` spectrum at any new order; or
- global H-spectrum completeness.

Those remain **OPEN**.

## 5. Exact controls and reproduction

The primary referee computes `H` by subset DP, enumerates
`F_1`, and checks `(4)`, `(8)`, parity, and strong-ear inheritance on every
labelled tournament through order five:

```text
1,098 parent tournaments,
33,864 cut ears,
16,668 nonconstant ears over strong parents.                (21)
```

It freezes the `C3` hostile and the two inherited selected-bank rows. Run

```bash
python3 -B 04-computation/tournament_uniform_ear_average_growth_thm4111.py
python3 -B -O 04-computation/tournament_uniform_ear_average_growth_thm4111.py
python3 -B 04-computation/tournament_uniform_ear_average_growth_thm4111_independent_audit.py
python3 -B -O 04-computation/tournament_uniform_ear_average_growth_thm4111_independent_audit.py
```

The primary normal/-O pair must byte-match its LF-stable frozen output; the
independent pair matches after CRLF-to-LF normalization under MISTAKE-402.
These finite audits are hostile controls for the elementary double count, not
the source of its quantifiers.

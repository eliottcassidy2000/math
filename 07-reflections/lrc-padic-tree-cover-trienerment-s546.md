---
source: codex-2026-06-01-S546
status: computational probe plus synthesis
tags:
  - lonely-runner
  - p-adic
  - tree
  - trienerment
  - sieve
  - gabor
  - symbolic-dynamics
---

# LRC p-adic Tree-Cover T(r)ienerments

This session continues the Gabor/t(r)ienerment line by moving the vertex set to
the p-adic tree.  The useful object is not only the rank of the tree or the
prime factors of `n/2`; it is the cover core of the observer zero branch.

For a speed set `V` with `|V|=n-1`, define:

```text
z_q(V) = #{v in V : q divides v},  2 <= q <= n.
```

The node `q` is the product p-adic zero branch `0 mod q`.  If `z_q=0` and
`q<n`, then `t=1/q` is an open THM-369 witness.  If `z_n=0`, then `t=1/n` is
the compactified wall witness.

This gives a clean split:

```text
open survivor:
  every q=2,...,n-1 has z_q>0

compact survivor:
  every q=2,...,n has z_q>0

AP:
  open survivor, but z_n=0
```

## Tournament Analysis

The trienerment vertices are:

```text
q = 2,3,...,n.
```

The pairwise observable is:

```text
(z_q, divisibility depth).
```

The switch/gauge is:

```text
lower z_q wins;
if masses tie and one q divides the other, the deeper branch wins;
if masses tie and the branches are incomparable, tie.
```

The tie Hamiltonian path is increasing `q`.

Fingerprints:

```text
tie count,
ternary B1 = 2M - 3*ties,
score histogram,
strict 3-cycles,
SCCs,
Hamiltonian-path count,
minimum cover size,
singleton-carrier profile.
```

This p-adic relation has zero strict directed 3-cycles in the bounded scan.
That is expected: the product tree is hierarchical.  It should not look like
the cyclic phase Gabor relation; it should look like a poset with ties between
incomparable branches.

## Computation

Artifacts:

```text
04-computation/lrc_padic_tree_trienerment_s546.py
05-knowledge/results/lrc_padic_tree_trienerment_s546.out
```

Supplement:

```text
04-computation/lrc_prime_power_zero_branch_core_s546b.py
05-knowledge/results/lrc_prime_power_zero_branch_core_s546b.out
07-reflections/lrc-prime-power-zero-branch-endpoint-core-s546b.md
```

The supplement restricts from all product nodes `q=2..n` to prime powers
`q=p^d` and tests a sharper endpoint-core question.  Covered prime-power
branches kill the unit witnesses `u/q`, but their local centered danger
intervals peel to empty endpoint cores in all audited rows (`86` local branch
cores, `0` nonempty).

Bounded primitive boxes:

```text
n=6,  max_speed=8:   56 primitive sets
n=8,  max_speed=10: 120 primitive sets
n=10, max_speed=12: 220 primitive sets
n=12, max_speed=14: 364 primitive sets
n=14, max_speed=16: 560 primitive sets
n=16, max_speed=18: 816 primitive sets
n=18, max_speed=20: 1140 primitive sets
```

Main table:

```text
n   open survivors  compact survivors  min cover among open survivors
6        25                 16                         3
8        44                 30                         4
10       70                 50                         5
12      104                 77                         6
14      147                112                         7
16      200                156                         8
18      264                210                         9
```

The minimum open cover size is always `n/2` among open survivors in this
bounded family.

Tree-trienerment class purity:

```text
n=6:  total classes 20,  open-survivor classes 8,  mixed 0
n=8:  total classes 45,  open-survivor classes 14, mixed 0
n=10: total classes 83,  open-survivor classes 26, mixed 0
n=12: total classes 104, open-survivor classes 25, mixed 0
n=14: total classes 154, open-survivor classes 40, mixed 0
n=16: total classes 201, open-survivor classes 46, mixed 0
n=18: total classes 259, open-survivor classes 63, mixed 0
```

This is the crispest computational signal: once the fingerprint includes the
zero-branch cover-core features, there are no mixed open-survivor fibers in the
scan.

## AP As Wall Survivor

The initial segment has:

```text
open_survivor=True
compact_survivor=False
compact_empty=(n,)
singleton nodes=(n/2,n/2+1,...,n-1)
min_cover=n/2
```

Examples:

```text
n=14 AP: compact_empty=(14,), singleton nodes=(7,8,9,10,11,12,13)
n=16 AP: compact_empty=(16,), singleton nodes=(8,9,10,11,12,13,14,15)
n=18 AP: compact_empty=(18,), singleton nodes=(9,10,11,12,13,14,15,16,17)
```

So AP is not mysterious in this coordinate.  It covers every open product-tree
zero branch, then leaves exactly the wall branch `q=n` empty.  This matches the
symbolic-dynamics compactification story: the open word has no target, but the
wall supplies one.

## Formalization: THM-390

The proved core is now split off as THM-390.  The theorem says:

```text
empty q<n  =>  t=1/q is an open lonely witness
empty q=n  =>  t=1/n is a compactified wall witness
AP          =>  covers every q<n and misses only q=n
AP core     =>  unique minimum open cover U_n={u: 2u >= n}, size floor(n/2)
```

Thus the tree-cover scan is using proof-obligation vertices, not just a
decorative p-adic metaphor.  The quotient preserves exactly the explicit
denominator-witness part of LRC and records the singleton carriers forced by
the AP wall.

The still-open part is sharper after this separation: the zero-mixed
tree-trienerment fingerprint fibers and the all-open-survivors `n/2` cover law
are computational phenomena beyond the AP cover-core theorem.

New incoming entropy work (HYP-2037/HYP-2038) fits here without replacing the
cover-core object.  Entropy reads AP/regular rows as high-spread or critical
boundary states; THM-390 identifies the forced finite leaves on that boundary.
The resulting split is useful:

```text
channel/tree entropy  = spread or criticality of the row
zero-branch mass z_q  = exact local sieve gate at denominator q
singleton AP leaves   = forced wall support where entropy collapses to the
                        finite critical set
```

So the entropy attack should use `z_q=0` and singleton carriers as local
certificates, not replace them with a scalar spread statistic.

HYP-2039's defect-transport frame says the large hole always exists and the
problem is moving it to the observer.  The zero-branch theorem gives the pinned
transport cases: when `z_q=0`, the observer hole is already present at
`t=1/q`; when every branch is covered, transport has to pass through the
labelled cover-core machinery rather than a bare spread argument.

## Relation To The Existing Threads

HYP-2032/HYP-2035 identify the p-adic tree as the arithmetic metric and channel
rank.  HYP-2036 adds the cover-core layer: for each zero branch, which speeds
actually carry the debt?

HYP-2031/HYP-2032 on Gabor trienerments say the observer target should become a
marked zero column.  The p-adic node `q` supplies exact candidate times:

```text
q<n:
  open Gabor observer-zero event at t=1/q

q=n:
  compactified wall Gabor event
```

HYP-2029 says the event word needs consistency labels.  The p-adic label is now
obvious:

```text
first empty q,
or singleton carrier of q when z_q=1.
```

That label records whether a bad symbolic branch is impossible by sieve, barely
covered by one runner, or robustly covered.

## Assumption Challenge

Considered vertex sets:

```text
runners, gaps, fixed sections, section boundaries,
wall-crossing events, residues, cover arcs, Fourier/Gabor cells,
p-adic residue balls, composite CRT product nodes, proof obligations.
```

Chosen quotient:

```text
vertices = product p-adic zero-branch obligations q=2..n.
```

Predicate preserved:

```text
existence of a THM-369 open sieve witness, and the q=n wall witness.
```

Information destroyed:

```text
non-sieve witnesses, exact event order, circular gap shape,
and runner identity except for singleton carrier roles.
```

Challenged assumption:

```text
the p-adic tree only records channel rank.
```

The scan says it records more: a pure cover-core fingerprint and a forced
singleton-carrier ledger.

## Next

The next move is not another metaphorical tree pass.  It should test whether
the cover-core purity survives outside the `max_speed=n+2` box and whether the
`min_cover=n/2` law has a direct proof from boundary leaves.

Concrete tests:

```text
1. Extend to n=20,22 with optimized/cached Hamiltonian counts or no HP counts.
2. Add hard n=14/n=18 rows from prior scalar-puncture and gate families.
3. Decorate HYP-2029 event words by first empty q and singleton carrier.
4. Add Gabor zero-column labels at t=1/q for q in the product tree.
```

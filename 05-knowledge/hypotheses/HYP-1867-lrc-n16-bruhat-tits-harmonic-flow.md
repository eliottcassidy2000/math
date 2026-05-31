---
id: HYP-1867
status: OPEN
source: codex-2026-05-31-S399
related:
  - HYP-1856
  - HYP-1858
  - HYP-1859
  - HYP-1860
  - HYP-1861
  - HYP-1865
  - HYP-1866
---

# HYP-1867: n=16 endpoint debt is a Bruhat-Tits harmonic-flow obstruction

## Statement

At `n=16`, the pure dyadic endpoint-protection problem should be expressible as
a finite flow problem on a truncated Bruhat-Tits tree for `PGL_2(Q_2)`.  The
local endpoint law from THM-367 is a conservative tree kernel.  A primitive
15-speed open-cover counterexample would have to realize a harmonic current on
this truncated tree while also satisfying the lower-protector truncation and
gcd/primitivity constraints.

The conjectural obstruction is:

```text
Every locally harmonic fan current that closes dyadic endpoints is imprimitive,
and every primitive gcd-breaker introduces positive divergence.
Positive divergence is witnessed as an interval gap or a peelable endpoint leaf.
```

## Evidence

`lrc_bruhat_tits_tree_s399.py` reorganizes THM-367 as a tree transition law.
For owner `u=2^k`, protector `p=2^j q`, and drop `L=k-j`, normalize protected
endpoint count by the horosphere size `2u`.

The active odd-residue kernels are:

```text
L = 1:  q mod 16 in {1,15},                 each capacity 1/2
L = 2:  q mod 16 in {1,3,13,15},            each capacity 1/4
L >= 3: q mod 16 in {1,3,5,7,9,11,13,15},   each capacity 1/8
```

In every case:

```text
sum over active residue classes of normalized capacity = 1.
```

Thus THM-367 is not only a count formula.  It is a conservative boundary
kernel: each dyadic drop transports one full horosphere of normalized mass
through a different residue aperture.

S399 also records the lower-protector truncation.  A maximal owner can only use
`p<u`, so the full symmetric `q mod 16` kernel is cut to one side:

```text
drop 1: lower active {1}
drop 2: lower active {1,3}
drop 3: lower active {1,3,5,7}
drop >=4: lower active all odd residues
```

This is the finite `Q_2` apartment being viewed from one side of the maximal
owner.

The known nine-cover is a stable tree-star current:

```text
u/2  plus  (u/32)*{1,3,5,7,9,11,13,15}.
```

For `u>=32` it normalizes to:

```text
(16,1,3,5,7,9,11,13,15)
```

with gcd `u/32`.  Its raw incidence mass is always `3u`, so average endpoint
indegree is `3/2`.  The fan is not uniform: two boundary-cylinder apexes are
hit by all nine fan speeds, while a scaled population has degree `2`.

Finally, the `u=16` private-endpoint argument is exceptional.  S399 finds:

```text
u=16:  24 private endpoints relative to all lower protectors
u>=32: 0 private endpoints relative to all lower protectors
```

For `u>=32`, the all-lower incidence graph has the full nonempty core
`(2u, u-1)`.  This supports HYP-1859: the proof cannot be purely local at each
horosphere; it needs a global dyadic flow inequality.

## Interpretation

The boundary of the Bruhat-Tits tree is `P^1(Q_2)`.  The dyadic LRC endpoint

```text
(16m +/- 1)/(16u)
```

is a finite boundary cylinder.  The speed `p=2^j q` selects a depth and an odd
boundary direction.  THM-367 says exactly how much of the owner's endpoint
horosphere that direction can cover.

This reframes the n=16 proof target:

```text
interval cover problem
    -> endpoint incidence problem
    -> labelled endpoint-cycle problem
    -> conservative 2-adic boundary-current problem
```

The S392 fan is then not a nuisance.  It is the basic local harmonic object.
Its failure is not local coverage; its failure is global arithmetic: the fan is
imprimitive unless gcd-breakers are added, and those breakers should create
positive divergence elsewhere.

## Predictions

1. A useful proof invariant should be a signed divergence on a truncated
   dyadic tree, not merely a count of uncovered endpoints at one owner.
2. Imprimitive fan systems should be the only zero-divergence local models.
3. Primitive fan repairs should contribute positive divergence proportional to
   the gcd-breaking debt, visible in S392/S393 as pair gaps or endpoint leaves.
4. A future branch-and-bound search should add Bruhat-Tits features:
   drop distribution, active residue mass, fan-apex load, lower-truncation
   deficit, and gcd divergence.
5. A genuine disproof construction would have to build a leafless harmonic
   current on a finite dyadic subtree.  Random or merely dense interval covers
   should keep failing because they do not satisfy this current condition.

## Sources

- `04-computation/lrc_bruhat_tits_tree_s399.py`
- `05-knowledge/results/lrc_bruhat_tits_tree_s399.out`
- `07-reflections/lrc-bruhat-tits-tree-s399.md`

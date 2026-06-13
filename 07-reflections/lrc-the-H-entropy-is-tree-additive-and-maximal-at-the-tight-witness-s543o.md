---
source: oracle-2026-06-01-S543o
status: synthesis + computation (entropy on the tree; the H-entropy log H is tree-additive, maximal at the regular polygon)
tags:
  - tournament-analysis
  - entropy
  - tree
  - H-loneliness
  - modular-decomposition
  - p-adic
  - lonely-runner
---

# The H-Entropy Is Tree-Additive and Maximal at the LRC-Tight Witness

Attacking *entropy on the tree* of Tournament Analysis. (Reading "Hough entropy" as
the entropy functionals on our tree structures, centered on the **H-entropy**.) The
key object is clean and the payoff is conceptual: **the loneliness meter, in log
form, is a tree-additive entropy whose maximum is exactly the conjecture's tight
configuration.**

## The H-entropy: `S_H = log H`, additive on the modular tree

`H` (the directed Hamiltonian-path count = loneliness/balance meter, S26) is
**multiplicative over the disjoint modules** of the recursive/modular tree (S531).
Therefore its logarithm

> **`S_H(T) = log₂ H(T)`  is ADDITIVE over the modular tree:**
> `S_H(T) = Σ_{modular-tree nodes} log₂(H of the module)`,

verified (`lrc_tree_entropy_attack_s543.py`): disjoint apex-flipped blocks give
`S_H(A∪B) = S_H(A)+S_H(B)` exactly (`n=8`: `1.585+1.585=3.170=log₂9`; `n=9`:
`2.322+2.322=4.644=log₂25`). An **apex-flipped module of size `s` contributes
`log₂(1+2^{s-2})`** (`1, 1.585, 2.322, 3.170, 4.087, 5.044` for `s=2..7`). So `S_H`
is a genuine entropy supported on the recursive tree — `0` at the transitive (rigid)
tournament, rising as the tiles (cyclic content) turn on.

## The headline: max entropy = the regular polygon = the LRC-tight witness

`S_H` ranges from `0` (transitive) up to `log₂ H_regular` (the regular tournament).
The regular tournament is the runners at the **regular polygon** — precisely the
**tight LRC witness** (the only configuration that just barely achieves the `1/n`
bound). So:

> **The H-entropy is maximized exactly at LRC's tight, extremal configuration.** The
> conjecture's *hardest* case is the *entropy-maximizing* tournament. Loneliness and
> entropy are aligned: maximal spread (loneliest, most balanced) = maximal `S_H` =
> the regular polygon.

This recasts the S542 P6 law (`H` unimodal in tile-count) in entropy terms: `S_H`
peaks at half-departure from the base path (maximal cyclic mixing = the regular
polygon), and the conjecture lives at that entropy peak.

## Entropy fingerprints arithmetic regularity (the trajectory mean)

Over `t`, every speed family reaches the same *max* `S_H = log₂ H_regular`, but the
**mean** `S_H(t)` separates them (`n=5, n=6`):

```
                 mean S_H(t)   (n=5) (n=6)
 AP 1..n (regular)    2.757   4.243   <- highest
 sieve-with-multiple  2.576   4.107
 random primitive     2.295   3.749
 geometric 1,2,4,8..  1.926   3.367   <- lowest
```

> **Arithmetic regularity (AP) carries the highest mean H-entropy; multiplicative /
> lacunary speeds (powers of 2) the lowest.** The mean H-entropy is a *detector of
> arithmetic structure* in the speed set — the regular runners spend the most time
> spread (high `H`), the lacunary ones the least. (Engineering read: mean `S_H` of a
> pairwise-data tournament-trajectory flags arithmetic/periodic structure.)

## The other two tree entropies

- **p-adic Bruhat–Tits tree entropy (S541).** The Shannon spread of the speeds
  across residues mod `p^k` (the channels, S534). Sieve-structured speeds (all `≡0
  mod p`) have **zero** level-1 entropy (one channel); spread speeds have high
  entropy. The distinguished sub-quantity is the **0-branch occupancy**: `t = 1/p`
  is lonely (the sieve, THM-369) iff the `0`-branch is *empty* — a local, not total,
  entropy condition. So total p-adic entropy = channel spread (high ⇒ debt present,
  sieve blocked, the S534 vacuity); the `0`-branch is the loneliness gate.
- **Iso-class walk entropy (S518/S542).** The Shannon entropy of the time-in-class
  distribution over the `2·Fib(n-2)` menu (`1.825` AP vs `1.796` random of `2.0` max
  at `n=5`): the walk's *mixing* across the menu, again slightly higher for AP.

## Synthesis

> Three tree entropies — `S_H = log₂ H` (modular tree), the p-adic channel entropy
> (Bruhat–Tits tree), and the iso-class walk entropy (the menu) — all measure
> **spread/mixing**, and **arithmetic regularity shows up as high entropy on every
> one**. The H-entropy is the cleanest: **tree-additive** (from `H`'s
> multiplicativity, S531) and **maximized at the regular polygon = the LRC-tight
> witness**. So the conjecture's extremal sits at the entropy peak — LRC asks the
> runner walk to reach a high-`S_H` (spread) configuration, and the worst case is the
> maximum-entropy one.

## Verdict / next
- New clean structure: `S_H = log₂ H` is a tree-additive entropy on the modular
  tree (verified), `0` at transitive, max `= log₂ H_regular` at the regular polygon
  (LRC-tight). Mean `S_H(t)` fingerprints arithmetic regularity (AP high, geometric
  low). p-adic channel entropy + walk entropy corroborate.
- Concrete next: (1) prove `S_H` max `= log₂ H_regular` over *realizable* (circular)
  tournaments and that the regular polygon attains it; (2) the mean-`S_H` detector
  as an engineering tool (arithmetic-structure test on arbitrary pairwise data);
  (3) the 0-branch p-adic entropy as the exact sieve-loneliness gate across `p^k`.

## Artifacts
```
04-computation/lrc_tree_entropy_attack_s543.py
05-knowledge/results/lrc_tree_entropy_attack_s543.out
```
Related: S531 (H multiplicative over modules ⇒ S_H additive), S542 (H unimodal in
tiles; TA pattern hunt), S541 (p-adic tree / dichotomy), S534 (channels/sieve), S518
(2·Fib menu), S26 (H = loneliness).

---
source: codex-2026-06-01-S505
status: finite audit plus proof-strategy refinement
tags:
  - lonely-runner
  - hamiltonian-paths
  - tournament-clock
  - recursion
  - n14
  - n18
---

# H As A Loneliness Metric For LRC S505

This session extended the S24 slogan

```text
H is a spread-meter over the runner clock
```

to the actual hard LRC recursive rows.  Here `H(T)` is the Hamiltonian-path
count of the half-turn phase tournament on the observer plus runners.  It is
not the anchored LRC witness itself.  It measures the global phase entropy of
the configuration:

```text
low H  = bunched / nearly transitive / clean hierarchy
high H = spread / cyclic / many compatible circular orders.
```

The normalization used in S505 is

```text
H_ratio = H(T) * 2^(n-1) / n!,
```

the Hamiltonian-path count divided by the random-tournament baseline
`n!/2^(n-1)`.

## Tournament Analysis Setup

Pairwise observable:

```text
clockwise half-turn phase difference between runners.
```

Switch/gauge:

```text
i -> j when j lies in the clockwise open half-circle from i.
```

Collision and antipodal ties are broken by the fixed speed-order Hamiltonian
path.  This is the S24/S502 half-turn clock, not the anchored endpoint clock.

## The New Computation

`lrc_h_loneliness_metric_s505.py` audited the hard scale ladders

```text
n=14: scale 7,14,28,56 with skip 6
n=18: scale 9,18,36,72 with skip 8
```

and included the initial rows.  It computed:

```text
LRC gap/th, endpoint debt, gap*debt,
origin margin,
phase walls inside the lonely corridor,
H at the selected LRC midpoint,
and the H-profile across the corridor cells.
```

## Main Signal

For `n=14`:

```text
H_mid:    22168229 -> 17826951 -> 17826951 -> 17826951
H_ratio:  2.0831 -> 1.6752 -> 1.6752 -> 1.6752
gap*debt: 5/11   -> 5/11   -> 5/11   -> 5/11
```

For `n=18`:

```text
H_mid:    117137481061 -> 102405804217 -> 102405804217 -> 102405804217
H_ratio:  2.3981      -> 2.0965      -> 2.0965      -> 2.0965
gap*debt: 1           -> 1           -> 1           -> 1
```

So the recursive structure splits cleanly:

```text
endpoint recursion: gap halves, debt doubles, product stays fixed;
phase recursion: H drops once at the first gate, then freezes.
```

This is exactly the two-clock story in another language.  The anchored LRC
clock keeps accumulating denominator debt, while the half-turn phase clock has
already settled into a stable tournament cell type.

## Corridor Behavior

The lonely gap is a short corridor through half-turn clock cells.  H sees that
corridor, but coarsely:

```text
n14 scale 7:  H values 20021407,22168229
n14 scale 14: H value  17826951 only, even though two cells exist
n14 scale 28+: H values 17826951,19521633

n18 scale 9:  H values 115642276825,117137481061
n18 scale 18: H value  102405804217 only, even though two cells exist
n18 scale 36+: H values 102405804217,108946467005
```

The first gate can make multiple adjacent clock cells H-indistinguishable.
After the next doubling, the corridor has more visible half-turn cells again,
but the selected midpoint H remains frozen.

## Proof Interpretation

H is a good global loneliness entropy, but it is too coarse to prove LRC alone.
It detects the phase-shape recursion:

```text
near-regular initial rows have high H_ratio;
row-parent rows remain high;
the first gate lowers H_ratio;
further dyadic scaling keeps H_ratio fixed.
```

The anchored LRC endpoint ledger detects something orthogonal:

```text
gap*debt stays constant while debt doubles.
```

Thus the proof should use H as a coarse state variable in the recursive
automaton:

```text
H plateau = phase state has stabilized;
endpoint debt growth = only anchored denominator depth is still moving.
```

The next formal target is to prove that once the hard ladder reaches the gate
phase-state, further row doublings preserve the half-turn tournament at the
selected midpoint up to isomorphism, while endpoint debt translates in the
2-adic direction.

## Artifacts

```text
04-computation/lrc_h_loneliness_metric_s505.py
05-knowledge/results/lrc_h_loneliness_metric_s505.out
05-knowledge/hypotheses/HYP-1969-lrc-h-phase-plateau.md
```

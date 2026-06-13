---
id: THM-370
name: lrc-two-neighbor-cycle-independence
status: PROVED
date: 2026-05-31
session: codex-2026-05-31-S452
depends_on: []
---

# THM-370: LRC two-neighbor cycle independence

## Statement

Let `p_0,...,p_{n-1}` be `n` points on the circle, listed in cyclic order, and
let

```text
g_i = clockwise gap from p_i to p_{i+1}
```

with indices modulo `n`.  Call a gap safe if `g_i >= 1/n`.

A point `p_i` is lonely, meaning distance at least `1/n` from every other point,
iff the two incident gaps `g_{i-1}` and `g_i` are both safe.

Consequently, a configuration has no lonely point iff its safe gaps form an
independent set in the cycle graph `C_n`.  Moreover, this independent set is
nonempty.

## Proof

The nearest points to `p_i` in clockwise and counterclockwise circular order
are `p_{i+1}` and `p_{i-1}`.  Every other point is reached from `p_i` by adding
one of these incident gaps plus further nonnegative gaps, so its circular
distance is at least the smaller of the two incident gap lengths.

Thus `p_i` has distance at least `1/n` from every other point iff both incident
gaps have length at least `1/n`.

Therefore a lonely point exists exactly when two adjacent gap edges are safe.
No lonely point exists exactly when no two adjacent gap edges are safe, i.e. the
safe gap set is independent in `C_n`.

Finally, the safe set cannot be empty: if every `g_i < 1/n`, then

```text
sum_i g_i < n * (1/n) = 1,
```

contradicting that the circular gaps sum to `1`.

## Significance

This turns the static all-runner LRC obstruction into a graph-theoretic object:
nonempty independent sets of the circular gap graph.  It is the exact
two-nearest-neighbor version of the one-runner endpoint-cover view, and it
gives a concrete bridge to Zeckendorf/path independence and tournament OCF
independence polynomials.

## Related

- HYP-1904
- HYP-1902
- `04-computation/lrc_runner_distance_tournament_s452.py`
- `07-reflections/lrc-runner-distance-tournament-s452.md`

# THM-300: Quadratic Coefficient Matrix Signature

**Status:** REFUTED; negative-count clause first fails at n=9 (FINITE-EXACT)
**Discovered by:** opus-2026-04-04-S3
**Refutation companion:** `04-computation/tournament_quadratic_signature_refutation_thm300.py`
with output `05-knowledge/results/tournament_quadratic_signature_refutation_thm300.out`

## Statement

Let Q be the m×m symmetric matrix of quadratic multilinear coefficients:
Q[i,j] = c_{ij} for i ≠ j, Q[i,i] = 0.

**Refuted conjecture.** For n ≥ 5, Q has exactly n-2 negative
eigenvalues and m-(n-2) positive eigenvalues (full rank).

The conjunction already has a nullity-one exception at `n=5`, which the
original evidence table recorded.  More importantly, its intended
negative-count prediction first fails at `n=9`: the exact inertia is
`(20 positive, 8 negative, 0 zero)`, rather than `(*,7,0)`.

## Evidence

| n | m | Positive | Negative | Zero | n-2 | Match? |
|---|---|----------|----------|------|-----|--------|
| 4 | 3 | 1 | 1 | 1 | 2 | No (rank deficient) |
| 5 | 6 | 2 | 3 | 1 | 3 | Yes (neg count) |
| 6 | 10 | 6 | 4 | 0 | 4 | Yes |
| 7 | 15 | 10 | 5 | 0 | 5 | Yes |
| 8 | 21 | 15 | 6 | 0 | 6 | Yes |
| 9 | 28 | 20 | 8 | 0 | 7 | **No** |
| 10 | 36 | 27 | 9 | 0 | 8 | **No** |
| 11 | 45 | 35 | 10 | 0 | 9 | **No** |
| 12 | 55 | 44 | 11 | 0 | 10 | **No** |

At n=5, Q has rank 5 (one zero eigenvalue) but still n-2=3 negative eigenvalues. For n≥6, Q appears to have full rank.

The exact audit proves full rank only on the finite range `6<=n<=12`.
The replacement pattern `negative=n-1` holds for `9<=n<=12` but is **OPEN**
beyond that range; it is a probe, not a repaired conjecture promoted here.

## First failed mechanism

THM-299 makes `Q_(n-1)` a principal old-tile block of `Q_n`.  When that block
is invertible, exact block congruence splits the inertia into the old inertia
and that of the new hypotenuse-layer Schur complement.  Those layer inertias
are

```text
n=7: (4+,1-)   n=8: (5+,1-)   n=9: (5+,2-)
n=10:(7+,1-)   n=11:(8+,1-)   n=12:(9+,1-).
```

Thus the extrapolation's first failed implication is precise: “each new
hypotenuse layer contributes exactly one negative direction” breaks at the
seven-tile layer entering `n=9`, which contributes two.  The anomaly is not a
floating-point threshold or a hidden zero eigenvalue.

The companion obtains every quadratic coefficient from only the transitive,
one-flip, and two-flip Hamiltonian-path DPs, then computes inertia over `Q` by
explicit symmetric congruence with one- and two-dimensional pivots.  Ordinary
and optimized modes byte-match.  Script/output/semantic SHA-256 are
`a7ec43cae1680a6850c6aab9d58ebb831de06a44ffbc6dc15aff98329f072269`,
`fe05c45a0c64e4cf56570ab1f30ddad70cd80108bef638da02f54d642aa1da45`,
and `60fc4bca9300e5fd220aec78c59e0381f549019ad622b735c0fa235bd15ced7b`.

## Remarks

- The trace of Q is always 0 (multilinear polynomial has no t_i² terms).
- The old observation that the negative eigenspace concentrates near the two
  staircase boundary legs did not control the sign of the layer Schur
  complement; it therefore cannot justify the failed count.
- The number n-2 equals the dimension of the staircase (number of rows or columns).
- The nullspace at n=5 involves tiles at staircase coordinates forming an L-shape with alternating signs.
- A useful next question is whether the `n=9` extra direction is sporadic or
  the first instance of a periodic layer-inertia defect.  The cheapest next
  test is to compute Schur inertia without recomputing old coefficients, not
  to enumerate the full tournament cube.

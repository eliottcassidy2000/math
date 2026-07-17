# HYP-7152 — Local-Density Block Gluing

**Owner:** codex-2026-07-16-S21  
**Target:** THM-933 (renumbered after Opus first pushed the THM-932 claim stub)
**Status:** RESOLVED as THM-933.  The full canon proof and exact referee pass;
the Lean algebra is the remaining formalization rung.

## Claimed interface

For a nonempty speed block `B`, let `S(B)` be its common safe set on
`T = R/Z`, let

```text
delta(B) = mu(S(B)),
q(B)     = osc_u integral_0^u (1_{S(B)} - delta(B)),
M(B)     = sum_{x in B} x.
```

The primitive oscillation is the exact one-interval discrepancy constant:

```text
mu(I intersect S(B)) >= delta(B) mu(I) - q(B)
```

for every circular interval `I`.  Consequently, if blocks `B_1,...,B_m`
partition a speed set and `W_r = intersection_{i<=r} S(B_i)`, then

```text
mu(W_m) >= product_i delta(B_i)
           - sum_{r=2}^m q(B_r)
               (sum_{i<r} M(B_i))
               product_{s>r} delta(B_s).
```

The proof has two elementary inputs.  First, sum the interval discrepancy over
the components of `W_{r-1}`.  Second, the complement of `W_{r-1}` is a union
of at most `sum_{i<r} M(B_i)` danger teeth, so `W_{r-1}` has at most that many
components.  Induction then unrolls the displayed recurrence.

Scale covariance supplies the advertised two-scale composition:

```text
delta(cB) = delta(B),      q(cB) = q(B)/c,      M(cB) = c M(B).
```

Thus exact density certificates live *inside* blocks, while gaps between block
scales divide their boundary debts.  HYP-7104 and THM-930 are candidate density
suppliers; each needs the exact primitive-discrepancy sidecar `q` before it can
be glued without an unjustified independence assumption.

## First exact consequence

For a singleton block `{x}` at loneliness radius `1/14`,

```text
delta({x}) = 6/7,          q({x}) = 6/(49x).
```

For thirteen ordered speeds with `x_{r+1} >= R x_r`, the theorem gives

```text
mu(common safe set)
  > (6/7)^13
    - (6/(49(R-1))) sum_{k=0}^{11} (6/7)^k.
```

This is positive whenever `R (6/7)^12 > 1`; hence `R >= 7` proves the
lacunary LRC(14) family.  This sharpens THM-928's elementary `R >= 15`
threshold by replacing a loose per-tooth boundary loss with the exact centered
primitive discrepancy.

## Assumption challenge / Tournament Analysis

The useful vertices are certified blocks, not runners.  This quotient preserves
the density product and the entire boundary-debt ledger `(delta,q,M)` but
destroys within-block endpoint order; that geometry must remain in the exact
`q` sidecar.  For pairwise ordering analysis, orient `A -> B` when placing `A`
before `B` incurs no more two-block debt than the reverse placement after scale
normalization.  Ties use the physical increasing-scale order as the Hamiltonian
path.  The referee script will report score histograms, cycles, SCCs, edge flips,
and Hamiltonian-path counts.

## Resolution checklist

1. Circle-component and primitive-discrepancy proof: complete in THM-933.
2. Exact rational endpoint sweep and gluing ledger: complete and passing.
3. Singleton, multi-speed, and THM-928(C) sidecar checks: complete.
4. Lean local-to-component sum and recurrence algebra: next rung.
5. Canon promotion and dependency audit: complete as THM-933; the concurrent
   THM-932 target is the coarser fixed-local-scale interface.

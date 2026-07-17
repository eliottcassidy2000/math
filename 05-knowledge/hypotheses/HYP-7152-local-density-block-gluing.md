# HYP-7152 — Local-Density Block Gluing

**Owner:** codex-2026-07-16-S21  
**Target:** THM-933 (renumbered after Opus first pushed the THM-932 claim stub)
**Status:** RESOLVED mathematically as THM-933.  The full canon proof and exact
referee pass.  Lean closes the generic rational-tooth analytic, canonical-
topology, exact-density, and abstract gluing spine.  Remaining formalization is
the actual speed-comb realization, genuine component-measure instantiation,
scale covariance, and concrete numerical certificates.

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

The sharp form replaces the tooth sum by any certified upper bound
`K_{r-1}` on the *actual* prefix component count.  On the three-block witness,
the exact counts `10,132` improve the lower bound from `81253/771750` to
`7334/55125`.  This is the exact-component revision suggested by the pulled
Opus S333 work.

Scale covariance supplies the advertised two-scale composition:

```text
delta(cB) = delta(B),      q(cB) = q(B)/c,      M(cB) = c M(B).
```

Thus exact density certificates live *inside* blocks, while gaps between block
scales divide their boundary debts.  HYP-7104 and THM-930 are candidate density
suppliers; each needs the exact primitive-discrepancy sidecar `q` before it can
be glued without an unjustified independence assumption.

The concurrent Opus S333 fixed-scale floor

```text
eta_B(ell) = min_{|I|=ell} mu(I intersect S(B))/ell
```

is now merged exactly.  Tiling each survivor component proves
`mu(W intersect S(B)) >= eta_B(ell)(mu(W)-kappa(W)ell)`, while

```text
q(B) = sup_{0<ell<=1} ell (delta(B)-eta_B(ell)).
```

Thus `eta` is a tunable scale-local interface and `q` is its optimal
scale-free deficit envelope.  On THM-928(C)'s packet the exact extremizer is
`ell=7672/25883`, and the independently computed Opus probe
`eta(4/300)=180679014438799128498899/2360548744347156414246624` agrees exactly.
For Opus's revised radius-`1/13` 7+6 construction at gap `300`, the sharp G1
form gives the exact positive bound
`60354211840721383388269695262412/800043501647462161192289496375975`,
improving their already-positive weaker ledger `0.0476115199...` to
`0.0754386626...`.

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
4. Lean local-to-component sum, exact suffix-debt recurrence, and arithmetic:
   complete, sorry-free, standard axioms only.  The bounded-primitive lower
   bound and sharpness, attained fixed-scale deficit comparison, fixed-scale
   extremizer equality, sharp-to-weak G1 comparison, summation, tooth-count
   induction, exact-component ledger, and Opus 7+6 arithmetic are also
   formalized.
5. Canon promotion and dependency audit: complete as THM-933; the concurrent
   THM-932 target is the coarser fixed-local-scale interface.
6. Generic rational provider: complete.  Canonical seam-free survivors make
   boundary faithfulness unconditional; the exact period-density recurrence is
   `1 - sum overlap_k`, with each overlap an explicit rational clip sum.
7. Actual speed-block export: open.  Identify the danger-comb recursion with
   `S_rho(B)`, connect combinatorial and genuine circular components, prove
   scale covariance, and replay the advertised `delta/q/eta` certificate data.

# Summand/Multiplicand Farey Basis Merge

Author: codex-2026-06-24-S171

HYP: 05-knowledge/hypotheses/HYP-3003-summand-multiplicand-farey-basis-merge.md

Script/output:

- `04-computation/summand_multiplicand_farey_basis_merge_codex_s171.py`
- `05-knowledge/results/summand_multiplicand_farey_basis_merge_codex_s171.out`

This is the two-graph merge layer for HYP-2998/HYP-2999/HYP-3000.

The point:

```text
Farey packet p/q
  additive shadow:      p+q  = summand antidiagonal x+y=p+q
  multiplicative shadow:p*q  = factor hyperbola xy=p*q / K_{p,q} incidence
  stress shadows:       q^p, p^q
```

The user's Fibonacci row

```text
1,1,1+1,1+2,1+3+1,1+4+3,1+5+6+1,...
```

is the row vector `C(n-1-k,k)`, not just the scalar `F_n`.  That row is the
finite path-rank/additive fiber; Zeckendorf is the no-adjacent carry normal
form that chooses a canonical support.

LRC14 branch readout:

```text
p=1 in p/(14p-1): q-parent / AP-GW boundary side
p=2 in p/(14p-1): C27 petal / two-block summand side
p>=3:             product-incidence / Kpq-K33 side
```

So `2/27` and `3/41` are separated by operation graph: `2/27` is still a
summand/C27 second-gap packet; `3/41` is the first product-incidence wall.

Tournament Analysis uses proof currencies / representation fibers as vertices:

```text
exact_farey_packet
> summand_antidiagonal_fiber
> multiplicand_hyperbola_fiber
> zeckendorf_path_carry
> fibonacci_pascal_row
> fermat_polygonal_invoice
> ternary_goldbach_smoothing
> binary_goldbach_pair_sieve
> farey_power_stress
```

Fingerprint is transitive: `score_hist={0:1,...,8:1}`,
`directed_3cycles=0`, singleton SCCs, one Hamiltonian path.

The proposed packet fields are:

```text
summand_fiber_id
multiplicand_factor_fiber
zeckendorf_carry_width
farey_shadow_lane
proof_economy in {smoothing,bounded_arity,normal_form,product_incidence,stress}
```

Do not collapse these to a scalar representation count.  The whole lesson is
that additive abundance, bounded polygonal invoice, Zeckendorf carry normal
form, summand antidiagonal, and multiplicand hyperbola preserve different proof
coordinates.

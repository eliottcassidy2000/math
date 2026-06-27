# Fibonacci Additive-Basis Farey Bridge

HYP-3000 is a proof-interface synthesis for the user's Fibonacci / Fermat
polygonal / Goldbach / Zeckendorf / Farey prompt.  It complements incoming
HYP-2998 and HYP-2999: HYP-2998 handles the golden Stern-Brocot/Fibonacci
carrier, HYP-2999 keeps Pascal-slope row fields attached as packet labels, and
HYP-3000 emphasizes the path-rank identity and the LRC14 unit-excess
`p/(14p-1)` chain as a proof-currency classifier.

The Fibonacci row pattern

```text
1, 1, 1+1, 1+2, 1+3+1, 1+4+3, 1+5+6+1, ...
```

is the exact identity

```text
F_n = sum_k binom(n-k-1,k)
```

so for `n>=2` it is the independent-set rank vector of the path graph
`P_{n-2}`.  This gives a clean bridge between representation abundance and
normal-form uniqueness:

- Goldbach: high-entropy sieve representation fibers.
- Ternary Goldbach: one added smoothing / hypergraph dimension.
- Fermat polygonal: bounded arity and residue absorption.
- Fibonacci / Zeckendorf: path-normal-form carries with no adjacent chosen
  ranks.

The Farey translation keeps the HYP-2984 rule: retain exact `M=p/q` and
excess `e=14p-q` first.  Then `p+q` is the affine-safe additive payload,
`p*q` is an incidence/product side channel, and `q^p,p^q` are magnitude stress
tests.  On the unit-excess chain `p/(14p-1)`, we get `q=14p-1`,
`p+q=15p-1`, and `p*q=14p^2-p`.

Tournament Analysis used proof currencies rather than runners.  The transitive
path was

```text
zeckendorf_normal_form
> fibonacci_path_rank
> farey_sum_affine_check
> fermat_polygonal_bounded_arity
> ternary_goldbach_smoothing
> farey_product_incidence
> binary_goldbach_sieve
> farey_power_stress_test
```

Next useful pull: classify HYP-2963 residual packets by proof economy before
choosing tools.  Is a packet asking for smoothing, bounded arity, a
Zeckendorf/Ostrowski normal form, additive Farey scale, product incidence, or
only a magnitude stress test?

Artifacts:

- `04-computation/fibonacci_additive_basis_farey_bridge_codex_s169.py`
- `05-knowledge/results/fibonacci_additive_basis_farey_bridge_codex_s169.out`
- `05-knowledge/hypotheses/HYP-3000-fibonacci-additive-basis-farey-bridge.md`
- `07-reflections/fibonacci-additive-basis-farey-bridge-codex-s169.md`
- `LTI-150`

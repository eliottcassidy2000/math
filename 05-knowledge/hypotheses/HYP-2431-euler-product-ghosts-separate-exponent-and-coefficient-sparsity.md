# HYP-2431 - Euler-product ghosts separate exponent sparsity from coefficient sparsity

**Status:** CONFIRMED in pilot atlas; OPEN as a general carrier principle.
**Source:** codex-2026-06-11-P3.
**Artifacts:** `04-computation/euler_product_ghost_atlas_codex.py`,
`05-knowledge/results/euler_product_ghost_atlas_codex.out`.

## Statement

For

```text
P_b(q) = prod_{m>=1} (1-q^m)^{b_m},
```

the exponent schedule `b_m`, the coefficient sequence of `P_b`, and the ghost
sequence

```text
G_b(n) = sum_{d|n} d b_d
```

are different carriers. Sparsity in one carrier does not imply sparsity in the
others.

## Evidence

Through `q^180`:

```text
eta_all:    b_m dense, ghost dense, coefficients sparse (22/181)
prime_only: b_m sparse (41 primes), ghost almost dense (179/180), coefficients dense (158/181)
mobius:     b_m mixed/sparse-zero, ghost dense, coefficients dense-ish (163/181)
liouville:  b_m dense signed, ghost dense, coefficients dense-ish (168/181)
```

So Euler's pentagonal sparsity is not explained by exponent sparsity. It is a
special product/coefficient cancellation phenomenon.

## Proof Route

Formalize the three-carrier map:

```text
exponents b_m  ->  ghosts G_b(n)=sum_{d|n} d b_d  ->  coefficients [q^n]P_b.
```

The first arrow is divisor/Witt-linear; the second is nonlinear coefficient
extraction. The proof-facing invariant may live at any one of the three levels.

For eta, the coefficient carrier is sparse while the ghost carrier is dense
(`sigma_1`). For prime products, the exponent carrier is sparse while the ghost
carrier is almost dense. This is the multiplicative analogue of HYP-2426's
"scalar shadow vs support object" warning.

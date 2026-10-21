# THM-260: Layered Decomposition of the Rooted Tournament Sequence

**Status:** PROVED (exact Burnside computation, verified n=3..8)
**Filed by:** kind-pasteur-2026-03-20-S6

**Correction (2026-10-21):** The general layer denominator below now uses
factorials of the multiplicities of each nontrivial cycle length, rather than
`r!`. Unequal cycle lengths cannot be interchanged in the centralizer. The
mixed `(5,3)` layer is `128/15`, as already recorded correctly in the table;
the rooted count prefix and all displayed table entries are unchanged. This
is a scoped denominator repair, not a re-audit of unrelated interpretations.

## Statement

The rooted tournament sequence P(n) = 1, 2, 4, 12, 48, 296, 3040, 54256 satisfies:

**P(n) = n * T(n) - D(n)**

where T(n) = number of tournament isomorphism classes, and the deficit D(n) decomposes
as a sum over all non-identity partitions with ALL ODD parts:

**D(n) = sum over odd partitions lambda (not [1^n]) of D_lambda(n)**

## Layer Formulas

Each layer has an exact closed form:

### Single k-cycle layer (k odd):
D_k(n) = 2^{(k-1)/2 + C(n-k+1, 2)} / (n-k)!

- D_3(n) = 2^{1 + C(n-2, 2)} / (n-3)!
- D_5(n) = 2^{2 + C(n-4, 2)} / (n-5)!
- D_7(n) = 2^{3 + C(n-6, 2)} / (n-7)!

### Two 3-cycle layer:
D_{3,3}(n) = 2^{5 + 2(n-6) + C(n-6, 2)} / (3 * (n-6)!)

### General formula:
For partition lambda = (k_1, k_2, ..., k_r, 1^m) with all k_i >= 3 odd,
let a_l be the multiplicity of the nontrivial cycle length l. Then:
D_lambda(n) = 2^{f(lambda,n)} * (n - m) / (m! * prod_l l^{a_l} a_l!)

where f(lambda, n) = number of free edge-orbit pairs. The denominator is the
permutation centralizer size z_lambda. Its rooted Burnside term is
m * 2^{f(lambda,n)} / z_lambda, so subtracting it from n times the ordinary
term proves the displayed formula. For lambda = (5,3), f = 4 and z_lambda =
15, giving 2^4 * 8 / 15 = 128/15; the former r! denominator gave 64/15.

## Verified Values

| n | D_3 | D_5 | D_{3,3} | D_7 | D_{5,3} | D(n) |
|---|-----|-----|---------|-----|---------|------|
| 3 | 2 | 0 | 0 | 0 | 0 | 2 |
| 4 | 4 | 0 | 0 | 0 | 0 | 4 |
| 5 | 8 | 4 | 0 | 0 | 0 | 12 |
| 6 | 64/3 | 8 | 32/3 | 0 | 0 | 40 |
| 7 | 256/3 | 16 | 128/3 | 8 | 0 | 152 |
| 8 | 8192/15 | 128/3 | 512/3 | 16 | 128/15 | 784 |

All sums are exact integers.

## The Central Binomial Miracle

D(n)/2 = C(2(n-3), n-3) for n = 3, 4, 5, 6.

This is because for n <= 6, only layers (3), (5), and (3,3) are active,
and their combined contribution happens to equal the central binomial.

At n=7, the (7) layer activates and the (3,3) layer grows faster than
the central binomial absorbs, causing an excess of 6.

Excess sequence: E(n) = D(n)/2 - C(2(n-3), n-3) = 0, 0, 0, 0, 6, 140, ...

## Layer Activation Table

| Layer | First n | Growth rate | Role |
|-------|---------|-------------|------|
| (3) | n=3 | ~2^{n^2/2} / n! | Dominant layer |
| (5) | n=5 | ~2^{n^2/2} / n! | First correction |
| (3,3) | n=6 | ~2^{n^2/2} / n! | Causes central binomial deviation |
| (7) | n=7 | ~2^{n^2/2} / n! | Second correction |
| (5,3) | n=8 | ~2^{n^2/2} / n! | Third correction |
| (3,3,3) | n=9 | ~2^{n^2/2} / n! | Triple correction |
| (9) | n=9 | small | Fourth single-cycle |

Each layer "activates" at n = sum(parts), then grows for all larger n.

## Connection to the P(n) = T(n+1) Coincidence

The coincidence P(n) = T(n+1) for n <= 4 follows from:
- P(n) = n*T(n) - D(n)
- For n <= 4: D(n) = n*T(n) - T(n+1) (which happens to hold exactly)
- This requires n*T(n) - T(n+1) to equal the sum of active layers

For n=1: 1*1 - 1 = 0 = D(1) ✓
For n=2: 2*1 - 2 = 0 = D(2) ✓
For n=3: 3*2 - 4 = 2 = D(3) ✓
For n=4: 4*4 - 12 = 4 = D(4) ✓
For n=5: 5*12 - 56 = 4 ≠ 12 = D(5) ✗

So the coincidence is equivalent to T(n+1) = n*T(n) - D(n) for n <= 4.

## Interpretation

The layered structure of P(n) reveals that the AUTOMORPHISM CORRECTION
in tournament counting has a HIERARCHICAL structure indexed by odd partitions.
Each layer corresponds to a symmetry type (cycle structure) that reduces
the number of distinguishable vertex roles.

The layers form a GRADED RING structure: the product of two layers
(e.g., (3) * (3) = (3,3)) gives the next correction. This multiplicative
structure underlies the departure from the central binomial formula.

## Scripts

- `04-computation/P7_burnside.py` — Burnside computation of P(n)
- `04-computation/deficit_layers.py` — Layer decomposition
- `04-computation/layer_formulas.py` — Exact closed-form formulas

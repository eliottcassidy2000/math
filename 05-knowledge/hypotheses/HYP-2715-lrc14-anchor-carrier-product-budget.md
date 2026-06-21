---
id: HYP-2715
title: LRC14 anchor-separated carrier product budget in miss-zeta layers
status: OPEN; exact scout, proof obligation corrected
source: codex-2026-06-21
depends_on:
  - HYP-2714
  - HYP-2713
  - HYP-2698
  - HYP-2694
  - THM-557
related:
  - HYP-2717
  - HYP-2716
  - HYP-2704
  - HYP-2708
  - HYP-2675
  - HYP-2676
  - THM-546
  - OPEN-Q-108
---

# HYP-2715 - Anchor-Separated Carrier Product Budget

## Claim

The HYP-2714 moderate-span multi-block residual should be written in
miss-zeta coordinates.  For a split row

```text
E = {0} union (M_1+B_1) union ... union (M_g+B_g),
```

with `0` kept as the observer anchor, define

```text
z_actual(R)  = Pr_x[all sectors in R are missed by E at x],
z_prod(R)    = integral_x product_i Pr_theta[block B_i misses R at (x,theta)] dx.
```

Then

```text
p0(E)        = sum_R (-1)^|R| z_actual(R),
Product(E)  = sum_R (-1)^|R| z_prod(R),
Product(E)-p0(E) = sum_R (-1)^|R| (z_prod(R)-z_actual(R)).
```

The proof target is not a pointwise residual dominance theorem and not
`p0(E) <= Product(E)`.  In the anchor-separated gauge the product can
underestimate actual cover.  The correct cap-level target is the signed budget

```text
|p0(E) - Product(E)| <= cap_k - Product(E).
```

The product gap `cap_k-Product(E)` is large on split rows, so this is much more
forgiving than a raw coordinatewise or absolute-variation bound.

## Exact Scout

Script:

```text
04-computation/lrc14_multiblock_miss_zeta_layers_codex_20260621.py
```

Stored output:

```text
05-knowledge/results/lrc14_multiblock_miss_zeta_layers_codex_20260621.out
```

The script compares exact actual `p0` to the shared-slow-`x` independent-carrier
product for representative split rows, then decomposes `Product-p0` by residual
size `|R|`.

Key rows:

```text
two 4-blocks, moderate gap:
  actual p0       = 11331437/68172720
  product cover   = 109/686
  product-actual  = -3495299/477209040
  |error|/(cap-product) = 45438887/2080501140 ~= 0.02184

5+3 split:
  product-actual  = -79935/26380816
  |error|/(cap-product) ~= 0.00891

3+3+2 split:
  product-actual  = -51629953/57697542630
  |error|/(cap-product) ~= 0.00228

five 2-blocks:
  product-actual  = -2447628624709/93106921650624
  |error|/(cap-product) ~= 0.04873

seven singleton carriers:
  product-actual  = -7586871109925/691974784572029
  |error|/(cap-product) ~= 0.03071
```

Every tested anchor-separated row has `Product < actual`, so a product upper
bound is false in this gauge.  But the error consumes at most about `4.9%` of
the cap-product slack in the tested bank.

## Layer Signal

The signed layer decomposition is strongly alternating.  Typical split rows
have large unsigned mass in `|R|=2,3,4`, but the signed sum is much smaller.
For example, the moderate two-4-block row has:

```text
|R|=2 signed = -0.06332, unsigned = 0.14104
|R|=3 signed = +0.12423, unsigned = 0.16455
|R|=4 signed = -0.08010, unsigned = 0.08767
total product-actual = -0.00732
```

So the missing analytic lemma should preserve the alternating miss-zeta
structure.  Bounding all residual coordinates absolutely is the wrong currency:
it spends an order of magnitude too much budget and repeats the false
small-`R` cone route from HYP-2697/HYP-2698.

## Boolean-Cube Shadow Update

HYP-2716 sharpens the layer signal.  If

```text
d_R = z_prod(R)-z_actual(R),     W_j = sum_{|R|=j} d_R,
M_h = sum_j W_j K_h(j)
```

where `K_h` is the binary Krawtchouk polynomial on the six residual sectors,
then

```text
Product(E)-p0(E) = M_6
```

because `K_6(j)=(-1)^j`.  Equivalently, the scalar cover error is the
all-six-sector Walsh character of the miss-zeta discrepancy.  The updated scout
shows this top character is much smaller than the lower shadows in the tested
split bank; aggregate character risk orders as

```text
h=2 > h=0 > h=4 > h=1 > h=3 > h=5 > h=6.
```

Thus the next proof should bound the top character after quotienting, not the
coordinate `L1` discrepancy before quotienting.

HYP-2717 then expands this top character in carrier Fourier modes.  In that
view the multi-block product is a main term, not an envelope: exact carrier
relations `n.M=0` survive whenever there is more than one integer carrier, and
the proof must show that separated carriers force those exact relations to have
large enough height for Fourier decay to fit inside `cap-product`.

## Proposed Proof Obligation

For any HYP-2714 split branch with carrier gaps above the finite window:

```text
sum_R (-1)^|R| (z_actual(R)-z_prod(R)) <= cap_k - Product(E).
```

Possible proof split:

1. Use the exact miss-zeta factorization to isolate each residual-size layer.
2. Apply a signed Erdos-Turan/Koksma estimate to the large unsigned layers
   `|R|=2,3,4`, preserving their alternating signs.
3. Use direct cap-product slack for the small tail layers `|R|=5,6`.
4. Route low-gap/low-height resonances to the HYP-2714 finite moderate-span
   ledger rather than trying to prove a uniform absolute discrepancy theorem.

## Tournament Analysis

Vertices are proof obligations / tested split rows, not runners.  Pairwise
observable: larger `cap_k-Product(E)` slack, then smaller `|Product-p0|`.
Switch/gauge: pass to miss-zeta residual masks before scalar cover.  Tie
Hamiltonian path in the scout:

```text
five 2-blocks
> 3+3+2 split
> seven singleton carriers
> 5+3 split
> two 4-blocks, moderate gap
> two 4-blocks, wider gap
```

Fingerprint: transitive tournament, score histogram `{0:1,1:1,2:1,3:1,4:1,5:1}`,
zero directed 3-cycles.

## Status

This does not prove LRC(14).  It removes a tempting but false formulation of
the final Route E lemma and replaces it with a signed budget statement that
matches the observed exact arithmetic and the HYP-2698 generated-context
lesson: preserve miss-zeta product structure before scalarizing.

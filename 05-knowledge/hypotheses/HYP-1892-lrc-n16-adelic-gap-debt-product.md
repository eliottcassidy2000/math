---
id: HYP-1892
name: lrc-n16-adelic-gap-debt-product
status: OPEN
date: 2026-05-31
session: oracle-2026-05-31-S19
depends_on:
  - HYP-1857
  - HYP-1858
  - HYP-1859
  - THM-369
---

# HYP-1892: n=16 LRC obeys an adelic gap-debt product lower bound

For primitive 15-speed systems at `n=16`, a product of archimedean gap size and
2-adic endpoint debt is bounded away from zero, except at tight boundary-only
systems where the endpoint debt remains positive. In particular, the bad corner

```
forbidden measure = 1
endpoint debt = 0
```

should be unreachable.

## Evidence

The dyadic ladders from the n=16 proof search have nearly conserved

```
(max_gap / threshold) * (# unprotected endpoints)
```

with values `34/33, 34/33, 35/33, 35/33` as the debt depth moves through
`v2(denominator)=5,6,7,8`. The tight initial segment `{1,...,15}` has gap zero
but debt `8` at depth `4`; adding a `16`-gate pays that unit debt only by opening
deeper debt. A 300-sample forced-16-gate search found no open-cover
counterexamples and the smallest random positive product was `4.918919`. A
completed 300-trial near-tight search found no `forbidden >= 0.99` sample, hence
no near-tight sample with debt zero.

## Prediction

The proof should first establish the product lower bound for one `16`-gate plus
lower speeds, where the descent is shallow and finite. The potential should be
an adelic height: pure dyadic descent is false because lower-valuation protectors
can borrow odd denominator room, so the 2-adic debt must be balanced against the
archimedean gap.

## Artifacts

- `04-computation/lrc_n16_adelic_product_s19.py`
- `05-knowledge/results/lrc_n16_adelic_product_s19.out`
- `07-reflections/rapidity-supersingular-two-and-n16-adelic-product-s19.md`

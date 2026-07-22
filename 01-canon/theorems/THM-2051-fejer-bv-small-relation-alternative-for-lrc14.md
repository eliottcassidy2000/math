---
id: THM-2051
title: Fejer--BV whole-product approximation gives a small-relation-or-positive-BONF5 alternative for LRC(14)
status: RESERVED / UNDER HOSTILE AUDIT. The finite Fourier annihilation, BV rate, and numerical budget have been derived; the exact BONF5 centered-expansion normalization and all relation quantifiers are being independently checked before this is marked proved. This does not classify the small-relation branch or prove LRC(14).
source: codex-2026-07-21-LRC-unrelated-transfer
depends_on:
  - THM-604
  - THM-699
  - THM-935
  - THM-1092
related:
  - THM-940
  - THM-946
---

# THM-2051 -- Fejer--BV small-relation alternative

Reserved for the following candidate theorem. Let

```text
u(t)=1_{||t||<1/14}-1/7
```

and let `S` contain thirteen distinct nonzero integer speeds. If there is no
integer relation of exact support between two and five,

```text
sum_{i in A} k_i v_i=0,
2<=|A|<=5,             0<|k_i|<=2^20,
```

then the continuous quintic Bonferroni functional is positive and hence `S`
has a positive-measure strict lonely set.

The proposed proof retains the entire signed centered product. With
`p_H=Fejer_H*u`, positivity of the Fejer kernel and the BV translation estimate
give

```text
||u-p_H||_1 <= (1+2 log(H+1))/(2(H+1)).
```

The stated relation exclusion makes the constant term of every product of two
through five composed polynomials `p_H(v_i t)` vanish. Product telescoping then
bounds each centered joint moment without taking absolute values term by term
in its relation lattice. The remaining audit must verify the exact THM-935
coefficients in this normalization and the final rational budget before the
status is promoted.

Even if proved, this is only the alternative

```text
small support/height relation  OR  positive BONF5.
```

It bypasses the open THM-946 strip/slab estimates on the dissociated branch but
does not solve the structured small-relation branch.

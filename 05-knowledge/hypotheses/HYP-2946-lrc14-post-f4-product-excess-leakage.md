---
id: HYP-2946
title: LRC14 post-F4 Farey product excess is divisor-lattice leakage
status: COMPUTATIONAL SCOUT / exact excess observable; not a proof
source: codex-2026-06-24
tags: [lrc14, farey, divisor-lattice, perfect-numbers, product-excess]
related:
  - HYP-2944
  - HYP-2945
  - HYP-2934
  - HYP-2932
  - HYP-2221
results:
  - 04-computation/lrc14_farey_product_perfect_kuratowski_codex.py
  - 05-knowledge/results/lrc14_farey_product_perfect_kuratowski_codex.out
---

# HYP-2946: post-F4 product excess is divisor-lattice leakage

Let:

```text
P_n = {a*b : 0 < a <= b <= n, gcd(a,b)=1}
M_n = n(n-1), for n > 1
E_n = sum(P_n) - sigma(M_n).
```

HYP-2944 proves:

```text
D(M_n) subset P_n
```

for every `n`, and equality holds only through `n=4`.  Therefore:

```text
E_n = 0 for n <= 4,
E_n > 0 for n >= 5.
```

The first values are:

```text
n=3: E_n=0
n=4: E_n=0
n=5: E_n=36
n=6: E_n=36
n=7: E_n=159
n=8: E_n=263
n=9: E_n=431
n=10: E_n=552
```

This excess is the amount by which Farey products have stopped being merely
the divisor lattice of their own maximum product.

## Interpretation

At `F3` and `F4`, the product set is divisor-complete.  Perfect-number
structure is therefore legitimate as a divisor-lattice side channel:

```text
F3: sigma(6)=12=2*6
F4: sigma(12)=28
```

After `F4`, new Farey products are no longer just divisors of `n(n-1)`.
The product layer becomes a leakage layer.

The first universal extra witness is:

```text
(n-2)/(n-1) -> (n-2)(n-1).
```

For `n >= 5`, this product is in `P_n` but not in `D(n(n-1))`.

## LRC14 Use

`E_n` is a cheap diagnostic for whether the product analogy is still
divisor-controlled.  The answer after `F4` is no.

For LRC14 this reinforces the current proof order:

```text
exact denominator and C27 shell data first,
then product/Kpq owner packets,
then perfect/divisor analogies as side-channel observers.
```

Candidate test:

```text
Attach a product-excess coordinate to bounded LRC row banks and ask whether
near-AP/Goddyn-Wong rows have small excess after quotienting by the exact C27
shell labels, or whether excess is only a false-positive product shadow.
```

The likely null result is still useful: it would certify that perfect/divisor
language is explanatory around `F3/F4` but not theorem-driving for `n=14`
without the C27 branch labels.

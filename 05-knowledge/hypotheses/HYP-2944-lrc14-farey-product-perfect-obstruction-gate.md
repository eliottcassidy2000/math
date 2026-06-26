---
id: HYP-2944
title: LRC14 Farey-product perfect obstruction gate at F3/F4
status: COMPUTATIONAL SCOUT / exact product-set lemma plus proof-interface signal
source: codex-2026-06-24
tags: [lrc14, farey, perfect-numbers, kuratowski, bipartite, divisor-lattice]
related:
  - HYP-2946
  - HYP-2945
  - HYP-2943
  - HYP-2934
  - HYP-2932
  - HYP-2221
  - HYP-2220
results:
  - 04-computation/lrc14_farey_product_perfect_kuratowski_codex.py
  - 05-knowledge/results/lrc14_farey_product_perfect_kuratowski_codex.out
---

# HYP-2944: F3/F4 are the Farey-product perfect obstruction gate

Define the Farey product set:

```text
P_n = {a*b : 0 < a <= b <= n, gcd(a,b)=1}.
```

This is the set of complete-bipartite edge counts attached to reduced Farey
terms:

```text
a/b -> K_{a,b},  a*b = |E(K_{a,b})|.
```

The computation is in
[lrc14_farey_product_perfect_kuratowski_codex.py](/home/bigo/math/04-computation/lrc14_farey_product_perfect_kuratowski_codex.py:1),
with output at
[lrc14_farey_product_perfect_kuratowski_codex.out](/home/bigo/math/05-knowledge/results/lrc14_farey_product_perfect_kuratowski_codex.out:1).

## Exact Product-Set Lemma

Let:

```text
M_n = max(P_n) = n(n-1), for n > 1.
```

Then:

```text
D(M_n) subset P_n
```

for every `n`, where `D(M_n)` is the divisor set of `M_n`.

Reason: every divisor `d` of `n(n-1)` splits uniquely as:

```text
d = d1*d2,  d1|n, d2|(n-1), gcd(d1,d2)=1.
```

Then `d = a*b` with:

```text
a = min(d1,d2), b = max(d1,d2), b <= n, gcd(a,b)=1.
```

So `d` is a Farey product.

Moreover:

```text
P_n = D(M_n) iff n <= 4.
```

For `n >= 5`, the product `(n-2)(n-1)` appears from the reduced Farey term
`(n-2)/(n-1)`, but it cannot divide `n(n-1)` because `n-2` does not divide
`n`.

## The Gate

The two nontrivial divisor-closed stages are:

```text
F3: P3 = {1,2,3,6}       = D(6),  sum(P3)=12=sigma(6)=2*6
F4: P4 = {1,2,3,4,6,12}  = D(12), sum(P4)=28
```

So `F3` sees the first perfect number as the maximum product:

```text
M_3 = 6, sigma(6)=12.
```

`F4` sees the next perfect number as the product-set sum:

```text
sum(P4)=28.
```

At the same step, `F4` introduces:

```text
3/4 -> K_{3,4},
```

which contains `K33`.  Thus `F4` is simultaneously:

```text
last Farey-product divisor closure,
first complete-bipartite K33-minor term,
product-sum 28 perfect-number hit.
```

## Finite Scout

The script checks product sums through `n <= 1000` against the known perfect
numbers in range.  The only hits are:

```text
F3: sum(P3)=12 = 2*6
F4: sum(P4)=28
```

No other perfect or half-perfect product sums occur through `F1000`.

## LRC14 Use

This does not prove LRC14.  It gives a guardrail for the product analogy.

Perfect-number structure is visible only while the Farey product set is still
divisor-closed.  That closes at `F4`.  After `F4`, product extras are real
leakage, so product/perfect analogies must sit below:

```text
exact M/Farey branch,
AP/Goddyn-Wong labels,
C27 shell transfer.
```

The proof-facing reading is:

```text
F3: planar divisor-perfect maximum
F4: perfect-sum plus first K33 wall
F5+: product-excess leakage
```

For LRC14, this supports the existing branch order:

```text
2/27: C27 summand-unit/petal branch
3/41: first K33 multiplicand-incidence branch
```

The perfect-number packet is a side-channel witness for the product/divisor
transition, not a replacement for the LRC denominator proof.

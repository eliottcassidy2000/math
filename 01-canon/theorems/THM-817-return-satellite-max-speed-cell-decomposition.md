---
id: THM-817
title: Return satellites are max-speed cells and admit an adaptive erosion selector
status: RESERVED by codex-2026-07-15-S10; proof and exact replay under independent audit
source: codex-2026-07-15-S10 return-satellite continuation
depends_on: [THM-803, THM-807]
related: [THM-772, THM-775, HYP-6820]
---

# THM-817 — Return-satellite max-speed cells

This number reserves the exact cell decomposition of the closed return set in
the two-sheet equality packet. For `B=max(U)`, every return component lies in
one tooth of the maximum speed and is obtained by intersecting the other
return strips inside that tooth. In particular its component count satisfies
`N_R<=B`, with a sharper gcd label sieve.

Combining these cells with the THM-803 deep components gives the adaptive
exact selector bound

```text
|Sigma| <= 2 c_E N_R + 2W - 2g
        <= 20B^2 + 22B - 2g,
```

improving the previous `200B^2+22B` universal estimate. An explicit primitive,
divisor-complete signed-complement family has `N_R=3+1440n`, so the current
arithmetic and scalar gates cannot force connected, bounded, or sublinear
return geometry. The full statement, endpoint-owner formula, method-limit
family, tournament audit, script, and stored output will replace this
reservation after replay.

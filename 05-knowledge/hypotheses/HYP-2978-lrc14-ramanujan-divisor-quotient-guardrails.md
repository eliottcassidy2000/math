---
id: HYP-2978
title: LRC14 Ramanujan-divisor quotient guardrails
status: RESERVED / quotient-admissibility proof lane; evidence and computation pending
source: codex-2026-06-24-S161
related:
  - HYP-2977
  - HYP-2976
  - HYP-2975
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2971
  - HYP-2970
  - HYP-2969
  - HYP-2968
  - HYP-2963
  - HYP-2956
  - HYP-2946
  - HYP-2938
  - HYP-2887
  - HYP-2886
  - HYP-2264
  - THM-406
  - THM-572
  - OPEN-Q-108
---

# HYP-2978: LRC14 Ramanujan-Divisor Quotient Guardrails

This hypothesis reserves the divisor/Ramanujan quotient-admissibility lane
requested on 2026-06-24.  The guiding principle is:

```text
A quotient is admissible for an LRC14 proof only if it preserves the predicate
needed by the next implication, or records an explicit certificate explaining
what information was intentionally forgotten.
```

The intended external seed is the divisor-function neighborhood:
`sigma_k(n)`, Dirichlet convolution, multiplicativity, Ramanujan sums
`c_q(n)` as primitive-root power sums, and the bridge from divisibility data to
Fourier/cyclotomic packets.  The intended internal seed is the repeated repo
lesson that scalar quotients are useful only after labelled fibers are retained:
irreducible cores, unital designs, Faulhaber moment positivity, Pollock degree
jumps, unit-distance norm layers, tiling/solid analogies, and the current LRC14
dual stack.

Pending work:

1. Read the divisor-function page and one-hop related pages around multiplicative
   functions, Ramanujan sums, Dirichlet convolution, abundant/perfect numbers,
   and divisor summatory problems.
2. Build a small LRC14 audit comparing coarse divisor quotients against
   Ramanujan/cyclotomic packet quotients on AP/GW, K33, petal, covering, and
   HYP-2963-bank rows.
3. State an explicit quotient admissibility theorem target:

```text
Any quotient used to rule out an LRC14 counterexample must retain enough
divisor/cyclotomic phase data to distinguish AP/GW boundary atoms, positive
Toeplitz/Ramanujan exits, and K33/state-lift debts.
```

Expected falsifier shape: two rows with the same scalar divisor signature but
different LRC route.  If such pairs exist, the scalar quotient is demoted to a
feature and the retained quotient must include Ramanujan packet labels or exact
endpoint-owner/Farey data.

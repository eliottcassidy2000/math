---
id: HYP-8895
title: "Return semigroups record reachability, not mixed-sign noncancellation"
status: >
  PARTIAL / CORRECTED BY MISTAKE-234. Return lengths form an exact support
  semigroup and decide positive-coefficient constant-term reachability. They
  do not control cancellation for signed or complex coefficients.
source: boxeph-2026-07-21-S223; codex audit 2026-07-21
depends_on: []
related:
  - HYP-8878
  - HYP-8890
  - MISTAKE-234
script: 04-computation/one_dimensional_coprime_intervals_return_semigroup_boxeph_S223.py
output: 05-knowledge/results/one_dimensional_coprime_intervals_return_semigroup_boxeph_S223.out
---

# HYP-8895 — support reachability with a cancellation sidecar

For a finite Laurent support `S subset Z`, the return set

```text
R(S)={m>=0: 0 belongs to the m-fold sumset mS}
```

is an additive semigroup. For positive coefficients,
`CT(f^m)!=0 iff m in R(S)`. Pair periods, gcds, and Frobenius conductors are
therefore exact support data.

For mixed or complex coefficients, reachability is only necessary. With

```text
f(z)=z-z^(-1)+z^2-z^(-2),
```

the support is aperiodic and every length `m>=2` is reachable, but
`f(z^(-1))=-f(z)` forces every odd constant term to vanish. Any DvdK bypass
must retain channel phases or prove a separate noncancellation theorem.
HYP-8878's unique-minimum channel is one valid sufficient case; the general
complex saddle route HYP-8890 remains open.

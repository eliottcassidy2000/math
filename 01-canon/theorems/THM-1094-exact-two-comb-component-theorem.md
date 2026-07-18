---
id: THM-1094
title: Exact two-comb component theorem for the three-killer clustered stratum
status: ID RESERVED / PROOF BEING FROZEN — an exact 9,246,070-pair endpoint scan plus a uniform analytic tail proves the component inequality needed to repair THM-1061; the certificate and independent extremal replay are being written
source: codex-2026-07-18-S73 frontier audit
depends_on: [THM-1051, THM-1061]
related: [THM-1042, THM-1071, THM-1081, MISTAKE-163]
---

# THM-1094 — reserved exact two-comb component theorem

For every ten-speed core `P subset {1,...,12}` and ordered killers
`13 max(P)<k1<k2`, the exact core-safe remainder

```text
S(P) \ (D(k1) union D(k2))
```

has a component of length strictly greater than `1/(3k2)`.

The finite endpoint-order bank has 9,246,070 pairs and no failure.  The exact
extremal currently reduces to `3*k2*L=158/119` at `(k1,k2)=(153,159)`.
Pairs beyond the finite bank are covered by an explicit component-count and
mass inequality.  The full proof, hashes, and carrier audit will replace this
reservation before promotion.

If verified as stated, a third killer `k3>k2` cannot cover that component, so
the three-killer clustered stratum closes uniformly by measure alone.  This
does not address the four-killer all-scale bridge or global LRC(14).

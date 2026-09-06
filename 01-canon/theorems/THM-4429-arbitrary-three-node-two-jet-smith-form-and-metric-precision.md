---
id: THM-4429
title: "Arbitrary three-node two-jet Smith form and metric precision"
status: >
  PROVED by explicit determinantal-divisor witnesses, with independent exact
  Smith audits. Arbitrary three integer nodes and their complete p-adic
  precision law are closed. THM-4435 supplies all-node two-jet precision
  and refutes general metric-only four-node partitions; higher jets remain OPEN.
source: overnight-hexagon-sep05 research session
depends_on:
  - THM-4419-twojet-prime-wall-precision-and-dyadic-triple-smith-law
  - THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall
proof: 05-knowledge/results/three-node-smith-overnight-hexagon-sep05.md
script: 04-computation/three_node_smith_probe_overnight_hexagon_sep05.py
output: 05-knowledge/results/three_node_smith_probe_overnight_hexagon_sep05.out
independent_audit_script: 04-computation/three_node_smith_independent_overnight_hexagon_sep05.py
independent_audit_output: 05-knowledge/results/three_node_smith_independent_overnight_hexagon_sep05.out
---

# THM-4429 -- Arbitrary three-node two-jet Smith form and metric precision

**PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
[full proof and paired audit manifest](../../05-knowledge/results/three-node-smith-overnight-hexagon-sep05.md)
are part of this theorem. No external priority is claimed.

For arbitrary integer nodes x0<x1<x2, let a=x1-x0, b=x2-x0,
g=gcd(a,b), P=ab(b-a)>0, u=a/g, v=b/g and epsilon=gcd(3,g,u+v).
The map from integer polynomials of degree <6 to their values and first
Hasse derivatives at these nodes has Smith factors

```text
(1,1,D1,D2/D1,D3/D2,D4/D3),
D1=g*gcd(g,2), D2=gcd(g^4,6P),
D3=2*g^4*P*epsilon, D4=P^4.
```

Translation clears a unit two-by-two block. The remaining 69 minors reduce
to four explicit gcd ideals. Two degree-seven witnesses in the third ideal
leave exactly an exceptional factor two and, when appropriate, three;
a higher-degree witness decides whether that three persists. The determinant
alone does not determine this Smith list.

For a prime p, write the sorted pairwise difference valuations as (e,e,f),
f>=e (necessarily f>e at p=2). The determinantal valuations are

```text
d1=min(2e,e+v_p(2)), d2=min(4e,2e+f+v_p(6)),
d3=6e+f+v_p(2)+[p=3 and e=f>0], d4=8e+4f.
```

The full p-Smith list depends only on the distance tree, not higher unit
digits. Its largest exponent is

```text
L=2e+3f-v_p(2)-[p=3 and e=f>0].
```

Observations modulo p^(N+L) determine all six source coefficients modulo
p^N, and this uniform precision loss is sharp. Translations, nonunit scales,
permutations and small-prime mutations are independently audited. The degree
box and full observer are essential; no corresponding moving-module,
higher-jet, LRC(14), or JC(2) claim follows. The later
[THM-4435, four-node metric blindness and universal Hermite precision](THM-4435-four-node-metric-blindness-and-universal-hermite-precision.md)
extends the largest-factor formula to every node count but shows why the
full metric-only partition stops at three nodes in general.

**Separate audited four-node continuation:** the
[ternary double-pair theorem](../../05-knowledge/results/overnight2_20260906_smith_double_pair.md)
proves the full partition for occupancy `(2,2)` at prime three, including
arbitrary depth and unit lifts. Its weighted-minor certificate is independent
of this three-node proof; it preserves all competing intermediate ideals.
The other four-node shapes and higher jets remain separate questions.

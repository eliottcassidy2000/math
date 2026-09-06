---
id: THM-4419
title: "Two-jet prime-wall precision and dyadic triple Smith law"
status: >
  PROVED with independent exact audits. The full-residue s=p two-jet
  Smith partition and n<=p^2 consecutive extension independently reproduce
  concurrent synthesis work. New consequences give sharp p-adic precision,
  finite-level kernel counts, tensor-grid precision, and the all-depth
  dyadic three-node Smith partition. General multiscale clusters and
  higher jet orders remain OPEN.
source: synthesis-sep05 + independently concurrent synthesis_20260905 prime wall
depends_on:
  - THM-4080-confluent-two-jet-single-scale-smith-partition
  - THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall
proof: 05-knowledge/results/confluent-twojet-prime-wall-synthesis-sep05.md
script: 04-computation/confluent_twojet_prime_wall_synthesis_sep05.py
output: 05-knowledge/results/confluent_twojet_prime_wall_synthesis_sep05.out
independent_audit_script: 04-computation/confluent_twojet_prime_wall_synthesis_sep05_independent.py
independent_audit_output: 05-knowledge/results/confluent_twojet_prime_wall_synthesis_sep05_independent.out
---

# THM-4419 -- exact precision and the first dyadic three-node corner

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** The
[complete proof and evidence manifest](../../05-knowledge/results/confluent-twojet-prime-wall-synthesis-sep05.md)
credits [the concurrent full-residue theorem](../../05-knowledge/results/synthesis_20260905_wildcard_smith_boundary.md);
the common prime-boundary closure is one result, not two advances.

For p distinct residues u_i modulo p, nodes a+p^e u_i, e>=1, and
polynomials of degree below 2p, the matrix of values and first Hasse
derivatives has p-Smith exponents

```text
(0,0,e,...,(p-2)e,(p-1)e+1,
 (p+1)e,...,(2p-2)e,(2p-1)e-1).
```

Empty ranges are omitted. The modular derivative trace loses one direction;
dividing that trace by p recovers a unit Vandermonde minor. All intermediate
determinantal valuations gain exactly one, independent of e. Local CRT
therefore computes the full partition for consecutive nodes through n<=p^2.

Writing L=(2p-1)e-1, observations modulo p^(N+L) determine all source
coefficients modulo p^N, and this uniform loss is sharp. At level p^N
the kernel has p^(sum_i min(N,alpha_i)) elements. Independent rectangular
d-coordinate grids at depths e_1,...,e_d have all tensor sums as exponents
and sharp loss (2p-1)sum e_i-d. Coupled grids are outside this statement.

For the first two-scale family (0,2^e,2^(e+1)), the complete dyadic exponents
are, for every integer e>=0,

```text
e=0:  (0,0,0,0,2,2);
e=1:  (0,0,2,2,5,7);
e>=2: (0,0,e+1,2e+1,4e,5e+2).
```

After clearing the first node, all 69 minors of a fixed four-by-four
integer matrix give determinantal valuation envelopes
min(e+1,2e), min(3e+2,4e), 7e+2, 12e+4. Their consecutive differences
prove the formula at every depth. Adjoining a node need not preserve
earlier Smith factors; (0,8)->(0,8,16) is an exact hostile.

The later [THM-4429, arbitrary three-node Smith form](THM-4429-arbitrary-three-node-two-jet-smith-form-and-metric-precision.md)
closes every distinct integer triple and its metric-only p-adic precision
law; four-node and higher-jet continuations remain open.

Primary and independent SymPy/minor paths retain arbitrary residue lifts,
nonunit scales, small-prime hostiles, tensor and finite-kernel controls.
The main audit has 441,971 explicit checks. These precision results apply
to the declared integral degree box; they supply no automatic statement
for moving arithmetic modules, LRC(14), or the planar Jacobian conjecture.

## Audited odd-prime continuation, 2026-09-06

The [mixed-cluster theorem](../../05-knowledge/results/overnight_20260906_smith_mixed_cluster.md)
now proves the full partition for `p+1` nodes covering every residue with
one duplicated pair: odd prime `p`, outer depth `e>=1`, duplicate depth
`d>=1`, and arbitrary integral lifts. Competing minors select `min(e,d)`;
the old `2p`-entry Smith prefix survives adjoining exactly when `d>=e`.
The sharp precision loss is `(2p+1)e+3d`. Local CRT extends consecutive
nodes to `n<=p(p+1)` for odd primes. This is independently proof-audited
and supported by 1,821,396 exact checks; general cluster trees and higher
jets remain open. The original dyadic law is unchanged.

# A000568 correction by divisibility interactions

**Proved:** an exact signed correction calculus from a size-only upper envelope
of the Burnside sum; two-sided certified error bounds; an exact 3-cycle
resummation; a sharp four-term exactness threshold n<=23; and a proof that every
fixed correction depth eventually loses integer accuracy.

[Definitions and proofs](PROOF.md). Reproduce:

    python3 gas.py > results.json

All computations are rational. `results.json` stores exact signed errors,
lower/upper remainder certificates, and independent small-order controls.
The logarithms below are rounded displays only.

| n | plain upper error, log2 | 3-resummed upper error, log2 | resummed S3 error, log2 |
|---:|---:|---:|---:|
|20|62.8430|33.4830|exact|
|24|109.4279|65.9422|29.3690|
|32|246.6404|173.7002|74.4132|
|40|444.0148|340.8334|169.1189|

S0 is an upper bound; S3 a lower bound. Exact-residue recovery can consume the
intersection of these bounds, using Session277's existing uniqueness criterion.
The current C_j implementation enumerates partitions; no speed win is claimed.

The surprising cross-repo bridge is from the nested Brownian kernel in
THM-4346 (LRC sheet currents) to the arithmetic gcd kernel in tournament
Burnside symmetry. Their defect is positive on occupancy vectors but indefinite
as an operator. That distinction prevents an attractive but false spectral
analogy while yielding the useful combinatorial correction graph.

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
|64|not displayed here|1209.1966|810.2100|

S0 is an upper bound; S3 a lower bound. Exact-residue recovery can consume the
intersection of these bounds, using Session277's existing uniqueness criterion.
The current C_j implementation enumerates partitions; no speed win is claimed.

The surprising cross-repo bridge is from the nested Brownian kernel in
THM-4346 (LRC sheet currents) to the arithmetic gcd kernel in tournament
Burnside symmetry. Their defect is positive on occupancy vectors but indefinite
as an operator. That distinction prevents an attractive but false spectral
analogy while yielding the useful combinatorial correction graph.

## Certified interval, not just an error estimate

For each k<=3, combine the signed partial sum with BOTH bounds in(3), then
intersect all eight bounds and round inward to integer endpoints. The resulting
width W needs only `bit_length(W)` residue bits for a unique integer. At n=40
this requests142 bits of a621-bit answer, versus283 for the unresummed envelope.
This is a precision comparison, not a runtime comparison. The new script uses
the exactly enumerated answer's residue only as an independent uniqueness
**control**; it does not implement Session277's modular backend again or claim
its recovery run avoided exact partition enumeration.

At n=64, the1721-bit answer is confined to an interval needing738 residue bits,
versus1063 for the unresummed envelope and1544 in Session277's support-cut3
certificate. These are different certificate workloads; the new result does
not assert a matched runtime improvement. The full control run through64 took
about5.3 seconds in one local run; timing is descriptive only.

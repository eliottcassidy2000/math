---
id: THM-4456
title: "Sharp finite-length signed-root stability asymptotics"
status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED
depends_on:
  - THM-4454 / sharp-global-signed-root-duplication-stability
  - THM-4455 / three-atom-minimizing-sequence-rigidity
source: long-frontier-sep06 sharp dimension asymptotic
primary_script: 04-computation/long_frontier_sep06_dimension.py
primary_output: 05-knowledge/results/long_frontier_sep06_dimension.out
primary_script_sha256: badb2489aeb8e4a27a0e5b61c203323dfb988be693488637056a1f9f3c19e667
primary_output_sha256: d2e72609bb6b248bed6d721e75e0a2b0350316b046cd6bed5ad06ddae8973e30
report: 05-knowledge/results/long_frontier_sep06_dimension.md
audit: 05-knowledge/results/long_frontier_sep06_dimension_audit.md
hash_basis: raw LF repository bytes
---

# THM-4456 — sharp finite-length signed-root stability asymptotics

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
[complete proof](../../05-knowledge/results/long_frontier_sep06_dimension.md)
and [independent proof/source audit](../../05-knowledge/results/long_frontier_sep06_dimension_audit.md)
are part of this theorem. Its dependencies are
[THM-4454 / sharp global signed-root duplication stability](THM-4454-sharp-global-signed-root-duplication-stability.md)
and [THM-4455 / three-atom minimizing-sequence rigidity](THM-4455-three-atom-minimizing-sequence-rigidity.md).

Let I_N be the infimum of R=(J-c_*)/d2 from those theorems over finite
real lists of total length at most N, with p1=p2=1 and E>0. Then

    I_N=K3+C/N+o(1/N),
    C=-22-(16/3)sqrt(6)+10sqrt(2)+(40/3)sqrt(3)
     =2.1722010964723645447... .

This is an asymptotic theorem, with no claim of attainment or an exact
optimizer at a specified finite N. A three-equal-positive-atom family
with N-3 equal negative entries gives the sharp upper coefficient.

For any eligible sequence of length at most N tending to infinity,
let a>=b>=c be its three leading positive coordinates, ell their mean,
m the remaining square mass, and h=ell-c. The exact first-order
equality characterization is

    N(R-K3)->C iff Nm->(sqrt(3)-1)^2 and Nh->0.

Pad the tail to N-3 entries and put S=1-3ell. Such sharp sequences satisfy

    sum_tail (r_i-S/(N-3))^2=o(1/N),
    sum_tail |r_i-S/(N-3)|=o(1).

Consequently the actual used length and the number of negative entries
are both N-o(N). Positive tail mass tends to zero; negative magnitude
tends to sqrt(3)-1. This stronger first-moment conclusion uses the
sharp finite-length order and is not inferred from plain R->K3.

The mechanism is the length-uniform local expansion

    R-K3=A m+B h+O(m^(3/2)+h^2),
    A=-2-sqrt(6)/3+2sqrt(2)+(7/3)sqrt(3)>0,
    B=8-5sqrt(2)-4sqrt(3)+3sqrt(6)>0.

The ordered top-two distance creates the positive linear term B h.
THM-4455 places arbitrary approximate minimizers in this neighborhood;
local coercivity gives m+h=O(1/N). Signed Cauchy--Schwarz on the whole
tail gives (N-3)m>=(3ell-1)^2, proving the sharp lower coefficient
C=A*(sqrt(3)-1)^2. Its equality defect proves the tail conclusions.

The exact checker has117 always-active gates, with independent
atom-gradient and moment-expansion paths and twelve outward-rational
actual controls. Full independent proof/source audit and normal/optimized
replay passed. No finite parameter bank supplies an asymptotic quantifier.
No external priority, proof-assistant or actual Laurent-transport claim
is made.

---
id: THM-4427
title: "All cumulative signed-cycle gaps by transposition rigidity"
status: >
  PROVED with an independently audited FINITE-EXACT K8 base and analytic
  transposition rigidity at every Hamilton order n>=9. All cumulative
  D>=7 gaps and their equality classes are closed.
source: overnight-hexagon-sep05 research session
depends_on:
  - THM-4416-even-graph-cumulative-d5-d6-spectral-gap
  - THM-4078-even-graph-triangle-quotient-spectrum-and-boolean-noncommutation
proof: 05-knowledge/results/overnight_hexagon_sep05_d7_d8_gap.md
script: 04-computation/overnight_hexagon_sep05_hamilton_walsh.cpp
output: 05-knowledge/results/overnight_hexagon_sep05_hamilton_walsh8.out
independent_audit_script: 04-computation/overnight_hexagon_sep05_hamilton_direct.cpp
independent_audit_output: 05-knowledge/results/overnight_hexagon_sep05_hamilton_direct8.out
---

# THM-4427 -- All cumulative signed-cycle gaps by transposition rigidity

**PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
[full proof, controls, reproduction commands and manifests](../../05-knowledge/results/overnight_hexagon_sep05_d7_d8_gap.md)
are part of this theorem. The reservation is promoted after independent
mathematical and computational review.

For an edge signing H of K_n, modulo cut switching, let c_k count negative
unoriented simple k-cycles and let B,A be the balanced and antibalanced
classes. For every k>=8 and n>=k,

```text
c_k(H) >= (n-2)!/(n-k)! outside {B} for odd k, {B,A} for even k.
```

Equality consists exactly of single-negative-edge classes for odd k and
their union with globally negated single-edge classes for even k.
For every D>=7 and n>=D+1,

```text
min_(H not balanced) sum_(k=3)^(D+1) c_k(H)
  = A_(n,D) = sum_(k=3)^(D+1) (n-2)!/(n-k)!.
```

There are exactly binom(n,2) labelled equality classes, the single-edge
classes, and one relabelling orbit. In THM-4078's multiplicity-weighted
cycle operator the gap is 2A_(n,D). Together with THM-4083/4416 this closes
the entire previously open cumulative family D>=3.

Mechanism: a vertex swap changes H by a signed K_(2,r), with s=n-2-r.
Its exact Hamilton count is

```text
T_n(r)=2rs((n-2)(n-3)-2rs)(n-5)!.
```

For n>=9 and 2<=r<=n-4 this exceeds twice the single-edge count. Thus any
signing at or below that count has only r=0,1,n-3,n-2 pair disagreements;
a root-gauged graph classification leaves B,A and the single-edge classes
and their negatives. A two-type deletion argument transfers any Hamilton
base to all ambient orders. The K8 base is checked on all 2^21 switching
classes by literal parity, Walsh transform, and independent deleted-vertex
Hamilton paths. K9's complete 2^28 comparison is corroboration, not a
premise of the analytic n>=9 proof.

The negative C5 plus positive apex on K6 has c6=20<24; individual Hamilton
statements cannot be extended to every smaller length. Antibalance kills
all even layers, but paired odd layers pay its cumulative deficit strictly.
Booleanized adjacency, tournament H>=disc and LRC(14) remain outside scope.

**Additional independent audit, 2026-09-06:** the
[Gray-code bitset audit](../../05-knowledge/results/overnight2_20260906_signed_cycle_audit.md)
checks the complete analytic proof and rebuilds every K8 switching class
without a Walsh transform or a repository mathematics import. It recovers
both zero classes, minimum 720, and all 56 equality classes, with full
spectrum digest `e2dfba14125e7983`. This is additional evidence for the
same theorem, not a new mathematical claim or a K9 replay.

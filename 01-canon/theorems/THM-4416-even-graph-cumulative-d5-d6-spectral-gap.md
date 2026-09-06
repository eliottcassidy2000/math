---
id: THM-4416
title: "Even-graph cumulative D5/D6 spectral gap"
status: >
  PROVED with FINITE-EXACT local comparison lemmas and base cases, independently
  audited by direct cycle parity and integer Walsh transformation. For D=5
  at every n>=6 and D=6 at every n>=7, the cumulative negative-cycle minimum
  over nonbalanced switching classes is the single-edge value, attained only
  by single-edge classes. D>=7 remains OPEN.
source: synthesis-sep05 / recovered local-profile sidecar and deletion induction
depends_on:
  - THM-4078-even-graph-triangle-quotient-spectrum-and-boolean-noncommutation
  - THM-4083-even-graph-cumulative-d3-d4-spectral-gap
proof: 05-knowledge/results/even-graph-d5-d6-closure-synthesis-sep05.md
script: 04-computation/even_graph_d5_six_vertex_synthesis_sep05.py
output: 05-knowledge/results/even_graph_d5_six_vertex_synthesis_sep05.out
independent_audit_script: 04-computation/even_graph_d5_six_vertex_synthesis_sep05_independent.cpp
independent_audit_output: 05-knowledge/results/even_graph_d5_six_vertex_synthesis_sep05_independent.out
script_sha256: 48610079e151725cf3cbadadda1a205b131ac6ed194768cd210f6192a3d67eee
output_sha256: e4c2fd382669ebccc55a3b05e242fe9e5def36c096eba4daa4920b857cfb0035
independent_audit_script_sha256: ade5194ef6fda97d8c8cf8ecea4cb30bc120c62cab2b8a5f129ffc96a84f60d5
independent_audit_output_sha256: bf8e7499e7c1fbb470b0695789dbf5612e58b35b8e4a25c5f700ae049d1a5dca
hash_basis: raw LF bytes
---

# THM-4416 -- cumulative D5/D6 gaps

**PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [full proof and exact reproduction](../../05-knowledge/results/even-graph-d5-d6-closure-synthesis-sep05.md)
are part of this theorem. It supersedes the D5 open boundary in
[THM-4084](THM-4084-even-graph-matching-character-profile-and-d5-firewall.md)
and [THM-4200](THM-4200-even-graph-four-edge-d5-frustration-firewall.md).

Let c_k count negative unoriented simple k-cycles in a signed complete K_n,
and S_D=sum_(k=3)^(D+1)c_k. For D=5,n>=6 and D=6,n>=7,

```text
min_(nonbalanced classes) S_D = A_(n,D)
  = sum_(k=3)^(D+1) (n-2)!/(n-k)!.
```

Equality consists exactly of the binom(n,2) labelled single-edge switching
classes. In THM-4078's multiplicity-weighted cycle operator the Laplacian
gap is 2A_(n,D); its relabelling quotient has gap multiplicity one.

The decisive finite lemmas are 3c6>=2c4 on K6 and 2c7>=c4+c5 on K7.
Summing them over induced subsets gives c6>=c4 for n>=7 and
c7>=2c4+c5 for n>=8. Thus the higher layer pays the lower-layer surplus in
the exact vertex-deletion identity. Balanced deletions already characterize
single-edge and strict star classes by THM-4083. In the remaining b=0
case, an assumed S_D<=A_(n,D) contradicts, respectively,

```text
n A_(n-1,5)-(n-5)A_(n,5)-2binom(n,3)
  = (n-5)(n-4)(n-3)(3n-25)/3 > 0,

n A_(n-1,6)-(n-6)A_(n,6)-3binom(n,3)
  = (n-4)(n-3)(2n^3-42n^2+285n-620)/2 > 0,
```

for n>=9. In the second line the cubic at n=9+t is
2t^3+12t^2+15t+1. Independent exhaustive implementations verify the
n=6,7,8 bases, equality masks, and both local lemmas on 2,130,944 total
switching classes. Induction is a paper proof, not finite extrapolation.

The antibalanced signing still defeats individual even-layer lower bounds.
D>=7, Booleanized quotient adjacency, tournament H>=disc, and LRC(14)
are outside this theorem.

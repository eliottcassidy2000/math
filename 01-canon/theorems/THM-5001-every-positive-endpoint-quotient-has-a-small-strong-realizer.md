---
id: THM-5001
title: "Every positive integer is an endpoint-stabilizer quotient of a strong tournament of logarithmic order"
status: PROVED, independently audited. Endpoint spectrum only; global Chan-Pak completeness remains OPEN.
session: Session233-composition-20261109
source: 05-knowledge/results/session233_composition_20261109/endpoint/PROOF.md
audit: 05-knowledge/results/session233_composition_20261109/endpoint_audit/REPORT.md
---

# Complete endpoint spectrum

For every integer N>=1 there is a strong tournament U and a vertex x with

- exactly N directed Hamilton paths ending at x;
- trivial vertex stabilizer Aut(U)_x;
- at most max(5,ceil(log2(N+1))+3) vertices.

Thus the fibre over the endpoint orbit has exactly N path orbits. This closes
the endpoint-spectrum question, not the complete-path or complete-path-orbit
spectrum questions.

**Construction.** Take n=max(4,ceil(log2(N+1))+2) and the transitive tournament
on 0,...,n-1 with only 0->n-1 reversed. This strong rigid base has endpoint
vector (1,1,2,4,...,2^(n-3),1). Adjoin x, with incoming neighbors among the
interior vertices selected by the binary digits of N, and x beating all
other vertices. Deleting terminal x gives the required sum. The base is
strong and x has incoming and outgoing neighbors. Its score classes have
size at most two, and no automorphism can swap an oriented pair, proving
base rigidity and hence trivial stabilizer of x.

**Exact response sidecar.** If B_v(A) counts ordered spanning two-path covers
whose distinguished path starts at v (the other may be empty), then
B_x(A+x)=H(A)+sum_(x->v) B_v(A). This follows by removing the first vertex
of the distinguished path, with the singleton path handled separately.
The endpoint selector uses binary weights, but this response uses different
weights. It supplies no independent control of the full tuple (H,e,B).

[Proof, exact examples and inheritance](../../05-knowledge/results/session233_composition_20261109/endpoint/PROOF.md) · [Independent audit](../../05-knowledge/results/session233_composition_20261109/endpoint_audit/REPORT.md).

The near-transitive total-count family is inherited from THM-250-flip-formula.md;
the endpoint orbit interpretation is THM-4988-endpoint-721-boundaries-and-faithful-automorphism-graphs.md.
No external novelty or minimal-order claim is made.

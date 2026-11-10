---
id: THM-5090
title: "Universal induced odd-cycle carriers and the positive completion tax"
status: "PROVED; independently audited. Chan–Pak completeness remains OPEN."
session: Session237, 2026-11-09
source: 05-knowledge/results/session237_chanpak_synthesis_20261109/README.md
proof: 05-knowledge/results/session237_chanpak_synthesis_20261109/cycle_tilings/REPORT.md
audit: 05-knowledge/results/session237_chanpak_synthesis_20261109/cycle_tilings_audit/AUDIT.md
---

# THM-5090: Universal induced odd-cycle carriers and the positive completion tax

Every finite simple graph G with m vertices and e edges occurs induced in the
odd-directed-cycle conflict graph Omega(T) of a tournament T. The explicit
construction has e+5m+o vertices, where o is the number of odd-degree vertices
of G, and includes the corona G o K1 as an induced subgraph. Every completion
of its prescribed arcs satisfies the coefficientwise inequality

    I(Omega(T),x) >= (1+x)^m I(G,x/(1+x)) + m*x.

Thus H(T)>=3^m I(G,2/3)+2m. This is a carrier embedding, not an H-preserving
reduction. Its hereditary closure is all graphs, while full realization is
not universal (K3 would give the impossible H=7). Extra cycles and their
conflict neighborhoods are the required sidecar.

The proof uses the inherited Grinberg–Stanley odd-cycle formula, THM-002,
and exact arc-compatible edge-intersection construction. This does not
transfer THM-5040's even-graph domination spectrum to tournament paths.
The proof packet credits earlier induced C5/claw constructions and distinguishes
automorphism-even graphs from degree-even graphs.

[Full proof](../../05-knowledge/results/session237_chanpak_synthesis_20261109/cycle_tilings/REPORT.md) · [Independent audit](../../05-knowledge/results/session237_chanpak_synthesis_20261109/cycle_tilings_audit/AUDIT.md).

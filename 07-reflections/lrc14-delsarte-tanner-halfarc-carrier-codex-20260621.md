# LRC14 Delsarte/Tanner Half-Arc Carrier

**codex-2026-06-21.**  The Tanner graph prompt produced a useful guardrail:
the THM-534 Delsarte dual can be drawn as a bipartite check/atom graph, but
the resulting graph is not a sparse local code graph.  It is the triangular
incidence `r<=t` between factorial moments and missed-depth atoms, so it is
dense, has girth 4, and carries many 4-cycles.  Importing an LDPC expansion
argument would be the wrong move.

The positive residue is the sign split.  The undirected support has
automorphisms, but no automorphism orbit mixes positive and negative weighted
edges.  This is the exact sense in which the Doyle-Holt half-arc analogy
survives: orientation/sign is a real invariant, not a decoration that can be
forgotten and restored later.

The best next target is therefore not a Tanner theorem.  It is a parity
puncture/extend lemma for the Delsarte duals:

```text
K8  : even-only Krawtchouk support K0,K2,K4;
K9  : mixed low-degree support K0,K1,K2,K3;
K11 : mixed low-degree support K0,K1,K2.
```

That looks much closer to Li-style Delsarte parity and to the repo's half
tiling even/odd split than to a sparse Tanner graph.  If this can be spliced to
HYP-2737's row-slice odometer word, the final proof interface becomes cleaner:
generated word -> depth quotient -> parity/sign-rigid Delsarte cap.

The unit-distance prompt points the same way.  Unit-distance Delsarte/Hoffman
bounds are spectral and global; the local graph is valuable only when it
remembers the right side channel.  Here the side channel is the signed
Krawtchouk/binomial orientation, not raw adjacency.

# LRC14 Signed Tanner/Dessin Weak-Regularity Audit

**codex-2026-06-21-S73.**  The new prompt mixed Delsarte extremal functions,
Tanner graphs, unit-distance graphs, weakly regular graphs, the Doyle-Holt
half-arc split, and Belyi/modular-form language.  The exact audit says this is
useful, but only after changing the question.

The THM-534 Delsarte duals do have a Tanner-style carrier: moment rows `r` on
one side, missed-depth atoms `t` on the other, and edge `r--t` iff `r<=t`,
weighted by `y_r binom(t,r)`.  But it is a dense Ferrers carrier, not an LDPC
expander and not weakly regular.  Its degree quotients are equitable, while
common-neighbor counts form a chain rather than a two-value regularity law.

The surviving invariant is signed orientation.  Support automorphisms never
mix positive and negative edge orbits:

```text
K8 : automorphisms=6,   edge-orbits=15, mixed-sign=0
K9 : automorphisms=24,  edge-orbits=10, mixed-sign=0
K11: automorphisms=120, edge-orbits=6,  mixed-sign=0
```

That is the precise Doyle-Holt-style residue: edge support and signed arcs are
different data.  A natural bipartite rotation system also gives a dessin
passport, but the unsigned passport remembers branch degrees, not the
Delsarte inequality.  The signed passport is the proof-relevant object.

This changes the proof-order language after the HYP-2741 correction.  L7 is
not a gap-free theorem just because the discrepancy is exact: HYP-2739 closes
the residue numerator, and HYP-2741 fixes the finite-f1 rate, but r>=3
pairwise lift, base-size domination, and consec-maximality remain as real
obligations.  The graph analogies should feed the remaining generated-depth
and Delsarte parity layer:

```text
generated depth word -> signed Delsarte/quasicode parity -> aggregate cap ledger.
```

Do not replace that with a sparse Tanner theorem or a unit-distance regularity
claim.  The computed carrier explicitly lacks the local regularity those
proofs would require.

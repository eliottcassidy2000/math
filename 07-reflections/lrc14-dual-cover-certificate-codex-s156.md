# LRC14 Dual Cover Certificate

Date: 2026-06-24
Agent: codex-2026-06-24-S156
Related: HYP-2970, HYP-2969, HYP-2968, HYP-2966, HYP-2965, HYP-2252, HYP-2082, HYP-1890

This pass intentionally avoids extending the family classifier.  The new object
is the open-cover transition graph of the danger arcs.

The useful formula is the endpoint credit

```text
K((s,m),(r,n)) = 14*(r*m - s*n) + r + s.
```

Its sign is exactly the sign of the endpoint gap after multiplying by
`14rs`.  The zero case forces `r+s == 0 mod 14`, which recovers the AP/GW
boundary-pair skeleton from pure endpoint arithmetic.

The proof shape is graph-dual:

```text
strict counterexample
  <=> open endpoint graph has a unit-winding cycle
  <=> no potential Phi satisfies Phi(b)+epsilon<=Phi(a).
```

That is the Farkas certificate the older endpoint-protection notes wanted, but
stated on a concrete finite graph.  AP/GW should be closed-graph winding atoms
with zero credit and no open winding cycle.  Any genuine strict counterexample
must instead supply a positive-credit winding SCC, and that SCC is the object to
feed into K33/state-lift checks.

The next computation should be tiny: build this graph for named rows and run a
positive-winding cycle detector plus the dual potential solver.  If AP/GW are
the only closed-unit/open-empty atoms locally, HYP-2970 becomes a clean
certificate language for the global proof.

# LRC14 Circuit-Complexity Past-Work Synthesis

The useful circuit-complexity question is not "how small is the classifier for
the current finite bank?"  HYP-3115 already shows why that is dangerous: one
tuned literal, `apex7_error<=5`, isolates the top bounded-bank row.  That is a
signal, not a proof.

The better question is: what is the proof circuit's gate basis, and which
coordinates does each gate retain?  A past-work search gives a surprisingly
concrete answer.  HYP-2112 `Phi` is the numeric output gate.  HYP-2108 `P` is
the endpoint-cover sign gate.  HYP-2109's `L/M/R` automaton is the wall-crossing
state gate.  HYP-3023's magnitude cocycle is the route-purity gate.  HYP-3077's
Horn closure is the legality gate.  HYP-3082's protected branch graph is the
terminal no-naked-bridge gate.

That framing also keeps the broader invariant language honest.  Equinumerosity,
equidecomposability, and equidistribution are not interchangeable shadows:
fiber cardinality, scissors decomposition, and measure/root distribution lose
different coordinates.  A theorem-safe circuit can combine them only after
declaring the lost coordinate and the sidecar that repairs it.

The executable synthesis in HYP-3116 ranks proof carriers rather than runner
rows.  Its proof-payload tournament is acyclic and puts protected branch safety,
Horn sidecar closure, exact `Phi`, and magnitude-cocycle route purity ahead of
the smaller fitted classifier.  The gauge check is the real result: "smallest
circuit first" flips 38 pairwise decisions against the proof-payload gauge.
This says a tiny classifier is the wrong optimization target until uniformity
and destroyed-coordinate discharge are formal inputs.

The next useful computation is a HYP-2963 row ledger with columns for `Phi`,
`P`, endpoint owners, `L/M/R` state, magnitude cocycle, automatic word,
root/ear payload, Horn closure, protected-branch status, proof-depth stage, and
terminal exit.  The proof target is then sharp: every residual row must close
through one of these gates or emit named THM-572/F7 debt.  Anything that closes
only through a scalar or a finite threshold is an unsafe quotient.

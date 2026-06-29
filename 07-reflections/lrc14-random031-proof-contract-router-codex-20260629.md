# LRC14 random031 proof-contract router

Reserved for a proof-interface synthesis after HYP-3523 and HYP-3524.

The working claim is that random031 is now best treated as a finite contract
between terminal clauses and sidecars.  The stream and emitter audits have
made the residual small, but the remaining work is formal, not cosmetic: show
which clause consumes which certificate, which tail remains, and which
quotients are illegal because they forget owner labels, route sidecar `R`,
free-hole carry, or the vertical-halfturn guard.

Assumption challenge: the tournament vertices should be proof contracts,
sidecar obligations, and quotient-forbidden cuts.  Runners, arcs, raw witness
counts, and chamber volumes are secondary shadows unless the contract says
they reconstruct the terminal predicate.

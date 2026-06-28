# LRC14 Actual-Packet Sheaf Instantiation

HYP-3301 ended with a clear next hook but still on a toy matrix: instantiate
the first-obstruction / cusp-transfer language on actual packet rows.  The
natural bank was already in the repo: the curated HYP-2969 boundary-moment
ledger, which imports HYP-2963 packet labels and theorem exits.

The useful surprise is how small the real ambiguity is on that bank.

The coarse HYP-3301 packet

```text
(q-threshold bucket,
 six-unit contact profile,
 strict-safe-mass zero/nonzero,
 state-lift flag)
```

does not explode into many mixed fibers.  It leaves exactly one: seven
`q=14`, six-unit-boundary, positive-open, non-state-lift rows mixing
unit-petal-named exits with positive-Haar-open covering exits.

That is already a concrete first-obstruction statement.

The sharper surprise is the repair.  I expected the HYP-3310 `v2` layer to be
immediately necessary, but on this curated bank the nonunit residue word mod
`14` already kills the ambiguity completely.  The `v2` word alone is weaker.
So the first actual-packet obstruction here is not an unnamed qdiv>14 kernel;
it is a finite covering-layer residue-word sidecar.

This is exactly the kind of bank-local theorem-facing fact the repo wants:
honest, finite, and useful.  It also stays within the current guardrails.
HYP-3260 and HYP-3310 already warn that residue data is not globally enough:
same-residue height moves exist, and off-grid floor / endpoint-owner sidecars
still matter outside this bounded sample.

The next worthwhile extension is obvious now.  Push the same instantiation out
from the curated HYP-2969 bank to a larger HYP-2963 residual sample and ask
where residue-word exactness first breaks.  That breakpoint should name the
next real sidecar: `v2`/height, endpoint owner, or off-grid floor.

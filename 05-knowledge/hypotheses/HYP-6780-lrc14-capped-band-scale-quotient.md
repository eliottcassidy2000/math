# HYP-6780 -- The capped-envelope residual is scale-quotient, not a raw finite band

**Status:** reserved; exact scaling lemma proved, quotient classification open.

For every finite positive core `P` and integer `c >= 1`, the `1/14`-good set of
`cP` is the inverse image of the good set of `P` under the degree-`c` circle
cover. Consequently its measure is unchanged, its number of interval
components is multiplied by `c`, and the THM-755 cutoff satisfies

`v*(cP) = c v*(P)`.

Thus a per-core finite interval `(max(P), v*(P)]` does not imply a uniform raw
speed bound. The primitive covering ray

`V_c = {c, 2c, ..., 12c, 13c+1}`

with `13 | c` and `(2 | c or 7 | c)` is unbounded, satisfies the top-peel band
condition for every such `c`, and is already closed by THM-757 with
`M(V_c)=1/13`. This refutes only the *raw finite-enumeration interpretation*,
not LRC(14) or the capped-envelope theorem.

The open structural target is a finite or recursively controlled
classification after quotienting by core scale, normalized shape, and the
extra runner's residue/offset. Computation, a tournament audit over candidate
quotients, and corrections to the current frontier maps are in progress.

**Reserved by:** codex-2026-07-14-S1.

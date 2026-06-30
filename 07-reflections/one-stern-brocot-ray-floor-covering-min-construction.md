# One Stern-Brocot ray: the floor, the covering-min, and the construction

*mac-mini-2026-06-30-S52. A reflection on where the irregular covering-min actually lives.*

For weeks the LRC covering work has produced what looked like three unrelated rational quantities at each
level `n`:

- the **floor** `1/n` (the conjecture's bound),
- the **construction** `n/Phi_6(n)` (S44-S46's hexagonal/Eisenstein object, which S47 then showed is *not*
  the covering-min),
- the **primitive covering-min** `M_prim(n)` (this session's ILP: `2/13, 2/15, 4/33, 4/37, 3/31, 1/11, ...`),
  which looked hopelessly irregular.

They are the **same object at three depths.** Write the continued fraction `[0; n-1, k] = k/((n-1)k+1)`. As
`k` runs from `1` to `infinity` this traces a single Stern-Brocot ray from `1/n` up to `1/(n-1)`:

- `k=1`: `1/n` -- the floor.
- `k=a(n)`: the covering-min -- the **smallest depth a primitive covering set can reach**.
- `k=n`: `n/Phi_6(n)` -- the construction.
- `k -> inf`: `1/(n-1)` -- the top of the interval.

So the construction was never a rival theory; it is the **deep end** (`k=n`) of the ray whose **shallow end**
(`k=a(n)`) is the covering-min. The floor is the ray's foot. The "irregular" covering-min is just *which rung*
of one ladder the covering constraint can stand on -- and that rung index `a(n) = 2,2,4,4,3` (for `n=7..11`)
is the irregular degree of freedom.

*(Correction, MISTAKE-088: I first thought the covering-min "settled onto the top rung `1/(n-1)` for
`n>=12`". It does not -- that was a search-bound artifact. The construction sits at depth `k=n` with value
`n/Phi_6 < 1/(n-1)` for every `n`, so the covering-min is always strictly below `1/(n-1)`; for `n>=12` its
exact rung is open, needing a search out to speed `~n(n-1)`. The ray frame survives the correction; the clean
`n>=12` story did not.)*

The lesson keeps repeating in this project: **what looks irregular is a single structured object viewed in the
wrong coordinates.** The hexagon's "revenge" (HYP-3726) was the margin wearing five hats; the parity was the
bipartiteness of `C_n` (HYP-3729); and now the covering-min's irregularity is the achievable depth on one
Stern-Brocot ray. Each time, the continued fraction / triangle / cycle was already the right coordinate -- the
mathematics had drawn the ladder; we were reading the rungs out of order.

Two things the ray buys concretely. First, it says the covering-min and the construction should never have
been argued *against* each other (S47's "refutation" was really a *re-indexing*: depth `a(n)` vs depth `n`).
Second, it predicts the **achievability threshold** is the whole game: the covering-min is the lowest rung you
can build a covering set on, and proving the LRC floor is proving you can never build one on the `k=1` rung
(`1/n`) itself -- you always fall at least to `k=a(n)` (or, for `n>=12`, all the way to the top `1/(n-1)`).
The irregular `a(n)` is the fingerprint of *how far up the achievability forces you* before the divisibility
covering can be satisfied.

What's still open is exactly the achievability map `k -> (is [0;n-1,k] realizable by a primitive covering
set?)` -- non-monotone, number-theoretic, and (for the LRC) the thing that must be bounded away from `k=1`.
That map is the covering-min's true content, and it is a clean question now that the ray is the coordinate.

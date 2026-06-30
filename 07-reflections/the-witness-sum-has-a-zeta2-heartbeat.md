# The witness sum has a ζ(2) heartbeat

*klein-2026-06-30-S43. A reflection on HYP-3743.*

We set out to sum the lonely-runner witness hierarchy across primes — to combine the one-prime lower bounds
`M ≥ (r+1)/p` into a global statement. The sum turned out to be a *Farey object*. The dense core `{1,…,m}`
blocks the radius-`r` witness at a prime `p` exactly when the fractions `{d/j : j ≤ m, d ≤ r}` exhaust the
units of `Z/p`; so the reach is `2·F(m,r)+1`, twice the number of distinct fractions in the `m × r` grid, plus
one. The lonely runner's lower bound is bookkept by counting Farey fractions.

And the moment it became a fraction count, `ζ(2)` walked in. The density of distinct fractions in the grid is
`6/π² = 1/ζ(2)` — the coprimality density, the chance two integers share no factor. The collisions that thin
the grid below the full rectangle (`2/2 = 1/1`, `2/4 = 1/2`) are exactly the non-coprime pairs, and they cost
the reach a factor of `ζ(2)`. The witness sum doesn't approximate `ζ(2)`; it *is* a coprimality count.

This is the same `ζ(2)` the project met at the very beginning, in the floor bound `1/(2ζ(2))` (task #11, long
before any of this witness machinery existed). Back then `1/(2ζ(2))` was an empirical infimum of a normalized
covering radius, with no obvious mechanism. Now it reappears as half the witness-reach density — the `±` of
`2F+1` is the factor of two. Two completely different routes into the same problem, separated by months of
work, both deposit `ζ(2)` at the floor. When that happens it is not a coincidence; it is the problem telling
you where its constant lives. The covering radius of a lonely runner is, at bottom, a question about how often
fractions collide — and the answer to *that* is always `ζ(2)`.

There is a discipline in this that keeps paying off: when a bound resists a closed form, do not push harder on
the bound — ask what it is *counting*. The Farey rung (HYP-3732) counted Stern-Brocot depth; the witness sum
counts coprime fractions. Each time, the irregular analytic object dissolved into a clean arithmetic count, and
the count carried a famous constant (`Φ_6`, the golden ratio, now `ζ(2)`). The lonely runner keeps proving to
be number theory wearing a geometry costume. Follow the count, and the constant is waiting.

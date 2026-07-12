# Isolating the cancellation: a partial escape from the signed-cancellation wall

*klein-2026-07-11-S254, on THM-717.*

For nine independent decompositions across two routes (S222–S226, the "standing law"), the
verdict on LRC(14)-S3 was uniform: **no absolute / order-blind bound survives contact with the
problem's resonance objects — only signed methods work.** Every absolute Bonferroni / envelope /
count-cap died the same way: the true quantity survives by cancellation among terms of opposite
sign, and any bound that discards the signs overshoots by 10–150×.

THM-717 does not break this law — but it *quarantines* it. The k=9 base functional
`J = E[N(7−N)]` Abel-decomposes exactly into

    J = (6T₁ + 4T₂ + 2T₃)  −  (2T₅ + 4T₆),
        \_____ POS _____/     \___ BUNCH ___/

and the entire signed content lives in BUNCH. POS is a nonnegative combination of *monotone*
covering tails — no cancellation, so absolute Bonferroni applies to it freely. BUNCH is the
signed/cancelling part, but it is (a) small (≤ 2/7), (b) sharp only at the AP, and (c) a
"measure of near-origin bunching" that decays with the diameter. The same split works for the
k=8 deg-3 row, with BUNCH even smaller (≤ 0.027).

The lesson refines the standing law rather than overturning it:

> The signed cancellation that defeats absolute bounds is not spread uniformly through the
> functional — it concentrates in the **high-emptiness tail** (all runners bunched into ≤ 2 of
> the 7 sectors). Peel that tail off algebraically and the remaining bulk is cancellation-free.

Why the tail? `N(7−N)` is unimodal (peak at `N = 3,4`), so it *dips* at both ends. The `N = 0`
dip (over-covering) contributes zero and needs no correction. The `N = 5,6` dip (heavy bunching)
is where `N(7−N)` falls back down, and matching a monotone lower envelope there requires the
negative correction. That correction is exactly the bunching measure — the one genuinely
diameter-sensitive, AP-extremal, LRC-adjacent piece.

So the search for an absolute proof of the base was never hopeless — it was mis-scoped. The
right object to bound absolutely is POS (a coupled covering floor, Bonferroni-safe); the right
object to bound by a *separate*, diameter-aware argument is BUNCH (the near-origin bunching
`p₅ + 3p₆ ≤ 1/7`). Two clean sub-problems replace one contaminated one. This is the shape of the
escape: not a signed method that beats the cancellation, but an algebraic split that isolates it
into a term small enough to bound crudely.

Open question this raises: is BUNCH — "all `k` runners fit in an arc of length `2/7`" — itself a
shadow of the LRC covering problem one scale coarser? Its extremal is the AP, its mechanism is
three-gap bunching, and its value `1/7` is one sector. If so, the base's residual hardness is a
self-similar copy of the whole, and the moment ladder's last mile is genuinely irreducible to
elementary bounds — exactly what the isolated-tight-point picture (THM-708/709) predicts.

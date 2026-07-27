# kps-S135's first-untested-point law is the 2r rank-determination bound (and where it meets the nullcone rank)

**opus-2026-07-26. Provenance note, not truth source.** A one-mechanism sharpening of kps-S135's
"first-untested-point law" (minimal interpolants die at the first untested point; 24 breaks / 6
survivors), tested, with the over-strong form honestly released.

## The mechanism (✓ tested)

kps-S135 already recorded the decisive case verbatim — *"diag-sum order-3 recurrence, 3 points, broke
at m=8, rule: need 2r+2 terms."* That is the general law: **a linear-recurrence (rank-`r`) sequence is
not determined by fewer than `2r` terms** — any `(2r-1)`-prefix has a strictly-lower-order continuation
that agrees on the prefix and diverges at term `2r`. Demonstrated by Berlekamp–Massey
(`interpolation_rank_law_opus`): Fibonacci (`r=2`) locks to order 2 only at `k≥4=2r`; Tribonacci
(`r=3`) locks only at `k≥6=2r`; every shorter fit reports a lower-rank *imposter*. So **"the first
untested point" = the `≈2r` rank-determination threshold**: a minimal interpolant is a rank-`<r`
imposter, and it dies exactly where rank `r` first becomes visible. This is why the repo's
instrument-validation gate (S399) must "test past `2r`," and why minimal fits (the shortest formula)
are selected *for* early death (kps-S135 §2 — the shortest formula has the most equally-short
competitors).

## Where it meets — and does NOT meet — the nullcone rank (⌦ strong form released)

Tempting speculation: this "interpolation rank" is the same as the nullcone rank (LRC's spectral
`7⊗13`=2, JC's geometric dim-3; opus spectral-vs-geometric synthesis 2026-07-26). **Tested and
released.** kps-S135's break-instances are fits of *tournament/combinatorial sequences* (pure-blue
interleave, Paley `H`, diag-sum); their rank is the order of the true recurrence / degree of the closed
form — an **interpolation rank of a fitted sequence**. The nullcone rank is an invariant of the
*problem's obstruction*, not of any fitted sequence. The two coincide **only when the sequence being
fit is itself the nullcone spectrum** — e.g. the Farey/Ostrowski rung sequence `s = 55,69,83,…`
(`M=D/s`, the LRC spectrum): there the recurrence order *is* the spectral rank. Elsewhere they are
independent invariants. So the honest statement is:

> **WEAK (kept):** first-untested-point death = the `2r` rank-determination bound; validate past `2r`.
> **STRONG (released):** interpolation rank ≡ nullcone rank — false in general, true only on the
> spectrum-sequence.

## Takeaway

The first-untested-point law is not luck and not mysterious: it is the Berlekamp–Massey / interpolation
bound wearing a retrospective hat. The transferable rule is a number — **collect `2r` before trusting a
rank-`r` guess** — and the one place it fuses with the deep structure is when the sequence you are
fitting is the problem's own spectrum.

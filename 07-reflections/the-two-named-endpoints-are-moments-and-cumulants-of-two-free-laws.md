# The two named endpoints are moments and free cumulants — of two *different* free laws

*monad-explorer-2026-06-07 (deep-research, 13th session). Builds on THM-438 ADDENDUM-4
(two-point law free cumulants) and ADDENDUM-10 (the diagonal is A088368 = noncrossing lists).*

## The thing that clicked

The cycle-rank triangle `t(k,m)` of THM-438 has, after twelve sessions, exactly **two**
sequences on it that OEIS recognises (ADD-10):

- the **diagonal** `t(k,k) = A088368(k) = Σ_{π∈NC(k)} ∏_B |B|! ~ e·k!`  (the "wild" end), and
- the **signed row sum** `Σ_m (−1)^m t(k,m) = (−1)^k C_k`  (Catalan, the "tame" end).

Everything between them is uncatalogued. ADD-10 phrased this as "the Catalan is a cancellation,
not a count: both endpoints named, the path anonymous." True — but it left the two endpoints
looking like *unrelated* coincidences (one a Callan partition statistic, one a Catalan number).

They are not unrelated. **Both are free-probability objects.** And once you see that, the
"anonymous path" stops being a mystery and becomes a *theorem about why no path could be named.*

## Two canonical sequences, two different laws

Free probability attaches to every law `μ` exactly two sequences: its **moments** `m_k` and its
**free cumulants** `κ_k`, tied by Möbius inversion on the non-crossing lattice,
`m_k = Σ_{π∈NC(k)} ∏_B κ_{|B|}`, equivalently `M(z)=C(zM(z))`.

The triangle shows both — but **each endpoint belongs to a different law**:

| endpoint | free-prob role | the law | source |
|---|---|---|---|
| diagonal `t(k,k)=A088368(k)` | **moments** | the *factorial law* `κ_n=n!` | ADD-11 (this session) |
| signed sum `(−1)^kC_k` | **free cumulants** | the *two-point law* `½(δ_a+δ_{−a})` | ADD-4 |

The diagonal identity `Σ_{NC(k)}∏|B|!` is *literally* the free moment–cumulant formula with
`κ_n=n!`. That is not an analogy — substitute and the functional equation `M=F(zM)`,
`F=Σ n!w^n`, is Callan's A088368 g.f. on the nose, which **proves** the closed form (no
bijection, no number theory). The signed sum is *literally* a free cumulant of the arcsine /
two-point spectrum every doubly-regular tournament carries.

## Why the path is anonymous — now with a reason

If the two endpoints were the moment **and** the cumulant of *one* law, the connecting transit
would be that law's moment↔cumulant Möbius transform — a structured, nameable object. They are
not. The diagonal is moments of the *factorial* law; the signed sum is cumulants of the
*two-point* law. There is **no single law** whose `(moments, free cumulants)` pair are these two
sequences. So the deterministic Möbius transit the triangle performs from one endpoint to the
other is not the moment-cumulant map of any law at all — and there is no reason for it to pass
through anything catalogued. **The structurelessness of the path is the shadow of the endpoints
belonging to two laws, not one.** That is the cleanest statement I can give of "the Catalan is a
cancellation, not a count."

## A second face of an old wall

The same identification dissolves a separate puzzle. ADD-6 could not find a finite catalytic /
Tutte equation for the bivariate `U(x,y)`; the obstruction was that `A088368~e·k!` makes `U`
Gevrey-1 / **resurgent**, not algebraic. The factorial-law reading *names that wildness*: the
law's free-cumulant generating function is `R(z)=Σ_{n≥1} n! z^{n-1}` — the canonical divergent
factorial series, the textbook resurgent object. ADD-4's free probability and ADD-6's resurgence
were never two findings; they are **one object seen twice**: the diagonal is the moments of a law
whose free cumulants diverge factorially, and that divergence is exactly the analytic wall ADD-6
hit. A constant (`e`) and a wall (resurgence) turn out to be the same fact — the factorial growth
of `n!` — wearing the hat of a moment growth-rate in one place and a Gevrey class in the other.

## What the project keeps doing

This is the project's signature move again: a quantity that looked like a brute over-count (ADD-3
called the diagonal "the all-pairings overcount, uncatalogued") turns out, two corrections later
(MISTAKE-060/061/062/063, then ADD-10's naming, then this), to be a **named object with a
generating-function meaning** — here, free moments of a specific law. The mathematics does not
tolerate a permanently anonymous endpoint; pushed hard enough, both ends of the cancellation name
themselves, and naming them reveals the cancellation's *grammar*: moments of one law, cumulants
of another, with a transit between that belongs to neither. The wild end (`e`, resurgence,
factorial cumulants) and the tame end (Catalan, two-point spectrum) are the two halves of free
probability's single duality, split across two laws — which is why the triangle had to spend its
entire interior being structureless to get from one to the other.

*Open thread for the next explorer: the factorial law `κ_n=n!` has moments `A088368~e·k!`,
suspiciously close to the exponential law's `∫₀^∞ x^k e^{−x}dx = k!`. Is `κ_n=n!` the free
cumulant sequence of a *named* distribution? If so, its analytic structure may hand over the
off-diagonal columns `P_m` (currently OEIS-negative) for free.*

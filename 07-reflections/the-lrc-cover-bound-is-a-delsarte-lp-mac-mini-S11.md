# The LRC cover bound is a Delsarte linear program

**Session:** mac-mini-2026-06-21-S11. The user pointed at MDS codes, arcs in projective geometry,
and a 56-coincidence (size-3 challenger shapes = 56 = tournaments on 6 vertices). kps had just
opened HYP-2723 reading the carrier error as the weight enumerator of a relation code. Following
both into the coding-theory literature collapsed three separate framings of the arithmetic crux into
one classical object.

## The unification

The three live readings of "consec maximizes the Z/7 cover" are one Delsarte linear program:

- **kps HYP-2723** — the relation code `Λ(E) = {n : Σ nᵢeᵢ = 0}`, a `[k, k−1, d]` code whose
  minimum support `d` is small for consec (anti-MDS, many short relations) and large for Sidon/arcs
  (MDS, general position). This is the **Delsarte scheme**.
- **my HYP-2724** — `measS7 = P(N=0)`, a functional of the residue-depth law `π_E`, read through the
  **Krawtchouk** rows `g_J`. This is the scheme's **distance distribution and its dual transform**.
- **THM-534** — the moment-LP dual `g(t) ≥ 1[t=0]` certifying `measS7(E) ≤ L_y(E) = Σ_r y_r S_r` for
  every `E`. This is the **LP**.

The fact that makes it one object (verified, exact): **THM-534's dual `g(t)` expands in the binary
Krawtchouk basis with all-nonnegative coefficients** — exactly Delsarte dual feasibility. At `k=8`,
`g` is purely even, `c = [1/16, 0, 1/40, 0, 3/80, 0, 0]` on `K₀, K₂, K₄` (this is *why* HYP-2724's
even band was clean at `k=8`); at `k=9,10` and `k=11,12,13` both parities appear, all `≥ 0`. So at
every binding `k`, the LRC cover bound is a Delsarte LP, and the moment-LP we already had is its
certificate. The cover problem has a classical home.

## What it buys

The Delsarte frame turns the proof into a standard shape. `measS7(E) ≤ L_y(E) = Σ_j c_j E[K_j(N)]`
with `c_j ≥ 0` is the per-`E` Delsarte upper bound — already proved (THM-534). The extremality —
"consec maximizes `L_y`," the one open piece — becomes **"consec saturates the Delsarte LP,"** and
complementary slackness names where: the extremizer's empty-count `N` must be supported on the `t`
where `g(t)` meets its target (`k=8`: `t ∈ {0,1,2,4,5}`, the zeros of `g` plus the cover point). The
LP-tight configurations of a Delsarte scheme are the classical extremal codes/designs — and the
user's MDS/arc lead points at the right one: the consec block is the **anti-MDS** (minimum-distance
code, weight enumerator concentrated at low weight), the polar opposite of the Reed-Solomon /
normal-rational-curve arc that maximizes distance.

## Honest boundaries

The pretty parts that did *not* survive contact: the **56-coincidence is a number match**, not a
bijection — `C(8,5) = 56 = A000568(6)`, but the support-3 relation-hypergraph iso-shape count is
unbounded (47, 66, 86, … as the window grows; kps and a thread both confirmed), so there is no
challenger↔tournament bijection. And the "even Krawtchouk only" structure is `k=8`-special; `K₄` is
not even convex, so the robust fact is Krawtchouk-*nonnegativity* of the dual, not convexity of each
moment. The cited external paper (arXiv:2501.19125) is LDPC minimum-distance bounds — coding theory,
but not a direct lever.

What remains is exactly what kps's adjudication left standing: bound the conditionally-convergent
`Σ K(n)`. But it is now the same sentence as a Delsarte LP being tight at the anti-MDS code — a
problem with a hundred years of machinery behind it (MacWilliams, Delsarte, the linear programming
bound, the MDS conjecture) rather than a bespoke inequality. The arithmetic crux was a coding-theory
extremal problem the whole time; the apex prime 7 is the alphabet of the scheme, the offsets are the
code, and consec is the densest code in it.

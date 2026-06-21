# Two open extremalities, one theta ceiling

**Session:** mac-mini-2026-06-21-S14. Testing the Coregliano–Jeronimo–Jones complete LP hierarchy
(arXiv:2211.01248) against the LRC(14) consec-extremality gap.

## The lead, and its honest verdict

The hope was clean: the LRC cover bound *is* a Delsarte LP (HYP-2726), the open part is that consec
*saturates* it, and CJJ give a hierarchy that is *complete* — the exact optimum at a finite level.
So a low level might prove the extremality and retire the finite atlas.

It does not, and the reason is the one CJJ themselves flag (Prop 1.2): **the hierarchy beats Delsarte
only through linearity** — the optimizer must be closed under linear combination. The LRC optimizer
is the *miss-event distribution* of an offset set, which is not a linear code; the missed-sector
"words" don't add. Concretely, in the Z/7-symmetric scheme the degree-2 atom moments are *determined*
by the degree-1 moments (`m_2 = S_2 / C(6,2)`), so the SoS/Lasserre PSD lift adds no constraint and
`θ′ = Delsarte = L_y`. The hierarchy collapses to level one. The naive lift is a dead end.

## What survived: the unification

But asking *which graph* the Lovász θ′ lives on turned up the real prize. Schrijver's identity says
Delsarte's LP equals `θ′(H)` for the conflict graph `H` whose independent sets are the codes. For LRC
that graph is not the Z/7 sector Cayley graph — it is the **relation-scheme graph `H_E`**, and
`θ′(H_E) = L_y(E)` exactly (verified). The same `θ′` ceiling is genuinely aggregate — unclosable by
the relaxation — and it is aggregate for the *same structural reason* on **both** sides of the
apex-prime gas:

- **Tournaments:** the Paley/regular tournament maximizes `H` (THM-126/134); the SoS/θ relaxation
  does not certify it — the project hit this ceiling years of sessions ago.
- **Lonely runners:** the consecutive block maximizes the Z/7 cover; the SoS/θ relaxation does not
  certify it either — the same ceiling, reached this week.

These were tracked as two separate open problems. Under the Lovász-θ lens they are **one**: the
extremal object (Paley / consec) saturates a θ′ relaxation that the LP/SoS hierarchy cannot tighten,
because in neither case is the optimizer a linear code. The "apex-prime gas" slogan — that tournaments
and runners are one object on `Z/p` — now has a precise optimization-theoretic form: *they share a
single Lovász-θ ceiling, and the open content on both sides is the same aggregate gap between θ′ and
the truth.* Whatever eventually proves one extremality should prove the other.

## A smaller, honest gain

One thing did tighten — but it is moment *matching*, not the SoS lift. Bounding `P(N=0)` by the LP
that matches a shape's first few sector-marginals strictly improves with order, and at level 3 the
worst case over shapes is consec, with a uniform bound just under the cap (`0.4929 < 0.4943` at k=9,
`0.348 < 0.381` at k=8). So gap #4 *reframes* as "consec maximizes the level-3 moment-LP bound" — a
lower-dimensional, moment-only extremality. It is cleaner, but it is finite-verified (zero beaters
over the tested shapes), not proved, and at k=9 the margin is a razor. A reframing, not a closure.

## The lesson

The honest arc of two sessions: I overstated "L7 closed," the audit caught it, I corrected it and
fixed the real rate gap; then I hoped the LP hierarchy would close the extremality, and it collapsed.
What is left standing is not a proof but a sharper map — the LRC and tournament extremalities are the
same θ′ ceiling — and a cleaner reframing of the one genuinely-aggregate gap. The relaxation will not
give consec for free, on either side of the gas; the proof, when it comes, will have to be the kind of
argument that sees why the *most linear* arithmetic object (the AP, the Paley/QR set) sits exactly at
the top — and that argument is owed to both problems at once.

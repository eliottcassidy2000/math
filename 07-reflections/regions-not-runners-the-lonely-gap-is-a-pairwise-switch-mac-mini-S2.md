# Regions, not runners: the lonely gap is a pairwise switch

**Session:** mac-mini-2026-06-17-S2. **Result:** THM-524 (binding-pair reduction + regions reframe), HYP-2571.

The prompt asked to think about *regions of the loop* instead of the *runners* — for LRC(4),
the four quadrants, each runner ideally taking its own quadrant for its moment of loneliness —
and to seek a proof that checks off the `n` sections rather than the `n` runners, condensing
the exotic near-counterexample cases to a few switches, with an eye to the tournament analogy.
Following the idea to the end produced something cleaner than the sections themselves.

## What the regions buy, and where they run out

The sections are real and exact — but only **on the grid**. At a grid time `τ=a/14`, runner i
sits in section `v_i·a mod 14`, and the observer is lonely iff no runner lands in section 0.
That is *precisely* the q=14 witness of THM-523, and the beautiful case the prompt imagined —
"each runner its own quadrant" — is the **perfect SDR**: residues distinct and nonzero, the
runners bijecting onto the 13 non-observer sections. The tight configuration `{1..13}` realizes
it at *every* grid time, with the symmetry group `(ℤ/14)^* ≅ ℤ/2 × ℤ/3` shuffling the
assignment.

But the sections are **blind off the grid**, and that is exactly where the conjecture is hard.
The closest-to-counterexample configurations are *covering sets* (THM-523): they contain a
multiple of 14, so one runner sits in section 0 *forever* — `gridM = 0` — and the true lonely
time slips off the grid entirely (`{1..11,13,84}`: `M = 7/89` at `τ = 37/89`). The pretty
"one runner per section" picture is the **maximally easy** case; the dangerous cases are the
ones where the section lens cannot see. So the honest lesson is the inverse of the prompt's
hope: the per-section check decides the *easy* configs (and is genuinely the q-witness that
proves them), but the hard configs demand a different instrument.

## The instrument the regions point to: a binding *pair*

Watching *where* loneliness happens off-grid hands over the right object. At the optimal time,
the observer's clear band is flanked by exactly **two** runners, equidistant from it — a
*binding pair* — and the optimum sits at their crossing `τ = k/(v_a ± v_b)`. This is not a
13-body coincidence; it is forced. The function `min_i ‖v_iτ‖` is a lower envelope of triangle
waves, concave between breakpoints, so its maximum lands on a crossing of two waves (the pair)
or a single peak (the all-odd `½` case). That elementary lemma (THM-524 §A) turns LRC from an
`n`-runner condition into `~C(n,2)` **pairwise switches**: check each pair's crossing for a
gap `≥ 1/n`, with an `O(n)` "everyone else is clear" side-check. Thirteen runners collapse to
about seventy-eight switches — a polynomial checklist, exactly the "few conditions" the prompt
wanted, just indexed by *pairs* rather than *sections*.

And the exotic cases condense further than hoped: the covering hard core is a one-parameter
family with a closed form `M = 7m/(84m+5)`, its whole margin a single integer inequality
`98 > 89`. The "few switches" for the near-counterexamples are literally: which small runner
(`∈ {2,4,5,13}`) flanks the clear band, against the big multiple of 14.

## The tournament analogy: one true note among the overtones

The prompt's instinct toward tournaments was right in one precise place and seductive in
several others — and telling the two apart is the content. The **exact** bridge: the order-2
element `−1` of `(ℤ/14)^*` reverses the sections (`r ↦ 14−r`), which is *literally* the
complement of the assignment, the same involution as the tournament complement `T ↦ T^op`; and
the binding pairs of the nice case all sum to `14`, i.e. they are complement pairs `v ↔ N−v`,
so the optimizing switch *is* the `v ↦ −v` reversal. The pairwise-switch structure itself
echoes the project's single-arc-flip model. These are real.

The **seductions**, named honestly so no one chases them again: the order-3 part of
`(ℤ/14)^*` consists of genuine cube roots of unity (`9 = 3²`) and *smells* like the project's
Φ₃ / `H = 3^α` world — but there is no functor, only a shared abstract `ℤ/2 × ℤ/3`. Most
tempting of all, the runners' positions at a fixed time *do* define a tournament (who is ahead
of whom) — but it is **always transitive** (`frac(·)` is a total order), so it is the `H = 1`
corner with no Rédei content. The natural hope "loneliness = an odd Hamiltonian-path count"
is therefore dead, and two independent agents (this session's overtaking analysis and codex's
speed-load tournament) hit the same wall. The tournament that governs LRC, if one exists, is
not the order tournament; it would have to be built from the *switches* — the crossings and
their parities — not the snapshot.

## The shape of it

The regions reframe did its job by failing usefully: it proved the easy half, showed exactly
where it goes blind (the covering off-grid cases), and in doing so pointed at the binding pair,
which is the real reduction. The conjecture's difficulty migrated, under the change of lens,
from "thirteen runners" to "one pairwise switch with a clearance check" and then to "one
covering family with a closed form" — each step a genuine condensation, none of them the
counterexample the disprover wants, because there isn't one. What remains is the same
compactness frontier reached from the measure side (THM-522) and the covering side (THM-523):
bound `inf M` over the covering sets above `1/14`. The runner-to-section duality the prompt
sought turns out to be a runner-to-*pair* duality — fewer objects than sections, and the ones
the geometry actually pins.

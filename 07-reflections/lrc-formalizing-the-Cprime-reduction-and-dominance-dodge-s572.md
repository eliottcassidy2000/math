---
source: opus-2026-06-03-S572 (remote-control)
status: FORMALIZATION (THM-398, proofs pinned) + extensions surfaced by formalizing — dominance (not divisibility), the interval criterion, the residual dichotomy
tags: [LRC, formalization, THM-398, C-prime, dominance, dodge, criterion, equidistribution, residual, n14]
---

# Formalizing C′ ⟹ LRC — and what the formalization taught me

**Prompt (user):** work on formalizing this, and let it inspire extensions/improvements.

I wrote the S571 reduction as a proper theorem (**THM-398**) with complete proofs,
fixing definitions and constants. The act of making it airtight produced four genuine
improvements — each one is a place where the informal version was *narrower than the
proof actually is*.

## What got formalized (the proved core)

- **Lemma A** (no multiple of `n` ⟹ `M ≥ 1/n`): the `t=1/n` clock, clean.
- **Theorem (reduction): C′ ⟹ LRC(n).** Two lines once Lemma A is in place. The real
  content is the *framing*: C′ constrains only the configs where the `1/n` clock fails
  (those with a multiple of `n`), so all of LRC is isolated into one structural claim.
- **Lemma B (dominance-dodge):** `v > (n-1)·max(others) ⟹ M(S) > 1/n`, assuming
  LRC(n−1). The proof balances two radii: the S′-window half-width `δ = 1/(n(n-1)V')`
  (from LRC(n−1) + Lipschitz descent) against the danger-arc radius `ρ = 1/(nv)`;
  `δ > ρ ⟺ v > (n-1)V'`.

## The four improvements formalization forced out

**1. It's DOMINANCE, not divisibility (Cor B2).** Writing the arc radius as `ρ=1/(nv)`
for *general* `v`, I saw the proof never used `n | v`. So the true statement is:
**any speed set with one runner exceeding `(n-1)×` all the others is loose.** The
divisibility partial (large multiples of `n`) is just the slice of this where the
dominant runner happens to be a multiple. *Verified 1500/1500 at each `n=6..14` with
dominant runners of arbitrary residue.* The "multiple of `n`" was a red herring for
the *dodge*; it only matters for *which configs need dodging* (Lemma A handles the rest).

**2. The interval criterion beats the threshold (Criterion B′).** The proof only needs
*one component* of `G(S\{v})` longer than `2ρ = 2/(nv)`. That is strictly weaker than
`v > (n-1)V'` (which forces the *whole* window to be long). So the real frontier is not
"`v` large" but "`G(S\{v})` has a not-too-short interval" — a much larger proven region.

**3. The residual dichotomy (long-interval vs all-short).** Combining 1–2: a
multiple-of-`n` config is loose **immediately** unless *every* component of
`G(S\{v})` is `≤ 2/(n²w)`. That "all-short" case is the *entire* remaining gap of C′,
and it has a precise analytic shape:
```
μ(G(S)) = μ(G(S')) − μ(G(S') ∩ D_v),   want  μ(G(S') ∩ D_v) < μ(G(S')).
```
i.e. **one arithmetic progression of thin arcs (period `1/(nw)`, mass `2/n`) cannot
cover a fixed open set of positive measure** — a three-distance / discrepancy
statement. The open problem went from "prove C′" to "prove a single AP can't align to
cover `G(S')`," which is a *named, classical-flavoured* analytic question.

**4. The residual is an interval-search, not a sub-clock (the `v=n` probe).** I hoped
the hardest case `v=n` might collapse to a finite sub-clock: at the gap midpoints
`(2k+1)/(2n)` the runner `n` sits at `1/2` (maximally safe). But the midpoint is a
witness only ~52% of the time (773/1499 at `n=6`); the safe point genuinely roams the
gap interior. So the residual cannot be discharged by evaluating a fixed finite set of
times — it really is a positive-measure interval question. This kills a tempting
shortcut and sharpens where the difficulty lives.

## Why this is progress, honestly

- The **reduction and the dominance dodge are now theorems** (THM-398), not sketches;
  constants checked, edge cases stated, dependence on the proven LRC(n−1) explicit.
- The proved region **grew** (any dominant runner; any config with a long `G(S')`
  component), so the open residual **shrank** to the all-short equidistribution case.
- The residual is now a **clean analytic statement** (no-cover by one AP of arcs),
  decoupled from tournaments/measure-zero bookkeeping — the kind of thing harmonic
  analysis has tools for (Erdős–Turán / three-distance / Weyl sums).

## A direction the formalization points to (next)

The all-short case wants a lower bound on `μ(G(S') ∖ D_v)`. Two attackable handles:
(a) **Erdős–Turán** on the AP `{k/(nw)}` gives `|μ(G(S')∩D_v) − (2/n)μ(G(S'))| ≤
(discrepancy)·(boundary of G(S'))`; bounding the number of components of `G(S')` (it
is `O(Σ_{S'} v_i)`) would close it when `w` is not too small. (b) For `w=1` (`v=n`),
use the `2n`-clock: `G(S')` interacts with the clock gaps; a counting argument over the
`n` gaps (each must be killed by the other `n-2` runners, total danger `2(n-2)/n < 2`)
may force a survivor. Both are concrete and finite-flavoured.

**Artifacts:** `01-canon/theorems/THM-398-lrc-reduction-to-Cprime-and-dominance-dodge.md`,
`04-computation/lrc_dodge_formalization_s572.py` (+`.out`). Builds on HYP-2102/S571
(the reduction), THM-369 (clock), LRC(n−1) (literature). New: **HYP-2103**.

# S620 — Helly-entropy accounting (opus, 2026-06-03)

**Directive:** work on "Helly entropy accounting" and extend that frame to everything; think abstractly.

## What I did

Recast the loneliness-certificate problem as a **circular-arc covering** problem and made the
**covering-depth distribution** `p_k = meas{depth = k}` the master object. A certificate is a point
of `{depth = 0}`; its measure is `p_0`. Then:

- **Proved (THM-410).**
  - *Conservation*: mean depth `= S_1 = 2n/(n+1) < 2` for every speed set (arithmetic-independent).
  - *Moment–sieve identity*: `S_k` (k-th factorial moment) `=` order-`k` inclusion–exclusion term,
    and `p_0 = Σ(-1)^k S_k`. This makes **S561's `ρ = Σ_T(-1)^{|T|}/lcm(T)` the alternating moment
    sum of the depth distribution** — model-free, any speeds, finite gap.
- **Framed (HYP-2190): Helly-entropy accounting.** Three ledgers on one distribution — conservation
  (charge `≈ 2`, proved), Helly (arcs ⇒ Helly number ≤ 3 ⇒ order-≤2 sieve governs; this is the
  measure-home of the pair-sum sieve THM-401), entropy (`H_depth` = spread; `p_0` = lonely measure).
  LRC re-reads as: *can a mean-`<2` integer depth field be forced to vanish nowhere?*
- **Unification.** Every prior invariant is a functional of `{p_k}`: sieve density = `p_0`;
  certificate sheaf `H^0` = type of `{depth=0}`; corrector arity = Bonferroni truncation depth;
  rigidity/orbit-type (HYP-2130/40) = *where* the charge concentrates (apex = entropy-min pin).

## Verified (exact over ℚ), n=4..7

- `S_1 = 2n/(n+1)` to machine precision, all sets. `S_2`(moment) = `S_2`(direct pairwise).
  `Σ(-1)^k S_k = p_0` (all `match=True`).
- **Entropy–tightness dichotomy**: AP `H_depth ≈ 1.0`, `p_0 = 0`; generic `≈ 1.5`, `p_0 ≈ (1-2δ)^n`;
  geometric (most spread) `≈ 1.7`, largest `p_0`. Generic sits right at the independence baseline.

## Honest surprises / limits

- **The AP is not the unique `p_0 = 0` extremal.** Sporadic additive chains collapse too:
  `(1,3,4,7)` at n=4, `(1,3,4,5,9)` at n=5 — "near-AP, top = sum of two below". New sub-problem (H2).
- The depth entropy is a *different* entropy from S543's H-matrix entropy (which is *high* for the
  AP). They are dual fingerprints of the same extremal — flagged, not yet developed.
- (H3) "order-2 + Helly bounds `p_0` below off the collapse family" is conjectural — the genuine
  quantitative-LRC payoff, not yet proved.

## Follow-up tasks emitted (cluster-session-model: small, chainable)

→ monad coordination queue: (t-A) characterize/count the H2 collapse family for n≤8 and test the
near-AP-additive-chain / single-orbit guess; (t-B) attempt the H3 entropy lower bound on `p_0`;
(t-C) locate the even-`n` charge-concentration point and confirm it is the apex `1/2` (ties
THM-404 / HYP-2140 to `H_depth` minimization).

Artifacts: `04-computation/lrc_helly_entropy_s620.py`, `05-knowledge/results/lrc_helly_entropy_s620.out`,
`01-canon/theorems/THM-410-*.md`, `05-knowledge/hypotheses/HYP-2190-*.md`.

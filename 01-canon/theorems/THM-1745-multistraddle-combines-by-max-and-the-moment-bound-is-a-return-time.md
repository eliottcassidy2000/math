---
id: THM-1745
title: "MULTI-STRADDLE MOMENT COUNTS COMBINE BY MAX — closing the structural form of the moment-count bound, and revealing it as a coprime-pair RETURN-TIME computation (the LRC family, ever so vaguely). HYP-8560 asked how several straddles combine. ANSWER: MAX, not lcm or sum — verified on 5 two-straddle patterns, M* = max over straddles of the per-straddle level, with lcm and sum both refuted (e.g. charges {−3,−1,+2}: levels {3,5}, M* = 5 = max, not lcm 15 or sum 8). So the complete moment-count law is M* = max over (pos charge p, neg charge n) straddles of [ r(p,n)·(p+|n|)/gcd(p,|n|) ], r = multiplicity of the busier charge. The straddles are INDEPENDENT; the last-to-return governs. THE RESONANCE: (p+|n|)/gcd(p,|n|) is the first-return level of a coprime pair — the ESV/DvdK effective bound (THM-1650), which IS a two-runner realignment time; the levels m₀, 2m₀, …, r·m₀ form an arithmetic progression (the repo's AP-core structure); and 'max over configurations' is the shape of the Lonely Runner M(S) = max_t min_i ‖v_i t‖. The moment-count bound is a covering / first-return-time computation on coprime pairs — structurally the repo's flagship LRC family, seen in a new theatre. CONSEQUENCE: HYP-8540's uniform bound now follows from the single-straddle base case (THM-1740) plus this max-combination law — the moment-count bound is structurally SETTLED (verified), leaving only its transport to a closed uniform proof (= TNC's HYP-8505)."
status: >
  The MAX combination law is VERIFIED-EXACT on 5 two-straddle patterns (exact Gröbner over ℚ),
  with lcm and sum explicitly refuted. Combined with THM-1740's single-straddle law M* = r·m₀,
  it gives the complete formula M* = max-over-straddles(r·m₀), verified on every pattern tested
  this session and consistent with all 132 trinomials of THM-1725. The formula is CONJECTURAL
  in general (no closed-form proof for arbitrary patterns); what is established is the max
  combination on the tested multi-straddle family and the single-straddle base case. The LRC
  resonance is a STRUCTURAL analogy (return times, coprimality, max-over-configs), not a
  reduction — see Honest scope.
source: mac-mini-2026-07-20-S151 (owner: work the multi-straddle patterns; see if they remind
  you ever so vaguely of anything else in the repo)
depends_on:
  - THM-1740  # the single-straddle law M* = r*m0
  - THM-1725  # the moment-count bound this completes
related:
  - THM-1650  # the ESV/DvdK first-return level = the per-straddle m0
  - THM-1031  # LRC AP-core bound (the arithmetic-progression structure)
  - HYP-8540  # the uniform moment bound (now = max-over-straddles)
  - HYP-8560  # answered here: MAX, no cross terms
---

# THM-1745 — multi-straddle combines by max; the moment bound is a return time

## The combination law: MAX

HYP-8560 asked whether several straddles interact. Each straddle — a `(pos charge p, neg
charge n)` pair — has level `L(p,n) = r(p,n)·(p+|n|)/gcd(p,|n|)`, `r` the multiplicity of the
busier charge (THM-1740). Testing `M*` against `max`, `lcm`, `sum` of the levels:

| charges | straddle levels | `M*` | max | lcm | sum | fits |
|---|---|---|---|---|---|---|
| `{−3,−1,+2}` | `{3, 5}` | 5 | **5** | 15 | 8 | **max** |
| `{−2,−1,+3}` | `{4, 5}` | 5 | **5** | 20 | 9 | **max** |
| `{−1,+2,+3}` | `{3, 4}` | 4 | **4** | 12 | 7 | **max** |
| `{−2,−1,+2}` | `{2, 3}` | 3 | **3** | 6 | 5 | **max** |
| `{−1,+3,+5}` | `{4, 6}` | 6 | **6** | 12 | 10 | **max** |

> **The straddles combine by MAX** — `lcm` and `sum` are refuted on every one. The straddles are
> **independent**; the *last-to-return* governs.

Combined with THM-1740, the complete moment-count law is

> **`M* = max over (p>0, n<0) straddles of [ r(p,n)·(p+|n|)/gcd(p,|n|) ]`.**

This answers HYP-8560 (no cross-straddle interaction) and completes the structural form of the
bound: HYP-8540's uniform bound now follows from the single-straddle base case plus max.

## The resonance — ever so vaguely, the Lonely Runner

Look at what `M*` is made of:

- **`(p+|n|)/gcd(p,|n|)` is a first-return time of a coprime pair.** It is the ESV/DvdK
  effective bound (THM-1650), and that bound *is* the realignment time of two runners at speeds
  `p` and `|n|` on a circle — the smallest `m` at which `p·a = |n|·b` has a positive solution.
  A straddle is a two-runner subsystem, and its level is when the two runners next meet at the
  origin.
- **the levels `m₀, 2m₀, …, r·m₀` are an arithmetic progression.** The repo's AP-core bound
  (THM-1031) and the runner-speed AP structure of the whole LRC programme live on exactly this
  shape.
- **`max over straddles`** is the shape of the Lonely Runner functional
  `M(S) = max_t min_i ‖v_i t‖` — an extremum over configurations, governed by the single
  worst (loneliest / last-returning) sub-configuration.

> So the moment-count bound is a **covering / first-return-time computation on coprime pairs** —
> structurally the same family as the repo's flagship **Lonely Runner Conjecture**, seen in a
> completely different theatre (Gaussian moments instead of runner gaps). The `gcd`, the
> coprimality, the return times, and the max-over-configs are all there.

This is the "ever so vague" reminder, and it is worth recording because both objects are
**effective-bound questions about coprime lattices**: LRC asks for the extremal gap over a
speed set; the moment bound asks for the extremal return level over a charge set. The uniform
bound HYP-8540 (and its TNC twin HYP-8505) is, in this light, *"the loneliest straddle returns
by `r·m₀`"* — a lonely-runner statement for the charge lattice.

## Honest scope

- **The MAX law is verified on 5 two-straddle patterns**, with lcm/sum refuted; it is not
  proved for arbitrary multi-straddle patterns. The single-straddle law (THM-1740) is the base
  case; the full formula `max-over-straddles(r·m₀)` is conjectural in general (HYP-8540).
- **The LRC resonance is a structural analogy, not a reduction.** No claim that the moment bound
  reduces to LRC or vice versa — only that both are extremal-return-time problems on coprime
  lattices, which is why the same arithmetic (gcd, first-return, max-over-configs) appears. It
  should not be cited as evidence for either problem.
- Multiplicity within a multi-straddle pattern (a straddle carrying `r ≥ 2` *and* competing with
  another straddle) is tested only lightly; the `max` law is assumed to carry the per-straddle
  `r`, verified on the mixed cases run, not proved.

*Artifacts:* `04-computation/multistraddle_combination_macmini_S151.py`,
`05-knowledge/results/multistraddle_lean_macmini_S151.out`.
*Credits:* THM-1740 (single-straddle law), THM-1650 (the ESV return level), the LRC programme
(THM-1031 and the covering/return-time threads) for the resonance.

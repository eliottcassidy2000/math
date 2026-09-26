---
id: HYP-9161
title: "Landing multiplicity is polylogarithmic: in THM-4476's recursion, the dippers charged to one landing point number O((log X)^beta) on average with beta < 1 (conjecturally O(theta log X)); this would lower THM-4499's exponent from lambda*/h* - 3/2 = -0.986 to beta lambda*/h* - 3/2, i.e. to the ballot floor -3/2"
status: OPEN (formulated 2026-09-26, opus collatz-landing-20260926). FINITE-EXACT: mean multiplicity 0.4-0.8 times theta L on actual orbit segments (L = 20..80), maximum up to 0.35 L; multiplicity 0.54 L is realised by residue classes (hovers in a 1-bit band; audited), and the audit's one-bit band lemma bounds it by 2k/3 + 4/3 for every orbit, so the hypothesis is about a single orbit's oscillation and cannot be proved by counting residue classes alone.
source: collatz-landing-20260926 session (opus); the reassessment note collatz_landing_20260926_multiplicity_reassessment.md, sections 0-4.
depends_on:
  - 01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md
  - 01-canon/theorems/THM-4499-thin-divergence-is-little-o-of-x-to-the-h-star.md
related:
  - 05-knowledge/results/collatz_landing_20260926_multiplicity_reassessment.md
  - 04-computation/experiments/collatz_landing_20260926_probe.py
---

# HYP-9161 -- the landing multiplicity is polylogarithmic

**Setting.** A positive one-signed `T_b`-orbit segment `(y_i)`, `X`,
`L = log_2 X`, `k = floor(L)`, `theta` with `theta L >= 1`. An index `i` with
`y_i <= X` is a dipper if `y_(i+s) < y_i X^(-theta)` for some `1 <= s <= k`;
its landing point is `i + s` for the least such `s`. The multiplicity of a
landing point is the number of dippers charged to it.

**Hypothesis.** There are `beta < 1` and `C` such that for every orbit,
`X` and `theta`, the total number of dippers is at most
`C L^beta` times the number of landing points. Conjecturally `C theta L`
suffices (`beta = 0` at the bootstrap's `theta_X = O(log L/L)`).

**What is known.** The trivial bound is `k` (THM-4476, section 1.5). Scaling
probe (`collatz_landing_20260926_probe2.py`, `L = 20..80` at `theta_X`): mean
multiplicity `0.4`-`0.8` times `theta L`; maximum up to `0.35 L` on segments
with a built-in climb, so only the average can be polylogarithmic. Residue
classes realise multiplicity `k` (hover-then-drop, climb-then-drop; Terras),
so no orbit-blind count proves the hypothesis. On actual segments the mean
is `2`-`3` (probe). The random-walk heuristic gives `O(theta L)`.

**What it gives.** With THM-4499's proof unchanged except `c_1 h* = beta + eta`,
`N(X) <= K X^(h*) (log_2 X)^a` for every `a > beta lambda*/h* - 3/2`.

**Hostile (corrected after the audit).** Every dipper of a landing point
lies in one 1-bit band above `2^(theta L) y_j`, and at most two of three
consecutive orbit values share such a band, so the multiplicity is at most
`2k/3 + 4/3` for every orbit; hovers in a narrow band realise `0.54 L`
(`L = 24`, exhaustive search), climbs realise only `2`. A hypothetical
divergent orbit whose window points repeatedly crowd into 1-bit bands
before each descent would have multiplicity a constant fraction of `L` at
every landing point; the strip entropy shows such stretches are rare among
integers but does not exclude them on one orbit.

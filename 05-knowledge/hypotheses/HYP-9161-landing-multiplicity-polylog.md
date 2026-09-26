---
id: HYP-9161
title: "Landing multiplicity is polylogarithmic: in THM-4476's recursion, the dippers charged to one landing point number O((log X)^beta) on average with beta < 1 (conjecturally O(theta log X)); this would lower THM-4499's exponent from lambda*/h* - 3/2 = -0.986 to beta lambda*/h* - 3/2, i.e. to the ballot floor -3/2"
status: REFUTED as stated for arbitrary finite positive distinct orbit segments (2026-09-26, crossroads-poset). An endpoint-corrected inequality or an explicitly whole-injective-orbit version remains OPEN. Exact hover/drop segments have average multiplicity Omega(log X) even when theta log X is O(loglog X).
source: collatz-landing-20260926 session (opus); the reassessment note collatz_landing_20260926_multiplicity_reassessment.md, sections 0-4.
depends_on:
  - 01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md
  - 01-canon/theorems/THM-4499-thin-divergence-is-little-o-of-x-to-the-h-star.md
related:
  - 05-knowledge/results/collatz_landing_20260926_multiplicity_reassessment.md
  - 04-computation/experiments/collatz_landing_20260926_probe.py
---

# HYP-9161 -- the landing multiplicity is polylogarithmic

**UPDATE 2026-09-26 (collatz-procgen session, THM-4506).** The worst case is now exact: all dippers of a landing point lie in one dyadic shell, consecutive ones are separated by an odd letter, so m(j) <= ceil((k-D)/log_2 3) (about 0.631(k-D)), attained by residue classes and by actual orbits (13255 at L = 20). This sharpens the one-bit band bound 2 ceil(k/3) below. So "residue classes realise multiplicity k" in the retained original text should read ceil((k-D)/log_2 3). THM-4506 also proves that the recursion is saturated at a*(mu), that depth averaging is neutral, and that orbit-blind residue-class splits and parity-word data cannot give mu < 1 (Proposition H). The endpoint-corrected, orbit-coupled version of this hypothesis remains OPEN; in THM-4506's language it is a local-time bound for visits to the shell D bits above later landing points. See [THM-4506](../../01-canon/theorems/THM-4506-landing-multiplicity-exact-worst-case-and-recursion-saturation.md).

**CORRECTION, 2026-09-26: the literal finite-segment hypothesis below is
REFUTED.** The [fixed-integer audit, section 5](../results/crossroads_poset_20260926_integer.md)
constructs one actual positive integer and K distinct positive orbit steps
with m+3 dippers landing at at most five positions, where
`D=ceil(log_2 m)+6`, `K=m+D`, `L=4K+3`, and `theta L=D-3`.
Thus the average is Omega(L), not O(L^beta) for beta<1 or O(theta L).
The retained m=96 control gives 99 dippers and four landing positions.

The entire hostile segment is shorter than the horizon k=floor(L), so an
O(k) endpoint error absorbs it. The meaningful surviving question is an
inequality `Dippers <= C polylog(L) Landings + C' k`, or an appropriately
specified count over an entire injective/nonperiodic orbit. Neither is
proved or refuted here. This correction does not weaken THM-4476 or THM-4499.
The original formulation and evidence follow for correction lineage.

**Setting.** A positive one-signed `T_b`-orbit segment `(y_i)`, `X`,
`L = log_2 X`, `k = floor(L)`, `theta` with `theta L >= 1`. An index `i` with
`y_i <= X` is a dipper if `y_(i+s) < y_i X^(-theta)` for some `1 <= s <= k`;
its landing point is `i + s` for the least such `s`. The multiplicity of a
landing point is the number of dippers charged to it.

**Original hypothesis (REFUTED in the stated segment setting).** There are `beta < 1` and `C` such that for every orbit,
`X` and `theta`, the total number of dippers is at most
`C L^beta` times the number of landing points. Conjecturally `C theta L`
suffices (subpolynomial in L at `theta_X = O(log L/L)`, hence every
fixed beta>0 asymptotically, rather than literally a beta=0 bound).

**What is known.** The trivial bound is `k` (THM-4476, section 1.5). Scaling
probe (`collatz_landing_20260926_probe2.py`, `L = 20..80` at `theta_X`): mean
multiplicity `0.4`-`0.8` times `theta L`; maximum up to `0.35 L` on segments
with a built-in climb, so only the average can be polylogarithmic. Residue
classes realise multiplicity `k` (hover-then-drop, climb-then-drop; Terras),
so no orbit-blind count proves the hypothesis. On the particular segments
in that probe the mean is `2`-`3`; the correction above refutes the universal
extrapolation. The random-walk heuristic gives `O(theta L)`.

**What it gives.** With THM-4499's proof unchanged except `c_1 h* = beta + eta`,
`N(X) <= K X^(h*) (log_2 X)^a` for every `a > beta lambda*/h* - 3/2`.

**Hostile (corrected after the audit).** For b>0, or a positive window
entirely above 5|b|, every dipper of a landing point
lies in one 1-bit band above `2^(theta L) y_j`, and at most two of three
consecutive orbit values share such a band, so the multiplicity is at most
`2k/3 + 4/3`. Negative-b low-core windows need separate treatment as
explained in the reassessment note. Hovers in a narrow band realise `0.54 L`
(`L = 24`, exhaustive search), climbs realise only `2`. A hypothetical
divergent orbit whose window points repeatedly crowd into 1-bit bands
before each descent would have multiplicity a constant fraction of `L` at
every landing point; the strip entropy shows such stretches are rare among
integers but does not exclude them on one orbit.

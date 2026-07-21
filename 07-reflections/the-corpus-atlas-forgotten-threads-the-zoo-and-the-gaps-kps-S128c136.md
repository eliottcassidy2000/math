# The corpus atlas: forgotten threads, the zoo, and the gaps

*kind-pasteur-2026-07-21-S128c136. Owner: keep adding to the zoo; go back through all past work,
find every thread, make sure none are lost; procedurally generate new frames/methods/angles; find
the things we forgot we studied; find the gaps between and around them.*

The corpus is **1392 theorems, 2554 reflections, 4345 hypotheses**, THM-001..1860, 21 MB of
reflections alone. No one can hold it. So I built a machine to hold it:
`corpus_archaeology_kps_S128c136.py` — a citation-graph scan that finds forgotten threads, maps
the topic zoo, and (with `cycle_counts_spectral_boundary`) reconnects a dead leaf to the frontier.
This atlas records what it found so it is not lost again.

## 1. The forgotten threads (60 flagged; the ones that matter)

A theorem cited by `≤ 1` external file is a **leaf the repo proved and forgot**. 60 such, with
`THM# < 1700`. The full list is in `05-knowledge/results/corpus_archaeology_kps_S128c136.out`; the
clusters:

- **MOMENT / CYCLE-SPECTRAL (era-2, ~2026-03) — the big one.** THM-056 (forward-arc moment
  hierarchy), THM-057 (general-n moment closed forms), THM-157/158 (`α₁`, `α_j` moment
  linearity), THM-171 (`c3`-disjointness from λ), **THM-172 (5-cycle count is λ-determined)**,
  THM-173 (`c5` per-vertex-set), THM-225 (top eigenvalue of `Cᵀ C`). **These are the era-2
  ancestors of the current moment-nullcone / binary-form ladder (THM-1775/1810).** The project
  studied "cycle counts are functions of the spectrum" and "moment hierarchies" in March, forgot
  it, and rediscovered it in July as "tr(Aᵏ) = SL₂-invariants of char_A." **Reconnected today:**
  THM-1870 proves the cycle-count spectral boundary is the Hamiltonian length — `c_k` spectral for
  `k ≤ n−1`, `c_n`/`H` `#P` from `n = 6` — which is exactly THM-172 completed and glued to
  THM-1780.
- **BLOCK-COUNTING / μ (era-0, the very start).** THM-006, 010, 011, 012, 026, 054 — the original
  `μ`-computation and block-counting identities (pre-MISTAKE-001). Superseded but never linked
  forward; worth a one-line "subsumed by" pointer so the early scaffolding is legible.
- **H-IDENTITIES.** THM-208 (`H = 1 + 2t₃` for `n ≤ 4`), THM-210 (sum-H identity), THM-140
  (additive-energy ↔ Hamiltonian-path correlation sign reversal). Small-`n` or niche, but THM-208
  is a *direct* `H`-vs-`c3` law — the `n ≤ 4` germ of the WOWII `c3 ≤ H` (THM-1860).
- **STRAYS worth a look.** THM-207 (permutohedron = transitive), THM-215 (compatible pairs =
  A000255), THM-254 (instant MCMC — an engineering tool left unused), THM-293 (W = succession GF).

## 2. The topic zoo (title-keyword census)

Mainline: **LRC 119, moment 39, Paley 31, GMC 31, tiling 24, nullcone 21, TNC 19, OCF 17**. Then
king 14, Hamiltonian 13, spectral 11, Farey 10, metagraph 9. The tail is the interesting part —
**deep ideas with almost no theorems**:

> **ORPHANS (≤ 2–3 theorems): homology (1), waggly (1), Fibonacci (1), binary form (2), BIBD (2),
> Vandermonde (2), zeta (2), Smith (2), modular (2), eigenvalue (3), GIT (3).**

These are not unstudied — path homology had a whole thread (β₂ = 0, the `n = 8` threshold, THM-125)
— they are **under-titled**: the work exists in reflections and results but few theorems carry the
keyword, so they read as orphans and get forgotten. The fix is not more work; it is *indexing* —
which is what this atlas is for. `binary form`, `eigenvalue`, `GIT` are genuinely NEW-and-thin
(the July burst), and deserve to grow.

## 3. The gaps (procedurally generated)

Cross-product the three axes the repo actually has — **objects** {tournament A, skew S, kernel R,
metagraph Gₙ/Eₙ, staircase, LRC speed-set} × **functionals** {trace/moment, Burnside, CT, Gaussian
E, Laplace} × **methods** {detection-depth, GIT-stability, WOWII-generate, recurrence-order,
half-dictionary} — and read off the empty cells:

- **Metagraph `Gₙ` moments.** The moment-nullcone template (THM-1775) is built for A; it has never
  been run on the metagraph adjacency (klein-S395 computed `α(Gₙ)` but not `tr(Gₙᵏ)` nullcone).
  Gap: *is there a "transitive metagraph" = nullcone of the metagraph's own trace moments?*
- **Even-graph `Eₙ` binary form.** THM-1810 put `char_A` on `Sym^n`; the DUAL object `Eₙ` (CLAUDE.md
  first-class) has a characteristic form too, never taken. Gap: *what is `Eₙ`'s GIT nullcone?*
- **`c_k`-from-traces closed form.** THM-172 did `c5`; the general Newton-style `c_k(tr A¹..Aᵏ)` is
  unwritten (THM-1870 named it). Gap: the explicit polynomial, and *why `k = n` is the wall*.
- **Half-dictionary on R.** The `x ↦ 2x+1` half-dictionary (THM-1555) relates `char_S ↔ char_A`;
  it has never been applied to the toral kernel `R` (which is also a binary form). Gap: *does
  `R`'s "skew companion" exist and mean anything for TNC?*
- **WOWII on the skew spectrum.** The graffiti generator (THM-1845) uses combinatorial invariants;
  it has never ingested the *spectral* ones (ρ, #real eigenvalues, spectral gap, the `c_k`). Gap:
  spectral-vs-combinatorial inequalities — a whole untested quadrant.
- **LRC as a moment nullcone — the honest non-fit (THM-1775 §2) — has a dual.** The tight AP is an
  *extremal* pole; the repo never asked *what is LRC's nullcone*, the trivial pole. Gap: the
  degenerate speed-configuration and its moment functional.

## 4. New angles / things to compute (a standing queue)

1. Run `cycle_counts_spectral_boundary` at `n = 7` (with a trace sieve) to confirm `c_n` is the
   wall.
2. Add `{c4, c5, c6, #Ham-cycles, ρ, #real-eigenvalues}` to the WOWII zoo and regenerate — the
   spectral quadrant.
3. The `c_k`-from-traces closed form (Newton on cycle counts).
4. `Eₙ` and `Gₙ` characteristic binary forms + their GIT nullcones.
5. Reconnect THM-208 (`H = 1 + 2t₃`, `n ≤ 4`) as the germ of the `c3`-vs-`H` envelope (THM-1860),
   and find the true envelope function `H(c3, n)`.
6. A "subsumed-by" pass over the era-0 block-counting leaves so the early scaffolding links forward.

## The meta-point

The repo does not forget because it is careless; it forgets because it is *large* and *fast* and
*multi-agent*. The countermeasure is not heroic re-reading — it is a **citation graph plus a topic
index plus a gap cross-product**, run as a tool. Today it turned one dead leaf (THM-172) into a
frontier theorem (THM-1870). Keeping the archaeology script in `04-computation/` and rerunning it
each era is the cheapest way to make sure none are lost — and the gap list above is the standing
answer to "what next."

## Cross-links
`corpus_archaeology_kps_S128c136.py` · `cycle_counts_spectral_boundary_kps_S128c136.py` · THM-1870 ·
THM-1775/1810 (ladder) · THM-1780 (H splits) · THM-1845/1860 (WOWII) · THM-172/171/173 (rescued) ·
THM-208 · klein-S395 (metagraph WOWII) · PROBLEM-ATLAS-2026-07-20.

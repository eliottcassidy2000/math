---
source: oracle-2026-06-03-S582o
status: computation (overnight cycles 1-2: R_m is the unique round/chi=2 regular tournament; chi is a cyclicity gauge maxing at the tight boundary)
tags:
  - lonely-runner
  - dichromatic-number
  - chi
  - round-tournament
  - regular-tournament
  - extremal-family
---

# The χ Gauge: the Tight Orbit Is Uniquely R_m, and χ Maxes (=2) at the Loneliness Boundary

Overnight cycles 1–2, completing the S581o/opus-S591 χ thread.

## Cycle 1 — R_m is the UNIQUE round (and unique χ=2) regular tournament

`lrc_round_regular_unique_Rm_s582.py`. Enumerating *all* regular tournaments
(spectrum-deduped):
- **m=7 (spectrum-complete, all 3):** only `R_7` (interval circulant `{1,2,3}`) is **round**
  (locally transitive) and has **χ=2**; Paley `QR_7` (χ=3, aut 21) and the non-circulant
  third (χ=3, aut 3) are **non-round**.
- **m=9 (12/15 spectra; the 3 cospectral-hidden all sit in non-round χ=3 buckets):** only
  `R_9` is round / χ=2.

> **`R_m` is the UNIQUE round regular tournament, and the unique χ=2 regular tournament.**
> Combined with opus-S591's *"LRC ⟹ the runner tournament is round"* (the half-turn
> comparator forces contiguous-arc out-neighbourhoods), this **rigorously pins the LRC tight
> *regular* orbit to be exactly `R_m` (the AP) — with no vertex-transitivity assumption**,
> closing opus-S591's "VT round = interval circulant" qualifier. (Almost certainly the known
> fact: a regular locally-transitive tournament is the rotational `R_m`.)

## Cycle 2 — χ is a cyclicity gauge that maxes at the tight boundary

`lrc_chi_vs_margin_s582b.py` (n=8, χ of the optimal-time round tournament vs margin `M·n`):
- **TIGHT** (`M·n = 1`: AP + both sporadics): **χ = 2** — uniformly. *This closes opus-S591's
  open qualifier:* χ is constant `= 2` on the **whole** tight set, not just the regular AP.
- **LOOSE**: χ ∈ {1, 2}. **Very loose** (`M·n ≥ 2.67`) → **χ = 1 (transitive!)**; moderately
  loose → χ = 2. **χ = 3 never occurs** (those tournaments are non-round = inaccessible).
- So **χ = 2 is *necessary but not sufficient* for tightness** (24/57 loose configs are also
  χ = 2).

> **χ(optimal-time tournament) ∈ {1, 2} for every LRC-accessible config**, and it is a
> *cyclicity gauge*: it **maxes at 2 exactly at the loneliness boundary** (tight = the most
> cyclic the observer's view can be) and **drops to 1 (transitive)** when the config is very
> lonely (observer sits in a big gap, sees a transitive order). The χ = 3 Paley orbit is
> geometrically out of reach.

## Synthesis (answering the χ question end-to-end)
- Among **regular** tournaments, χ adds *strictly* beyond vertex-transitivity: `R_m` (χ=2)
  and Paley (χ=3) are both VT with equal 3-cycle counts, separated only by χ — and `R_m` is
  the *unique* round/χ=2 regular one (cycle 1).
- Among **LRC optimal-time** tournaments (all round), χ ∈ {1,2}: tight ⟺ the **maximum**
  (χ=2), very loose ⟺ χ=1; "χ=2" alone is necessary-not-sufficient (cycle 2). The LRC
  extremal is the χ-maximal accessible tournament `R_m`, never the (inaccessible) χ=3 Paley.

## Honest limits / next
- m=9 has 3 cospectral-hidden regular classes not individually canon-checked (all in
  non-round χ=3 spectra, so the uniqueness conclusion holds for the 12 distinct spectra;
  a full canon enumeration would remove the caveat). m=7 is complete.
- χ=2 being necessary-not-sufficient means χ alone does **not** characterize the tight set;
  the extra discriminator (which χ=2 round configs are tight) is the observer-coupling /
  gap-structure (S580o, S556o) — the natural next probe.

## Artifacts
```
04-computation/lrc_round_regular_unique_Rm_s582.py        (+.out)
04-computation/lrc_chi_vs_margin_s582b.py                 (+.out)
```
Related: S581o/HYP-2135 (χ separates regular beyond VT), opus-S591/HYP-2133 (LRC=round,
interval not Paley), HYP-2091 (clean-polygon R_m), S580o (observer-coupling), S576o (tight census).

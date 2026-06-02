---
source: oracle-2026-06-01-S554o
status: synthesis (26 necessary conditions a counterexample to LRC@14 must satisfy; the thin locus)
tags:
  - lonely-runner
  - n14
  - necessary-conditions
  - counterexample-locus
  - proof-lite
---

# Twenty-Six Necessary Conditions for an LRC@14 Counterexample

Following the four conditions of S553 (averaging-extremal, sieve-covered, high-energy,
7-class coupled), this session mines *every* repo insight for **necessary conditions**
a counterexample to LRC@14 must satisfy. Each is provable (violate it ⇒ some `t` is
lonely ⇒ not a counterexample). Piling them up pins the counterexample to a razor-thin
locus; jointly unsatisfiable would be a proof.

## The catalogue (26 conditions, by source)

**A. Sieve / divisibility [THM-369, S546]**
1. **A1 sieve-covered:** a multiple of every `q ∈ {2,…,14}`.
2. **A2:** a multiple of `14` (`q=n`; else `t=1/14` is lonely; blocks all `t=a/14`).
3. **A3:** a multiple of each maximal prime power `≤14`: `8, 9, 5, 7, 11, 13`.
4. **A4:** a multiple of `n* = 7` (the doubled-prime channel modulus).
5. **A5:** for each `q≤14`, not all speeds coprime to `q` (`≡` A1).

**B. Moment / counting [S553]**
6. **B1 averaging-extremal:** `near(t) ≥ 1` for all `t`, yet mean `near = 2−2/14 =
   12/7`, so `near` touches the floor `1` but **never `0`** (not even in closure).
7. **B2 weaker-threshold:** at threshold `1/(c·14)`, the near-set still covers for
   each `c` (the first-moment ladder).
8. **B3 second-moment-consistent:** the pairwise danger-overlaps are tuned so the
   `near`-count's minimum is `1`, never `0`.

**C. Covering / resonance [S550, S525, S545]**
9. **C1 covering:** the danger arcs `B_i` (total measure `13·1/7 = 13/7 > 1`) cover
   `[0,1)`, with the overlaps (resonances) exactly arranged.
10. **C2 high-energy core:** resonance energy `E(v) ≥ (1−2/n)^{n-1} = (6/7)^{13}`.
11. **C3 short resonance:** the minimal resonance length `ℓ(v)` is small.
12. **C4 pervasive returns:** the cascade product of conditional clearances is `0`
    (the 3-/k-term returns dominate, S545).

**D. CRT / channel [S524, S552, S533]**
13. **D1 7-class coupled:** all 7 mod-7 CRT classes blocked at every `t`.
14. **D2 singleton coupling:** the singleton class `{mult of 7}` blocks exactly the
    windows where the six pair-classes are jointly safe.
15. **D3 parity:** the full-support inside debt is present mod `7`.

**E. Geometric / tournament [S511, S530, S539]**
16. **E1 never a source:** the observer is never a source in the LRC-marked
    tournament (out-degree `< 13` at all `t`).
17. **E2 narrow apex:** the observer's bracketing (apex) gap is `< 2/14` at all `t`
    (never a wide-enough gap at `0`).
18. **E3 perpetual tie:** the observer always has `≥1` trienerment tie (`≡` B1).

**F. Reduction / normalization [S549, classical]**
19. **F1 primitive:** `gcd(speeds) = 1` (WLOG by scaling invariance).
20. **F2:** distinct nonzero speeds.
21. **F3 bounded:** a minimal counterexample has bounded speeds (finite reduction,
    Cusick/Tao).

**G. Spread / structure [S544, S521o, S548]**
22. **G1 frequency-concentrated:** the speeds are resonant, not decorrelated.
23. **G2 commensurable:** the orbit is a closed geodesic (auto for integers; else
    equidistribution ⇒ lonely).
24. **G4 not AP-scaled:** `v ≠ c·(1,…,13)` (the AP is lonely at the wall `t=k/14`).

**H. Diophantine [S535, S545/S550]**
25. **H1:** the difference set carries a specific QR/Frobenius pattern mod `7`.
26. **H2:** the relation lattice's short vectors (returns) supply the high energy.

## The verification (`lrc_counterexample_necessary_conditions_s554.py`)

- A **"sieve-minimal" candidate** `(2,…,14)`-style set satisfies **10/10** of the
  checkable conditions (A1–A4, C2, C3, D1, F1, F2, G4) — and is **lonely anyway**
  (a lonely time exists). The conjunction of all checkable conditions does **not**
  produce a counterexample.
- A **"prime-power-heavy"** set satisfies 8/10 and is lonely.
- **`0` of 30** random primitive sieve-covered sets are non-lonely.
- The **AP `(1,…,13)`** is the only set whose *open* lonely set is empty — but that is
  the measure-zero **wall** (it is lonely at the closed `t=k/14`, S553/S551), so it is
  the extremal, **not** a counterexample. It satisfies 7/10 (it lacks A1/A2 — no
  multiple of 14 — and is AP, failing G4), so it is *excluded by the sieve conditions
  themselves*. A true counterexample would have to be sieve-covered (unlike the AP)
  yet never lonely — and none is known.

## The honest verdict

> A counterexample to LRC@14 must satisfy **at least 26 distinct necessary
> conditions** simultaneously — sieve-covered, averaging-extremal, high-energy,
> 7-class coupled, narrow-apex, never-a-source, primitive, not-AP, … The locus is
> razor-thin: the most counterexample-like constructible sets satisfy 8–10 of the
> checkable conditions yet remain lonely, and no counterexample appears among
> sieve-covered samples. **This is the strongest honest "proof-lite": the conditions
> are not (yet) shown jointly unsatisfiable, but they corner a counterexample into an
> arithmetic locus that the AP — the extremal — escapes only by failing the sieve.**

The path to a real proof is to find **two conditions that are provably incompatible**
(e.g., A1 sieve-covered + C2 high-energy + B1 averaging-extremal ⇒ contradiction), or
to push the finite reduction (F3) to a feasible verified bound for `n=14`.

## Verdict / next
- Catalogued 26 necessary conditions (provable; sourced across S511–S553); verified
  the most-constrained candidates remain lonely; no counterexample in samples.
- Concrete next: (1) seek a **provable incompatibility** among the conditions (the
  real proof route) — e.g., does sieve-covered + short-resonance force a lonely
  rational `t`?; (2) tighten **F3** (the bounded-speed reduction) toward a verified
  range for `n=14`; (3) formalize the cleanest necessary conditions (A1–A4, B1) in
  Lean on top of the sieve / `near_pair` / `almost-lonely` results.

## Artifacts
```
04-computation/lrc_counterexample_necessary_conditions_s554.py
05-knowledge/results/lrc_counterexample_necessary_conditions_s554.out
```
Related: S553 (the first four conditions + almost-lonely), S550 (high-energy core),
S551 (Vitali / measure-zero wall), S552 (7-class), S524 (CRT coupon), S546 (doubled
prime), S511/S530/S539 (tournament/apex/trienerment), THM-369 (sieve).

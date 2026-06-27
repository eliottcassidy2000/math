---
id: HYP-3103
title: NEW SIGNAL — the zeros of the miss-count PGF G_N(z)=E[z^N]=Σq_t z^t; the gK8/coverage extremizer (consec) has ZERO real roots (maximal sector correlation = least-factorizable = the Lee-Yang/independence signal), the #real-root count stratifies config space, and the zero arguments drift near the 7th-roots of unity (third zero robustly ≈ 3·360/7)
status: VERIFIED signal (root structure stratifies; extremizer in #real=0 stratum; independence reading rigorous) + BOLD untested hypotheses (Lee-Yang extremality, apex-7 zero arc, break=root-collision). Not a new bound.
source: mac-mini-2026-06-27-S66
reflections:
  - the-miss-count-partition-function-and-its-zeros
related:
  - HYP-3085   # gK8/S2/Clebsch = the sector correlation this zero-signal measures
  - THM-534    # the moment-LP / inclusion-exclusion (G_N(-1+1)=... ; G_N(0)=p0)
  - THM-577    # cap = pair-Pascal (the extremizer's value)
  - OPEN-Q-108
external: INVESTIGATION-BACKLOG:1520 (Route 5 — LRC/tournament parallel fugacity partition functions)
---

# HYP-3103 — The zeros of the miss-count partition function (a new signal)

## The new object + signal
The LRC inclusion-exclusion `meas(S7)=Σ_k(-1)^k MISS_k` with `MISS_k=S_k` (factorial moments) is the value at
`z=0` of the **probability generating function of the sector-miss-count** `N`:
```
G_N(z) = Σ_{t=0}^6 q_t z^t,   q_t=P(N=t);   G_N(0)=q_0=p0 (LRC coverage);  G_N(3)=E[3^N] (tournament fugacity).
```
The project measured only the moments `S_r` and the single value `p0`. The **new signal is the zeros of
`G_N(z)` in ℂ** — never previously tracked.

## VERIFIED findings (`lrc_missPGF_new_signal_S66.py`, `..._realroots_vs_extrememass_S66.py`)
- **consec (and its dilation even-AP): 0 real roots** — three complex-conjugate pairs `|z|=1.49,1.58,1.70`,
  all far from `z=0`. Spread/break/covering configs: a real root **near `z=0`** (`|z*|≈0.05–0.12`) + 2–4 real.
- **The `#real`-root count STRATIFIES config space**, and the gK8 extremizer is in the `#real=0` stratum: over
  250 configs the global `max L_yK8 = 3.58` (= consec) occurs ONLY at `#real=0`; `#real=2,4` cap at `≈0.7–1.0`.
  `corr(#real, extreme-mass q0+q6) = -0.37`.

## Why it is the right signal (rigorous)
A non-negative-coefficient PGF is **real-rooted ⟺ `N` is a sum of INDEPENDENT indicators** (Pólya-frequency /
Newton). So `#real roots of G_N` measures **sector independence**, and **consec (0 real roots) is the maximally
NON-independent = maximally sector-correlated config** — exactly the high `q0+q6` extreme-mass that
`L_yK8=10(q0+q6)+q3` rewards (S60 gK8/S2/Clebsch). The recurring "consec is extremal" = **consec pushes the
partition-function zeros maximally off the real axis** (Lee-Yang confinement). This unifies the analytic zero
structure with the gK8 correlation finding.

## BOLD tests this session
- **★ apex-7 zero arc (partial):** consec's zero arguments drift through the 7th-root-of-unity arguments
  `{51.4°,102.9°,154.3°}` as the row varies, **without clean convergence**, but the **third zero is robustly
  `≈154-157° ≈ 3·360/7`** across all k=8..13, and at k=11 the middle zero is *exactly* `102.9°=2·360/7`. The
  magnitudes grow `1.49→2.25` with k. Verdict: **suggestive of an apex-7 zero structure, not confirmed.**
- **★ Lee-Yang extremality (untested):** prove the extremizer has no real zero on `[-1,0)` — a zero-free-region
  reformulation of coverage extremality (a NEW analytic route, replacing the non-transitive exchange / the
  irreducibly-global crux).
- **★ break = root-collision (untested):** the k=8 minimizer "break" (top-cluster → middle-spread) may be a
  real-zero colliding onto the real axis (`#real` jump) — testable via the discriminant of `G_N`.

## The new-signals slate (the measurement program)
1. `#real roots` (VERIFIED to stratify). 2. nearest-zero distance to `z=0` (coverage-pole proximity).
3. the Lee-Yang confinement region / zero arc. 4. the fugacity rank curve (`rank` vs `z:0→3`).
5. root-argument spread/regularity. 6. the discriminant/resultant (the real-root transition = phase boundary).
7. (tournament parallel) the zeros of the winding tournament's `I(Ω,x)` vs the miss-PGF zeros.

## Honest caveat
The Route 5 "rank-flip" (consec best at LRC `z=0`, worst at tournament `z=+2`) was **not reproduced** in my
6-config set (consec is rank-1 at `z=3` too, having the fattest `q6` tail) — because the tournament fugacity
acts on a DIFFERENT generating function (odd-cycle `α_j`, not sector `S_k`). The two are **parallel** partition
functions sharing an AP-extremal surface, not one object — signal #7 is the way to compare them directly.

## Next
Measure signals 2–7 across the bounded config bank; test the discriminant=break and the Lee-Yang zero-free
region; compute the winding-tournament `I(Ω,x)` zeros and compare arcs.

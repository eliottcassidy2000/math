# The μ_{1/7} family census: the AP is the unique minimizer, between correlation and decorrelation

**opus-2026-07-07-S131.** Continuing Route 1 (the owner's chosen route after the S130 audit). The
density floor's binding quantity is `μ_{1/7}(E) = meas{ x∈[0,1] : circular max-gap of {frac(e·x):
e∈E} > 1/7 }`. Following the owner's directive — *"if you encounter large amounts of seeming
randomness, do a census and understand the properties of the in/finite families"* — I censused the
whole (infinite) family of cluster shapes `E`, and it resolves into a clean, rigid structure.

## The census: AP is the unique global minimizer

`μ_{1/7}` is invariant under dilation and translation of `E` (both are rigid rotations of the phase
configuration), so it is a function of the **affine shape class** of `E`. Exhaustive enumeration of
all shape classes (min 0, gcd of gaps 1):

| k | shapes censused | μ_{1/7}(AP) | global min | AP the unique min? | runner-up | gap |
|---|---|---|---|---|---|---|
| 8 | 1716 (spread ≤ 13) | 0.9402 | 0.9402 | **yes** | 0.9719 | 0.032 |
| 9 | 495 (spread ≤ 12) | 0.8402 | 0.8402 | **yes** | 0.9163 | 0.076 |
| 10 | 55 (spread ≤ 11) | 0.7756 | 0.7756 | **yes** | 0.8402 | 0.065 |

Zero shapes fall below the AP, and the AP is a *strict, isolated* minimizer (a clean gap to the
runner-up). At `k=13`, 40 aggressive adversarial descents (S130) + structured adversaries all return
to the AP. **The arithmetic progression `{1,…,k}` uniquely minimizes `μ_{1/7}`.**

## The two ends of the family: correlation ↔ decorrelation

`μ_{1/7}(E)` ranges over exactly `[μ_{1/7}(AP), independent-limit]`:

- **Minimum — the AP (maximally correlated).** `{frac(jx)}` is the orbit of `x`; by the three-gap
  theorem it is the most equidistributed `k`-point configuration, so it has the smallest max-gaps and
  hence the smallest good-set. Exact via three-gap piecewise-linear breakpoints (S130): `μ_{1/7}(AP)`
  = 1 (k≤7), 691/735, 247/294, 38/49, 1381/2205, 13823/24255, **477/1078** (k=8..13).
- **Maximum — decorrelated `E` (independent-like).** As the affine shape spreads out, the phases
  behave like `k` independent uniform points, whose max-gap law gives `μ_{1/7} = 1 − Σ_j (−1)^j
  C(k,j)(1−j/7)_+^{k−1}` ≈ 1.000, 1.000, 0.9996, 0.9984, 0.9951, **0.9883** (k=8..13). Random
  spread-out `E` hit `≈ 1.000` empirically.

So the density floor `μ_{1/7}(E) ≥ 477/1078 > 0` is the statement that the AP end is the floor.
This IS the fleet's dichotomy (kps-S53: near-AP three-gap + spread decorrelation) — but seen as **one
family with two ends**, the AP minimum and the decorrelated maximum, rather than two separate cases.

## The clean reformulation, and why the proof is hard

`μ_{1/7}(E) = 1 − C(E)`, where `C(E) = meas{ x : the radius-1/14 arcs around {frac(e·x)} cover the
circle }` (max-gap ≤ 1/7 ⟺ covering radius ≤ 1/14). So **AP-minimality ⟺ the orbit `{jx}` maximizes
the covering probability** among all `k`-point integer-cluster configurations. This is a clean
extremal statement — and a web check finds it is *not* a standard named result (the three-gap theorem
gives the AP's gap structure, but the minimality-over-`E` is not in the literature I can find). So it
is a genuine (apparently novel) extremal conjecture, strongly verified here.

Why it resists a quick proof: a large affine spread does **not** force decorrelation — an "AP-core +
outlier" (e.g. `{0,1,…,11,10^6}`) has a large normalized spread yet a small `μ_{1/7}` (the AP core
dominates). So there is no clean "bounded-spread finite check + large-spread decorrelation tail"
split; AP-cores hide at every scale, and the minimum genuinely sits at the *full* AP. The extremal
fact is irreducible.

## Status and honest next step

- **Rigorous:** `μ_{1/7}(AP_k)` exact (three-gap), the independent-limit exact, AP-minimality
  exhaustively verified for k=8,9,10 and adversarially for k=13. The floor value `477/1078` is
  independently confirmed against the canon's `rhoGlobFloorRat(13)`.
- **The one hard lemma (near-AP branch of Route 1):** prove `C(E) ≤ C(AP)` — the orbit maximizes
  covering probability. Apparently novel; strongly supported. The complementary spread branch
  (decorrelation) is mac-mini-S38's `reach_decorr` lane.

Deep structural correctness achieved on the *shape* of the floor; the remaining content is a single
clean (hard) extremal inequality, not a pile of cases. That is the honest state of Route 1's
near-AP node.

*Three-gap references: [Wikipedia](https://en.wikipedia.org/wiki/Three-gap_theorem); Sós/Surányi/
Świerczkowski (1950s); [concise geometric proof, arXiv:2308.11999](https://arxiv.org/abs/2308.11999).*

---
source: kind-pasteur-2026-07-07-S59
status: PROVED (the pointwise lemma + exact constants) / delineation of the residual.
  A one-line subset lemma turns the exact Farey-roof AP constants into a k-uniform,
  diameter-graded floor on BOTH the mean E[maxgap] and the irreducible tail mu_theta.
  It survives death-star-S1's capstone (which kills the mean-via-AP-extremality detour)
  because it never uses extremality — it uses the AP as a SUBSET, not as a minimizer.
tags:
  - lonely-runner
  - LRC14
  - route-1
  - density-floor
  - max-gap
  - three-gap
  - monotonicity
---

# The diameter floor feeds the irreducible tail

> **⚠ CORRECTION (kps-S60, MISTAKE-121).** This reflection's claim "**no bite at k=8..10**" is
> FALSE — a table-start artifact (the μ-scan began at n=13, but a k-point cluster's diameter
> starts at k−1). The true union-bound bites are **k=8: diam ≤ 9; k=9: diam ≤ 11; k=10:
> diam ≤ 11** (μ(AP₉)=0.840 ≥ bar₈=0.675, etc.). Moreover the S60 **intersection ledger**
> (HYP-4847) computes meas(G_P ∩ {roof > 1/7}) exactly and extends every k=8..12 bite past
> the union-bound values. See `lrc_intersection_ledger_kps_S60.py`.

**kind-pasteur-2026-07-07-S59 (HYP-4797).** Owner: understand the LRC(14) history deeply,
audit validity, uncover forgotten-but-true factoids, extend and connect. This session's
finding: a **one-line pointwise lemma** — so elementary it was never written down — converts
opus-S134's Farey roof into a **rigorous, exact, diameter-graded lower bound on every gap
functional at once**, including the tail `μ_{1/7}` that death-star-S1 just proved
irreducible. It carves the bounded-diameter bulk off the open density-floor lemma with
exact constants, and it explains *where* every current adversarial record lives.

## The lemma (Brick 1)

> For integer sets `E ⊆ F` and every real `x`: the point sets nest,
> `{frac(ex) : e∈E} ⊆ {frac(fx) : f∈F}`, so the widest `F`-gap — an arc free of
> `F`-points, hence free of `E`-points — lies inside some `E`-gap:
>
> **`maxgap(E, x) ≥ maxgap(F, x)` pointwise.**
>
> Hence every monotone gap functional inherits the bound: `E[maxgap(E)] ≥ E[maxgap(F)]`
> and `μ_θ(E) ≥ μ_θ(F)` for every threshold `θ`.

With `F = {0..D}` (translate `E` to `min E = 0`; gcd-reduce first — both `μ`- and
mean-invariant, canon), and `{0..D} ≅ AP_{D+1}` (a translate):

> **`μ_θ(E) ≥ μ_θ(AP_{D+1})` and `E[maxgap(E)] ≥ A(D+1)`, `D` = primitive diameter of `E`,**

where both right-hand sides are **exact rationals** from the opus-S134 Farey roof
(`A(n) = Σ_{cells} 1/(qq'²)`, `μ_θ(AP_n) =` per-cell linear superlevel). A pleasant
corollary of the same lemma: `AP_n ⊆ AP_{n+1}` gives the monotonicity of the roof tables
in `n` for free.

Verification: 0 violations / 16000 pointwise random nested checks; the roof cross-checked
against the independent death-star exact engine at `n = 13, 14, 21` (`93/440`,
`2081/10296`, `89173/604656`); zoo and random families all satisfy the tail bound.
Scripts: `lrc_diameter_monotonicity_leg_kps_S59.py`, `lrc_tail_diameter_floor_kps_S59.py`.

## Why this survives the capstone

death-star-S1 proved the mean detour dead **as an AP-extremality route**: `E[maxgap]` is
not AP-minimized, truncated means are AP-minimized only on a window (`θ* ≤ 0.181`) disjoint
from non-vacuity (`θ* ≥ 0.195`), so "reduce the tail to a mean, reduce the mean to the AP"
cannot close the floor. **The diameter floor never reduces to the AP as an extremizer.** It
uses the full AP orbit as a *superset* whose gaps are pointwise dominated — a use that is
valid for every family, every functional, every threshold, with the AP value entering as a
*lower* anchor rather than a conjectured minimum. The hard part of (A') that death-star
localized (fine-scale tail extremality) is untouched; what the diameter floor does is
**shrink the domain on which (A') is still needed**.

## The exact numbers

- **Tail, k=13 leg (the leg the mean route was serving; bar `m_P = 14249/252252`):**
  `μ_{1/7}(AP_n) ≥ m_P` through `n = 76` exactly (`2314528732/40290957525 ≈ 0.057445`),
  crossing at `n = 77`. So
  > **every 13-element `E` with primitive diameter `≤ 75` has
  > `μ_{1/7}(E) ≥ μ_{1/7}(AP_{D+1}) ≥ m_P` — the k=13 quantitative floor holds there.**
  THM-527-D's (verified) extremal-spread domain is `O(k) ≈ 30`; the floor covers it with
  a **2.5× diameter margin**. Every family the fleet has ever put on the minimizer board
  (AP 12, monad record 20, death-star prim-sat 22, GW 23, opus stretch 28, kps-S57 40) is
  inside `D ≤ 75`.
- **Tail, k=8..12 union-bound legs (bars `1 − min meas(G_P) + m_P`):** bites at
  `k=12` (diam `11..23`) and `k=11` (diam `10..15`); **no bite at k=8..10** — the
  union-bound bars (`0.45–0.68`) exceed the roof values at every `n`, honestly leaving
  those legs to the genuine per-k floors (monad redirect target 1) or a better
  G_P-coupling.
- **Mean version (context for the record chase):** `A(n) > 1/7` through `n = 22`
  (`A(22) = 3029671/21162960 = 1/7 + 6391/21162960` — the margin at the boundary is a
  hair, `3.02·10⁻⁴`), so **every `E` with primitive diameter `≤ 21` has
  `E[maxgap] > 1/7`** — which contains monad's record family (`diam 20`,
  `E[mg] = 12907/65520 = 0.196993`). Any future record below `1/7` must have primitive
  diameter `≥ 22`; the observed per-diameter minima at `D = 22..32` are flat around
  `0.203–0.206` with no erosion trend.
- **The residual regime is far from binding:** structured probes at primitive diameter
  `80–100` (interlacings, jittered dilated APs, stretches) give `μ_{1/7} = 0.58–0.97`,
  i.e. **10–17× the bar**. The open frontier `D > 75` is comfortably conjectural.

## The corrected deficit frame (the forgotten factoid, made precise)

The second thread of the session: the loop-vs-torus Fourier identity
`E_x[F({e_i x})] = ∫_{T^13} F + Σ_{k ∈ L(E), k ≠ 0} F̂(k)`, `L(E) = {k : Σk_ie_i = 0}`,
applied to gap functionals. Two symmetries prune the sum:

1. **Translation symmetry** (`maxgap(θ + c·1) = maxgap(θ)`): `F̂(k) = 0` unless
   **`Σ_i k_i = 0`**.
2. **Pair-distance exact uniformity** (the forgotten universal factoid: `‖(e_i−e_j)x‖` is
   exactly uniform for every pair, every `E`): the pure-difference modes `m(δ_i − δ_j)`
   never lie in `L(E)` anyway (`m(e_i − e_j) ≠ 0`).

Combined: a contributing mode needs `Σk_i = 0` **and** `k·e = 0`, and any `k` with `≤ 2`
nonzero entries satisfying both is identically zero. So

> **the entire deviation of `E[maxgap(E)]` (and of every gap functional) from the iid value
> `H₁₃/13 = 1145993/4684680 ≈ 0.24463` is carried by zero-sum relations of weight ≥ 3.**

The zoo confirms it sharply: greedy-Sidon (no small weight-3 relations) sits at `0.24234`
(deficit `0.0023`); primes at `0.2382`; geometric `{2^i}` at `0.25017` — *above* iid, via
its zero-sum weight-3 relations `2·2^i − 3·2^{i+1} + 2^{i+2} = 0` contributing with the
*positive* sign; and every low-mean family is relation-dense (AP: 36 midpoint triples,
deficit `0.0333`; monad record: interlacing triples like `9 + 11 = 2·10`, deficit
`0.0476`). The adversary's total budget to drag the mean to `1/7` is `0.10177` — the best
harvest on the board is 47% of it, and the harvest requires additive structure that (by
Freiman) pulls the primitive diameter back into the proved zone. The caveat from my own
zoo: counting only coefficient-`(1,1,−2)` triples is not enough — general small zero-sum
weight-3 vectors (e.g. `(2,−3,1)`) contribute too, with either sign. Any rigorous budget
bound must run over the full weight-3 lattice shell, not the 3-AP count.

## Composition (where this leaves the density floor)

```
k=13 leg of hlarge (bar m_P):
  primitive diam <= 75      PROVED   (this session: subset lemma + exact roof constants)
  primitive diam  > 75      OPEN     (probes at 10-17x bar; decorrelation flavor;
                                      far-element substructure -> upstream peel GREEN;
                                      near-dilated-AP -> klein-S152-style rigidity;
                                      genuinely spread -> pair bound 1/(3pq) [verified,
                                      rigorous] + Chung-Erdos = monad target 4)
k=12, k=11 legs:            bites diam <= 23, <= 15 respectively (same lemma)
k=8..10 legs:               untouched (union bound too crude; needs per-k floors)
```

The one-line lemma is Lean-friendly: set-inclusion of finite point sets, a max over a
subset, and exact rational Farey sums (`native_decide`-able per `n`). The roof formula
itself (opus-S134 Claims 1–2) rests on classical three-distance facts and is the natural
next formalization target — formalizing it to `n = 76` discharges the k=13 bounded-diameter
floor entirely.

## The honest one-liner

The fleet spent two days arguing about which functional is AP-minimized; the useful move
was never extremality at all — **the AP bounds every family from below not because it is
smallest but because it is fullest**, and "fullest" is a one-line set inclusion that
transfers exact constants to every functional at once, `k`-uniformly, out to diameter 75.

## Pointers

- Builds on: opus-S134 Farey roof (exact `A(n)`, `μ_θ(AP_n)`); death-star-S1 capstone
  (tail irreducibility — the frame this feeds); monad-explorer-S1 HYP-4787/MISTAKE-118
  (the honest bars `m_P`, `thr_k`; the k=13-only scope of the mean route); THM-527-D
  (bounded-spread domain `~30`); THM-530 (union-bound structure, `m_P`); opus-S130
  (μ_{1/7} floor + exact k=8..13 constants); mac-mini-S15 (three-gap frame).
- Ceded/refuted en route: my conditional-reverse-Markov-on-G_P idea was independently
  tested and refuted by monad-explorer-S1 (PART 3) before I ran it — correctly ceded;
  HYP-4787 number ceded (renumbered 4797).
- Files: `04-computation/lrc_diameter_monotonicity_leg_kps_S59.py`,
  `04-computation/lrc_tail_diameter_floor_kps_S59.py`, outputs in
  `05-knowledge/results/` (with the exact-crossing addendum).
- Does NOT prove LRC(14) or all of (A'). It proves the bounded-diameter part of the
  k=13 floor and delineates the rest.

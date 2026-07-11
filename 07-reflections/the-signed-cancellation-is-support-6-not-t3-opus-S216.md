---
source: opus-2026-07-11-S216
status: SYNTHESIS — the "cross-triple signed cancellation" that closes LRC(14) is the support-6
  relation-lattice Minkowski count (THM-538 + HYP-2608a), NOT the t≥3 character-layer sum (THM-684,
  superseded). Every route converges here. Backed by 4 concept-cluster scouts + the transfer THM-685.
tags:
  - lrc14
  - signed-cancellation
  - THM-684
  - THM-685
  - THM-538
  - HYP-2608
  - measure-floor
  - convergence
---

# The signed cancellation is support-6, not t≥3

**opus-2026-07-11-S216.** Asked to close "the t≥3 signed cancellation" (THM-684) — the last analytic
ingredient of the OffLine/Fourier route — I ran four concept-cluster scouts across the repo and cross-read
the transfer, density-floor, and covering-reformulation canon. The verdict is decisive and worth recording
in one place, because the fleet has circled this for months from six directions and the picture is now
crisp.

## 1. The t≥3 character-layer cancellation is superseded — do not chase it

`LM/Q = (b/Q)¹³ + Σ_t layer_t` (THM-684, exact). Death-star's −5/7 self-similarity
`A(ψ)=(1−2β)conj(ĉ(ψ))` is exact (`g²=(1−2β)1_B+β²`; verified 5e-15) and collapses one Cauchy–Schwarz
level, so the t≥3 layers *are* signed ĉ-weighted t2-transforms. But **every agent who tried to close t≥3
by bounding that recursion hit the identical wall**: the layer series is order-one, alternating, and
**non-truncatable** on relation-stacked families (THM-685(B), PROVED-exact — even generic GEN has
`layer₃=−0.073` overshooting the total deficit `−0.038`, with `R≥4` pulling back `+0.034`). Absolute
per-triple bounds run 55× over budget and don't decay. **"The cancellation is the theorem" is true, but the
character-layer sum is the wrong object to bound.**

The resolution — **THM-685(A), the Kronecker transfer**, PROVED and formalized sorry-free in
`LRCMeasureTransfer.lean`:

> `|LM(q) − q·μ(S)| ≤ K(S) ≤ Σ_l v_l`,  μ(S) = meas{α : frac(v_l α) ∈ [1/14,13/14] ∀l}.

You never sum the layers. You compute the exact Kronecker-line measure μ(S) in one breakpoint sweep, and a
floor `μ(S) ≥ μ₀` gives a live multiplier at every `q > Σv/μ₀` — certification a-priori complete. **The
entire remaining analytic content moved from "t≥3 signed cancellation" to the measure floor μ(S) > 0**
(`SafeMeasureFloor` in `LRCResidualMeasureFloor.lean`: the single open hypothesis of `lrc14_of_measureFloor`).

## 2. The measure floor is one extremal lemma: **AP minimizes μ**

Raw safe measure, computed exactly this session
(`04-computation/lrc14_ap_minimizes_mu_koksma_route_opus_S216.py`):

| family | μ(S) exact | μ − (6/7)¹³ |
|---|---|---|
| AP {1..13} and dilate 3·{1..13} | **0.000000** | −0.1348 |
| near-AP (one detune) | 0.033 | −0.102 |
| mild spread | 0.125 | −0.010 |
| GEN dissociated | 0.117 | −0.018 |
| wide dissociated (max 401) | **0.1354** | +0.0006 |

The AP is the exact μ-minimizer (an isolated tight point, μ=0 — matching the discrete `LM/Q=0` for the AP);
spread raises μ monotonically to the decorrelated value `(6/7)¹³=0.1348`. **"AP/consecutive minimizes μ"** is
exactly THM-530/THM-657/THM-527's single open lemma (HYP-2602, HYP-2607, HYP-2608), the k=8,9,10 wide-spread
bound `span(E)>B(k) ⟹ meas ≤ cap_k`. Six independent angles (moment-LP dual THM-534, Beurling U4 THM-537,
covering reformulation THM-657, G_P union bound THM-530, tent/PZ THM-660/661, transfer THM-685) reduce
LRC(14) to it. All say the same words: **verified exhaustively, zero exceedances, NOT proved.**

## 3. Where the real signed cancellation lives: a lattice reciprocal-product sum

Why do compact census + pairwise decorrelation not finish it? Because the deficit is a **signed sum over the
offset-relation lattice** (HYP-2606):

> `meas(S7(E)) = M7(k) + Σ_{0≠n∈Λ°(E)} K(n)`,   `K(n) = D7(n mod 7) / ∏_j n_j`   (coset factorization, HYP-2646),

with `D7` a finite mod-7 character coefficient (`|Re D7| ≤ 0.1431`, nonzero on all 46656 cosets). The naive
absolute envelope `Σ|K(n)| = Σ 1/∏|n_j|` **diverges** — the barely-covers wall (arc length 13/7 > 1; pair
masses `1/49` exact; Bonferroni useless). The signed sum is only **conditionally convergent**: the mod-7
character `D7` oscillates and the reciprocal-product tail must be summed *with signs*. **The cancellation is
not a luxury; it is the theorem.**

⚠️ Honest correction (THM-538, kps-2026-06-19): the clean "annihilates support ≤ 5" holds for the
*active-coordinate* sum `Q(n)`, but is **FALSE for the zero-padded kernel `K(n)`** in the measure — the zero
coords carry `(1−|T|/7)^z` factors that break it, so short (support 2–5) relations **do** contribute (the AP
is support-3-dominated). The correct ruling quantity is the **support-6 relation *density* R6** (HYP-2645),
not a hard floor: wide-spread admissible sets have R6 small (corr ≈ 0), which de-risks the count.

This is why my validation splits the way it does:
- **support-2 (pairs):** `|A(a,b)−(6/7)²| ≤ (24/7)/max(a,b)` (Koksma, THM-686(C), PROVED — the continuum twin
  of my formalized t2 bound LEM-022) holds (0 violations); the pairwise deviations shrink with spread.
- **higher support:** the conditionally-convergent tail `Σ D7(n mod 7)/∏ n_j`, ruled by R6 and closable by a
  **Minkowski / successive-minima count** on Λ°(E) — the lattice coupling is what makes the signed sum
  converge where the free sum diverges (HYP-2608a; "the single highest-value next step" per OPEN-QUESTIONS).

**The sharp structural connection.** `Σ_{n∈Λ°(E)} D7(n mod 7)/∏_j n_j` is a **character-twisted reciprocal-
product sum over a lattice** — precisely the higher-rank analog of the t2 hyperbola sum I *already formalized*
in LEM-022: `Σ_{h≠0} 1/(cdist h · cdist wh)`, bounded by death-star's `harmonic_ratio_sum_mul_le`
(`S·P ≤ 20(log q)²`) via the rank-1 ratio lattice `P(w)`. My t2 result **is** the rank-1 (two-coordinate)
case of the exact object that closes LRC(14); the last gap is the same sum at rank ≤ 5 with the D7 twist. So
the user's instinct — a **cross-relation signed cancellation** is the crux — is correct, and it is *the same
species* of sum as the one already in Lean, one rank up. Its right form is the lattice reciprocal-product
count, not death-star's t≥3 character layers (bypassed by the transfer). Both are images of THM-681's W₀ (the
connected cascade's non-vanishing part), and only the lattice/Minkowski form converges.

## 4. What this means for the formalization

`LRC14Statement ⟸ LRCUpTo13 (cited) + SafeMeasureFloor`, kernel-pure, every other branch discharged with
foundational axioms only. The single remaining mathematical obligation is the measure floor, i.e.
AP-minimizes-μ, i.e. the support-6 Minkowski count. The t2 half of the Fourier route (LEM-022, five-agent
composition, S215) is done and feeds the support-2 level; it is real infrastructure but **not on the
critical path to the last gap** — the last gap is a lattice-point count, not a character-sum bound.

**Actionable handoff:** stop investing in t≥3 character-layer cancellation bounds. The live target is the
support-6 successive-minima count (THM-538 + HYP-2608a) — execute the Minkowski `|K(n)| ≤ c₁⁶/(λ₁···λ₆)`
count on Λ°(E), which converts the 0-exceedances-over-40k into a gap-free proof of the k=8,9,10 rows and
thence LRC(14).

→ THM-684 (superseded t≥3), THM-685 (transfer, the bypass), THM-538 (support-6 floor — the key localization),
HYP-2608a/2602/2607 (the open lemma), THM-530/657/527 (AP minimizes μ), LEM-022 (my t2, support-2 tool),
THM-686(C) (Koksma per-pair), THM-661 (density-floor legs), `LRCResidualMeasureFloor.lean` (the single Lean obligation).

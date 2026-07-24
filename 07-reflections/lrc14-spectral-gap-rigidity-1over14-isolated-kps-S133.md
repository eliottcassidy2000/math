---
source: kind-pasteur-2026-07-23-S133 (Opus 4.8)
status: RIGIDITY RESULT (empirical, apparently NOVEL) + a clean conjecture that IMPLIES LRC(14). Pushing the
  "handle the near-extremizers" crux from kps-S132. Found a SPECTRAL GAP in the lonely-runner gap value:
  no 13-speed config has gap in the open interval (1/14, 1/13). Confirmed exhaustively on small support and
  against large near-dilate perturbations; not found in the literature. Corrects an over-claim (extremizer
  is NOT unique — infinite tight family). Reframes LRC(14) as spectral-gap + variational-bulk + tight residual.
tags: [lrc, lonely-runner, rigidity, spectral-gap, tight-instances, extremizer, variational, Gamma-k, program]
related: [kps-S132, THM-518, THM-515, "Lonely Runner turns 60" survey, arXiv 2211.08749]
---

# LRC(14): the gap value 1/14 is ISOLATED — a spectral gap up to 1/13 — and why it's the finish line

**kps-2026-07-23-S133.** Owner: push the rigidity theorem (the S132 crux — near-extremizers block the
variational route). I found the right form of it, and it is a clean, apparently-new structural statement.

## 1. Convention (literature-confirmed)
For `k` distinct positive-integer speeds, `Γ(k) = inf_S sup_τ min_v ‖vτ‖`; LRC says `Γ(k) = 1/(k+1)`.
Our case: `k=13`, `Γ(13)=1/14`, canonical extremiser `{1,…,13}` (gap `=1/14` at `τ=1/14`). (arXiv:2211.08749.)

## 2. The empirical SPECTRAL GAP (the rigidity result)
> **No 13-speed configuration has gap in the open interval `(1/14, 1/13)`.**
> Equivalently: `gap(S) > 1/14  ⟹  gap(S) ≥ 1/13 = Γ(12)`.

Evidence (exact rational gaps, `g=max_{a/q}min_v‖v·a/q‖`):
- **Exhaustive** over all 13-subsets of `{1,…,16}` and `{1,…,17}` (2940 configs): the ONLY gap `=1/14` is
  `{1,…,13}`; **zero** configs in `(1/14,1/13)`; the next value up is exactly `1/13` (7 configs).
- **Large near-dilates** — the natural attack — all snap to exactly `1/13`, never the band: `10·{1..13}` with
  the top speed `→131` or `→129`, `14·{1..13}+1`, `7·{1..13}·…→92`, `2·{1..13}→27/→25`, etc. → all `=1/13`.
- A hand-built "band" candidate `{2,…,14}` (which at `τ=1/27` gives `2/27≈0.0741 ∈` band) actually has a far
  better time `τ=1/16 ⟹ gap=1/8`. Every attempt to sit in the band is defeated by a better `τ`.

**Literature check:** the recent survey/papers discuss *tight instances* (gap `=1/14`) but I found **no
spectral-gap / value-isolation result** — this appears novel.

## 3. Correction: rigidity is VALUE-isolation, NOT config-uniqueness
My first-pass "gap `=1/14 ⟺` dilate of `{1,…,13}`" is **FALSE**: there is an *infinite family of tight
instances* (gap exactly `1/14`), e.g. with largest speed `2^{k-1}`-scale — they just have large speeds and so
missed the small-support scan. So the extremiser is NOT unique. The correct rigidity object is the **isolation
of the value `1/14`** (§2), which survives: all tight instances sit *exactly* at `1/14`, with a clean buffer
`1/13 − 1/14 = 1/182` before any non-tight config.

## 4. The Spectral-Gap Conjecture and why it IS the finish line
> **SGC(k):** the `k`-speed gap spectrum omits `(Γ(k), Γ(k−1)) = (1/(k+1), 1/k)`.
> I.e. `gap(S) < Γ(k−1) ⟹ gap(S) = Γ(k)`.

**SGC(13) ⟹ LRC(14):** any non-tight 13-config has `gap ≥ 1/13 ≥ 1/14`; tight configs have `gap = 1/14`.
So `gap ≥ 1/14` for all `S` — the conjecture. (SGC is strictly stronger: it also gives the isolation.)

**Why this is the productive decomposition (with S132's variational bound):**
- **Bulk (non-tight, `gap ≥ 1/13`):** the S132 variational bound `gap ≥ V(S)` with a Fejér concentrator +
  config-adaptive centring needs only **loss `< 1/182`** to certify `V(S) ≥ 1/14` here (since
  `gap − loss ≥ 1/13 − 1/182 = 1/14`). The `1/182` buffer is EXACTLY what the spectral gap provides —
  without it the variational route dies on near-extremisers (S132). *This is why isolation matters.*
- **Residual (tight instances, `gap = 1/14`):** an infinite but **well-structured** family (literature:
  lacunary, `2^{k}`-scaled, "accelerate specific speeds"). Handle by exact/algebraic classification, not the
  lossy bound. This is the THM-518 resonance world, now with a finite-type target.

So LRC(14) = [prove SGC's bulk half: non-tight ⟹ `gap ≥ 1/13`, via the variational bound] +
[classify & verify the tight-instance family]. The spectral gap is what glues them.

## 5. Why SGC might be TRUE — the "effective runners are quantised" heuristic
`Γ(k)=1/(k+1)` is a decreasing step function of `k`. SGC says the 13-speed spectrum jumps from `Γ(13)=1/14`
straight to `Γ(12)=1/13` — **as if a config has an integer number of "effective" runners.** Heuristic: at the
optimal `τ*`, if a speed is *not* binding (`‖vτ*‖ > gap`), the config behaves like `≤12` genuine constraints,
forcing `gap ≥ Γ(12)=1/13`; only when all 13 bind in the rigid tight pattern do you reach `1/14`. Making
"effective runner count is an integer" rigorous (a semicontinuity/reduction lemma for the gap spectrum) would
prove SGC — and is a cleaner target than LRC's raw inequality.

## 6. Honest status + next steps
- SGC(13) is **empirical** (exhaustive small support + large near-dilates); it must be tested at **large
  speeds** (esp. perturbations of the `2^{k}`-scale tight instances I couldn't construct exactly) before trust.
- **Attack the bulk half** `non-tight ⟹ gap ≥ 1/13` with the S132 variational engine at high Fejér degree +
  config-adaptive centring; chart where loss dips toward `1/182` — that locus = the tight-instance neighbourhood.
- **Classify the tight instances** (pull the exact `2^{k}`-family from the survey; connect to THM-518 resonances).
- **Prove/refute SGC** via the effective-runner reduction lemma (§5).

Files: `/tmp/{rigidity,search,spectral,neardilate,direct}.py`. Builds on kps-S132 (variational bound), THM-518.

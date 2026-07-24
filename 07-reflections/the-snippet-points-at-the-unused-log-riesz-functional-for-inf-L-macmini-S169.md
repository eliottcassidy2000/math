---
source: mac-mini-2026-07-23-S169 (Opus 4.8)
status: INSIGHT / LEAD (not a theorem). Using the external artanh snippet's STRUCTURE as a lens on the
  repo's #1 open problem (LRC(14): inf_S L(S) > 0). Headline: the repo's entire Riesz analysis never takes
  the LOG of the Riesz product; the snippet's functional (2 artanh = Σ_{k odd}t^k/k) is exactly that unused
  log/entropy functional. Two subclaims are machine-verified (freq-blindness of ∫log R; the ∫M log R
  per-mode = arctanh-family). The payoff claim (an entropy certificate breaks the AP-core stall) is an
  UNPROVEN lead with a concrete next experiment. Soundness of ∫M log R as a certificate is NOT established.
tags: [lrc, riesz, entropy, log-energy, artanh, snippet, inf-L, thm-518, bedert, lead]
related: [THM-515, THM-518, THM-2000, HYP-1992, HYP-1866]
---

# The snippet points at the one Riesz functional the repo never formed: log/entropy

**mac-mini-2026-07-23-S169.** Owner directive: stop hunting the snippet's source; use its STRUCTURE to find
insights for repo problems. Three parallel repo-mining passes (rapidity/adelic, LRC-Riesz, THM-2000) + the
snippet's confirmed form converge on one sharp, actionable gap.

## 0. The confirmed snippet form (for reference)
`(2457/6592)·log(8847357/2974400) − log(1285/896) > 1/25`, `2·artanh(t)=log((1+t)/(1−t))=Σ_{k odd}t^k/k`.

## 1. The gap (exhaustive search, Riesz mining pass)
The repo's LRC(14) Riesz analysis (THM-515/518, `lrc14_riesz_*` scripts) stays entirely at the **linear**
pairing `∫M·R/∫R` and the **multiplicative** product `R=∏_m(1+a_m cos 2πmτ)≥0`. It **never** forms
`∫log R`, the entropy `∫R log R`, `∫M log R`, or the expansion `log(1−p cos)`. And `2·artanh` /
`log((1+t)/(1−t))` / `Σ_{k odd}t^k/k` occur in the repo ONLY in tournament/OCF/hyperbolic threads
(`TANGENTS.md:2579`, `SESSION-LOG.md:6`), NEVER on the LRC covering side. The snippet is external evidence
that a certified-artanh **log-energy** bound wins a WIDER-GAP result (`>1/25`) — in the regime the repo's
linear route stalls (THM-518: ratio **1.0096 ≥ 1** on the `{1..13}\{j}∪{14m}` AP-cores).

## 2. Why the NAIVE log fails — and the fix (both machine-verified)
`(1/2π)∫log(1+a cosθ)dθ = log((1+√(1−a²))/2)`, so `∫log R = Σ_m g(a_m)` depends ONLY on amplitudes, **not
frequencies** — frequency-blind, useless as a looseness certificate. That is presumably why it was skipped.

The FIX is the frequency-sensitive weighting. `log(1+a cosθ)` has Fourier coefficients
`2(−1)^{k+1}/k · ρ^k`, `ρ=(1−√(1−a²))/a` — a **harmonic×geometric** series = the artanh family
(`2·artanh(ρ)=log((1+ρ)/(1−ρ))`, the snippet's function OF ρ). Hence
```
∫ M·log R  =  Σ_v Σ_k  1̂_danger(kv) · (2(−1)^{k+1}/k) ρ^k      — FREQUENCY-SENSITIVE per-speed
```
carries the speed set's additive structure. The repo's linear `∫M·R` shares only the first-order (`k=1`)
term; the log/entropy functional is the concave all-orders resummation whose per-mode weights ARE the
snippet's arctanh coefficients. (Verified k=1,2,3 and the freq-blindness numerically this session.)

## 3. A genuine unification (not a coincidence)
At `a=0.6 ⇒ ρ=1/3 ⇒ 2·artanh(1/3)=log 2`. That is the SAME `artanh(1/3)=½log2` that (i) is the THM-2000
hexagonal mass `M(6,2)=2log2`, and (ii) opus-S2 used in the `t=1/3` sandwich to certify `M(6,2)>M(4,3)`
float-free. So the tournament/Abel–Dini arctanh and the (unused) LRC log-Riesz arctanh are the **same
function** `2·artanh(ρ)`, `ρ=(1−√(1−a²))/a`. THM-252's rapidity `arctanh(λ)=½ln((m−k)/k)` and HYP-1992's
LRC rapidity `r_i=atanh(1−2g_i)` are the same object again — the repo already has this function on the
tournament side; the snippet says to carry it to the covering side.

## 4. Why it might break the AP-core stall (LEAD, unproven)
Bedert's Riesz GAIN scales as `dim₂²/n³` (additive dimension). For the `{1..13}` AP-cores `dim₂~log N≈2–3`,
so the gain `≈0.006` is worthless (THM-518 stall at 1.0096). But the snippet's optimizer weight numerator is
`2457 = 3·Σ_{k=1}^{13}k² = 3·819` (klein-S404, verified: uniquely n=13 among consecutive sets) — the **L²
energy** `Σv²`, which is LARGE (819) exactly on the AP-cores, NOT small like `dim₂`. So an
energy/entropy-weighted log functional is **not obstructed by the low-additive-dimension issue that kills
Bedert's Riesz gain** on precisely the extremizers that carry `inf L`. That is the mechanism by which the
snippet's log-energy structure could reach where the linear/product route provably cannot.

## 5. Honest caveats
- Soundness of `∫M log R < c ⟹ loose` is NOT established (log R changes sign; the `M≥1`-on-danger argument
  that makes `∫MR/∫R<1` sound does not transfer verbatim). This is a candidate functional, not a certificate.
- The snippet's source problem is unidentified; `n=13` is klein's circumstantial fingerprint, strong but
  not proof. The value of this reflection does NOT depend on the snippet being LRC — it is the internal
  observation "the repo never tried the log/entropy Riesz functional, and here is its exact arctanh form."
- THM-2000 refinement (mining pass): the snippet is NOT a THM-2000 mass ordering — those reduce to single
  `log 2` inequalities (only `M(6,2)>M(4,3)` is cross-axis & certifiable, already done). The snippet has two
  INDEPENDENT logs. It uses THM-2000 §3.1's Abel–Dini telescoping machinery, nothing more.

## 6. Concrete next experiment
Compute a concave/relative-entropy Riesz functional (candidates: `∫M log R`, `∫R log R`, or the
KL-divergence of `R/∫R` from the safe indicator) on the `{1..13}\{j}∪{14m}` cores where the linear ratio
stalls at 1.0096, vs. loose sets, using the existing `R̂` enumeration (`lrc14_riesz_optimize_macmini_0615s5.py`).
Question: does the entropy functional separate loose from tight where the linear ratio cannot? If yes, build
the artanh-certified (Lean-portable, via opus-S2's `log_lower/log_upper`) version and target `>1/(2n−1)=1/25`.

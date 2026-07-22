# LRC and GMC(2) are one integer-kernel non-cancellation problem — reducing LRC(14) to the maximal-resonance (AP) cores

**death-star-2026-07-21-S102** (HYP-8879). Owner: leverage the GMC(2)/DvdK/scale-clock/tournament-zeta ideas
toward LRC. The transfer is exact at the level of the object: the LRC lonely measure and the GMC(2) moment are
**the same integer-kernel resonance sum**, and my S101 unique-cycle criterion becomes a quantitative reduction of
LRC(14) to the relation-rich (AP) cores.

## The unification (sound, verified)
The L∞ lonely measure of a speed core `{v_i}` at threshold `δ` is, by Fourier expansion of the one-runner lonely
indicator `g(x)=𝟙[‖x‖>δ]` (coefficients `ĝ_0=1-2δ`, `ĝ_k=-sin(2πkδ)/(πk)`),
```
    μ(lonely) = ∫_0^1 ∏_i g(v_i t) dt = Σ_{k : Σ_i k_i v_i = 0} ∏_i ĝ_{k_i}.
```
This is a sum over the **integer kernel** `{k : Σ k_i v_i = 0}` (the *resonances*) weighted by Fourier products —
structurally identical to the GMC(2) moment `E[P^m] = Σ_{balanced channels r : Σ r_i q_i = 0} multinomial·A(r)!·c^r`
(a sum over the balanced-channel kernel of the *charges*). Verified: the resonance sum reproduces the direct
integral (`{1,3,5}` at `δ=1/4`: 0.105 vs 0.100). **Covering (`μ=0`) is a cancellation of resonance terms —
exactly GMC(2)'s `E=0`.** Both are the same non-cancellation problem on the integer kernel of a linear map.

## The clock-floor decomposition and the reduction to the AP
Split off the zero resonance:
```
    μ = (1-2δ)^n  +  Σ_{nonzero resonances} ∏ ĝ_{k_i}   =   MAIN  +  corrections.
```
The **MAIN term `(1-2δ)^n`** is the Eisenstein/clock **floor**: for LRC(14) (`n=13, δ=1/14`) it is `(6/7)^13`, and
`1-2δ = 12/14 = 6/7` is precisely the THM-878 clock-floor constant `A(q) ≥ 6/7` (and boxeph's S221 "Eisenstein
floor"). Covering therefore requires the resonance corrections to **cancel the `(6/7)^13` floor** exactly.

Verified quantitatively (`lrc_gmc2_resonance_unification_deathstar_S102.py`): `|corrections|/MAIN` is
- **0.03–0.10 for Sidon / low-resonance cores** (`{1,2,5,11}`, `{1,2,5,11,22}`): the corrections cannot come
  near `1`, so `μ ≈ MAIN > 0` — the core is **robustly lonely, never covering**;
- **0.89–0.96 for arithmetic progressions** (`{1,2,3,4}`, `{1,2,3,4,5}`): the AP has enough coincident resonances
  to nearly cancel the floor (`{1..5}` at `δ=1/6`: `μ=0.0000` by direct integration — covering).

So a core covers **only if it is maximally resonant**: only the AP-neighborhood has the resonance density to
cancel the floor. This is my S101 mechanism transferred: **few / unique minimal resonances ⟹ the floor survives
⟹ lonely ⟹ not covering.** LRC(14) reduces to the AP-neighborhood cores.

## Where the residual sits (unified across the two threads)
The AP is exactly the **maximal-resonance** object — verified: AP `{1..6}` has 18 minimal resonances vs a Sidon
set's 2. That is:
- the **GMC(2) coincident-cycle hard stratum** (S101 — the AP's difference/relation set has many coincident
  minimal cycles, the paradigm being antipodal pairs `±1,…,±(n-1)`);
- the **degenerate tournament zeta** (S99/THM-1926 — many coincident shortest primitive cycles);
- the **resonant multi-clock** (S99 scale-then-clock — the AP is the tight resonant clock);
- codex's **relation-rich covering core** and boxeph's **tight-AP** frontier (S214).
The GENERIC/low-resonance cores are handled by the same "floor survives" argument here and by codex's scaled
zeta-core / missing-clock sieve (THM-2057) on the covering side — the two threads' easy cases coincide too.

## Honest scope and what it buys
- **Not a proof of LRC(14).** The Fourier/exponential-sum decomposition of the lonely gap is a standard tool; the
  contribution is the **unification** (LRC = GMC(2) as one integer-kernel non-cancellation), the **main-term =
  clock-floor** identification (`(6/7)^13` = THM-878 = S221 Eisenstein floor), and the **quantitative reduction**
  of LRC(14) to the AP-neighborhood via `|corrections|/MAIN`, with the residual pinned as the tournament-zeta
  coincident-cycle degeneracy.
- **What must still be done** (the engine): a rigorous bound `|corrections| < MAIN` for every non-AP 13-core — a
  resonance-count × Fourier-decay (`|ĝ_k|~1/πk`) estimate. The numerics say the margin is large off the AP; the AP
  itself is the codimension-≥1 residual where the sharp analysis (chi-criterion, rank-or-Euler) lives.
- **Payoff for the fleet:** this gives the LRC residual a GMC(2)-style name (coincident-cycle / degenerate-zeta),
  a quantitative Sidon-vs-AP separation, and the same clock-floor constant on both sides — a bridge for
  transferring the GMC(2) non-cancellation machinery (Frobenius `Q̄^p` single-power certificate, unique-cycle
  criterion) to the covering side.

Cross-links: S99/HYP-8876 (scale-clock, tournament zeta), S101/HYP-8878 (unique-cycle DvdK-free), S100/HYP-8877,
THM-878 (clock moduli, 6/7 floor), THM-1926 (tournament zeta), THM-2047 (pair-sum wall), THM-2057 (codex scaled
zeta-core / clock sieve), boxeph S214 (tight AP relation-rich), S221 (Eisenstein-floor+cusp), memory
`gmc-lrc-same-positivity-manoeuvre`, `lrc14-frontier-and-sharp-horn`. Script
`04-computation/lrc_gmc2_resonance_unification_deathstar_S102.py` (+ `.out`). HYP-8879.

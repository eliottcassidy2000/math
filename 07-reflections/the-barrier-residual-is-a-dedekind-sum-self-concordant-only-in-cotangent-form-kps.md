# The barrier residual is a Dedekind sum — self-concordant only in its cotangent (log) form, and reciprocity is the lever

*kind-pasteur-2026-06-30. A deep study of the far-element residual `Δ_w` I had loosely called a "self-concordant barrier residual." Three questions settled rigorously (`lrc14_barrier_residual_dedekind_selfconcordance_kps.py`), plus a survey of the project's "finite periodic objects." The residual is a Dedekind sum in three interchangeable disguises; "self-concordant" is literal only in one of them; and that one comes with reciprocity.*

> **Concurrent convergence (added at close-out, 2026-06-30).** The "concrete next step" below — period-max via Dedekind `O(log)` reciprocity descent — was executed *first and more completely* by opus **HYP-3770**: the covering-min margin `= −12·s(n,Φ₆)/n²` descends by reciprocity in `O(log)` (2 steps, since `Φ₆≡1 mod n`), verified to `n=10³⁰` and `lcm(1..80)~10³⁴`, with `s(n,Φ₆) → −1/12 = ζ(−1)`. Its **honest negative** (only the *construction* rung is a single Dedekind sum; the extremal covering-min rung `a(n)` is not) is the same "reciprocity computes sums, not the max" limit I hit. **My complement:** the hypercontractivity/variance route (below) bounds the *sup over the whole period* via sub-Gaussian moments — exactly the extremal-rung gap opus's per-rung reciprocity leaves open. klein/mac-mini **HYP-3768** proved the `B₂/E₂` margin = Dedekind sum. I **defer to them** on the covering-min Dedekind/`−1/12`; what is mine here is the sub-Gaussian **sup** bound and the sector-`p0` reading.

## What the residual actually is

THM-563 pins the single-far residual exactly: `Δ_w·w = Σ_j Σ_t ±S_j(frac(w·t))`, with `S_j` the centered **sawtooth antiderivative** — a generalized Dedekind sum, periodic with period `7·lcm(B)`, `|S_j| ≤ 3/49`. The sawtooth `((x))` is the first Bernoulli `B₁`; its antiderivative is `B₂`. So the residual is a **`B₂` / weight-2 object**, and the survey of the project's finite-periodic objects confirms this is the *one recurring object* everywhere:

- the cycle-graph resistance `G_P = 2t(P−t) = −2P²B₂(t/P)+P²/3` (the L7 discrepancy, HYP-2745),
- `ζ(−1) = −1/12 = B₂/2` (the regularized killer-sum, the infinite-AP limit),
- the doublet `= 2ζ(3)` Mordell–Tornheim (a depth-2 rung, THM-578),
- the Fejér kernel `F₇` (the level-7 `E₂` boundary),

are all shadows of the second Bernoulli / non-holomorphic Eisenstein `E₂`. The far-element barrier residual is the **single-far, degree-1 Fourier shadow** of that object — i.e. a **Dedekind sum**.

## The three disguises are one number

Rigorously verified (agreement to `10⁻¹⁶`):

> `s(h,k) = Σ_i ((i/k))((hi/k))` (Rademacher/sawtooth, `B₂`, L1) `= (1/4k) Σ_i cot(πi/k)cot(πhi/k)` (cotangent).

The **sector `Δ_w` is the sawtooth form**; the **cotangent form is what a log-barrier produces** (because `d/dx[−log(2 sin πx)] = −π cot πx`, so summing the log-barrier gradient over the orbit `{hi/k}` *is* the cotangent Dedekind sum). Same signed number, two barriers. This is the number-theoretic content of "attacker = defender at the certificate level": the signed cancellation the conjecture needs is a Dedekind sum, and the cycle-resistance reflection already flagged it as *the signed degree-1 shadow of the absolute `B₂`/L1 discrepancy*.

## Self-concordance: literal only in the cotangent (log) form

The correction to my earlier reflection. Self-concordance (`|F‴| ≤ 2(F″)^{3/2}`) needs `C³` smoothness and `F″>0`. Tested:

| barrier | `F″` | ratio `|F‴|/(F″)^{3/2}` | self-concordant? |
|---|---|---|---|
| `−log x` | `1/x²` | **2.0000** | yes (canonical) |
| `−log sin πx` (arc) | `π²/sin²` | `2|cos πx| ≤ 2` | yes |
| tent `|x−½|` (sawtooth) | **0** | ∞ | **no** |

So the **sector/tent residual is a Rademacher `B₂` sawtooth — not self-concordant** (piecewise-linear, `F″=0`, `F‴` a delta at the kink). The interior-point / central-path picture is a *literal* interior-point method **only if the functional is recast with the log barrier** `−log(dist to danger)` (equivalently the arc barrier `−log sin`). Then, and only then, do the interior-point guarantees (smooth central path, Newton) actually apply.

## What the log form buys: reciprocity = a Euclidean descent

The reason to pay the recasting cost is **Dedekind reciprocity**, which the sawtooth-period-max route cannot use:

> `s(h,k) + s(k,h) = −¼ + (1/12)(h/k + k/h + 1/hk)`  ⟹  `s(h,k) = R(h,k) − s(k mod h, h)`,

a **Euclidean/continued-fraction descent** computing `s(h,k)` in `O(log k)` steps (verified: `s(89,233)` in 10 steps, exact). THM-563's period-max is a per-base maximum over one period of length `7·lcm(B)` — and `lcm(B)` can be `10⁵–10⁷`, so the **general-bounded-base closure it leaves "in progress" is an enumeration over enormous periods**. Reciprocity replaces that enumeration with an `O(log)` descent per endpoint sum — making the all-`B⊆[0,14]` finite check *tractable*, and giving the bound a **structural continued-fraction form** rather than a brute period-max.

## The honest caveat, and why it's the right caveat

The classical Dedekind sum is **not uniformly bounded**: `s(1,k) = (k−1)(k−2)/(12k) ~ k/12`. So reciprocity does not hand you a small constant for free — the bound is **Diophantine**, `|s(h,k)|` controlled by the *partial quotients* of `h/k` (small ⟺ badly-approximable `h/k`; large ⟺ `h/k ≈ a/small-denominator`, the AP-resonant case). (The LRC residual itself stays bounded because `Δ_w·w` is a *finite* sum of `|S_j|≤3/49` sawtooth terms, one per arc endpoint — THM-563's period-max `≈1` — so the `k/12` growth is a property of a *single full-orbit* Dedekind sum, not of `Δ_w`; reciprocity is the tool for the *sub-sums*.)

This caveat is exactly right, because it is the **same continued-fraction object as the ζ(2) floor**: the floor came from Farey/coprime density summed by Mertens (`the-resonance-killing-game`, `zeta2-governs-the-lonely-runner-floor`), and the residual bound comes from the continued fraction of `h/k`. **The far-element barrier residual (Dedekind reciprocity) and the resonance floor (ζ(2) Farey density) are two ends of one continued-fraction / three-distance structure** — the residual is the Dedekind sum, the floor is the Farey sum, and the Euclidean algorithm is the bridge. The AP-resonant hard case (`h/k` near a low-denominator rational, large partial quotient) is precisely the resonance-killing extremal (the AP), unifying the two.

## Net

- **Pinned:** the far-element residual is a **Dedekind sum** = the single-far shadow of `B₂`/`E₂`; sawtooth (sector) and cotangent (log) forms are the same signed number.
- **Corrected:** "self-concordant" is literal only for the **log-barrier** form; the sector/tent form is not self-concordant. Interior-point becomes literal only after recasting with the log barrier.
- **New lever:** the log/cotangent form inherits **Dedekind reciprocity** — an `O(log)` Euclidean descent that makes the general-bounded-base closure tractable and gives the bound a continued-fraction form.
- **Unification:** the residual (Dedekind/reciprocity) and the floor (ζ(2)/Farey) are two ends of one continued-fraction structure; the AP-resonant hard case is shared.
- **Concrete next step — DONE, honestly (`lrc14_periodmax_via_variance_pairwise_kps.py`).** Reciprocity computes *sums*, and the period-max is a genuine *max*, so it is **not** an `O(log)` descent. But the honest payoff is real: the **variance** `E_w[(w·Δ_w)²]` **decouples into pairwise tent-correlations**, each averaged over `lcm(den θ_e, den θ_{e'}) ≤ (7·max B)²` — *never* the full `7·lcm(B)` period. Verified exact (pairwise `E[g²]` == full-enumeration `E[g²]`), and the cost is bounded by **max(B), not lcm(B)** — computed for `B={0..7,9,11,13}` (full period 1,261,260, infeasible to max-enumerate) in 33 s. And the period-max is a **faithful multiple of the rms**: `period-max / rms ∈ [3.48, 4.90]` across six bases (roughly constant, `~sqrt(2 log m)` for `m = #endpoints`, the sub-Gaussian heuristic). So: `period-max ≲ C·rms`, `C ≈ 5` empirically, and `rms` is enumeration-free. To make it rigorous, prove the sub-Gaussian sup bound `period-max ≤ C√(log P)·rms` via **hypercontractivity**. **Pushed (`lrc14_hypercontractivity_subgaussian_kps.py`):** `g(w)` IS sub-Gaussian — normalized moments `‖g‖_{2k}/‖g‖_2` sit *at* the Gaussian values `[(2k-1)!!]^{1/2k}` (`1.33/1.60/1.82/2.00` vs `1.32/1.57/1.79/1.98`), kurtosis `→3` as bases grow. The sup obeys the extreme-value law `max/rms ≈ 0.92·√(2 ln P)`, and since `P ≤ 7·lcm(1..14) = 2,522,520` for every `B⊆[0,14]`, this gives a **uniform `period-max ≤ 5.43·rms`**. So the general-bounded-base closure becomes a **per-base variance computation of cost `O(max B)`** (not `O(lcm B)`) times a fixed sub-Gaussian factor. Rigorous skeleton: `E[g^{2k}] = Σ_{resonances}∏ĝ₀(nᵢ)`, `ĝ₀(n)~1/n²` (tent decay) ⟹ convergent, and the non-pairing resonances are bounded by the **additive energy of the endpoint frequencies `{θ_e}={j/(7e):e∈B}`** (bounded for bounded `B`); the pairings give the Gaussian `(2k-1)!!(E[g²])^k`. The exact max, if ever needed, is an Ostrowski/three-distance continued-fraction sup — the same CF object as the ζ(2) floor. **And the `1/12` is not incidental:** `∫((x))² = 1/12 = B₂/2 = −ζ(−1)` is the sawtooth variance — the rms denominator of the whole hypercontractive bound lives in this Bernoulli seed.

— Related: THM-563 (single-far periodicity), THM-578 (doublet/Mordell–Tornheim), HYP-2745 (`G_P = B₂`/cycle resistance), `the-lrc-discrepancy-is-a-cycle-resistance.md`, `zeta2-governs-the-lonely-runner-floor.md`, `the-resonance-killing-game-and-the-zeta-duality-of-the-lonely-runner.md`, companion `the-far-runner-is-the-log-barrier-lrc14-as-an-interior-point-method-kps.md`, OPEN-Q-108. Artifacts: `04-computation/lrc14_barrier_residual_dedekind_selfconcordance_kps.py`, `05-knowledge/results/lrc14_barrier_residual_dedekind_selfconcordance_kps.out`.

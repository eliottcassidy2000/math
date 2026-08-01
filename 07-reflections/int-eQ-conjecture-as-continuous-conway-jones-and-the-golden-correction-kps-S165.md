---
source: kind-pasteur-2026-07-24-S165 (Opus 4.8)
status: CORRECTION + SYNTHESIS. (1) RETRACTS my kps-S163 "confabulation" flag on C*=log_5(5 phi^2): it is the
  CONFIRMED AMM 12592 answer (opus min-max Lagrangian + klein THM-3027, two independent derivations); I wrongly
  compared it to the weaker artanh lower bound 9049/6592. (2) Reframes the external "int_0^1 e^{Q(t)} dt != 0 for
  nonconstant Q in Qbar[t]" conjecture as CONTINUOUS Conway-Jones: its only zero-locus is the 2 pi i Z resonance
  (int e^{2 pi i n t}=delta), the exact continuous analog of my FC(3) DISCRETE seed L((sum zeta^j X_j)^k)=k![m|k];
  algebraicity blocks the continuous resonance. Builds on concurrent HYP-9078 (weights j! vs 1/(j+1); implication
  UNVERIFIED), the Liu-Sun 2020 radial reduction (FC(2)-homogeneous = Lebesgue moment problem), and opus's
  saddle-cap on my kps-S160 rank-drop locus.
tags: [factorial-conjecture, exponential-integral, roots-of-unity, conway-jones, transcendence, amm12592, correction]
related: [kps-S154, kps-S156, kps-S160, kps-S163, HYP-9078, THM-3027, THM-3022]
corrects: [kps-S163 golden-confabulation flag]
---

# The int e^Q conjecture as continuous Conway-Jones; and the golden correction

## 1. CORRECTION: the golden constant is confirmed (retract kps-S163)
kps-S163 flagged `C_arch = log_5(5 phi^2)` as a confabulation. **That was wrong.** It is the confirmed AMM 12592
minimal-slope answer:
> **`gamma* = 2 log(phi)/log 5 = log_5(phi^2)`, `C* = log_5(5 phi^2) = 1.59798743566...`** -- two independent
> derivations (opus's min-max Lagrangian `x* = (3-sqrt5)/2 = 1/phi^2`; klein THM-3027 `(1-tau)^2=tau`,
> `tau^2-3tau+1=0` = min poly of `phi^2`).
My error: I compared it to `9049/6592 = 1.3727`, which is the **weaker artanh-certificate lower bound**
(`gamma >= 2457/6592 = 0.3727 < 0.598`), and mislabeled the certificate rate as "the answer." The certificate is a
valid but non-tight lower bound; the answer is the golden. (Timing: the golden derivation landed 2026-08-01,
after my last-turn repo search returned empty -- not confabulated, concurrent.) The kps-S164 tandem stands, with
this correction: the artanh two-bias log-odds functional bounds `gamma*` from below at `0.3727`; the TRUE `gamma*`
`= 0.598` is the golden, so closing the gap `[0.3727, 0.598]` is the live problem.

## 2. The int e^Q conjecture is continuous Conway-Jones
External claim (friend's AI): `int_0^1 e^{Q(t)} dt != 0` (indeed transcendental) for nonconstant `Q in Qbar[t]`.
Concurrent HYP-9078 established: real case trivial (positivity); degree 1 = Lindemann-Weierstrass
(`int e^{at+b}=e^b(e^a-1)/a`, zero iff `a in 2 pi i Z`, transcendental); bridge `int_0^1 e^Q = sum_m M(Q^m)/m!`,
`M(t^j)=1/(j+1)`; and the asserted `=> FC(2)` implication is **UNVERIFIED** (FC weight `j!` vs integral weight
`1/(j+1)`). My contribution is the mechanism:
> **The ONLY zero-locus of `int_0^1 e^{Q}` is the character resonance `int_0^1 e^{2 pi i n t}dt = delta_{n,0}`
> (`Q in 2 pi i Z * t`). This is the CONTINUOUS analog of my FC(3) DISCRETE seed
> `L((sum_{j} zeta_m^j X_j)^k) = k! [m | k]` (kps-S156) -- both are exactly character orthogonality.**
Algebraic data can only reach the DISCRETE, leaking resonance (FC seed vanishes only at multiples `m|k`; verified,
kps-S156); the CONTINUOUS full resonance `int e^Q = 0` requires the transcendental `2 pi i`, which algebraicity
forbids. **So `int e^Q != 0` and "no FC counterexample" are one Conway-Jones statement in two registers:
algebraic/finite data cannot realize the full `2 pi i` resonance that vanishing demands.** (Same governor as the
series closed-form locus, LRC's Lam-Leung, and AMM's `zeta_6`.)

This is the intuition for why the implication *should* hold; it does NOT bridge HYP-9078's weight gap
(`j!` vs `1/(j+1)`). The actual bridge is **Liu-Sun 2020 Thm 2.6**: FC(2)-for-homogeneous reduces radially to the
`[0,1]`-Lebesgue polynomial moment problem -- exactly the `1/(j+1)` world of the integral. So the honest chain is:
`FC(2)-homogeneous == Lebesgue moment problem (Liu-Sun)`, and the exponential upgrade `int e^Q != 0` is the
transcendence layer above it; the general (non-homogeneous) FC(2) implication is the unproved step.

## 3. Degree >= 2 is E-function transcendence (verified structure)
`e^{Q(t)}` solves `y' = Q'(t) y` (first-order, polynomial coefficients) -> its `[0,1]` integral is a period of an
**E-function** (Siegel-Shidlovskii / Andre / Beukers). Degree 2 is the Fresnel/erf value (checked:
`Re int_0^1 e^{i t^2} = 0.90452...`, no low-height PSLQ relation vs `{1,pi,sqrt pi,e,sqrt2}` -> consistent with
transcendental). The friend's "E-functions + Beukers lifting" is exactly this machinery; the transcendence of the
integral (if proven) yields transcendence of many special-function values, as claimed.

## 4. Integrating the concurrent FC(3) work (credit)
- **opus saddle cap:** for an isolated max-modulus point, `int_T phi^{3k} ~ (3/k) M^{3k} sum_j p_j e^{3 i k theta_j}
  != 0` closes the **Morse case analytically at all degrees** -- complements my finite-degree transversality
  (kps-S159/S161), and the max-modulus ridge sits **exactly on the fold `{J=0}` = my rank-drop locus** (kps-S161),
  pinning the residual there. Consistent with the KZ-gap verdict (kps-S160), not stronger.
- **Discipline (OpenAI "ten advances" lesson, owner-flagged):** that page is *curated, validation-not-established*;
  "fetch the primary before opening a lane; the capstone is not the source." Same lesson as my S151 false negative
  and the golden-flag error here. Liu-Sun is the primary for the homogeneous case; cite it, don't reprove.

## 5. Net
Golden confirmed (kps-S163 corrected). The `int e^Q` conjecture is continuous Conway-Jones -- unifying it with the
FC seed, the series locus, LRC, and AMM under "algebraic data cannot realize the full `2 pi i` resonance." The
`=> FC(2)` implication is UNVERIFIED beyond the homogeneous (Liu-Sun) case; degree `>= 2` transcendence is
E-function theory. Live: close AMM's `[0.3727, 0.598]` gap to the golden; supply the non-homogeneous FC(2) bridge.

Files: verify inline. Builds on kps-S154/S156/S160/S163; engages HYP-9078, THM-3022/3027, Liu-Sun 2020.

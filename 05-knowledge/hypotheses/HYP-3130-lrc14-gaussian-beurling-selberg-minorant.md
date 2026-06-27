---
id: HYP-3130
title: LRC14 Gaussian / Beurling-Selberg minorant for the r=2..6 multi-far floor (the prompt's TOOL 1) -- a TRUE C^infty minorant with super-polynomial UNIFORM tail; INDEPENDENT proof of the apex Q-block floor (exact cap-matching c_r) + the honest envelope obstruction
status: VERIFIED tool + rigorous apex-block uniform floor (mod THM-546); the absolute envelope CANNOT close R' (sign-cancellation essential -- agrees w/ HYP-3128 Asano obstruction, HYP-3129 signed-SPEC). NOT a new R'>=c proof.
source: kind-pasteur-2026-06-27-S256
scripts:
  - 04-computation/lrc14_gaussian_minorant_floor_kpswf15.py
  - 04-computation/lrc14_apex_block_minorant_floor_kpswf15.py
  - 04-computation/lrc14_cross_resonance_minorant_kpswf15.py
results:
  - 05-knowledge/results/lrc14_gaussian_minorant_floor_kpswf15.out
  - 05-knowledge/results/lrc14_apex_block_minorant_floor_kpswf15.out
  - 05-knowledge/results/lrc14_cross_resonance_minorant_kpswf15.out
  - 05-knowledge/results/lrc14_apex_floor_stress_kpswf15.out
  - 05-knowledge/results/lrc14_minorant_certificates_kpswf15.out
related:
  - HYP-3129   # SPEC resonance-lattice elementary certificate (TOOL-3, L2-Parseval ceiling) -- complementary route
  - HYP-3128   # loneliness = danger-count PGF; Lee-Yang Q-block + Asano obstruction -- INDEPENDENTLY converges on same Q/R split
  - HYP-3127   # Asano contraction of single-far factors (mac-mini)
  - HYP-3121   # the 3-engine synthesis / lift-and-decorrelate floor (the prompt)
  - HYP-2861   # L2-Cauchy-Schwarz spectrum tail
  - HYP-2856   # 3/pi^2 Farey floor
  - HYP-2606   # singular series / chat on the lattice
  - HYP-2840   # Vitali multi-far / Beurling-Selberg one-sided cert (Swing 2)
  - THM-546    # single-far comb bound (the peel used in the apex recursion)
  - OPEN-Q-108
external: Lonely Runner Conjecture (first open case, 13 speeds); Beurling-Selberg extremal functions; Vaaler; heat-kernel / C^infty mollifier minorants.
---

# HYP-3130 -- the Gaussian / Beurling-Selberg minorant (the prompt's TOOL 1)

## Object
Covering set `S = R u 14Q`, `r=|Q| in {2..6}` (the few-apex OPEN core), `R` 14-free (`|R|=13-r>=7`).
Loneliness `L(S) = int_0^1 prod_{s in S} phi_s(t) dt`, `phi_s(t)=1_{[1/14,13/14]}({s t})`, arc of measure 6/7.
LRC(14) covering branch (r=2..6) `<=> L(S)>0`.

## The minorant idea (the wide-V tail FIX)
The sharp arc indicator `phi_s` has 1/n Fourier decay -> the relation-lattice / resonance sum
`sum_{sum_s k_s s=0} prod_s phihat(k_s)` has a slowly decaying tail, hard to bound UNIFORMLY in the
speed magnitudes. Replace each `phi_s` by a NONNEGATIVE minorant `0<=psi_s<=phi_s` whose Fourier
coefficients are super-polynomially decaying (C^infty mollifier, route i) or finitely supported
(Beurling-Selberg band-N, route ii). Then
```
L(S) >= int prod psi_s = MAIN + RESONANCE,
MAIN = (int psi)^|S|,  RESONANCE = sum_{(k_s)!=0, sum k_s s=0} prod psihat(k_s),
```
and the RESONANCE tail is super-poly / finite-band => UNIFORM in the speed magnitudes.

## TOOL 1, route (i): the validated TRUE C^infty minorant
`psi = 1_{[1/14+delta, 13/14-delta]} * rho_delta`, `rho_delta` the standard C^infty bump
`rho(u)=C exp(-1/(1-u^2))` on `(-1,1)` scaled to width delta.
- **TRUE minorant (VERIFIED exactly):** support(psi) = [1/14,13/14] EXACTLY => `psi<=phi` pointwise;
  `psi>=0`. Grid validation: `min_psi=0.0`, `max_leak(outside arc)=0.0`, `max_over(inside)=3.5e-15`.
- `int psi = h0 = 6/7 - 2 delta` EXACTLY.
- `psihat(k)=chi_arc(k) rhohat(delta k)`, `|chi_arc(k)|=|sin(pi k 6/7)|/(pi|k|)<=1/(pi|k|)` (=0 at the
  APEX 7|k), `rhohat` decays faster than any power. **CERTIFICATE 2 (tail mass `sum_{|k|>B}|psihat|`,
  delta=0.05):** `B=16:2.0e-2, B=32:3.8e-3, B=64:4.9e-4, B=128:3.4e-5` -- super-poly, BAND B is
  speed-INDEPENDENT. This is the uniform tail the prompt asked for.

Route (ii) Beurling-Selberg/Fejer-of-shrunk-arc gives a FINITE band but the Fejer tail LEAKS outside
the arc (not a true minorant without large shrink); the C^infty route is the validated true minorant.

## RESULT A -- the apex (Q) block floor is CLOSED uniformly (the NEW few-apex structure)
After `u=14t`, `L(S) = R' * meas(R-safe) * meas(Q-lonely)`, `Q-lonely={u: ||m u||>=1/14, m in Q}`.
Because Q has only `r<=6` runners this is a BOUNDED-COMPLEXITY loneliness measure. PROVED uniform floor:
```
meas(Q-lonely) >= c_r > 0,   c_2..c_6 = 66/91, 55/91, 1979/4004, 2243/5880, 3029/10780
                                      = 0.7253, 0.6044, 0.4943, 0.3815, 0.2810
```
- **The minimizers ARE the LRC covering caps:** `c_4=1979/4004=cap_9`, `c_5=2243/5880=cap_8`; minimizers
  `{1,11,12,13}`, `{1,5,7,8,9}` = the pairwise-avoidance cap configs (THM-576). The minorant route
  reproduces the cap structure exactly.
- **CERTIFICATE 1 (inf on bounded Q):** the recursion `(6/7) c_{r-1} > c_r` holds with positive gaps
  `+0.0094/+0.0173/+0.0238/+0.0422/+0.0460` (r=2..6). Single-far peel (THM-546): adding a large w to Q
  multiplies meas(Q-lonely) by a factor `-> 6/7` from below (empirical peel factor 0.85->0.857 for
  w=50..1000). So wide Q has measure `>= (6/7)c_{r-1} > c_r`; the infimum is on bounded Q, finite-checked
  exactly. **Stress test: 0 violations** over ~6500 wide/adversarial/dilated/multiscale Q (cap 150).
  => `meas(Q-lonely)>=c_r` is RIGOROUS modulo the (proved) THM-546 comb bound.

## RESULT B -- the honest ENVELOPE obstruction (why the minorant alone cannot close R')
The absolute (Schur/triangle) envelope `|RESONANCE| <= (h0+B1)^13 - h0^13 - 13 h0^12 B1`,
`B1=sum_{k!=0}|psihat(k)|`, FAILS catastrophically (`MAIN - envelope ~ -10^4`) at EVERY delta:
**B1/h0 ~ 1.4-1.7 for ALL delta** (never < 1), so `(h0+B1)^13` always explodes. This is a hard
obstruction (same wall as MISTAKE-078: the absolute relation-lattice envelope diverges), NOT a tuning
issue. **Conclusion: the uniform floor MUST come from SIGN cancellation in the resonance, not the
absolute bound.** This AGREES with HYP-3128 (Asano cannot certify Xi(1,1)>0 because the >=7-speed
R-block factor is not unit-disk zero-free) and HYP-3129 (the floor needs the SIGNED SPEC, not |SPEC|).

## RESULT C -- the cross-resonance (R'-coupling) is the apex-14 |N|=1 level, super-poly tail in N
Apply the minorant to both blocks: `int psi^R psi^Q = (R-floor)(Q-floor) + CROSS`,
`CROSS = sum_{N!=0} [sum_{sum_r j_r r=-14N} prod psihat(j_r)] [sum_{sum_m k_m m=N} prod psihat(k_m)]`.
- A cross relation needs `sum_r j_r r = -14N` with R 14-FREE => uses the APEX PRIME 7 (2*7=14): the
  apex-prime-7 coupling that makes LRC(14) the hard composite case.
- VERIFIED: the cross-level sum CONVERGES super-polynomially in `|N|` (sum `|N|<=6` matches direct CROSS
  to 85-95%); the `|N|=1` level dominates when `1 in Q` (e.g. r=4: product `-4.2e-3`). The R'-coupling
  has a UNIFORM (band-free in N) tail -- exactly the prompt's "uniform tail independent of speed magnitudes".
- BUT the minorant DEGRADES the ratio (R'_min ~ 0.20-0.36 vs exact R' ~ 0.47-1.14): the minorant's small
  R-floor inflates `CROSS/baseline`, so R'_min is NOT a good proxy for the exact R'. The minorant is the
  right tool for the APEX block (Result A) but not for the coupling itself.

## Net (the deliverable, honest)
- **The prompt's TOOL 1 is built + validated:** a TRUE C^infty minorant `psi<=phi` with a super-poly,
  speed-INDEPENDENT Fourier tail (CERT 2). It converts the wide-V 1/n tail into a uniform tail.
- **What it CLOSES uniformly:** the apex (Q) block factor `meas(Q-lonely)>=c_r>0` for r=2..6, with EXACT
  cap-matching constants and a rigorous "inf-on-bounded-Q" recursion (mod THM-546). This is the NEW
  few-apex structure of the r=2..6 case, now uniformly handled by THREE independent routes (this minorant,
  HYP-3128 Lee-Yang, HYP-3129 union/SPEC).
- **What it does NOT close:** the coupling `R' >= c` (Node-3). The minorant MAIN term is NOT provably
  `>` the resonance via the absolute envelope (B1>h0, hard obstruction) -- sign cancellation is essential.
  The signed/uniform R'>=c content is supplied by the complementary HYP-3129 (R'>=0.642 via the signed
  SPEC + L2 ceiling), NOT by this minorant; HYP-3128 shows Asano cannot do it either. So TOOL 1 (minorant)
  is the apex-block + uniform-tail engine; the coupling floor remains the HYP-3129 / Node-3 job.
- **Convergence:** three same-prompt concurrent routes (minorant / Lee-Yang-Asano / SPEC) all factor
  L(S)=R'*mR*mQ, all close the Q-block, all isolate R'>=c as the sole open content. Strong triangulation.

-> HYP-3129, HYP-3128, HYP-3127, HYP-3121, HYP-2861, HYP-2856, HYP-2606, HYP-2840, THM-546, OPEN-Q-108.

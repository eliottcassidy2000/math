---
id: HYP-3794
title: THE OPEN CORE IS FINITE -- the Malikiosis-Santos-Schymura velocity bound (unconditional, 2024) caps a LRC14 counterexample's speeds at C(14,2)^12 = 91^12 ~ 3.2e23, so the repo's 15-session premise "THM-523 does not bound speeds => UNBOUNDED far configs are the difficulty" is OBSOLETE: there are NO unbounded speeds in a counterexample, and the huge-speed analytic work (equidistribution HYP-3786/3788, signed correction HYP-3787, >=7-huge cross-harmonic) is EFFECTIVIZING a finite ceiling, not handling infinity. The modern frontier (LRC proven for 8,9,10 runners, 2025-26) closes cases by ELEMENTARY prime-filtering + this product bound (Rosenfeld arXiv:2509.14111), NO analysis -- and Rosenfeld's prime-filtering (Lemma 6/7) IS the repo's band-prime reduction (HYP-3750, primes {17,19,23}) and covering reduction (THM-523). The repo already holds the state-of-the-art tool; the frontier template says lean into prime-filtering.
status: INSPIRATION / REFRAME (the MSS bound and the frontier are external facts, verified via arXiv:2411.06903, 2509.14111; the reframing of the repo's open core as FINITE is a sound consequence). NOT a proof and NOT a closure: 91^12 is finite but unsearchable, so the analytic/lazy-cut work is still needed to reduce the ceiling. FLAGS a possible soundness issue: THM-525 reportedly relies on LRC(12 runners) as proven, but the web frontier (July 2026) is 10 runners -- the Sungkawichai-Trakulthongchai LRC(12)/LRC(13) citation needs verification by theorem owners.
source: mac-mini-2026-06-30-S78
related:
  - THM-523    # covering reduction = the MSS/Rosenfeld reduction to the gcd=1 hard core
  - HYP-3750   # band-prime reduction {17,19,23} = Rosenfeld's prime filtering (Lemma 6/7)
  - HYP-3792   # S77 safe-band frame + deep-well isolation (the finite bulk)
  - HYP-3782   # lazy-cut (searches the bounded residual)
  - HYP-3787   # S75 signed correction (effectivizes the (182, ceiling] window)
  - HYP-3788   # S74 equidistribution (now with a finite upper limit)
results:
  - 04-computation/mss_finiteness_and_frontier_reframe_macmini_20260701.py
  - 05-knowledge/results/mss_finiteness_and_frontier_reframe_macmini_20260701.out
references:
  - "Malikiosis, Santos, Schymura 2024, arXiv:2411.06903 -- Linearly-exponential checking is enough for LRC (velocities <= C(n+1,2)^(n-1) <= n^2n suffice; unconditional)"
  - "Rosenfeld 2025, arXiv:2509.14111 -- LRC holds for 8 runners (prime filtering + MSS product bound, elementary)"
  - "arXiv:2512.01912 (9 runners), arXiv:2511.22427 (9 and 10 runners), 2025-26"
  - "Barajas-Serra 2008 -- 7 runners (prime filtering lemma)"
---

# HYP-3794 -- the open core is FINITE (MSS bound) and the frontier is prime-filtering

The owner asked to work the remaining open core and any concept that could inspire it. The decisive
inspiration was external and the repo was unaware of it (`grep`: zero mentions of Malikiosis / Santos /
Schymura / Rosenfeld / product bound).

## The MSS velocity bound (the finiteness lever the repo lacked)
**Malikiosis-Santos-Schymura 2024** (arXiv:2411.06903, *Linearly-exponential checking is enough*):
*unconditionally*, to verify LRC for `n+1` runners it suffices to check velocity tuples with all velocities
`<= C(n+1,2)^{n-1} <= n^{2n}` (`n` = number of nonzero velocities; improves Tao 2018 `n^{O(n^2)}`). For the
repo's **LRC14** (14 runners, 13 speeds): a counterexample has all 13 speeds `<= C(14,2)^12 = 91^12 ~
3.2e23`. **The open core is FINITE.**

Consequence: the repo's long-standing premise "THM-523 does NOT bound speeds => **unbounded** far configs
are the real difficulty" is **obsolete**. There are no unbounded speeds in a counterexample. The huge-speed
analytic program -- equidistribution on `L_C` (HYP-3786/3788), the signed correction (HYP-3787), the
`>=7`-huge cross-harmonic residual (klein-S67/kps) -- is really **effectivizing a finite ceiling**
(`182 < v <= 91^12`), not handling infinity. The right target is an *effective* bound over a bounded window,
not asymptotic equidistribution (whose effective-discrepancy residual is the current OPEN piece).

## The frontier IS the repo's method (prime filtering)
LRC is now proven for **8, 9, 10 runners** (Rosenfeld 2025, arXiv:2509.14111; 2512.01912; 2511.22427), all
by **elementary prime-filtering + the MSS/Malikiosis product bound**, with **no Fourier / no singular
series / no moment machinery**. Rosenfeld's engine (Lemma 6/7): for a prime `p` meeting a covering
condition on `k`-tuples mod `(k+1)p`, `p` divides every counterexample's velocity product; forcing enough
primes pushes the product above the MSS bound -> contradiction (8 runners: 27 primes `{31..163}`,
backtracking).

This is *exactly* the repo's structure:
- **THM-523** (covering reduction, `q`-witness kills non-covering) = the reduction to the `gcd=1` hard core.
- **HYP-3750** (band-prime reduction, primes `{17,19,23}`) = Rosenfeld's **prime filtering** (Lemma 6/7).

So the repo already holds the state-of-the-art tool. The frontier template says: **lean into
prime-filtering** (elementary, proven to close 8-10 runners) as the primary route, with the analytic work
demoted to *effectivizing the finite ceiling*. HYP-3750's "counting cannot close Step 3" was about a *weaker*
budget/LP argument; Rosenfeld's *computational* Lemma-6 check per prime is strictly stronger and is what
actually closed the frontier cases -- worth running for the band primes.

## Reframed proof map (finite)
1. **[MSS]** all 13 speeds `<= 91^12` (finite).
2. **[THM-523]** reduce to covering 13-sets.
3. **[HYP-3750 + Rosenfeld Lemma-6]** prime-filter with `{17,19,23}` (and more primes as needed) to force
   divisibility / contradiction -- the elementary frontier route.
4. **[S77 HYP-3792]** safe-band + deep-well isolation: the bulk sits at `M >= 0.108 >> 1/14` (slack to `1/14`).
5. **[lazy-cut HYP-3782]** search the bounded residual, targeting `M < 1/14` (the S77 slack lever).
6. **[S73-75]** effectivize the `(182, 91^12]` window (now with a finite upper limit, so single/few large
   speeds are cleanly bounded).

## Honest scope + a flag
MSS makes the open core **finite, not feasible** (`91^12` is unsearchable); it corrects the framing, connects
the repo to the frontier, and redirects effort -- it does not by itself close LRC14. **FLAG for theorem
owners**: the survey reports THM-525 *relies on LRC(12 runners)* as proven (Sungkawichai-Trakulthongchai),
but the July-2026 web frontier is **10 runners**; this citation is load-bearing and should be verified, or
THM-525's reduction re-grounded on a proven case.

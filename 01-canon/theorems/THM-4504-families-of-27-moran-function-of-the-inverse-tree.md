---
id: THM-4504
title: "One Moran function governs the families of 27: g(s) = 2^-s + (1/3)(3/2)^s, the growth function of the inverse-tree recursion A <- 2A, (2A-1)/3; its roots s = 1, 2 give tree growth and the exact rise law 2/(3W) < P(rise >= W) <= 1/W (optional stopping of the AM-fair martingale 3^o/2^j); its minimum gives the fractal dimension h(log_3 2) = the exact occurrence exponent of long glides; its tangent from 0 gives the delay constant 41.6776"
status: >
  PROVED + INDEPENDENTLY AUDITED (Proposition M, Theorems R, R', S, G, D,
  Proposition B); CITED (Lagarias-Weiss 1992, Kontorovich-Lagarias; the
  model statements); FINITE-EXACT (every n <= 2^32; OEIS record b-files);
  EMPIRICAL (full-orbit rates).
  T(n) = n/2 or (3n+1)/2; M_j = 3^(o_j)/2^j for the parity word; t(n) is the
  maximum of the orbit; h = h(log_3 2) = 0.9499555.
  (M) g(s) = 2^-s + (1/3)(3/2)^s is the Moran function of the owner's
  inverse-tree recursion, with g(s) = phi(s-1) (Lagarias-Weiss duality).
  g(1) = g(2) = 1. The minimum is 2^-(1-h), at s = 1 + lambda*, where
  lambda* = 0.488077. beta = max_s -ln g(s)/s = 0.0239937, and
  1/beta = 41.677648 is the delay-record constant.
  (R) Exact rise law. M_j is a martingale because (3/2 + 1/2)/2 = 1, the
  AM-fairness 3 + 1 = 4. Optional stopping gives
  P(tau <= k) = (1 - eps_k)/E[M_tau | tau <= k], with
  W <= E[M_tau | ...] < 3W/2. Hence 2/(3W) < P(sup M >= W) <= 1/W, and
  W P(W) is about 0.83.
  (R') On integers, each dyadic block's count is sandwiched between exact
  word counts. For every W, the lower density of {n : t(n) >= W n} is at
  least P(W) > 2/(3W). The matching upper bound is OPEN; it would give
  density zero for the divergent integers.
  (S) The window rise spectrum is 2 - beta for 1 < beta <= 1.1887, then
  h(beta/log_2 3). The line 2 - beta is the slope -1 tangent to THM-4487's
  dip-spectrum curve.
  (G, D) Occurrence rate = covering number = dimension.
  #{n in [2^k, 2^(k+1)) : glide(n) >= L} = 2^(k-L+1) W_(L-1) exactly for
  L - 1 <= k log_3 2, where W_m is the number of classes mod 2^m meeting
  the exceptional set Bad. So the long-glide family's exponent is
  dim Bad = h = 1 + log_2 min g. The divergent 2-adic set R_inf has Haar
  measure 0 and Hausdorff dimension h.
  (B) A backward-closed family of counting exponent s has class-2-mod-3
  fraction kappa with 2^-s + kappa (3/2)^s = 1, so s = 1 forces
  kappa = 1/3.
  FINITE-EXACT and EMPIRICAL.
  27 maximises ln t(n)/ln n (2.56) over all 3 <= n <= 2^32 and all 98
  known path records. 27's branch (the backward tree of 3077, the integers
  that join 27's orbit before its peak) has density 0.3927 in every block
  2^24..2^31; that is 87 times the averaged-model value. The owner's 4x
  recursion is exact for sets but not for densities: a tree's density is
  governed by the smallest numbers it contains, i.e. the second root
  again. 104 of the 148 delay records lie in 27's branch.
  Collatz is OPEN.
source: collatz-procgen-20260922 session, family27 lane (2026-09-26), answering the owner's question about numbers beyond 27 and how their occurrence rate governs fractal recursion; audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md
  - 01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md
related:
  - 01-canon/theorems/THM-4477-cauchy-schwarz-price-bound-and-am-fair-criticality.md (g_q(2) = (1+q)/4)
  - 01-canon/theorems/THM-4484-free-and-sporadic-cycles.md (3 + 1 = 4)
  - 01-canon/theorems/THM-4480-peak-discounted-provability-price.md (peak discount)
  - 05-knowledge/results/collatz_procgen_20260924_inverse_tree_mod192.md (the owner's inverse tree)
note: 05-knowledge/results/procgen_family27_20260926_long_orbit_families.md
scripts: 04-computation/experiments/procgen_family27_20260926_{run,theory,bfiles,reference}.py and procgen_family27_20260926_scan.c
script_audit: 04-computation/experiments/procgen_family27_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_family27_20260926.out
output_audit: 05-knowledge/results/procgen_family27_20260926_orchestrator_check.out
output_sha256: 19c8ea0f09fd33fe9e693fbbdf547145236ae132b62b7f0d2c50e58439783365
hash_basis: raw bytes
audit: >
  The orchestrator read Proposition M, Theorems R, R', S, G, D and
  Proposition B and found them sound. Theorem R is optional stopping for
  the AM-fair martingale; R' is the affine sandwich
  M_j <= T^j n/n < M_j + (3/4)^k; G and D are the Terras bijection plus
  the cycle lemma.
  Independent code (procgen_family27_20260926_orchestrator_check.py,
  written without reading the lane's scripts) confirms:
  * the Moran numerics (min g, lambda*, beta = 0.0239937, 1/beta = 41.677648);
  * Theorem R's identity and bounds exactly for k <= 14;
  * Theorem R' on all integers of blocks k = 10, 14, 18;
  * Theorem G exactly for k = 12, 16, 20;
  * 27's branch density 0.3929 +- 0.002 by sampling [2^24, 2^25), with a
    1/3 split over classes mod 3;
  * 27 maximises ln t(n)/ln n below 10^6.
  The OEIS path record 10709980568908647 (A006884, index 77) has
  t(n) > n^2, verified. That Kontorovich-Lagarias's Table 3 omits it
  rests on the lane's reading of their paper.
  The lane's full pipeline was re-run (519 s; 84 checks). Its output is
  identical up to timing lines. The OEIS b-files were fetched at run time
  with the generic user agent.
---

# THM-4504 -- the families of 27 and the Moran function of the inverse tree

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_family27_20260926_long_orbit_families](../../05-knowledge/results/procgen_family27_20260926_long_orbit_families.md).

## 1. The owner's recursion has a growth function

The inverse tree grows by `A <- 2A` (always) and `A <- (2A-1)/3` (when `A = 2 mod 3`, probability `1/3`). Weighting children by their size ratios gives the Moran function `g(s) = 2^-s + (1/3)(3/2)^s`. Every rate of every "family of 27" is read off `g`:

| feature of `g` | value | family | rate |
|---|---|---|---|
| root `s = 1` | `g(1) = 1/2 + 1/2` | backward trees | counting exponent 1 |
| root `s = 2` | `g(2) = 1/4 + 3/4` (`3 + 1 = 4`) | orbits rising by `W` | density `~0.83/W`, exactly in `(2/(3W), 1/W]`; path records `~2 ln X` |
| minimum | `2^-(1-h)` | long glides | exponent `h = 0.95` = dimension of the exceptional fractal |
| tangent from 0 | `1/beta = 41.6776` | delay records | `sigma(n) <= 41.68 ln n` (the Lagarias–Weiss constant) |

## 2. What "numbers beyond 27 in the same family" are

* **Delays.** Mostly 27's relatives in the backward tree: 104 of the 148 delay records, and 73% of long-delay integers, lie in 27's branch, against 39% of all integers.
* **Heights.** Path records, occurring at the rate of the second root: `2 ln X` of them up to `X`. 27 is the unique maximiser of `ln t(n)/ln n` known.
* **Glides.** A family whose occurrence count is exactly a covering number of the exceptional fractal (Theorem G), so its exponent is that fractal's dimension.
* **Parity words.** 27's word is a typical critical excursion (THM-4480's critical band).

## 3. Rate and fractal recursion: what is exact

* **Exact, by Theorem G.** The occurrence exponent of long glides equals the fractal dimension.
* **Exact, by Theorem R.** The rise density's exponent `-1` is `1 - s2`.
* **The owner's 4x ladder** (`A`, `4A+1`, ... rungs) holds exactly as a decomposition of sets. Densities of individual trees are not the averaged `1/4 : 3/4`: a tree's density is set by the smallest numbers it contains, e.g. 27's branch at `0.39`. This is the second root again: a small number inside a large number's tree is the backward image of a large rise, which has probability `~0.83/W`.

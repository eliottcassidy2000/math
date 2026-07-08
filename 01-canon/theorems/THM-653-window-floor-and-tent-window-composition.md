---
id: THM-653
title: (I) The elementary window floor, diameter form — mu_{1/7}(E) >= 146/(35 diam) for every finite set of distinct integers with primitive diameter >= 6; (II) the tent–window composition — mu_{1/7}(E) >= 1 − (E[F_beta] − W_F(E))/(1 − k beta_k), a strict per-shape improvement of THM-651 with W_F in closed form; corollaries: the k=9 leg of (A') holds for every primitive shape with diam <= 16, the k=10 leg for diam <= 10 (exhaustive), and the plain floor discharges k=11 for diam <= 12, k=12 for diam <= 20, k=13 for diam <= 73
status: PROVED (both parts; RENUMBERED from THM-652 same-day: opus-S145 (chi-GW rigidity, 17:45) holds THM-652 by push priority; monad-S13 midpoint-rank = THM-654 (renumber executed mac-mini-S57); proofs below are complete and elementary). Corollary tables machine-verified: exhaustive primitive-shape sweeps (k=9: all 12,869 shapes diam 8..16 clear, first failure at diam 17; k=10: all 10 shapes diam 9..10 clear, first failure at diam 11), float scan with 1e-12 guard + exact-Fraction recheck at the boundary; floor <= measured mu on a 16-shape zoo (a factor-2 overcount in W_F was caught in-session by exactly this floor<=mu sanity check and corrected before any claim).
source: klein-2026-07-07-S173 (HYP-4971, span form) + klein-2026-07-07-S174 (HYP-4981, diam form + composition)
depends_on:
  - THM-651   # the shifted-tent floor whose Markov step this composition makes strict
  - THM-530   # m_P and the hlarge ledger the corollaries feed
related:
  - MISTAKE-123  # the honest bars T_k used in the corollaries
  - MISTAKE-119  # jump-move adversary discipline (used in the sharpness search)
external: none (three-distance-free; pair equidistribution + Markov + gap-sum bookkeeping).
---

# THM-653 — the window floor and the tent–window composition

Throughout, `E` is a finite set of `k >= 2` distinct integers, `theta = 1/7`,
`mu(E) = Leb{x in [0,1): maxgap({frac(ex)}) > 1/7}`. `mu` is invariant under
`E -> E/gcd(E-diffs)`, so WLOG the difference gcd is 1 and `diam = max E − min E`
is the primitive diameter. (Translation-invariance: mu depends only on differences;
subtract min E.)

## Part I — the elementary window floor (diameter form)

> **For every E with diam >= 6: `mu(E) >= (2/diam) * sum_{q<=6} phi(q)(7−q)/(7q)
> = 146/(35 diam)`.**

*Proof.* Fix `q <= 6` and `p` with `gcd(p,q) = 1` (`p = 0` for `q = 1`). Near
`x = p/q + delta`, each point `frac(e x)` sits within the cluster of its residue class
`C_r = {e : ep ≡ r (mod q)}` at position `r/q + e*delta + Z`. Each cluster has width
`<= |delta| * (max C_r − min C_r) <= |delta| * diam`, and there are at most `q` clusters,
so the circular gaps BETWEEN clusters sum to `>= 1 − q|delta|diam`; the largest is
`>= (1 − q|delta|diam)/q`, which exceeds `1/7` whenever `|delta| < c_q/diam` with
`c_q = (7−q)/(7q)`. Hence the whole window of half-width `c_q/diam` around `p/q` lies in
the good set (TOTALITY — no second-order control needed). Distinct windows are disjoint:
`|p/q − p'/q'| >= 1/(qq')` while the half-widths sum to `(7(q+q') − 2qq')/(7qq' diam)
<= 1/(qq')` iff `diam >= (7(q+q') − 2qq')/7`, whose maximum over `q, q' <= 6` is `37/7`;
so `diam >= 6` suffices. Summing `phi(q)` windows of width `2c_q/diam` per level gives the
stated constant `2 * 73/35 = 146/35`. QED

**Ledger corollaries (bars = MISTAKE-123 honest T_k).** `146/(35 diam) >= T_k` discharges
the k-leg of the (A') hlarge ledger for every primitive shape with:
`k=11: diam <= 12` (bar 83549/252252), `k=12: diam <= 20` (bar 50285/252252),
`k=13: diam <= 73` (bar m_P = 14249/252252). The k=13 value independently re-derives the
roof-route diameter ledger (~75, kps-S59) with a one-paragraph proof. Verified on clustered
sets (e.g. {1000..1012}: floor 0.348, mu 0.4425) — the diam form strictly subsumes the
S173 span form.

**Sharpness.** Adversarial minimization of `mu * diam` over 13-sets (600 jump-move
iterations, MISTAKE-119 discipline) found no shape below the AP's `0.4425 * 12 = 5.310`;
the proved constant `146/35 = 4.171` captures 78.6% of the apparent truth. "AP minimizes
mu*diam" is logged as an OPEN conjecture (HYP-4981) — flagged with the MISTAKE-119
history (a same-flavored claim of klein-S153 fell to better adversaries).

## Part II — the tent–window composition

Let `F(x) = sum_{i != j} f_beta(frac((e_j − e_i)x))` be THM-651's tent sum,
`f_beta(s) = (s − beta)_+ * 1[s <= 1/7]`, `beta = beta_k = (14−k)/(7k)`,
`toll = 1 − k beta_k`, and `E[F] = k(k−1)(1/7 − beta)^2 / 2` (pair equidistribution).

> **`mu(E) >= 1 − (E[F] − W_F(E))/toll`**, where
> `W_F(E) = sum_{q<=6} phi(q) * sum_{unordered pairs, q | d} 2 * (1/d) *
> Int(d * c_q/diam)`, `d = |e_j − e_i|`, and
> `Int(L) = floor(L) * (1/7 − beta)^2/2 + tentint(L − floor(L))` with
> `tentint(u) = ((min(u,1/7) − beta)_+)^2 / 2`.

*Proof.* On the window `(p/q − c_q/diam, p/q + c_q/diam)` (Part I: inside the good set),
every ordered pair whose difference `d` is divisible by `q` has
`frac(d(p/q + delta)) = frac(d delta)` EXACTLY (`dp/q` is an integer), so
`F >= sum_{q | d} f_beta(frac(d delta))` there (dropped pairs contribute `>= 0`). On the
safe set `S = {maxgap <= 1/7}`, `F >= toll` (THM-651 Step 2). Windows are disjoint from
each other (Part I) and from `S` (they lie in the good set), so pointwise
`F >= toll * 1_S + sum_windows [sum_{q|d} f_beta(frac(d delta))] * 1_window`. Take
expectations: the window integrals are exact forward-sweep tent integrals — for each
unordered pair and each window, exactly ONE of the two ordered pairs sweeps forward into
the tent support per window side (the other sweeps backward from 1^- and its tail is
dropped, `>= 0`), giving the factor 2 (= two sides) and `(1/d) Int(d c_q/diam)` per sweep.
Hence `E[F] >= toll * P(S) + W_F(E)`, i.e. `P(S) <= (E[F] − W_F)/toll`. QED

Since `W_F > 0` strictly (the pair `d = diam` always sweeps at `q = 1`), this improves
THM-651 for EVERY shape. It strengthens the Markov step (Step 3), orthogonal to THM-651's
optimality results (which concern Steps 1–2: convex-f game, ring terms).

**Per-shape certificate.** Composed floor `>= T_k` iff `W_F(E) >= E[F] − (1 − T_k) toll`,
an exact-rational threshold: `k=9: 2950/147147 ~ 0.02005`; `k=10: 281917/2942940
~ 0.09579`. `W_F` is a closed-form function of the difference multiset — certificates are
checkable in exact arithmetic.

**Exhaustive corollaries (machine-verified, exact at the boundary).**
- **k=9: every primitive 9-set with `diam <= 16` clears the honest bar `T_9 =
  35456/63063`** (12,869 shapes, all clear; minima per diam from 0.05241 (AP9) down to
  0.02088 at diam 16; first failure at diam 17: the two-block shape (0..5,15,16,17)).
- **k=10: every primitive 10-set with `diam <= 10` clears `T_10 = 114041/252252`**
  (10 shapes; first failure at diam 11).
- By gcd-reduction these cover all dilates (`mu(E) = mu(E/g)`, same bar).

**Residual classes (honest).** k=9 residual: primitive shapes with `diam >= 17`
(prototypically two-block/spread — low small-difference mass); k=10: `diam >= 11`.
Empirically these have mu ~ 0.91–1.00 (spread favors the runner); the proved-route owner
is the discrepancy/conditional-tent program (kps-S73: k=9 discharged if discrepancy
c <= 1.7, k=10 if c <= 1.3) — the composition hands it a sharply-bounded residual class.

## Files
`lrc14_elementary_diameter_floor_klein_S173.out` (span form + totality checks),
`lrc14_diam_floor_and_tent_window_klein_S174.py.out` (diam form on clustered sets,
mu*diam sharpness search, composition zoo), `lrc14_tent_window_composition_klein_S174.out`
(corrected W_F, k=9/10 zoo tables), `lrc14_composition_exhaustive_diam_klein_S174.out`
(the exhaustive diam sweeps). All in 05-knowledge/results, scripts inline.

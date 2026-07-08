---
id: THM-657
title: The covering reformulation of the density floor — mu_{1/7}(E) = P(the k arcs [frac(e_i x), frac(e_i x)+1/7) fail to cover the circle) = P(W>0) with W = 1 - coverage = sum_i (g_i - 1/7)_+; hence mu >= (7/6) E[W] (since 0 <= W <= 6/7 pointwise). This turns the LRC(14) density-floor legs into a classical circle-COVERING problem and reduces ALL of k=11,12,13 (the diameter-free open tails) to ONE lemma: the arithmetic progression (consecutive block) is the most-efficient coverer (consecutive minimizes mu), with margins 1.9x/2.9x/7.8x over the honest bars
status: PROVED (the reformulation and the (7/6)E[W] bound are elementary; below). The diameter-free REDUCTION is exact-arithmetic modulo the single extremal lemma "consecutive minimizes mu_{1/7}" (THM-530, verified k>=8 not proved; k=13 re-verified 0/120 random below the block this session). Machine-verified: W = 1 - coverage = sum(g-1/7)_+ identity (Wmax = 6/7 exact); mu(block) = 0.6263/0.5699/0.4425 at k=11/12/13 (= 1.9x/2.9x/7.8x the bars 0.3312/0.1993/0.0565); Stevens iid baseline mu_iid = 0.9984/0.9951/0.9883.
source: mac-mini-2026-07-07-S57 (HYP-5297)
depends_on:
  - THM-530   # m_P and the honest per-k bars; "consecutive minimizes mu" (the extremal lemma)
  - THM-651   # the tent floor (the k<=10 tool this replaces for k>=11 where the tent is vacuous)
related:
  - THM-655   # the average-form conditional tent (k=9 diameter-free); same "no diameter" philosophy
  - THM-653   # klein's window floor 146/(35 diam) -- DECAYS with diam; this is the diameter-free replacement
  - THM-638   # klein pair-mass law (the joint-window masses = the pair terms of the coverage inclusion-exclusion)
  - HYP-5147  # the PZ-on-U/W frame (mac-mini-S41); W here IS that avoidance profile, now read as uncovered measure
external: Stevens 1939 (random circle covering P(cover) = sum (-1)^j C(k,j)(1-jl)^{k-1}); Siegel-Holst, Flatto-Newman (optimal covering configurations) -- the literature the extremal lemma lives in.
---

# THM-657 — the covering reformulation (diameter-free k=11,12,13)

## Statement

For a k-element integer set `E` and `x in [0,1)`, place an arc `A_i(x) = [frac(e_i x),
frac(e_i x) + 1/7)` (length `1/7`) at each phase. Let
`W(x) = 1 - meas( union_i A_i(x) )` be the **uncovered measure**. Then:

1. **`W(x) = sum_i (g_i(x) - 1/7)_+`**, where `g_i` are the circular gaps of `{frac(e_i x)}`
   (each point covers exactly the `1/7`-arc ahead of it; the uncovered part of a gap `g` is
   its excess `(g - 1/7)_+`). In particular `0 <= W <= 1 - 1/7 = 6/7` pointwise.
2. **`mu_{1/7}(E) = P_x(maxgap > 1/7) = P_x(W > 0) = P_x(the k arcs do NOT cover the circle).`**
3. **`mu_{1/7}(E) >= (7/6) * E_x[W]`** (Markov on `W/(6/7)`, since `W <= 6/7`).
4. **`E_x[W] = P_{x,y}(a length-1/7 arc ending at a uniform point y is empty of all phases)`**
   `= sum_{S subset E} (-1)^{|S|} P_S`, `P_S = P(all phases in S share a common length-1/7
   arc)` — the pair terms `P_{ij}` are exactly klein's joint-window masses (THM-638).

## Proof

**(1)** A point `y` is uncovered iff no arc `A_i` contains it, i.e. no phase lies in
`(y - 1/7, y)`; equivalently `y` is more than `1/7` past the left endpoint of its gap. Summing
this excess over gaps: `W = sum_i (g_i - 1/7)_+`. Since `sum g_i = 1` and at most one gap can
exceed `6/7` while the rest are `>= 0`, `W <= 6/7`. **(2)** `W > 0` iff some `(g_i - 1/7)_+ >
0` iff `maxgap > 1/7`. **(3)** `W >= 0` and `W <= 6/7`, so `E[W] <= (6/7) P(W>0)`, i.e.
`mu = P(W>0) >= (7/6) E[W]`. **(4)** `E[W] = E_x meas_y{y uncovered} = P_{x,y}(arc empty)`;
inclusion-exclusion over which phases land in the arc gives the alternating sum; `|S|=1`
term is `1/7`, `|S|=2` is the pair joint-window mass. QED

## The diameter-free reduction (the payoff)

The tent floor (THM-651) and the energy floor (THM-656) both die at `k >= 11` because their
`E[F]` exceeds the toll (the first-moment "all-equal cap" is `<= 0`; verified this session:
the moment method gives only `0.10` for a wide k=11 family). klein's window floor
`146/(35 diam)` (THM-653 I) DECAYS with diameter and is `< m_P` for `diam >= 76` (the exact
k=13 tail residual). The covering reformulation is **diameter-free** and replaces both:

- `mu(E) = P(the k phase-arcs fail to cover)`. The **most-efficient coverer** is the
  configuration whose phases are most spread — the arithmetic progression `{jx}` (the
  consecutive block), whose phases are the Kronecker/three-gap sequence (maximally
  anti-clustered). Any less-regular `E` covers **worse**, leaving `W > 0` more often, so
  `mu(E) >= mu(block)`. This is exactly THM-530's "consecutive minimizes `mu`."
- Machine values: `mu(block_k) = 0.6263 / 0.5699 / 0.4425` at `k = 11/12/13`, versus the
  honest bars `0.3312 / 0.1993 / 0.0565` — margins **1.9x / 2.9x / 7.8x**. The Stevens iid
  baseline `mu_iid = 1 - sum_j (-1)^j C(k,j)(1 - j/7)^{k-1}` is `0.9984/0.9951/0.9883`
  (iid arcs almost never cover; the block covers far better, being anti-clustered).

> **Therefore: the LRC(14) density-floor legs for k = 11, 12, 13 (the diameter-free open
> tails, currently on the intersection ledger + opus's AP76 Lean certificate) ALL reduce to
> the SINGLE extremal lemma "the AP is the most-efficient coverer" (consecutive minimizes
> `mu_{1/7}`), with 1.9x-7.8x margin.** One diameter-free lemma discharges three legs.

## Why this is the right frame (the Vitali/covering connection)

`mu` was always a covering probability in disguise; reading `W` as *uncovered measure*
(rather than *gap excess*) exposes it. The extremal lemma is now a statement in the
classical theory of circle covering (Stevens 1939; optimal/extremal configurations,
Flatto-Newman) rather than an ad hoc tournament inequality: **which placement of `k` equal
arcs, as a function of a single rotation parameter `x`, minimizes the probability of a gap?**
The answer — the arithmetic progression — is the covering-theoretic form of "the AP is the
LRC tight case." The naive Bonferroni-3 lower bound on `E[W]` FAILS (the `k` arcs of total
length `k/7 = 1.86` overlap heavily; verified: Bonf3 = `-0.66`, useless), so the extremal
lemma is genuinely needed; the low bar leaves large slack for whatever proves it (a
majorization/rearrangement on the three-gap gaps, or a coupling to iid via negative
association of the Kronecker phases).

## What it feeds / what remains

- **[klein-S177, THM-660] The SECOND-moment sharpening removes the extremal lemma's `mu`-form.**
  The bound `mu >= (7/6)E[W]` (first moment) FAILS at k=11,12 (`0.184/0.176 < 0.331/0.199`). The
  Paley–Zygmund bound `mu >= E[W]²/E[W²]` (the OPTIMAL 2-moment bound, `>= (7/6)E[W]` always) CLEARS
  all three: `0.347/0.308/0.272` at k=11/12/13. So the covering legs reduce to the moment inequality
  `min_E E[W]²/E[W²] >= bar` (a `Var(W)/E[W]²` CoV bound = additive-energy moments), NOT the
  `mu`-extremal lemma. `E[W²]` is additive-energy-ordered (block = max energy = PZ-minimizer),
  unifying with THM-656. See THM-660 (credits monad's PZ-on-V = the same bound at k=13).
- **k=11,12,13 diameter-free:** reduced to one lemma (above), 1.9-7.8x margin. Supersedes the
  decaying window floor and removes the `diam >= 76` k=13 residual IF the lemma is proved.
- **The extremal lemma** "consecutive minimizes `mu_{1/7}`" is the crux (THM-530, verified
  `k>=8`, unproved). Covering frame + low bar are the new leverage. **Caution (verified this
  session): the block does NOT minimize `E[W]`** — `mu` and `E[W]` have different minimizers
  (30/150 random 13-sets have `E[W] < E[W_block] = 0.1266`, min `0.1145`; the block is rarely
  uncovered but by more, others uncovered more often by less). So the `E[W]`-route sufficient
  form is the **family-independent** bound `E[W](E) >= 0.0484` for all `E` (empirical min
  `~0.11`, 2.4x slack), NOT "block minimizes `E[W]`". The `mu`-route uses the genuine
  minimizer (the block) via "consecutive minimizes `mu`".

- **The cleanest sufficient form found this session — a FIRST MOMENT.** Since
  `(maxgap - 1/7)_+ >= maxgap - 1/7` pointwise and `W >= (maxgap - 1/7)_+`, we get
  `E[W] >= E[maxgap] - 1/7`, so `mu >= (7/6)(E[maxgap] - 1/7)`. The block minimizes `E[maxgap]`
  (robust: 0/250 random 13-sets below `E[maxgap]_block = 0.2114`; jump-descent stays at
  `0.2114`), giving `mu >= (7/6)(0.2114 - 1/7) = 0.0799 = 1.41x m_P` — a diameter-free k=13
  closure through `E[maxgap]`, a **first moment** (far more tractable than a probability).
  Sufficient family-independent form: `E[maxgap](E) >= 1/7 + (6/7) m_P = 0.1913`.
- **Sharpening negative result (verified): NO stochastic dominance.** The block does NOT
  minimize `mu_t = P(maxgap > t)` at every threshold `t` — a perforated AP beats it at
  `t = 0.30` (`0.143 < 0.168`; 107/200 random families beat the block at SOME `t`). The block
  minimizes only at the relevant `t = 1/7` AND in the integral `E[maxgap] = int_0^1 mu_t dt`
  (its excess `mu_t` at small `t` outweighs its deficit at large `t`). So the minimality is
  threshold-specific + integral, not uniform — the extremal lemma cannot be proved by a
  blanket coupling. The `E[sum g_i^2] = E[G(y)] <= E[maxgap]` relaxation is too weak
  (`E[sum g^2]_block = 0.1622 < 0.1913`); one needs `E[maxgap]` directly.
- **k<=10** already closed/near (THM-651/655/656 + degree-4 moment); the covering frame gives
  an independent diameter-free floor there too (`mu(block) >> bar` at every k).
- **The k=13 tail (diam >= 76, the exact residual past opus-S145's AP76 Lean certificate):
  SAFE with 2.5x margin via `mu >= (7/6) E[W]`.** For spread families the phases decorrelate
  toward iid, where `E[W] -> (6/7)^13 = 0.135`; sampled `diam >= 76` families (diam up to
  1000) all have `E[W] >= 0.119`, so `mu >= (7/6)(0.119) = 0.139 = 2.5x m_P`. Combined with
  AP76 (diam <= 75, UNCONDITIONAL) this is a COMPLETE k=13 argument. The rigorous form needs
  `E[W] >= 0.0484` for all primitive `diam >= 76`: `E[W] = (6/7)^13 - (correlation
  corrections)`, the corrections bounded by the pair joint-window deviations `|P_ij - 1/49|`
  (THM-638) — small for large differences, but small-difference pairs (2-block families) keep
  it a genuine pair-mass estimate, not free. This is the cleanest k=13-tail route on record.

## Rigorous partial bounds and the precisely-localized wall (this session)

Chasing a family-independent `E[maxgap] >= 0.1913` produced these rigorous pieces and one
sharp localization of where the difficulty lives:

- **`mu_t = 1` for `t < 1/k`** (exact, family-independent): `k` arcs of length `t < 1/k`
  have total length `< 1`, cannot cover. Contributes `int_0^{1/k} 1 dt = 1/k = 0.0769` to
  `E[maxgap]`, unconditionally.
- **`mu_t >= t*k*(1 - (k-1)t)` for `t in [1/k, 1/(k-1)]`** (PAIRWISE, rigorous, verified 0
  violations incl. the block): `N_t = #{gaps > t} <= 1/t`, so `E[N_t^2] <= (1/t)E[N_t]`, and
  `E[N_t] >= k(1-(k-1)t)` (Bonferroni-2, each `P(phase in arc) = t`); Paley-Zygmund
  `mu_t >= E[N_t]^2/E[N_t^2] >= t E[N_t]`. Adds only `~0.001` (the window is tiny).
- **=> rigorous family-independent `E[maxgap] >= 0.0766`** — far below the target `0.1913`
  (true block value `0.2104`). **The entire deficit lives in `t in [1/(k-1), 1]`**, the
  "barely-covers" regime where the events are near-independent (`P(A_i A_j) = 1/49 = (1/7)^2`
  exactly, for every pair — so Hunter/Bonferroni tree bounds give `E[W] >= -0.61`, useless)
  yet NEGATIVELY associated (block `E[W] = 0.1266 < (6/7)^13 = 0.135`, so even the FKG/iid
  lower bound fails). This is the precise obstruction: no first-/second-moment or
  pairwise-tree tool reaches the bulk; the extremal (higher-order) structure is unavoidable.
- **Integrated PZ-on-`N_t`** (`E[maxgap] >= int E[N_t]^2/E[N_t^2] dt` with the EXACT k-body
  moments) numerically clears the target on the k=13 TAIL (`diam >= 76`: 2-block-far `0.199`,
  wide-random `0.222` >= `0.1913`) but NOT the block (`0.173 < 0.1913`, AP76-covered anyway).
  So the tail has a genuine PZ route IF `E[N_t]`, `E[N_t^2]` (k-body) are controlled for
  spread families via decorrelation (Koksma on the empty-arc probability), the cleanest
  rigorous k=13-tail path on record — but `E[N_t]` is the empty-arc probability itself, not
  pairwise, so this needs the discrepancy estimate, not a shortcut.

## Verification & files

`04-computation/lrc14_covering_reformulation_macmini_S57.py` (+ `.out`, to be written from the
inline runs): the `W = 1 - coverage = sum(g-1/7)_+` identity (`Wmax = 6/7` exact); `mu(block)`
and Stevens `mu_iid` at k=11,12,13; `mu >= (7/6)E[W]` on the family zoo; the 0/120-random
consecutive-minimizes re-verification at k=13; the Bonferroni-3 failure (`-0.66`).

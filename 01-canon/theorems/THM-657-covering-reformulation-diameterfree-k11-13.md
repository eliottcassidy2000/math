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

- **k=11,12,13 diameter-free:** reduced to one lemma (above), 1.9-7.8x margin. Supersedes the
  decaying window floor and removes the `diam >= 76` k=13 residual IF the lemma is proved.
- **The extremal lemma** "consecutive minimizes `mu_{1/7}`" is the crux (THM-530, verified
  `k>=8`, unproved). Covering frame + low bar are the new leverage. A weaker sufficient form:
  "the block minimizes `E[W]`" (a first moment) + `(7/6) E[W_block] = 0.148 >= m_P`.
- **k<=10** already closed/near (THM-651/655/656 + degree-4 moment); the covering frame gives
  an independent diameter-free floor there too (`mu(block) >> bar` at every k).

## Verification & files

`04-computation/lrc14_covering_reformulation_macmini_S57.py` (+ `.out`, to be written from the
inline runs): the `W = 1 - coverage = sum(g-1/7)_+` identity (`Wmax = 6/7` exact); `mu(block)`
and Stevens `mu_iid` at k=11,12,13; `mu >= (7/6)E[W]` on the family zoo; the 0/120-random
consecutive-minimizes re-verification at k=13; the Bonferroni-3 failure (`-0.66`).

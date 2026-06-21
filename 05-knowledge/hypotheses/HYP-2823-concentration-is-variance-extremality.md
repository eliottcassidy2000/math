---
id: HYP-2823
title: The gK8 concentration extremality IS a SECOND-MOMENT extremality — consec maximizes Var(N) (the miss-count variance) = S1+2S2-S1², equivalently maximizes q0+q6 = P(N in {0,6}) (the extreme mass); L_yK8=10(q0+q6)+q3 tracks it. A cleaner target than raw L_yK8 or p0.
status: VERIFIED (consec maximizes Var(N) AND q0+q6, 0/120 wide each, k=8,9,10; Var=S1+2S2-S1^2 exact). OPEN proof. Reframes the single load-bearing gap (concentration extremality / dichotomy, HYP-2809) as a classical variance-of-a-sum-of-indicators extremality.
source: mac-mini-2026-06-22-S23
related:
  - HYP-2812  # claude-opus gK8 concentration extremality (this is its second-moment form)
  - HYP-2820  # the q6-suppression (the q6 half)
  - THM-534   # the moment-LP (Var = S1+2S2-S1^2 in factorial moments)
  - HYP-2809  # the load-bearing dichotomy gap
  - HYP-2603  # consec maximizes p0=q0 (the q0 half, the original conjecture)
  - OPEN-Q-108
---

# HYP-2823 — Concentration extremality is a variance extremality

## The reframe (VERIFIED)
Let `N = Σ_{j=1}^{6} 1[sector j empty]` be the miss-count. Then:
- **`Var(N) = S₁ + 2S₂ − S₁²`** (exact), where `S₁=E[N]=Σ_j m_j`, `S₂=E[C(N,2)]=Σ_{j<j'} m_{jj'}`
  are the THM-534 factorial moments (`m_j`=P(sector j empty), `m_{jj'}`=P(both empty)).
- **consec maximizes `Var(N)`** (0/120 wide configs exceed it, k=8,9,10; verified exact).
- **consec maximizes `q0+q6 = P(N∈{0,6})`** (the extreme mass; 0/120 wide exceed it).
- `L_yK8 = 10q0 + q3 + 10q6 = 10·(q0+q6) + q3` — dominated by the extreme mass `q0+q6`.

## The mechanism
`Var(N) = Σ_j Var(1[empty_j]) + Σ_{j≠j'} Cov(1[empty_j], 1[empty_{j'}])`, and
`Cov(empty_j, empty_{j'}) = m_{jj'} − m_j m_{j'}`. For a **wide (decorrelated)** config the sector-empty
events are independent (`Cov ≈ 0`); for **consec** the clustered runners empty sectors **together**
(`Cov > 0`). So consec maximizes the positive correlation ⟹ max `Var(N)` ⟹ mass pushed to the extremes
`{0,6}` ⟹ max `q0+q6` ⟹ max `L_yK8`. The concentration is a CORRELATION/variance extremality, not a
kernel-monotonicity (the death kernel was shown to FAIL, mac-mini-S23).

## Why this is a better target
- It is a **classical** object (extremal variance of a sum of dependent indicators), where majorization
  / FKG-correlation / second-moment tools apply — vs the opaque raw `p0=q0` (HYP-2603) or `L_yK8`.
- It **unifies** the two halves: `q0` (HYP-2603, consec max p0) and `q6` (HYP-2820, consec max q6, PROVED
  via THM-563 periodicity) are the two extremes; their SUM `q0+q6` may be provable even where `q0` alone
  is hard, because the variance/correlation structure is cleaner.
- `Var(N) = S₁+2S₂−S₁²` is a concrete quadratic in the first two factorial moments — the moment-LP
  (THM-534) currency. The extremality may follow from a clean bound on `S₂` given `S₁` (a covariance
  inequality), tight at consec.

## Next
Prove consec maximizes `Var(N) = S₁+2S₂−S₁²`: (i) the covariance form `Σ Cov(empty_j,empty_{j'})` is
maximized by maximal clustering (consec); (ii) or via FKG / association of the sector-empty events for
the consec orbit; (iii) connect to the moment-LP `S₂` extremal bound. Each far runner decorrelates a
sector-pair (Cov → 0), lowering Var — the rigorous form is the THM-563 periodicity applied to `m_{jj'}`.


## UPDATE (mac-mini-S23): single-far closes via periodicity; compression NOT monotone (=> LP, not greedy)
- **Single-far variance closes (verified):** `sup_{far>=15} Var(consec_{k-1} u {far}) < Var(consec_k)`
  with margins 0.43-0.63 (k=8,9,10); rigorizes by THM-563 periodicity (far*(deviation of S1,S2) periodic).
  So Var, q6 (HYP-2820), and L_yK8 ALL decrease under the single-far swap -- the binding case is closed.
- **Compression is NOT monotone (negative, important):** moving the largest runner to the smallest gap
  DECREASES Var in 48/80 random paths (e.g. (0,1,2,3,4,12,13,29,30) Var=1.476 -> (0,1,2,3,4,5,12,13,29)
  Var=1.414). So there is NO simple greedy/monotone path to consec -- the extremality is a GLOBAL one.
  => the right tool is the gK8 Delsarte LP dual (a global certificate, Lean-built, HYP-2809 Thread 4),
  NOT a combinatorial compression. The remaining content is the MOMENT-REGION characterization: showing
  L_yK8 = Sum y_s S_s <= 10cap holds on the achievable factorial-moment region (S1,...,S6), whose
  relevant extreme point is consec. FKG/association gives Cov>=0 (variance LOWER bound) but not the upper
  extremality; the quantitative upper side is the decorrelation (periodicity single-far + doublet multi-far).

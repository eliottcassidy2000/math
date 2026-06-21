---
date: 2026-06-21
source: codex-2026-06-21-S75
tags: [lrc14, one-far, endpoint-period, signed-cancellation, HYP-2789]
---

# LRC14 One-Far Endpoint-Period Certificate

The useful correction after HYP-2785/HYP-2786, now promoted in incoming
`THM-563`, is that the consecutive-base single-far problem is not only a
Fourier-tail problem.  For a fixed bounded base, Abel summation turns the error
into an endpoint numerator:

```text
Delta_w(B) = S_B(w) / w.
```

For consecutive `B`, every endpoint is rational, so `S_B(w)` is periodic in
`w`.  The periods are not small (`420..17640`), which explains why the mod-7
table fails, but one complete endpoint-period scan is still finite and exact.

The strong signal is the numerator bound:

```text
max positive S_B(w) <= 1289/980 < 1.316.
```

Since the wide branch has `w>=15`, this gives `Delta_w <= maxS/15`, which is
below `0.453` of the cap margin in every consecutive binding row.  The true
maxima are even smaller and match the low-head scout rows from HYP-2786.

This narrows the remaining single-far proof under `THM-563`: attach a
period-max/slack ledger to the non-consecutive bounded bases that survive the
HYP-2788 single-perturbation reduction.  The exact consecutive-base branch now
has a finite certificate, so further effort should not chase a loose BV or
infinite-tail bound there.

After incoming `THM-563` and HYP-2788, this is best viewed as an addendum and
regression certificate under the wide-region architecture.  It extends the
stored period table through k=11,12 and records the row-pressure tournament;
the main proof route is now `HYP-2788` near-cap/genuine-wide reduction plus
finite `THM-563` period-max checks.  The latest incoming S6 period-max output
pushes that finite check through top dangerous k=8/k=9 rows, and the continuous
period-max addendum closes dilated consecutive bases.  The next compute target
is any residual non-consecutive bounded-base slack/certification, not the
consecutive table.

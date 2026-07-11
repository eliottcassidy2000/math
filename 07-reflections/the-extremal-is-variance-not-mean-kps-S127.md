# The extremal is a variance, and the k=8,9 boundary is one crossover

*kind-pasteur-2026-07-11-S127. Owner: "keep chipping at the extremal." The chip: the consec-extremality
that both proof routes bottom out on is a **variance-maximization**, and the exact reason it holds at k=8,9
(where the ladder needs it) and fails at k≥10 (where the ladder doesn't) is a single, legible crossover.*

---

## What everything reduces to

Both routes to LRC(14) — the B5/clean-ruler side and the moment-ladder side — end at one lemma: the
arithmetic progression `{1,…,k}` minimizes `J = E[N(7−N)]` (N = number of empty sectors) over all integer
cores. mac-mini's THM-710 collapses the ladder to needing this only at **k = 8 and k = 9**; opus reformulated
it as consec maximizing a coverage-variance. I went looking for what actually carries the extremality.

## It is not the mean — it is the variance

Write `J = 49/4 − Var(N) − (E[N] − 7/2)²`. Minimizing J is maximizing `Var(N) + (E[N]−7/2)²`. The natural
reading — both terms peak at consec — is wrong, and the correction is the whole point:

- **consec does *not* minimize `E[N]`.** Exhaustively over `[1..14]`, the best average coverer is
  `{2,4,5,6,7,8,10,12,14}` (`E[N] = 1.340 < 1.446`). consec is not the tightest coverer on average.
- **consec *exactly* maximizes `Var(N)`** — confirmed adversarially over large-value and dilated-AP cores,
  and exhaustively for every `k = 8, 9, 10, 11`. This is the fundamental extremal.

So the `(E[N]−7/2)²` term actually works *against* consec (a lower-`E[N]` core beats it there), and consec's
J-minimality is bought entirely by its variance. The extremal is bimodality: consec covers all seven sectors
often *and* fails deeply often, and it is that spread — not good average coverage — that minimizes J.

## Why k=8,9 and not k≥10 — one crossover

The cleanest thing the decomposition buys is an explanation of the base's shape. `J(E) − J(consec)` is
`consec`'s variance lead minus the competitor's mean-term gain:

- **k=9:** consec's `Var` lead over the runner-up is `0.594`; the runner-up's `(E[N]−7/2)²` gain is `0.447`.
  Variance wins — J-min holds.
- **k=10:** the low-`E[N]` core's mean gain is `0.462`; consec's variance lead is only `0.435`. The mean term
  wins — **J-min flips.**

`Var`-maximization holds at every k; J-minimization holds only while consec's variance lead outruns the
mean-deviation gain, and that race is lost exactly at k=10. And k=10 is precisely where mac-mini's
eigen-transfer lets the ladder *inherit* rather than check — so the functional stops being consec-extremal
exactly one step after the last row that needs it to be. The base is `{k=8} + {k=9}` not by convention but
because that is the last k where variance dominates.

## What this is, and is not

It is a reformulation with an exact mechanism, not a proof. `Var(N)`-maximization is itself LRC-hard — it says
the arithmetic progression maximizes the variance of sector-occupancy, which is the three-distance rigidity
of `{jx}` wearing a statistical costume. But it sharpens the target three ways: the quantity to prove extremal
is `Var(N)` (clean, single, robust across all k), not the composite J; the mean term is a *contaminant*, not a
co-conspirator; and the k=8,9-only base is explained rather than stipulated. Whoever proves this should aim at
the variance, and only afterward spend the k-specific margin (`+0.147` at k=9) to recover the J-bound the
ladder needs.

The recurring shape holds once more: the object that looked like a two-part maximum was one quantity —
variance — plus a term pulling the other way, and seeing that is what tells you where the real extremal lives.

*Files: HYP-6015; the k=9 exhaustive/adversarial data (`lrc14_k9_base_exhaustive_kps_S127`, cont.32).
Sharpens opus-S221/S222 (coverage-variance / longest-AP), builds on mac-mini THM-710/711. Extends
[[the-anchor-selects-the-invariant-coarse-freiman-kps-S127]].*

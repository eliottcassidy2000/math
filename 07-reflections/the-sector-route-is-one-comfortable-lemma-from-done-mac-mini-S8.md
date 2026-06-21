# The LRC(14) sector route is assembled — one comfortable lemma from done

**Session:** mac-mini-2026-06-20-S8 (a long overnight push toward the proof, with kps, codex, and
opus running concurrently). The result is not a proof, but a sharp *localization*: the entire sector
route now sits on a single well-posed analytic lemma with comfortable margin.

## The route, as it now stands

The sector route closes if `measS7(E) ≤ cap_k` for every cluster offset set `E` — equivalently, if
the consecutive block `consec_k` maximizes the Z/7 cover `measS7` (verified globally: 0 of 4000
random gcd-1 shapes of span ≤ 100 beat it) and `measS7(consec_k) ≤ cap_k` (tabulated). Splitting on
the gap structure of `E`, four regimes are now **closed**, and the budget is comfortable in all of
them:

- **Bounded small span (≤14):** finite exact check; consec is the max. `k=12,13` are bounded-only
  (their iid surjection rate already exceeds the cap, so no wide danger exists).
- **Single far point:** the *universal wide supremum* (kps-S23: the max cover is monotone
  decreasing in the number of far points, so one far point is the worst case). It is closed by a
  comb bound `|Δ_w| ≤ 2c₁(E')/(7w)` with `c₁(bounded base) = O(m)`, giving a per-base cutoff
  `W* ≤ 48` and a finite window with zero violations and margins ≥ 0.12. The decoupled limit is
  `measS7(B) + (1/7)·m₁(B)`, worst base `consec_{k−1}`, with margins 0.13–0.18 — comfortable, *not*
  the 0.01 the floor gate showed; the direct cover route avoids the tight floor.
- **Multiple far points / fully dissociated:** strictly easier (cover decreases with each far
  point), bottoming at the iid surjection rate `7!S(k,7)/7^k ≪ cap`; opus-S7 gives explicit
  thresholds via a Beurling constant.

## The one remaining regime

What is left is the **moderate-span balanced** case: gcd-1 shapes of span between ~15 and a few
hundred, with no gap large enough to peel, and not a dilation of a small shape. Here the single-far
comb peel is too weak (the relevant constant grows like the span, not the cardinality), and the
finite family — while genuinely finite by a Hadamard bound — is far too large to enumerate. This is
the multi-block regime, and the tool for it is the **carrier-product** bound: well-separated
sub-clusters share only the slow coordinate, so their colorings decouple, and (codex's Route E,
"splitting strictly costs cover") the single block dominates. What remains is to make the finite-
separation carrier error rigorous — a multi-dimensional Weyl / Erdős–Turán estimate.

And here is the encouraging part, the thing this session actually nailed down: **the budget is
comfortable.** The worst balanced moderate-span shape found (over thousands sampled) covers only
0.197, against a cap of 0.381 — a slack of 0.184. The whole problem has exactly **one** genuinely
tight point, and it is `consec_k` itself (slack 0.023, a single tabulated rational). Everywhere
else there is room to spare. So the remaining lemma does not need to be sharp; a crude multi-block
Weyl bound with a budget of 0.18 will close it.

## What changed

A month ago "consec maximizes the cover" was an opaque, irreducibly-aggregate wall — every
sub-structure (per-block, per-band, monotone, majorization, convexity, synchronization,
interference) was refuted. It still is irreducible *as a single global statement*. But the gap
structure cuts the global statement into regimes, four of which fall to peeling and decoupling with
comfortable margins, and the wall has retreated to one regime with one tool and a wide budget. The
lonely runner at fourteen is not proved. But for the first time the open content is a single
analytic lemma — "the finite-separation multi-block carrier error is below 0.18" — and not the
conjecture itself. That is the difference between a frontier and a finish line in sight.

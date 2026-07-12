# The k=13 Farey window is non-empty — but it holds no divisor-complete family

*kind-pasteur-2026-07-11-S127 cont.51. Owner: "spend another similar session of progress." Testing whether
the empty M-window (1/14, 2/27) is free of divisor-complete families — the same hour opus-S247 refuted the
window being empty at all — the two results fit together exactly and restore the two-bucket closing.*

---

## The two results, same hour

opus-S246 had proposed the LRC(14) crux as a **Farey-window rigidity**: the M-spectrum window `(1/14, 2/27)`
is empty, so `M < 2/27 ⟹` the dilated interval `{1..13}`. opus-S247 **refuted** it at k=13: the family
`{1,…,11,13,36}` (= `{1..13} ∖ {12}` plus `36`) has `M = 3/41 = 0.0732 ∈ (1/14, 2/27)`. The window is
**non-empty**; the second-rung gap that exists for k=12 (mod-13, a field) closes up for k=13 (mod-14,
composite, zero-divisor at 7).

I was testing the *same* window from the other side — restricted to **divisor-complete** families:

> **No divisor-complete family lies in `(1/14, 2/27)`.** Over the structural cores and an adversarial
> hill-climb (Vmax ≤ 32, near-tight seeds), the DC M-floor is **1/12 = 0.0833** (at `{1,2,3,4,10..18}`,
> matching my cont.41 exhaustive census), well above `2/27 = 0.0741`. Zero DC families in the window.

## Why they fit — and restore the closing

The resolution is immediate once you check opus-S247's witness: `{1,…,11,13,36}` **is not
divisor-complete** — it has no multiple of `14`. So it lives in the **non-DC bucket**, where THM-366's
two-bucket dispatch gives `t = 1/14 ⟹ M ≥ 1/14` directly (every speed avoids `0 mod 14`). It is lonely for
the elementary reason, not a hard case at all.

So the window's non-emptiness (opus-S247) and its DC-emptiness (this session) combine cleanly:

- The families filling `(1/14, 2/27)` are all **non-DC near-AP** families (like `{1..11,13,36}`) — dispatched
  by THM-366 (`M ≥ 1/d ≥ 1/14`).
- No **DC** family is anywhere near the tight bound — the DC floor is `1/12`, a `1/84` margin.

Therefore LRC(14) closes as **[non-DC: THM-366, `M ≥ 1/14`, including every window-filler] + [DC: `M ≥ 1/12`]**,
and opus-S247's non-empty window does **not** obstruct it. The correction breaks the *single-statement*
Farey-window rigidity (which tried to prove LRC(14) in one shot), but not the **two-bucket route** — because
the window-fillers, the very families that refute the rigidity, are non-DC and already handled. The loneliness
tightness lives entirely in the non-DC bucket; the DC hard core is loose.

## What this sharpens

opus-S247's own conclusion was "the crux is *AP minimizes M*, not window-empty." This session sharpens *where*
that crux bites: not on the DC families (they are loose, floor `1/12`), but on the **non-DC near-AP** families
— which are exactly THM-366's domain. So the closing does **not** need a new k=13 rigidity theorem; it needs
the two pieces the fleet already has in flight: THM-366 (proved) for non-DC, and the DC `M ≥ 1/12` floor
(boxeph-S20's finite check through Vmax ≤ 30, plus my cont.50 dilation transport carrying the near-dilate
orbit). The `14 = 2·7` composite difficulty that killed the window-empty form lives in the non-DC bucket,
where the *elementary* dispatch already resolves it.

## Honest scope

The DC-empty-window is verified (adversarial hill-climb Vmax ≤ 32 + cont.41 exhaustive Vmax ≤ 22 + the
structural cores), not proved uniformly — the DC floor `1/12` is a bounded-structure fact, carried to all
scales by dilation (MISTAKE-140 / my cont.50 transport), with the large-structure DC families looser still
(cont.49: `M ≈ 0.35`). The full `DC ⟹ M ≥ 1/12` is boxeph-S20's finite check plus the large-diameter
looseness. What this session adds is the timely complement to opus-S247: the non-empty window is a non-DC
phenomenon, so the correction costs the single-statement route but not the two-bucket one.

*Files: lrc14_dc_empty_window_kps_S127.py/out. HYP-6175. Complements opus-S247 (window non-empty at k=13,
HYP-6165) and opus-S246 (E3/Farey); rests on THM-366 (two-bucket dispatch), boxeph-S20 (DC finite check),
[[dilation-preserves-looseness-formalized-the-near-dilate-slice-reduces-to-cores-kps-S127]] (transport). Also
resolved two HYP-number collisions this session (kps-cont49 keeps 6140; opus-S247 keeps 6165; klein-S265 → 6170).*

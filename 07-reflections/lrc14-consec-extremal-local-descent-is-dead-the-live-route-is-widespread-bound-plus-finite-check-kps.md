# Pruning the proof tree for "consec is extremal" (LRC(14)): every local monotone descent is dead — the live route is a non-monotone wide-spread bound + the (already-done) finite check

**Source:** kind-pasteur-2026-06-18-S8. Dispatch: long creative session at an LRC(14) proof,
push/pull often, draw on other agents' work. The thread has converged (THM-530/534/536/537) to
ONE crux on three equivalent functionals: **consec `{0,…,k−1}` is extremal** for k=8..12 —
minimizes `μ_{1/7}(E)` / its lower bound `EWLB_A(E)` (THM-530/531), minimizes the moment
`S_1(E)` and maximizes the dual `L_y(E)` (THM-534), maximizes `meas(S7)` (HYP-2604). Verified
exhaustively, zero counterexamples; the gap is *symbolic*. This session **tested which proof
tactics can close it** — and pruned the natural ones.

## What I tested (exact rational engines) and found

I built an exact `EWLB_A` engine (staggered-window sweep; reproduces `EWLB(consec_8)=407/588`)
and an exact single-sector `S_1` engine, and stress-tested the two natural "consec is extremal"
proof tactics:

1. **Gap-contraction monotonicity (THM-530-D's flagged route): DEAD.** The move "shrink a
   sorted-offset gap by 1" (toward consec) **increases** the functional in **164/386** cases for
   `EWLB`, and **173/431** for `S_1`. Example: `(0,2,4,5,6,7,8,10)` has `EWLB=0.754`, but
   contracting → `(0,2,3,4,5,6,7,9)` *raises* it to `0.838`. So although consec is the global
   minimum, **there is no monotone single-gap descent path to it.**
2. **Stranger pull-in (large-gap / block-merge): DEAD.** For `E={0..6}∪{w}`, `EWLB(w)` oscillates
   (`w=7:0.692, 8:0.785, 9:0.777, 10:0.808, 14:0.765, 21:0.779, …`) — **not monotone** as the
   stranger comes in. The local dips sit at **resonant `w` (multiples of 7, and 1/7-arc
   resonances)**, not at smaller spread.
3. **Consec genuinely minimizes `S_1`** (0/150 below) — the cleaner moment is also extremized at
   consec — **but `S_1` is non-monotone under contraction too** (173/431). So even the simplest
   moment offers no monotone handle.

## The takeaway — the extremality is irreducibly global, on EVERY functional

mac-mini's Angle B (THM-536) found `meas(S7)`-extremality is "irreducibly aggregate" (no
per-IE-block / span-monotone / translation-symmetric proof). This session **extends that to
`EWLB`, `μ_{1/7}`'s lower bound, and the moment `S_1`**: across all of them, consec is the global
minimum but **no local monotone move (gap-contraction, stranger pull-in) descends to it**. The
class of proofs by *rearrangement / monotone reduction / one-step compression* is therefore
**ruled out** — confirming and sharpening that the closing argument must be global. (This does
not refute THM-530-D's *stratum-minimum* monotonicity claim — a statement about
`min{μ : maxgap=G}` — but it shows that claim cannot be realized by the obvious per-config
contraction tactic, which is what one would reach for.)

## The live route (the positive redirect)

The same data points the way. **Large-spread `EWLB` stays comfortably ABOVE consec**: for
`{0..6}∪{w}`, `EWLB(w)∈[0.76,0.83]` for all `w≥8`, vs `0.692` at consec — oscillating but never
dipping to the consec value. So the difficulty is *not* at large spread; it is concentrated in
the **bounded-spread region** (where the finite exhaustive check is already done, k=8..11). The
viable proof structure is therefore:

> **(a) a NON-MONOTONE uniform wide-spread lower bound** — "spread `> B(k)` ⟹ `EWLB(E)` (resp.
> `L_y(E)`) is on the safe side of consec/cap" — proved *not* by monotone descent but by a
> one-shot decoupling/fluctuation estimate (the orbit of a wide-spread `E` is near-equidistributed,
> pushing the functional toward its independent value, which the data shows is `≈0.78 > 0.692`);
> **plus (b) the bounded-spread finite check** (already exhaustive, zero counterexamples).

This **replaces the dead contraction lemma** as the precise remaining target (HYP-2608): the
wide-spread bound is a Weyl-equidistribution / discrepancy estimate, not a rearrangement. The
repo currently has *no* wide-spread decoupling bound (the Explore map flagged this gap), so it is
the clean next objective — and it is a standard-flavored analytic estimate, unlike the
combinatorial monotonicity that just died.

## The 7-adic resonance reading

The non-monotone landscape's local minima sit at **`w ≡ 0 mod 7` and 1/7-arc resonances** — the
difficulty is the **apex-prime-7 resonance structure** (the same `7` of `14=2·7`, the
`s(7|t)=0` vanishing, THM-503), not the spread. A wide-spread bound must therefore beat the
*resonant* wide configs specifically; the generic (Sidon-like) wide configs are far above consec.
This is why the bound is non-monotone: resonances create dips, but (the data says) never below
consec.

## Status / honesty

This session is **route-validity pruning + a redirect**, not a proof. PROVED/VERIFIED (exact):
`EWLB(consec_8)=407/588`; consec minimizes `EWLB` and `S_1` (sample); the gap-contraction and
stranger-pull-in tactics are **non-monotone** (164/386, 173/431, oscillating) — so the
monotone-descent class of proofs is dead; large-spread `EWLB` stays above consec. CONJECTURAL
(HYP-2608): the wide-spread lower bound + finite check completes "consec extremal." This neither
proves nor disproves LRC(14); it prunes the proof tree and sharpens the live target. Cross-links:
THM-530/531 (EWLB, μ_{1/7}, AP-invariance), THM-534 (moment-LP, the S_1/L_y functionals),
THM-536 (Sturmian "irreducibly aggregate"), THM-503 (the apex-prime-7), HYP-2602/2607 (the crux),
HYP-2608 (the redirect). Scripts: `04-computation/lrc14_{ewlb_contraction_test,largegap_stranger_test}_kps.py` (+.out).

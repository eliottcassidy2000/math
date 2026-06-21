# consec-max is the wrong wall: it is false for k≥12 and the proof never needed it

**Session:** mac-mini-2026-06-21-S20 (high-throughput trawl + drill, concurrent with kps/codex).
Builds on the S20 creative trawl and the 5-route LAYER-3 drill (workflows wiarhvx3i, wkkl2szpf).

## The thing everyone was proving

For many sessions the LRC(14) sector route has had one acknowledged open core: **gap #4 / LAYER 3 —
"the consecutive block `consec={0,…,k-1}` maximizes the Z/7 cover `measS7` over the full-residue
stratum."** kps, codex, and I have all attacked it: Delsarte LP, CJJ hierarchy, Lovász θ′, dilation
symmetry, conductance bottleneck, windows/binding-speed, Huffer-Shepp, Markoff, quasimodular. Every
route hit the same "genuinely aggregate" wall. The S20 drill localized it beautifully — and then
ROUTE 4 knocked the premise out.

## consec is not the maximizer

At **k=12** the shape `E* = [0,1,2,3,4,5,6,7,8,9,10,12]` — same balanced full-residue profile as
consec — **strictly beats** consec on `measS7` (consec ≈ 0.624, E* ≈ 0.645; verified two ways,
span-robust). The full boundary: consec-maximizes-`measS7` holds at k=8,9,10,11, **fails at k=12**
(isolated), recovers at 13,14,15, and **fails for all k≥16** with the beater-count growing
(12, 79, 233, …; these are the two-cluster / additive-energy shapes of HYP-2738). So the universal
statement "consec maximizes measS7" is simply **false**. Years of session-effort were aimed at a
theorem that isn't true.

## …and it never mattered

The LRC(14) reduction does **not** require consec to be the maximizer. It requires
`max_E measS7(E) ≤ cap_k`. And the cap holds with **growing** margin exactly where consec stops
winning:

| k  | measS7 (max) | cap_k        | margin |
|----|--------------|--------------|--------|
| 10 | 0.5045 (consec) | 55/91 = 0.6044 | **0.100** (binding) |
| 11 | 0.5815 (consec) | 66/91 = 0.7253 | 0.144 |
| 12 | 0.6452 (**E\***)  | 6/7  = 0.8571 | 0.212 |

The binding row is **k=10**, where consec *is* the max and sits 0.10 below cap (verified exhaustively
over the stratum). For k≥12 the cap has ≥0.21 of slack — a *crude* upper bound clears it; the fact
that some non-consec shape is the true maximizer there is irrelevant. consec-max was a **proof
strategy**, sufficient-if-true at the bounded rows, and it is neither necessary nor (beyond k=11)
true.

## The corrected architecture

The real decomposition was there all along: `measS7 = P2 + R`, with `P2` the decorrelated plateau
(Delsarte/moment-LP, `measS7 ≤ L_y ≤ cap`, proved and Lean-formalized, audit-verified over all 11432
k=8 shapes) and `R` the resonance correction (`≤ (20/7)/p`, sharpened this session, HYP-2750). Under
that lens:

1. **Bounded binding rows (k ≤ 11):** prove `max measS7 ≤ cap_k`. The Delsarte dual `L_y ≤ cap`
   (a single Krawtchouk-nonnegative certificate per row) bounds *all* shapes at once — it does **not**
   need consec to be the argmax. This is the place to spend rigor, and it is finite.
2. **Large rows (k ≥ 12):** the cap slack (≥0.21) admits a crude `measS7 ≤ cap` bound. No aggregate
   Schur argument, no consec-max.

So the "genuinely aggregate wall" was an artifact of asking for the wrong thing — the *unique
extremal identity of the maximizer* — when the proof only ever needed an *upper bound that the cap
clears*. The Delsarte LP supplies exactly that.

## The lesson

The clean-looking extremal statement ("the most-linear object wins") was seductive and the whole team
optimized toward it. It is a real phenomenon at small k and a false one at large k, and in neither
regime is it what the proof consumes. When a sub-goal resists every tool *and* keeps being restated as
"genuinely aggregate," it is worth re-deriving what the parent actually requires — here, a one-sided
cap with slack, not a two-sided extremal characterization. The honest move (the same one the L7 audit
taught) is to verify the load-bearing claim, not the elegant one. The LRC(14) sector route is closer
than the consec-max framing made it look: a finite Delsarte obligation at k≤11 plus a slack bound
above, not an infinite Schur-maximization.

(Caveat reconciled: the drill's WIN-decomposition code and my breakpoint code give k=12 measS7 values
differing by ~0.006. An independent audit (5 methods, incl. boundary-free prime grids) confirmed **my
breakpoint values are correct**; the WIN code silently drops the j=0/sector-0 strip (length 1/7,
mass-zero for consec at k≤11, positive from k=12) and breaks dilation-invariance. E*>consec and all
cap-margins stand.)

## AUDITED CORRECTION (same session) — the binding row is k=8, the load-bearing point is k=9 Delsarte

An adversarial audit of this very reframe corrected the title's "k10" (and my HYP-2778) on two counts;
recorded honestly per the L7-audit discipline:

1. **The binding cap row is k=8, not k=10.** Cap-margins (cap − max measS7), exhaustive + span-robust:
   +0.054 (k=8), +0.078 (k=9), +0.100 (k=10), +0.144 (k=11), +0.212 (k=12), +0.318 (k=13). k=8 is
   tightest. My "k=10" came from conflating this with the audit's earlier cap-vs-plateau margin.
2. **The genuinely load-bearing obligation is the Delsarte leg, tightest at k=9** — margin
   cap_9 − L_y(consec_9) = 10441/7567560 ≈ **0.00138**, razor-thin. That, not k=10 and not consec-max,
   is the delicate point.
3. **The k-range is finite {8,…,13}** for LRC(14) (|P|+|E|=13, so k≤13; k≤7 pigeonhole; k=13 trivial).
   The "k≥16 beater-growth" was about the abstract measS7 problem, not LRC(14).

So the corrected architecture stands and is sound (no cap violation, all cap values re-derived), but
its weakest links are the **k=9 Delsarte razor margin** and the **unbounded-span sup** (the multi-carrier
joint-decorrelation, HYP-2684/2694 — the real open analytic residual), not the consec-max wall. See
HYP-2781 for the full audited map.

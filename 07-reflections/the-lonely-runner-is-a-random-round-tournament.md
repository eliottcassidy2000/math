# The Lonely Runner is a question about random round tournaments

*kind-pasteur-2026-06-22-S36. Synthesizing the maturation of the LRC(14) understanding with the
project's core object — and naming the abstract thing both are shadows of.*

## What we have been sharpening toward

Trace the LRC(14) arc by the canon theorems it deposited:

- **THM-501–515** — the lonely *measure* `L(S)` as a singular series / theta function. First object.
- **THM-523** — `L = 0` is necessary but not sufficient; the real quantity is the *gap*
  `M(S) = max_τ min_v ‖vτ‖`, and the `τ = 1/q` witnesses reduce LRC(14) to **covering sets**.
- **THM-527** — the obstruction is a gap **width**, not a gap measure.
- **THM-528/531** — **scale-decouple** the small part `G_P` from the cluster; reduce to a
  three-distance extremal statement.
- **THM-532–538** — the **seven-sector cover**; `meas(S7) = P(N_E = 0)`, `N_E` = number of empty
  sectors; the moment-LP dual `L_y(E) = E[g(N_E)]` with integer-root Bonferroni `g`.
- **THM-546–557** — far-element **decorrelation**; the floor is a quasi-independence ratio.

Two sessions ago I added the floor side (the spectrum sum `SPEC = Σ ĉ(n)ĝ(n) = ` deviation from
independence, low-Farey dominated, complement-symmetric). This session, three things converged.

## The abstract object, named

**The empty-sector count `N_E(x) ∈ {0,…,6}`** is what the whole reduction is about. The cap side
(THM-534) is its `g`-moment `E[g(N_E)]`; the floor side (the spectrum) is its decorrelation from
`G_P`. And `N_E(x)` is not a bare integer — by **HYP-2605** it is a *coarsening of the score
distribution of a tournament*:

> At phase `x`, the cluster offsets give `k` points `p_i(x) = frac(e_i x)` on the circle. The
> difference-winding map `i → j ⟺ frac((e_i − e_j)x) ∈ (0, ½)` is (a.e.) a **round/local
> tournament `T(x)`**. Then (R3, verified k ≤ 10)
> **`μ_{1/7}(E) = P_x[ T(x) has a scale-1/7 local sink ]`**,
> and (R4) `c₃(T(x)) = C(k,3) − Σ_i C(s_i, 2)` — the project's 3-cycle count, exactly.

So the Lonely Runner gap problem, stripped to its hard core, is a question about a **random round
tournament `T(x)`** drawn by spinning the phase `x`: *does a positive fraction of phases leave the
tournament with a scale-1/7 local sink (an empty 1/7-arc)?* The lonely point exists iff the winding
tournament is "sink-prone" on a positive-measure set of phases.

This is the same object the project has studied all along, viewed through the other window:
- the project asks for `T`'s **Hamiltonian-path count** `H(T) = I(Ω(T), 2)` (Rédei, OCF);
- the LRC asks for `T(x)`'s **local-sink probability** `μ_{1/7}`.

And **both are extremized by the regular object.** HYP-2605: the AP cluster `{0,…,k−1}` maximizes
`E_x[H(T(x))]` (exact global max at k=7). THM-534: the AP (consec) maximizes the sector moment `L_y`
(for k ≤ 11). Just as Paley / self-complementary tournaments maximize `H`. The same principle —
**regularity extremizes the inclusion-exclusion count of non-conflicting configurations** — runs
through both domains, and now we can see why it is the *same* principle: the LRC cluster's winding
tournament and the project's tournament are the same kind of animal, and the AP is its "Paley."

## Why the AP wins (this session's mechanism)

The remaining crux is "consec maximizes `L_y = E[g(N_E)]`." I pinned the order it lives in: `g` is
large at the **extremes** `g(0)=g(6)=1` and tiny in the middle `g(3)=1/10`, so `L_y` rewards a
**bimodal** `N_E`. Computed: consec's `N_E` puts mass 0.347 on `{0,6}`; spread clusters put 0.03–0.05
there and 0.65–0.68 in the middle. **The AP wins because its three-distance rigidity makes coverage
all-or-nothing** — as `x` moves, the AP orbit either spreads to hit every sector (`N=0`) or collapses
onto a low-denominator grid leaving a long empty run (`N` large), with little partial-coverage time.
A floppy (non-AP) cluster dwells in the middle states `g` punishes. In tournament language: the AP's
winding ensemble is the most "either-strongly-connected-or-has-a-clean-sink," the least "muddled" —
the bimodality is the tournament analogue of the score sequence being extremal.

(Honest scope: consec-max-`L_y` is true only k ≤ 11; k=12,13 hold by cap margin, not extremality —
the canon correction. And `meas(S7)` co-moves with `E[H]` only as a strong trend, 36% discordant
pairs — so the tournament dictionary is an exact event-description, R2/R3/R4, but cyclicity is a
*correlate* of sector-fill, not its equal. The bridge is structural, not a numeric identity.)

## The pointer beyond

If `N_E` is the score-coarsening of a random round tournament and the AP is its Paley, then the right
tool for the crux is the one the project already has for `H`-extremality: a **rigidity/regularity
argument on the score distribution**, not a term-by-term moment inequality (which provably fails —
the dual terms trade off). Concretely: the AP's `T(x)` has the most *rigid* score sequence (three
distinct gap lengths ⟹ few distinct scores), and rigidity is exactly what forces bimodal `N_E`. The
conjecture to chase: **`N_E(consec)` dominates every `N_E(E)` in the "extreme-favoring" order, because
the consecutive winding tournament is score-rigid (three-distance) where every other cluster's is
score-floppy** — the LRC's last lemma is a statement about which random tournament ensemble is most
concentrated on its strongly-connected / clean-sink extremes. The Lonely Runner has been a tournament
problem all along; we were measuring its Hamiltonian content from the wrong side.

— Related: [[lrc14-thread]], HYP-2605 (the winding-tournament dictionary), THM-534 (the moment dual),
HYP-2607 (consec-max-`L_y`, scoped k≤11), HYP-2869 (the bimodality mechanism), HYP-2868 (relative
independence), `eight-roads-one-crux-the-empty-sector-moment-functional.md`; THM-001/002 (Rédei, OCF),
the SC-maximizer dichotomy, `everything-is-the-triangle.md`.

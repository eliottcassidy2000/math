# Two routes, one ladder — the B5 obligation is the seven-sector residue, one scale finer

*kind-pasteur-2026-07-11-S127. Owner: "work on finishing it up." I went to the Lean finish line to see what
"finishing" actually requires, and found that the obligation the machine is waiting on is the same object I
had been analyzing at the coarse scale — and that it composes with my own peel into a single proof.*

---

## What the machine is actually waiting for

The Lean development (`LRC14CompletionAudit`) is blunt about the state: LRC(14) is fully formalized,
kernel-sound, sorry-free, foundational-axioms-only, **modulo** the LRC(≤13) citation and **one** obligation,
`hB5`:

> every residual covering family has some modulus `q` with `B5(v,q) > 0`,

where `B5 = S₀ − S₁ + S₂ − S₃ + S₄ − S₅` and `S_d = Σ_p C(bandCount(v,q,p), d)`. Here `bandCount(v,q,p)` is
the number of runners whose `(v_i p mod q)` misses the safe band `[q/14, 13q/14]` — i.e. the number of
runners *not* lonely at the rational time `p/q`. `THM-671` gives `B5 ≤ liveCount(q)`, so `B5 > 0` exhibits a
live multiplier, hence a loneliness witness.

Stare at that for a moment. `bandCount(v,q,p)` is the **discrete empty-count `N`** — the same `N` from the
seven-sector residue — but at the fine scale `1/14` and over a discrete ruler `q` instead of the continuous
circle. And `B5 = Σ (−1)^d S_d` is the **alternating factorial-moment Bonferroni** estimate of the
lonely-multiplier count — which is *exactly* the moment-ladder majorant of `Φ = p₀ + (1/3)p₁` that mac-mini
(THM-703) and opus (S220) built at the coarse scale. The Lean obligation and the analytic residue the whole
fleet has been circling are **the same mathematics at two resolutions**: coarse `1/7` continuous (the
recursion's home) and fine `1/14` discrete (the certificate's home).

## The binding case is the same, and the ruler is a pair-sum

If they are one ladder, the extremal (hardest) family should be the same. It is. Adversarially minimizing
`max_q B5` over residual covering families drives straight to the **near-AP** `{1,…,12, 26}`, with floor
`B5 = 2` — and the winning modulus is `q = 27 = 1 + 26`, a **pair-sum ruler**, exactly the shape the audit
names. General bounded residuals sit far from the floor (`B5 ≥ 12`); only the near-AP shapes are tight. The
same AP-extremality that governs the coarse residue (opus-S222's longest-AP monotonicity) governs the fine
certificate. One extremal family, two scales.

The depth is legible too: for the binding near-AP, `B_3 = 0` — depth-3 Bonferroni *fails* to certify — while
`B_4 = B_5 = liveCount` is already exact. Depth-5 works because the surviving multipliers' bandCounts
concentrate on `{0,…,4}`; the required depth is just the max bandCount at the good ruler. That is why the
Lean obligation is stated at depth 5 and not 3.

## The finishing insight: the two routes compose

I had been treating the seven-sector recursion (my THM-701) and the B5 certificate (the Lean obligation) as
*alternative* routes to LRC(14). They are not. They are **complementary halves of one proof**:

- **THM-701's peel** handles the **unbounded** direction — a family with a far element loses it, and the
  recursion descends to bounded cores.
- **`B5 > 0`** discharges the **bounded** residual — the base case the recursion bottoms out on.

And the honest scope-check makes the division of labor precise. The naive explicit ruler `q = M+1` discharges
`{1,…,k−1, M}` only while `M` is bounded; push `M` out to 120 or 200 and `B5` goes *negative* (`−108`,
`−1620`). For a moment that looks like a counterexample to `hB5`. It is not: those families have a far
element `M`, so they are **peeled by THM-701** and never reach the B5 base case. Restricted to genuinely
bounded residuals (`max ≤ 3·second-max`), `B5 > 0` holds with room (floor 12 over 390 families). The peel and
the certificate cover exactly each other's blind spots.

So the a-priori proof of the single remaining Lean obligation has a clean shape:

> `hB5` ⟸ [THM-701: peel far elements down to bounded cores] + [`B5 > 0` on the bounded near-AP binding
> case, via a pair-sum ruler `q`].

Both halves are the same alternating-moment ladder; the first runs it as a recursion in the unbounded
direction, the second as a Bonferroni certificate in the bounded one.

## Why this is the shape of the finish

The recurring lesson of this whole arc — the wall was a decorrelation, the certificate did not factor, the
anchor selected the invariant, E2 was only a shadow — has been that the object keeps being **one structure
wearing different costumes**. Here is the cleanest instance yet: the measure-theoretic recursion I proved and
the discrete arithmetic certificate the machine is waiting on are literally the same alternating factorial-
moment ladder, at `1/7` and at `1/14`. The proof was never going to be two separate arguments meeting in the
middle. It was always going to be one ladder, climbed from both ends — the far elements peeled off the top,
the bounded near-AP certified at the bottom, and the whole thing is the Bonferroni inclusion–exclusion of
"how many runners are not yet lonely."

*Files: `04-computation/lrc14_B5_moment_ladder_kps_S127.py` (+`.out`), `lrc14_B5_adversarial_floor_kps_S127.py`
(+`.out`). HYP-5995. Unifies THM-701 (peel), `hB5` (Lean finish line, `LRC14CompletionAudit`), THM-671
(`B5 ≤ liveCount`, klein-S210), THM-703/opus-S220/S222/S223 (the moment ladder). Extends
[[the-anchor-selects-the-invariant-coarse-freiman-kps-S127]].*

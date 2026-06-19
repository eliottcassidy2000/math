# Refine the test function to buy back the slack

**Source:** mac-mini-2026-06-18-S7. Dispatch: improve the LRC(14) proof routes, creatively,
outside the box. Canon: THM-533 (the finer-cover certificate), plus an 8-angle creative
workflow. Builds on THM-532 (seven-sector relation-height split) and codex's HYP-2603/2604.

## The gap, and where it came from

THM-532 reduced the LRC(14) residual to `meas(S7(E)) ≤ cap_k` — bound the measure of `x`
where the cluster orbit hits all seven fixed `1/7`-sectors. The relation-height split made it
clean: a tiny main term plus a correction graded by a short-relation weight. But the crude
certificate `corr ≤ C·W` *missed*, by a hair: `C*·W(consec_8) = 0.384` against a budget of
`0.357`. The proof stalled three thousandths short.

The reason was not the method. It was the **test set**. The seven fixed sectors are the
*crudest* finite witness that a configuration fails to be a `1/7`-net: a set can hit all
seven sectors and still have a gap of `1/7` straddling a sector boundary. So `meas(S7)` is a
*lossy* over-estimate of the net measure — at the extremal arithmetic progression it sits at
`0.327`, while the true net measure is `0.060`. The whole `0.267` between them is slack that
`S7` throws on the floor, leaving only `0.054` against the cap, and the crude bound needs more
than that.

## The fix: a finer ruler

Replace the seven fixed sectors by `L` equally-spaced arcs of the same length `1/7`. A net
still hits every one of them, so `N ⊆ S_L ⊆ S_7`, and `meas(S_L)` decreases monotonically
toward the true net measure as `L` grows. Each refinement hands back slack: at `L=14` the
extremal AP drops from `0.327` to `0.196` (slack `0.185`), at `L=42` to `0.107` (slack
`0.275`). And the crude certificate, evaluated honestly (true main term, a bank that includes
the perforated shapes that beat the AP-ratio at `L=7`), now **closes** at `L=14` for every
dangerous row, with a comfortable factor-of-two margin. The knife-edge became a corridor.

Two things made the refinement free rather than costly:

- The **weight is the same**. `W(E) = Σ_{triples} 1/height` depends only on the offsets, not
  on the cover; only the Fourier constant `C_L` changes, and it *shrinks* with `L` (finer arcs
  have smaller coefficients). So refining the cover improves the bound without touching the
  hard combinatorial object.
- The **extremal collapsed to a one-liner**. `S7` needed "the AP maximises `meas(S7)`," a
  measure inequality that resisted proof and even failed at `k=12,13`. The finer-cover
  certificate needs only `W(E) ≤ W(consec)`, and that is elementary: distinct sorted integers
  have `e_j − e_i ≥ j − i`, so every triple's height dominates the consecutive triple's, term
  by term. The inequality the measure route fought for is, on the weight, free.

## The pattern, fourth time running

This project's obstructions keep being correctly stated about the wrong object. Width not
measure (S1). Speed scale not invariant (S5). Gap not sector (S6). And now: the *coarseness of
the test set*, not the method. Each time the bound was real and a hair too weak, and each time
the move was not to push harder but to change what is being measured against — here, to ask
the net condition with a finer ruler so the lossy slack comes back. The Beurling–Selberg
extremal-function view (one of this session's angles) is the same idea taken to its limit: the
optimal smooth majorant of the net indicator is the `L → ∞` test function, with the best
possible main-term/tail trade. The finite `L=14` cover already lands inside the corridor; the
smooth limit is the clean rigorous home.

## What is actually left

One analytic bound: `corr_L(E) ≤ C_L · W(E)` with an explicit `C_L`. The Fourier ingredient is
in hand — the single `1/7`-arc has coefficient `|sin(πn/7)|/(π|n|)`, vanishing on `7 | n` (the
same `7` that has run through this problem from the start), and `max_n |\hat s(n)|·|n| =
0.310`. Assembling the support-3 sum and a geometric tail — the HYP-2601 calculation re-run
in the arc basis — is bookkeeping, not a wall, and the comfortable margin means it need not
even be sharp. The residual that was "one arithmetic progression and its neighbours, on a
knife-edge" is now "one Fourier estimate, with room to spare."

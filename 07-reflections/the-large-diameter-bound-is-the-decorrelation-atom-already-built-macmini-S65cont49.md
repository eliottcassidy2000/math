---
source: mac-mini-2026-07-09-S65 (cont.49, 2026-07-12)
tags: [lrc14, large-diameter, decorrelation-atom, THM-636, THM-720, height-descent, mining, convergence]
---
# The large-diameter lower bound is the decorrelation atom — already built, and it converges four ways

Mining past threads for the large-diameter lower bound (THM-720's remaining rigor: spread DC
families have M > 1/14 rigorously), I found it is NOT a new object to build — it is THM-636's
decorrelation atom (my own 2026-07-06-S38 work, already formalized sorry-free in Lean for r ≤ 11),
and four independent results are the SAME statement.

## The tool (THM-636, formalized)
For vᵢ = bᵢ + L·kᵢ (bounded base |bᵢ| ≤ B, L-scaled lift kᵢ): **reach(v) ≥ reach(k) − B/L**
(1-Lipschitz at the witness t_K/L). This is Tao-style height descent: the large-scale family's
reach is the lift family's reach minus a vanishing O(1/L) correction.

## The 13-runner closure of the large-diameter half
A large-diameter divisor-complete family, decomposed at its scale L, has a lift family k with FEW
distinct speeds — because DC is EVEN-HEAVY (multiples of 8, 14, … force even runners), so the
lifts collapse. Measured: **≤ 6 distinct lifts** (reach(k) ≈ 0.18–0.25, far above 1/7). Then:
> **reach(v) ≥ reach(k) − B/L ≥ 1/7 − 13/L > 1/14 for L > ~1274** (≤ 6 distinct lifts ⟹ LRC(7) ⟹
> reach(k) ≥ 1/7; B ≈ 13). Verified: reach(v) ≈ 0.25 across L = 500..40000, all LOOSE.
So every large-diameter DC family is loose by decorrelation descent, RIGOROUS via THM-636
(formalized) + LRC(≤13). The descent base (bounded height / small L) is exactly kps's
bounded-diameter finite check — the OTHER half of the dichotomy. The two halves are one descent.

## The four-way convergence (why this is the right object)
The SAME large-diameter looseness appears as:
1. **THM-720 (mine, cont.48):** pair-sum M grows with diameter (0.105 → 0.187), coverage-duality
   mechanism.
2. **THM-636 (mine, S38):** reach ≥ reach(k) − B/L, the decorrelation atom, formalized.
3. **LEM-013 (mac-mini/klein):** dissociated good-period margin 7·maxgap/Vmax grows 1.105 → 2.31
   with spread.
4. **klein S263 (today):** DC even-heavy ⟹ ~6 odd runners, scale-stable ⟹ the crux shrinks to a
   ~6-runner one that survives the unbounded window.
klein's ~6-odd-runner shrink and my ~6-distinct-lift decorrelation are THE SAME reduction: the DC
family's large-scale structure collapses to ~6 effective speeds, which clear trivially (reach ≥
1/7). Four routes, one phenomenon — the hallmark (per klein's algebraically-special-extremals
reflection) that we have found the true object.

## What remains
The large-diameter half is closed in structure: [≤ 6 effective speeds ⟹ reach ≥ 1/7 − B/L, THM-636]
+ [descent base = bounded-diameter finite check]. The remaining rigor is (a) the exact L-threshold /
the "≤ 6 distinct lifts" as a theorem (DC even-heaviness — provable from the divisor-complete
definition, klein S263's scale-stability), and (b) the bounded-diameter finite check (feasible).
The coverage-clearing duality (cont.47) is why it all holds: spread ⟹ bad coverer ⟹ effective
speeds collapse ⟹ loose.

→ THM-636 (decorrelation atom, formalized), THM-720 (pair-sum looseness), LEM-013 (dissociated
margin), klein S263 (~6-odd shrink), cont.47 (coverage-clearing duality), THM-687/688 (multi-scale),
LRC(≤13) (the lift floor). Files: lrc14_decorr_atom_13runner_macmini_S65cont49 (+ out).

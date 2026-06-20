# The new-speed deviation Δ_w: codex (measure-side) ≡ kps (sector-side) — EXACT BRIDGE

**kps + codex-s45 convergence, 2026-06-20.** Two independent LRC(14) routes are bounding the
SAME object, and localize the SAME extremal family.

## The exact identity
codex's "new-speed constant" `Delta_w^+` (HYP-2671) **= kps's plateau deviation `Δ_w`** (HYP-2653) exactly.
Both equal `1371/4319 · p1(E')` at the shared extremizer `E'=(0,1,2,4,8,12,16,20), w=24` (verified, exact Fraction).
Here `Δ_w = p0(E'∪{w}) − [p0(E') + (1/7)p1(E')]`, the deviation of the 7-sector cover measure from the far-element plateau.

## Two bounds on the same Delta_w
- **codex (relative):** `Δ_w ≤ p1(E')/3`, condition = shell-1-full AND max(E')>14. Extremizer = dyadic block
  `E_m={0,1,2,4,8,3m,4m,5m}, w=6m`, unique spike at m=4 (ratio 1371/4319 = 0.3174 < 1/3).
- **kps (corrected absolute tail, HYP-2653d):** the proof target is
  `sup_{max(E')>B} Δ_w(E',w) ≤ cap_k-Q(k-1)`, with `max(E')<=B`
  handled finitely.  The earlier `w|Δ_w| ≤ C(k)` line was the wrong proof
  currency: along the same dyadic family, `Δ_w` has a small nonzero floor, so
  `w|Δ_w|` grows with scale.  The dyadic family remains the diagnostic
  extremal family, but the theorem must bound `Δ_w` itself after cutoff.

## What each is good for
- codex's relative bound is a CONSTANT in w (works for all w>max in the shell-full regime) — drives the
  measure-side two-gate route. But it NEEDS shell-full: at a non-shell-full far core `[0,2,3,5,6,15], w=18`,
  `Δ_w/p1 = 0.349 > 1/3` (codex's bound would fail; the shell-full hypothesis excludes it).
- kps's corrected absolute tail needs no shell-full precondition — it drives
  the sector-side far-element-peel route by proving
  `sup_{max(E')>B} Δ_w <= cap_k-Q(k-1)` after a finite `max(E')<=B` check.
  The previous bounded-`C(k)` decay sentence is superseded by HYP-2653d.

## The shared extremal family (both routes)
**Dyadic core `{0,1,2,4,8,...}` + a scaled cluster `s·{3,4,5}` + resonant `w` (multiple of s, or =6m).**
This is where `×w` collapses a scaled cluster's breakpoints, killing the Koksma cancellation — the unique
hard case for BOTH the relative and absolute bounds. Proving either bound = certifying this dyadic-resonant
family + showing every other config is strictly dominated.

**Correction note:** read "absolute bound" above in the HYP-2653d sense:
uniform `Delta_w` tail cap, not bounded `w*Delta_w`.

→ HYP-2653, HYP-2653c, HYP-2653d, HYP-2671 (codex), HYP-2644/2645, OPEN-Q-108.

# The new-speed deviation Δ_w: codex (measure-side) ≡ kps (sector-side) — EXACT BRIDGE

**kps + codex-s45 convergence, 2026-06-20.** Two independent LRC(14) routes are bounding the
SAME object, and localize the SAME extremal family.

## The exact identity
codex's "new-speed constant" `Delta_w^+` (HYP-2671) **= kps's plateau deviation `Δ_w`** (HYP-2653) exactly.
Both equal `1371/4319 · p1(E')` at the shared extremizer `E'=(0,1,2,4,8,12,16,20), w=24` (verified, exact Fraction).
Here `Δ_w = p0(E'∪{w}) − [p0(E') + (1/7)p1(E')]`, the deviation of the 7-sector cover measure from the far-element plateau.

## Two bounds on the same Δ_w
- **codex (relative):** `Δ_w ≤ p1(E')/3`, condition = shell-1-full AND max(E')>14. Extremizer = dyadic block
  `E_m={0,1,2,4,8,3m,4m,5m}, w=6m`, unique spike at m=4 (ratio 1371/4319 = 0.3174 < 1/3).
- **kps (absolute):** `w|Δ_w| ≤ C(k) ≈ c·(#scales) ≤ c·k`, UNCONDITIONAL. sup ≈ 2.86 (k=9, 8-elt core),
  extremizer = SAME family `{0,1,2,4,8}+s·{3,4,5}, w=s·10` (e.g. {0,1,2,4,8,84,112,140}, w=280).

## What each is good for
- codex's relative bound is a CONSTANT in w (works for all w>max in the shell-full regime) — drives the
  measure-side two-gate route. But it NEEDS shell-full: at a non-shell-full far core `[0,2,3,5,6,15], w=18`,
  `Δ_w/p1 = 0.349 > 1/3` (codex's bound would fail; the shell-full hypothesis excludes it).
- kps's absolute bound DECAYS (`|Δ_w| ≤ C(k)/w → 0`), needs no shell-full precondition — drives the
  sector-side far-element-peel route. The bounded-C(k) ⟹ peel w≥C(k)/margin ⟹ p0(E)≤Q(k−1)+margin<cap_k.

## The shared extremal family (both routes)
**Dyadic core `{0,1,2,4,8,...}` + a scaled cluster `s·{3,4,5}` + resonant `w` (multiple of s, or =6m).**
This is where `×w` collapses a scaled cluster's breakpoints, killing the Koksma cancellation — the unique
hard case for BOTH the relative and absolute bounds. Proving either bound = certifying this dyadic-resonant
family + showing every other config is strictly dominated.

→ HYP-2653, HYP-2653c, HYP-2671 (codex), HYP-2644/2645, OPEN-Q-108.

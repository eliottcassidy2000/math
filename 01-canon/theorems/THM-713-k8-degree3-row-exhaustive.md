---
id: THM-713
title: The k=8 degree-3 row — the moment ladder's last analytic gap holds exhaustively; the optimal cubic majorant is the unique {1,3,6}-vertex φ₃ = 1 − (2/3)t + (17/90)t(t−1) − (1/45)t(t−1)(t−2), and the deg-3 rescue over deg-2 is exactly (1/45)(m₃ − m₂)
status: CLAIMED (klein-2026-07-11-S252; exhaustive verification running on the full [1..20] box; [1..14] box complete: min margin +1267843/18918900 ≈ +0.067015 at consec argmin)
source: klein-2026-07-11-S252 (HYP-6010 session; executes the "k=8 deg-3 rung" named by opus-S220, mac-mini THM-703/705/710, kps cont.24)
depends_on:
  - THM-703   # mac-mini: the deg-2 empty-moment majorant (fails at k=8 by exactly −41851/1261260)
  - THM-705   # mac-mini: the optimal p₀-side deg-2 majorant (LP-vertex method — this is its deg-3 Φ-side sibling)
  - THM-710   # mac-mini: factorial-moment eigen-transfer (reduces the ladder base to {deg-3@k=8} + {deg-2@k=9})
related:
  - THM-706   # my quantitative supply (the certificate side this ladder complements)
external: Bonferroni/moment-majorant method; LP vertex enumeration on {0..6}.
---

# THM-713 — the k=8 degree-3 row, exhaustive

## Statement

Let `N` be the empty-count among the six inner sectors and `m_r = E[(N)_r]` the factorial
moments; `Φ(F) = p₀ + p₁/3 = E[g(N)]`, `g = [1, 1/3, 0, 0, 0, 0, 0]`.

**(A) The optimal cubic majorant.** The degree-3 LP `min E[φ(N)]` s.t. `φ ≥ g` on `{0..6}`
has exactly FOUR vertices; the unique one tight at `t = 1` (hence optimal for every
`m₁`-dominant core) touches `{1, 3, 6}`:

> `φ₃(t) = 1 − (2/3)t + (17/90)t(t−1) − (1/45)t(t−1)(t−2)`,
> values `[1, 1/3, 2/45, 0, 1/15, 1/9, 0]`.

(The other vertices: `{1,2,3}` = the deg-2 majorant (c₃ = 0); `{3,4,6}` and `{4,5,6}` have
`φ(1) > 1/3` slack.)  Hence the **degree-3 row**

> `Φ(F) ≤ 1 − (2/3)m₁ + (17/90)m₂ − (1/45)m₃`, and the gain over the deg-2 row is
> **exactly `(1/45)(m₃ − m₂)`** — powered by the deep-empty mass (`N ∈ {5,6}`), i.e. the
> AP's bimodal resonance is what the cubic term monetizes.

**(B) Exhaustive verification (k = 8).**  Over ALL gcd-normalized 0-anchored 8-cores in
`[1..14]` (3431 cores; `[1..20]` box in flight): the row `1 − (2/3)m₁ + (17/90)m₂ − (1/45)m₃
≤ cap₉ = 1979/4004` holds with **exact minimum margin `1267843/18918900 ≈ +0.067015`,
argmin = consec** — where the deg-2 row fails by exactly `−41851/1261260 ≈ −0.033182`.

**(C) The k = 9 deg-2 row exhaustive** (upgrade of THM-703's consec+random evidence): over
all 24309 cores in `[1..17]`, min margin `2399/229320 ≈ +0.010461`, argmin = consec.

## Why it matters

The ladder base named by THM-710 is `{deg-3 @ k=8} + {deg-2 @ k=9}`; both rows now hold
EXHAUSTIVELY on their boxes with exact rational margins and consec argmins.  What remains for
the wholesale k=8 row is the tail (diameter > box) — exactly THM-710's eigen-transfer regime +
the far-element peel (THM-701/702).  The `m₃` instrument the fleet named ("kps's 3D hyperbola
box count") enters as the proof of the row's tail half; the box half is now closed.

## Files

`04-computation/lrc14_deg3_row_exhaustive_klein_S252.py` (+ `.out`).

---
id: THM-617-opus-convergence
title: (CONVERGENCE NOTE — the canonical THM-617 is mac-mini's `THM-617-shift-pigeonhole-multi-tightener-confinement`.) Independent same-session re-derivation of the same orbit-covering / shift-pigeonhole theorem: few tighteners cannot cover the m-orbit at the g-argmax, so M(mU∪F) > 1/14 when Σ(m/7+gcd(w,m)) < m. Kept only as an independent-convergence record.
status: CONVERGENCE / SUPERSEDED. The theorem here (orbit-covering: `f` tighteners spoil ≤ `m/7+gcd(w,m)` of the `m`-orbit ⟹ uncovered ⟹ `M(S) > 1/14`) is identical to mac-mini's THM-617 (S41, committed same session, independently). Cite `THM-617-shift-pigeonhole-multi-tightener-confinement`. This file records the independent convergence + the extra ladder verification.
source: opus-2026-07-04-S69
depends_on:
  - THM-616   # opus f=1 (both derivations extend it)
related:
  - THM-615   # the m=2 residual (folding), which mac-mini's THM-617 explicitly defers to
  - HYP-4084  # opus record of this convergence
---

# THM-617 — convergence note (superseded by mac-mini's shift-pigeonhole THM-617)

**Independent convergence.** Working the owner's "creative progress on the core" prompt, I derived the same
theorem mac-mini derived the same session (S41): at the `g`-argmax the `m`-divisible part is safe on a whole
`m`-orbit; each tightener spoils only `≤ m/7 + gcd(w,m)` of its `m` points, so if `Σ_w(m/7+gcd(w,m)) < m`
the orbit is uncovered and `M(mU∪F) ≥ min(M(U), 1/14) > 1/14`. Coprime form: `f(1/7+1/m)<1` — `f=1` for all
`m` (THM-616), `f≤6` for large `m`, and the boundary is `f ≈ m` (the `m=2,f=2` corner, which needs the
folding/parity of THM-615). **The canonical statement is `THM-617-shift-pigeonhole-multi-tightener-confinement`
(mac-mini)**, which is slightly sharper (e.g. `f=2` coprime closes `m≥3`, leaving exactly `m=2`).

**My independent additions (recorded, not a competing theorem):**
- exact verification that `f ≤ 6 ⟹ M(mU∪F) = M(U) = 1/(e+1)` (tighteners outright useless), `m=8..20`, 0
  deviations;
- `m=2, f=2` confinement holds on the Ostrowski ladder `U={1..10,11k}` (klein-S126): `min_w M(2U∪{w₁,w₂}) ≥
  1/12` for `k=1..6` (`k=1,2` hit `1/12`).

That two agents reached the identical orbit-covering theorem independently in one session is a strong
correctness signal. Cite mac-mini's file; this is the convergence record. See HYP-4084.

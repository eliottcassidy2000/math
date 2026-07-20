---
id: THM-1585
title: "THE GAMMA BRIDGE'S DOMINATION STEP IS FALSE, SO GMC(2) IS NOT COMPLETE — and the reason is structural, not a slip. (0) CONTROL FIRST: I re-derived ψ_m independently (klein disclosed a no-op defect in their own code, so none of it is reused) via v = 1 + B(r)t·v + W(r)t²v² with B = b, W = r·a·c — confirming klein's ρ-freeness claim — and my implementation REPRODUCES klein's published numbers exactly on all three of their valid rows. The instrument agrees where theirs is valid. (I) THE CLAIM TESTED, verbatim: 'E_r[ψ_m] = Σc_k k! is dominated by its top term because consecutive terms have ratio c_{D−1}/(c_D·D) → 0'. Measured to m = 20: that ratio does NOT tend to 0 in ANY non-degenerate case — it GROWS LINEARLY in m for a=1,b=1,c=1 (0.5→5.0) and for a=1,b=3,c=1 (4.5→45.0), and is constant ≈0.45 elsewhere. The second term is up to 45× LARGER than the top term. (II) EQUIVALENTLY, the top term's share of E_r[ψ_m] falls to ZERO, not to one: 0.667→0.092 at b=1, and 0.182→0.0004 at b=3 — the toral part is 0.04% of the answer. The single case where the claim holds (b=0) is DEGENERATE: ψ_m is then a monomial in r, so 'the top term dominates' is vacuous. The mechanism was validated only on the case that cannot test it. (III) THE STRUCTURAL REASON: the top coefficient is the SAME UNIVERSAL SEQUENCE in every case, c_D(m=2k) = C(2k,k)/(2k) = 1, 3/2, 10/3, 35/4, 126/5, 77, …, independent of B and of a,c beyond their product — it is the leading symbol, and it is BOUNDED, while E_r[ψ_m] grows without bound with B. A bounded toral island inside an unbounded radial average cannot be dominant, and Gamma weights do not change that. This is mac-mini-S136's and boxeph-S168's ℓ¹-mass objection, confirmed quantitatively. (IV) CONSEQUENCE: TNC ⟹ NC2 is UNPROVED, hence death-star-S61g's 'GMC(2) is complete' is withdrawn-worthy — and the BODY of that same reflection (lines 66–68, 86) already says so, the headline contradicting its own document. NOTHING HERE REFUTES NC2 OR GMC(2): E_r[ψ_m] ≠ 0 in every case tested. What is refuted is the bridge — a hypothesis controlling only the top term cannot imply a statement about a sum in which the top term is 0.04%. (V) REPAIR DIRECTION, offered not claimed: in 4 of 5 non-degenerate cases every c_k shares ONE SIGN, giving Σc_k k! ≠ 0 with no asymptotics at all — strictly better than domination. But it is not general: a=1,b=1,c=−1 is sign-mixed at every m while still having E_r[ψ_m] ≠ 0. Whatever carries the SIGN-MIXED case is the real remaining content of GMC(2)"
status: >
  (0) CONTROL PASSED — independent implementation, agrees with klein's numbers on every
  valid row of theirs before any disagreement is reported.  This is the check that was
  missing in S128c113's broken fibre counter and it is run first here on purpose.
  (I), (II) VERIFIED-EXACT in rational arithmetic, m = 2..20, six cases.  These are
  measurements of the two quantities klein's sentence names, not reinterpretations of it.
  (III) VERIFIED (the c_D sequence is identical across all six cases and matches
  C(2k,k)/(2k) at every computed m); the reading of it as "bounded symbol inside an
  unbounded average" is an explanation, offered as such.
  (IV) This is a NEGATIVE about a proof, not about a theorem.  NC2 and GMC(2) remain
  OPEN and are consistent with all data here.  Filed as
  02-court/active/CASE-gamma-bridge-domination-step.md rather than by editing another
  agent's files.
  (V) OBSERVATION with an explicit counterexample to its own generality.  Not a theorem.
  DETECTION FLOOR: the {−1,0,1} stratum at M = 1 only, six (a,b,c) choices, m ≤ 20,
  exact rationals.  A wider stratum could in principle behave differently, though the
  universality of c_D in (III) argues against it.
source: kind-pasteur-2026-07-20-S128c120 (owner: work to finish the GMC(2) math, then formalize)
disputes:
  - klein-2026-07-20-S351   # the Gamma Bridge domination step (session-log level; no theorem file exists)
  - death-star-2026-07-20-S61g  # "GMC(2) is complete"
depends_on: []
related: [THM-1540, THM-1550, THM-1530, THM-1565, THM-1570, THM-1580, THM-1590]
court: 02-court/active/CASE-gamma-bridge-domination-step.md
script: 04-computation/gamma_bridge_domination_audit_kps_S128c120.py (+ .out)
---

# THM-1585 — the Gamma domination is false, and the toral coefficient is universal

Full argument, evidence table and relief sought are in
[the court case](../../02-court/active/CASE-gamma-bridge-domination-step.md). This file
records the mathematical content.

## The setting, re-derived independently

On the `{−1,0,1}` stratum at `M = 1`, with `P = Z̄·a(r) + b(r) + Z·c(r)`,
`Z = ρe^{iθ}`, `ρ = √r`, `u = e^{iθ}`: writing `R(r,u) = uP = ρa + bu + ρcu²` and
substituting `v := u/(t·g_{−1})` into the small-root equation `u = t·R(r,u)`,

> **`v = 1 + B(r)·t·v + W(r)·t²·v²`**,  `B = b`,  `W = g₁g_{−1} = ρc·ρa = r·a·c`,

which is `ρ`-free — klein's claim, independently confirmed — so `ψ_m := [tᵐ] log v` is a
polynomial in `r`, and `NC2 ⟺ E_r[ψ_m] = 0 ∀m` with `E_r[r^k] = k!`.

Note the stratum depends on `a` and `c` **only through the product `ac`**, which is
consistent with opus-S414's THM-1580 ("the nullcone sees `A, C` only via `h = sAC`")
reached independently on the same day.

**Control before dispute.** My implementation reproduces klein's published values on all
three of their valid rows (`3/2, 25/4, 331/6, 5937/8`; `1, 3, 20, 210`; `15/2, 1841/4,
659287/6, 469665345/8`). Only after that does anything below get reported.

## The measurement

| case | `RATIO = |c_{D−1}/(c_D·D)|`, `m = 2…20` | `SHARE = |c_D·D!|/|E_r[ψ_m]|` |
|---|---|---|
| `a=1, b=1, c=1` | 0.5, 1.0, 1.5, … **5.0** | 0.667 → **0.092** |
| `a=1, b=3, c=1` | 4.5, 9.0, 13.5, … **45.0** | 0.182 → **0.0004** |
| `a=1, b=1+r, c=1+2r` | ≈ 0.45 flat | ≈ 0.637 flat |
| `a=1+r, b=2, c=1` | 0.5 flat | ≈ 0.385 flat |
| `a=1, b=1, c=−1` | 0.5, 1.0, … 5.0 | 2.0 → 0.88 |
| `a=1, b=0, c=1` (degenerate) | 0 | 1.000 |

klein's claim is `RATIO → 0` and `SHARE → 1`. Neither happens anywhere except the
degenerate row, where `ψ_m` is a monomial and the statement is vacuous.

## Why — and it is structural

The top coefficient is **universal**:

> `c_D(m = 2k) = C(2k,k)/(2k)` = `1, 3/2, 10/3, 35/4, 126/5, 77, 1716/7, 6435/8, …`

identical in all six cases. It is the leading symbol, it does not see `B` at all, and it
is bounded. `E_r[ψ_m]`, by contrast, grows without bound as `B` grows — at `b = 3` it is
`1.7 × 10¹⁵` by `m = 20` while `c_D·D!` is `≈ 6.7 × 10¹¹`. So the toral quantity is a
bounded island inside an unbounded radial average. **Gamma weights cannot rescue a
bounded numerator.**

This is exactly mac-mini-S136's and boxeph-S168's objection — the lower degree levels
carry more `ℓ¹` mass than the Newton edge — now measured rather than feared.

## What this does and does not say

- It does **not** refute NC2 or GMC(2). `E_r[ψ_m] ≠ 0` in every case tested, at every `m`.
- It **does** break the inference `TNC ⟹ NC2` as published, because TNC constrains only
  `c_D`, and `c_D·D!` is `0.04%` of `E_r[ψ_m]` in the worst tested case. An implication
  whose hypothesis controls a vanishing fraction of its conclusion's object is not proved
  by that hypothesis, whether or not the conclusion is true.
- Therefore **GMC(2) is open**, and the toral-side advances (boxeph-S173 at `N=1` all `M`
  and `M=N=2`; klein-S357/THM-1590; opus-S414/THM-1580) are individually fine but do
  **not** compose into GMC(2) through this bridge.

## The repair direction

In four of the five non-degenerate cases **every coefficient `c_k` of `ψ_m` has the same
sign**, so `Σ c_k k! ≠ 0` follows immediately — no asymptotics, no domination, no ESV
saddle analysis. That is a much better mechanism than the one claimed, and it is *proved*
wherever it applies.

It is not general. `a=1, b=1, c=−1` is sign-mixed at every `m` (and `E_r[ψ_m] ≠ 0` there
anyway, alternating in sign). So:

> **Named next, and this is the real remaining content of GMC(2):** find the mechanism
> that carries the *sign-mixed* stratum. Sign-coherence handles the rest, and should be
> isolated and proved as a lemma in its own right — including a characterisation of which
> `(a,b,c)` give coherent `ψ_m`, which on this evidence is about the sign of `W = r·a·c`.

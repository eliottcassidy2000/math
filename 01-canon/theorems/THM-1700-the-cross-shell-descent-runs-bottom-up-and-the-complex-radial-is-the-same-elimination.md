---
id: THM-1700
title: "THE CROSS-SHELL DESCENT RUNS BOTTOM-UP, NOT TOP-DOWN — which is why HYP-8470's 'top shell one-sided' is not a real obstruction; and the COMPLEX RADIAL closes by the same integer-moment elimination. FIRST, a scope correction to my own S369: genuine Λ_s obeys the charge-radius LOCK (Z^pZ̄^q ↦ ρ^{p+q}u^{p−q}, charge and radius same parity), so charge balance in CT_u[Λ_s^m] forces EVEN total radius power — integer s-powers, NO √π. My S369 √π family was UNLOCKED. For a GENUINE two-sided P with one-sided top shell — witness P = αZ³ + βZ̄ + γZ (top +3 one-sided at h=3; h=1 straddles ±1) — the moments are E[P²]=2βγ, E[P⁴]=24αβ³+12β²γ², … with INTEGER coefficients, and exact elimination shows NO two-sided nullcone member: the correct one-sided locus is {β=0}∪{α=γ=0}, and (βγ)^k, (βα)^k both lie in the moment ideal. The mechanism: E[P²]=2βγ kills the bottom straddle FIRST, then higher moments force the top — a triangular descent from the low shell up, opposite to DvdK's top-down, so the one-sided top is never the obstruction. The β-dominant COMPLEX RADIAL (THM-1670) closes by the identical elimination: the real/complex branch split governs the analytic flat-term route but the algebraic route is blind to it and closes both."
status: >
  SCOPE CORRECTION to S369/THM-1695: PROVED (the lock is one line; my S369 √π family violates
  it, verified). Genuine Λ_s ⟹ CT_u[Λ_s^m] has only even radius powers ⟹ integer-moment radial
  functional, no √π.
  CROSS-SHELL WITNESS P = αZ³+βZ̄+γZ: CLOSED by exact Gröbner over ℚ (m ≤ 8), with the CORRECT
  one-sided predicate (βγ)^k,(βα)^k ∈ I — an interim run used the WRONG predicate (β^k only) and
  falsely reported an open case; corrected in §2.
  COMPLEX RADIAL: CLOSED by elimination for α=r(1+r) (non-monomial, two-sided), β=b₀+b₁r
  (Gröbner = ⟨1⟩); and for α=r, deg β ≤ 3 (THM-1660).
  The MECHANISM (bottom-up descent) is stated and verified on the witness; the GENERAL HYP-8470
  (many straddling shells, cancellation among charge-0 pairs) is NOT closed here.
  GMC(2) REMAINS OPEN.
source: klein-2026-07-20-S371 (owner: work piece 2 and work the complex radial and the cross-shell descent)
depends_on:
  - THM-1695  # my S369 two-pieces (Piece 1 = EMP, Piece 2 √π = unlocked) — corrected here
  - THM-1670  # the real/complex branch split; the complex radial is that β-dominant case
  - THM-1660  # the fixed-α charge-0 elimination
related:
  - THM-1645  # the polar bridge; the lock is stated there (s^{(a+b)/2}u^{a−b})
  - "death-star-S65: pinch locus γ = r·b'² < 0 ⟹ complex branch (the analytic side)"
  - "opus THM-1675 / mac-mini: the cross-shell / tuned-cancellation threads"
script: 04-computation/genuine_crossshell_klein_S371.py (+ .out)
---

# THM-1700 — the cross-shell descent runs bottom-up

## 0. Scope correction to my own S369 (THM-1695 Piece 2)

A genuine `Λ_s(u) = P(√s·u, √s/u)` obeys the **charge-radius lock**: `Z^pZ̄^q ↦ ρ^{p+q}u^{p−q}`,
so the charge `p−q` and the radius power `p+q` have the **same parity**. Consequently, in
`CT_u[Λ_s^m]` the charge-balance condition `Σ(p_i−q_i)=0` forces the total radius power
`Σ(p_i+q_i)` to be **even** — so `CT_u[Λ_s^m]` is a polynomial in `s = ρ²` (integer powers), and
the radial functional is the **pure exponential** `L(g)=∫g e^{−s}ds`, `L(s^l)=l!` — **no `√π`**.

My S369 Piece 2 family `Λ_s = √s·u + (a/u+b+cu)` has charges `±1` at `ρ^0` (an even shell) — a
lock violation, so it is **not a valid `Λ_s`**. The `√π` transcendence mechanism is a true fact
about the *abstract* (unlocked) nullcone but does **not** address GMC(2). This corrects THM-1695
§2. (Piece 1 = EMP is unaffected.)

## 1. The genuine cross-shell witness

HYP-8470's open configuration — *top shell one-sided, a lower shell straddles* — on a genuine
two-sided `P`:

```text
P = α·Z³ + β·Z̄ + γ·Z.
   charges  +3 (Z³, radius h=3, top shell, ONE-SIDED),  −1 (Z̄, h=1),  +1 (Z, h=1).
   the h=1 shell STRADDLES (±1); two-sided ⟺ β ≠ 0.
```

Exact moments (integer coefficients, odd moments vanish by charge balance):

```text
E[P²] = 2βγ
E[P⁴] = 24αβ³ + 12β²γ²
E[P⁶] = 720αβ⁴γ + 120β³γ³
E[P⁸] = 20160α²β⁶ + 20160αβ⁵γ² + 1680β⁴γ⁴
```

## 2. Closed — with the *correct* one-sided predicate

`P` is one-sided iff its charge support is single-signed:

```text
all positive: β = 0 (charge −1 absent);   all negative: α = γ = 0 (charges +3,+1 absent).
one-sided locus = V(β) ∪ V(α,γ) = V(⟨βα, βγ⟩).
```

So *no two-sided nullcone member* ⟺ `V(I) ⊆` one-sided locus ⟺ `(βγ)^k, (βα)^k ∈ I`. Both hold:
`βγ ∈ I` at `k=1` (it *is* `E[P²]/2`), `βα` at `k=3`. **Closed.**

> **Correction on record.** An interim run tested only `β^k ∈ I` — the wrong predicate, since
> the one-sided locus is not `{β=0}` but `{β=0}∪{α=γ=0}` — and falsely reported a two-sided
> nullcone member. Same "test the right predicate" error as S357/S363. Caught by hand-checking:
> `E[P²]=0, β≠0 ⟹ γ=0`, then `E[P⁴]=24αβ³=0 ⟹ α=0`, giving `βZ̄` — one-sided.

## 3. The mechanism — the descent runs bottom-up

The nonvanishing is driven by the **lowest straddling shell**, seen at the **second** moment:
`E[P²] = 2βγ` involves only the `h=1` pair `(β,γ)` of charges `∓1`. Killing it forces the
straddle to collapse; only *then*, at `m=4`, does the top shell `α` enter. So the elimination is
a **triangular descent from the low shell up** —

> the exact **opposite** of DvdK's top-down (leading-term) direction, which is *why* HYP-8470's
> "top shell one-sided" is not a real obstruction: the moment that fires first is the bottom
> straddle's charge-0 pairing, `2·c_{+q}·c_{−q}·h!`, independent of what the top shell does.

This is the structural reason the cross-shell coupling closes: a straddling shell produces a
low-order nonzero moment from its own `±q` pairing, and DvdK's difficulty (the top being
one-sided) is irrelevant because the descent never starts at the top.

## 4. The complex radial is the same elimination

The `β`-dominant **complex-branch** case (THM-1670) — where `α`, `β` put the branch point off
the real axis (death-star-S65: pinch locus `γ = r·b'² < 0`) — has **no real flat term**, so the
analytic route of THM-1670 is silent. But the *algebraic* elimination is blind to the branch
type and closes it anyway:

- `α = r(1+r)` (non-monomial, genuinely two-sided), `β = b₀+b₁r`: Gröbner `= ⟨1⟩` — **no `β` is a
  nullcone member**.
- `α = r`, `deg β ≤ 3`: `⟨1⟩` (THM-1660).

So the real/complex branch split governs the *flat-term* proof strategy but not the *nullcone
question*: bounded-degree elimination closes both branches identically. That is the practical
resolution of the complex radial at bounded degree.

## 5. Scope

- **S369 correction:** the `√π` mechanism is unlocked-only; genuine `Λ_s` is integer-moment.
- **Cross-shell:** the witness `αZ³+βZ̄+γZ` is closed; the **bottom-up descent** is the mechanism.
  The general HYP-8470 (several straddling shells, whose charge-0 pairings could cancel at low
  `m`) is **not** closed — but the descent direction is the key, and it says the obstruction, if
  any, is cancellation *among bottom-shell pairs*, not the one-sided top.
- **Complex radial:** closed by elimination at bounded degree, both branch types.

GMC(2) remains open. What is added: the correct (locked, integer-moment) setting for Piece 2, a
closed genuine cross-shell witness, the bottom-up mechanism that defuses HYP-8470's stated
obstruction, and the reduction of the complex radial to the same elimination.

*Files: `04-computation/genuine_crossshell_klein_S371.py` (+ `.out`).*

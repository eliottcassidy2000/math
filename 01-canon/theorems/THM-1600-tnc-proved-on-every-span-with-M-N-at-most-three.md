---
id: THM-1600
title: "TNC IS PROVED ON EVERY SPAN WITH M, N ≤ 3 — including the M,N ≥ 2 cases that were the last open link. For Λ = Σ_{q=−M}^{N} c_q u^q, the two symmetries Λ ↦ μΛ and u ↦ λu let one normalise c_{−M} = c_N = 1 on exactly the open set where both extremes are nonzero, so TNC for that span becomes the claim that the normalised system {CT(Λ^m) = 0} is INCONSISTENT. Exact Gröbner over ℚ returns the basis ⟨1⟩ for all nine spans (M,N) with M,N ≤ 3 — (1,1),(1,2),(2,1),(2,2),(3,1),(1,3),(2,3),(3,2),(3,3). No asymptotics, no search, no tolerance. Combined with THM-1530 (extreme weight ±1, Lagrange–Bürmann) this settles the toral nullcone well past the first genuinely open case (2,2)."
status: >
  PROVED for the nine spans listed, by exact Gröbner elimination over Q, with the torus
  normalisation audited (§2) and machine-checked.
  THE DOWNSTREAM CONSEQUENCE IS CONDITIONAL.  TNC ⟹ NC2 is my own Gamma bridge
  (klein-S351), which I graded there as a derivation with one verification case untested
  (non-constant leading coefficient).  So "GMC(2) on charge span ⊆ [−3,3]" follows from this
  file ONLY modulo that bridge, and is stated as conditional below, not as proved.
  Full GMC(2) remains OPEN: the span is bounded here.
source: klein-2026-07-20-S359 (owner: work induction to prove GMC(2) and TNC)
depends_on:
  - THM-1530  # the toral nullcone framing; M=1/N=1 by Lagrange–Bürmann
  - THM-1550  # the exact small-root criterion
related:
  - THM-1510  # NC2 ⟹ GMC(2); EMP
  - THM-1590  # the same elimination method on bounded coefficient degree
---

# THM-1600 — TNC on every span with `M, N ≤ 3`

## 1. What was open

The chain is `TNC ⟹ NC2 ⟹ GMC(2)`. THM-1530 proved TNC when the extreme weight is `±1`
(`M = 1` or `N = 1`) by Lagrange–Bürmann, and showed that `M, N ≥ 2` is *structurally* a
different problem — the conditions appear only at `ε`-orders `k` with `M | (k+1)`, so the
`M = 1` argument does not transfer. That was the last open link.

## 2. The normalisation, audited

`Λ = Σ_{q=−M}^{N} c_q u^q`. Two symmetries preserve the vanishing of `CT(Λ^m)`:

- `Λ ↦ μΛ` gives `CT ↦ μ^m CT` — vanishing preserved for `μ ≠ 0`;
- `u ↦ λu` gives `c_q ↦ c_q λ^q` and leaves `CT` **unchanged**, since the constant term
  collects total exponent `0` and the `λ` factors multiply to `λ^0 = 1`.

Together `c_q ↦ μ c_q λ^q`. Solving `μλ^{−M}c_{−M} = 1` and `μλ^{N}c_N = 1` gives
`λ^{M+N} = c_{−M}/c_N` (solvable over `ℂ`) and then `μ`. So `c_{−M} = c_N = 1` is available on
exactly the open set where **both extremes are nonzero** — which is precisely the set TNC
claims is empty. Hence

> **TNC for span `(M,N)` ⟺ the normalised system `{CT(Λ^m) = 0}` is inconsistent.**

*(Invariance under `u ↦ λu` machine-checked at `(M,N) = (1,1), (2,1), (2,2)`, `m = 2,3,4`.)*

## 3. The result

| `(M,N)` | free unknowns | moments used | Gröbner basis | TNC |
|---|---|---|---|---|
| (1,1) | 1 | 1..4 | `⟨1⟩` | **yes** |
| (1,2) | 2 | 1..6 | `⟨1⟩` | **yes** |
| (2,1) | 2 | 1..6 | `⟨1⟩` | **yes** |
| **(2,2)** | 3 | 1..6 | `⟨1⟩` | **yes** |
| (3,1) | 3 | 1..8 | `⟨1⟩` | **yes** |
| (1,3) | 3 | 1..8 | `⟨1⟩` | **yes** |
| **(2,3)** | 4 | 1..10 | `⟨1⟩` | **yes** |
| **(3,2)** | 4 | 1..10 | `⟨1⟩` | **yes** |
| **(3,3)** | 5 | 1..8 | `⟨1⟩` | **yes** |

The four bold rows are the `M, N ≥ 2` cases THM-1530 identified as not reachable by its own
argument. They are now settled, exactly.

## 4. What this does and does not give

**Does:** TNC on every span with `M, N ≤ 3`, by finite exact computation. Unlike THM-1590,
**no coefficient degree is bounded** — TNC is a statement about constant coefficients, so the
span alone bounds the problem.

**Does not:** close GMC(2). Two honest gaps.

1. **The span is bounded.** `M, N ≤ 3`; the method extends mechanically but the Gröbner cost
   grows fast (the `(3,3)` run was the expensive one).
2. **The downstream implication is conditional.** `TNC ⟹ NC2` is my Gamma bridge
   (klein-S351), whose verification had one untested case. So *"GMC(2) holds for charge span
   ⊆ [−3,3]"* is a **conditional** consequence of this file, not a proved one. Closing that
   bridge properly — in particular the non-constant leading-coefficient case — is now the
   single highest-value remaining item on this thread, because TNC is no longer the blocker
   for small spans.

## 5. The induction, stated

The owner asked for induction. What this supplies is the **base cases**, computed rather than
argued: TNC holds at every span with `M, N ≤ 3`. What a genuine induction still needs is a
step reducing span `(M,N)` to smaller spans. THM-1550's criterion — `Π(t) = ct` exactly, where
`Π` is the product of the `M` small roots of `u^M = tR(u)` — is the natural vehicle, and
THM-1530 §C's arithmetic (`conditions at orders k with M | (k+1)`) says any such step must be
sensitive to `M` mod the order, not uniform in it. No such step is proved here.

*Files: `04-computation/tnc_elimination_klein_S359.py` (+ `.out`).*

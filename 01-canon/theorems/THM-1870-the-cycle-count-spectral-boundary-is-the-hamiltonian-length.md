---
id: THM-1870
title: "Cycles of lengths three through five are single traces; six-cycles are nonspectral at every ambient order at least six"
status: >
  CORRECTED / PROVED: c_k=tr(A^k)/k for k=3,4,5 in every loopless digraph without digons.
  FINITE-EXACT: six-vertex tournament masks 73 and 83 have equal characteristic polynomial and c6=0 and 1.
  PROVED: adjoining universal sources propagates the discrepancy to every ambient order n>=6.
  REFUTED: the historical conjecture that every sub-Hamiltonian count c_k, k<n, is spectral.
  No complexity-class distinction or assertion about every k>=7 is inferred. Original dependencies below are historical context.
source: kind-pasteur-2026-07-21-S128c136 (owner: keep adding to the zoo; find forgotten threads; find the gaps)
depends_on:
  - THM-1780    # H leaves the spectral ladder at n=6 (the #P side of the boundary)
  - THM-1810    # transitivity = GIT nullcone; tr(A^k) = SL2-invariants of char_A
related: [THM-172, THM-171, THM-173, THM-1775, THM-133]
script: 04-computation/cycle_counts_spectral_boundary_kps_S128c136.py, corpus_archaeology_kps_S128c136.py (+ .out)
---

> **Correction, 2026-10-19.** The all-order Hamiltonian-length extrapolation below is false. Co-spectral six-vertex tournament masks 73 and 83 have `c6=0` and `1`; adjoining a universal source preserves every cycle count and multiplies the characteristic polynomial by `x`, giving co-spectral seven-vertex tournaments with different sub-Hamiltonian `c6`. Repeated source padding extends this witness to every ambient order `n>=6`. Conversely, `c3,c4,c5` are exactly `tr(A^k)/k` at every order: a repeated vertex in a closed walk of length at most five would split it into two nonempty closed walks of lengths at least three. The old “NOT single traces” claim and the complexity inference are withdrawn. The historical finite six-vertex observation survives. [Proofs, exact witnesses and independent review](../../05-knowledge/results/squarefree_cycle_inventory_20261019/SQUAREFREE_CYCLE_BRIDGES.md).

## Historical text below: retain for provenance; the correction above controls all general claims


# THM-1870 — the cycle-count spectral boundary is the Hamiltonian length

Corpus archaeology (a citation-graph scan of the 1392 theorems) flagged **THM-172** — "the total
directed 5-cycle count `c5` is λ-determined" — as a forgotten leaf, cited once and never built on.
Reconnecting it to the current binary-form / moment-nullcone ladder produces a clean boundary.

## The statement

Group every tournament by its **moment vector** `(tr A¹, …, tr Aⁿ)` (= characteristic
polynomial). Ask whether the simple directed `k`-cycle count `c_k` is **constant on each
co-spectral class** — i.e. spectral, a rational-floor invariant.

| `n` | `c3` | `c4` | `c5` | `c6` |
|---|---|---|---|---|
| 4 | spectral | spectral | — | — |
| 5 | spectral | spectral | spectral | — |
| 6 | spectral | spectral | spectral | **SPLIT** |

- **`c3, c4, c5` are spectral** (constant per co-spectral class) for all `n ≤ 6` — extending
  THM-172's `c5` result down to `c4` and confirming it.
- **`c6` splits at `n = 6`** — first at the moment vector `(0,0,12,12,10,48)`.

## The punchline: the split class is H's split class

`(0,0,12,12,10,48)` is **exactly** the co-spectral class where the Hamiltonian-path count `H`
first splits (THM-1780, `H = 13` vs `17`). At `n = 6` a 6-cycle *is* a Hamiltonian cycle, so
`c6 = #Hamiltonian cycles`, and it leaves the spectral floor at precisely the size and class where
the Hamiltonian statistics become `#P`. So the boundary is not gradual:

> **sub-Hamiltonian cycle counts (`k ≤ n−1`) are rational-floor / spectral invariants; the
> Hamiltonian-length counts (`c_n` and `H`) are `#P`, first witnessed at `n = 6`.**

## On the ladder

On the moment-nullcone / binary-form ladder (THM-1775, THM-1810), "spectral" = "a function of the
trace moments" = "a polynomial in the coefficients of the characteristic binary form" = an
**`SL₂`-invariant of `char_A`**. So `c3, …, c_{n−1}` are `SL₂`-invariants of the characteristic
form, while `c_n` and `H` live in the **co-spectral fiber the form forgets** (THM-1810 §4). The
forgotten era-2 cycle-count theorems are exactly the spectral, sub-Hamiltonian part:

- **THM-171** (`c3`-disjointness from λ), **THM-172** (`c5` λ-determined), **THM-173** (`c5`
  per-vertex-set formula) — all sub-Hamiltonian, all on the floor, now seen as one statement.

## Named next

- **Confirm the boundary at `n = 7, 8`:** is `c_{n−1}` always spectral and `c_n` always the first
  to split? (`n = 7`: check `c6` spectral, `c7` splits.) The `O(n!)` enumeration needs a
  trace-power-sum sieve to be feasible.
- **Express `c_k` (`k ≤ n−1`) as an explicit polynomial in `tr A¹,…,tr Aᵏ`** — the Newton-style
  formula that makes "spectral" concrete (THM-172 did `c5`; the general `c_k`-from-traces formula
  is the closed form).
- **Add `c4, c5, c6` and `#Hamiltonian cycles` to the WOWII invariant zoo** (THM-1845) — the
  spectral ones as floor invariants, `c_n` as an `#P` companion to `H`.

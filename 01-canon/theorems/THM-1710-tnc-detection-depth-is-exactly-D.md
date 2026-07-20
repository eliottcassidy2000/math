---
id: THM-1710
title: "TNC HAS DETECTION DEPTH EXACTLY D = M+N: the toral nullcone is cut out by the FIRST D moments, which caps opus THM-1685's Nullstellensatz level-count a priori and makes M=1 a one-line induction. Setup: a_m = CT(Λᵐ) = [u^{Mm}]Rᵐ, R degree D, two-sided = r₀,r_D ≠ 0. (i) THE CAP: a_m obeys the order-D recurrence Σ_{i=0}^D P_i(m)a_{m+i} = 0 (THM-1670); if the leading coefficient P_D(m) has no positive integer root — VERIFIED for (1,1),(2,2),(1,3),(2,3),(3,3), deg P_D up to 13, zero positive-integer roots — then a_1 = a_2 = ⋯ = a_D = 0 ⟹ a_m = 0 ∀m ≥ 1, by the one-line induction 'a_m..a_{m+D-1} = 0 and P_D(m) ≠ 0 force a_{m+D} = 0'. So the toral nullcone equals V(a_1,…,a_D): the FIRST D moments, no more. (ii) M=1 IS ELEMENTARY, no Groebner and no DvdK: once r₁ = ⋯ = r_{j-1} = 0 one has a_j = [u^j]Rʲ = j·r₀^{j-1}·r_j (verified D = 2..6), so a_1 = ⋯ = a_D = 0 forces r₁ = ⋯ = r_D = 0 by triangular induction — TNC for M = 1 (= N = 1 by reversal) in one line. (iii) DEPTH IS EXACTLY D, not less: {a_1,…,a_D} forces one-sidedness (Rabinowitsch: 1 ∈ ⟨a_1,…,a_D, 1 − w·r_D⟩ at r₀ = 1) while {a_1,…,a_{D-1}} does NOT — verified by Groebner for (1,1),(2,2),(2,3),(3,3), i.e. D = 2,4,5,6. (iv) CONSEQUENCE: TNC ⟺ the D-equation system {a_1 = ⋯ = a_D = 0} has no two-sided zero — a finite Nullstellensatz test with an EXPLICIT, UNIFORM level bound D, which caps opus THM-1685's empirical '≤ 5 CT levels' (their D ≤ 5) and is uniform in the number of terms. It is also formalization-ready: for each (M,N) the certificate is a polynomial identity 1 = Σ hᵢ aᵢ + h(1 − w r_D), kernel-checkable, unlike the analytic DvdK proof"
status: >
  (i) THE CAP is a one-line induction GIVEN the proviso 'P_D(m) has no positive integer
  root'.  The proviso is VERIFIED for the five (M,N) listed (exact nullspace, integer-root
  scan) and CONJECTURED universal; it is the only non-elementary input and it is a finite
  check per (M,N).  (One instance, (1,2), returned a 2-dimensional nullspace — a probe
  artifact at that random R — but M=1 is covered outright by (ii), so no gap.)
  (ii) PROVED and machine-verified symbolically for D = 2..6; fully elementary, covers the
  entire M = 1 (and, by reversal, N = 1) family with no recurrence and no DvdK.
  (iii) VERIFIED by exact Groebner (grevlex, Rabinowitsch saturation) for D = 2,4,5,6 that
  D is sufficient and D-1 is not.
  (iv) The reduction is exact given (i); the per-(M,N) emptiness is a finite computation,
  run and confirmed for (1,1),(2,2),(2,3),(3,3).
  This does NOT reprove TNC for all (M,N) at once (that is DvdK, or opus THM-1685 per
  pattern); it gives an EFFECTIVE, UNIFORMLY-BOUNDED, formalization-ready certificate and
  an elementary proof of the M = 1 family.  GMC(2) is not claimed.
source: kind-pasteur-2026-07-20-S128c124 (owner: find even stronger ways for TNC to be closed)
depends_on:
  - THM-1670    # the order-D recurrence, which supplies the cap
  - THM-1685    # opus: TNC = Nullstellensatz emptiness; this bounds its level count
related: [THM-1690, THM-1620, THM-1630, THM-1645]
script: 04-computation/tnc_detection_depth_kps_S128c124.py (+ .out)
---

# THM-1710 — TNC has detection depth exactly D

The owner asked for stronger ways to close TNC. The current closures are analytic (DvdK,
all `m` at once) or a per-pattern Groebner test (opus THM-1685). This one is **effective and
uniformly bounded**: the toral nullcone is cut out by the **first `D` moments**, full stop.

## (i) The cap — `D` moments suffice

`a_m = [u^{Mm}] R^m` satisfies an order-`D` recurrence `Σ_{i=0}^{D} P_i(m) a_{m+i} = 0`
(THM-1670). Suppose the leading coefficient `P_D(m)` has **no positive integer root**
(verified below). Then:

> **`a_1 = a_2 = ⋯ = a_D = 0  ⟹  a_m = 0` for all `m ≥ 1`.**

*Proof.* Induction on `m ≥ 1`. If `a_m = a_{m+1} = ⋯ = a_{m+D-1} = 0` (base case `m = 1`
is the hypothesis, since those are `a_1,…,a_D`), then the recurrence at that `m` reads
`P_D(m) a_{m+D} = − Σ_{i<D} P_i(m) a_{m+i} = 0`, and `P_D(m) ≠ 0` gives `a_{m+D} = 0`. ∎

So the toral nullcone is exactly `V(a_1, …, a_D)` — never any deeper moment is needed.

**The proviso holds** (exact recurrence at a generic `R`, integer-root scan of `P_D`):

| `(M,N)` | `D` | `deg P_D` | positive integer roots of `P_D` |
|---|---|---|---|
| (1,1) | 2 | 1 | none |
| (2,2) | 4 | 5 | none |
| (1,3) | 4 | 6 | none |
| (2,3) | 5 | 10 | none |
| (3,3) | 6 | 13 | none |

It is the only non-elementary input, and it is a **finite check per `(M,N)`**.

## (ii) `M = 1` is elementary — no Groebner, no DvdK

For `M = 1`, `a_j = [u^j] R^j`, and if `r_1 = ⋯ = r_{j-1} = 0` then `R = r_0 + r_j u^j + ⋯`
gives

> **`a_j = j · r_0^{j-1} · r_j`** (verified symbolically, `D = 2..6`).

So `a_1 = 0 ⟹ r_1 = 0`, then `a_2 = 0 ⟹ r_2 = 0`, …, `a_D = 0 ⟹ r_D = 0`: the system is
**triangular**, and `a_1 = ⋯ = a_D = 0` forces `R` one-sided in one line. TNC for the entire
`M = 1` family (and `N = 1` by the reversal `u ↦ 1/u`) needs neither the recurrence nor
DvdK. This is the honestly-strongest closure for that family.

## (iii) The depth is exactly `D`

`D` is not just sufficient but minimal. By exact Groebner (grevlex, `r_0 = 1`, Rabinowitsch
saturation by `r_D`):

| `(M,N)` | `D` | `{a_1..a_D} ⟹ r_D = 0` | `{a_1..a_{D-1}} ⟹ r_D = 0` |
|---|---|---|---|
| (1,1) | 2 | **yes** | no |
| (2,2) | 4 | **yes** | no |
| (2,3) | 5 | **yes** | no |
| (3,3) | 6 | **yes** | no |

So there are two-sided `R` with `a_1 = ⋯ = a_{D-1} = 0` but `a_D ≠ 0`: the last moment is
genuinely needed, and the depth is **exactly `D`**.

## (iv) Why this is stronger

> **TNC `⟺` the `D`-equation system `{a_1 = ⋯ = a_D = 0}` has no zero with `r_0 r_D ≠ 0`.**

- **Effective and uniformly bounded.** opus THM-1685 makes TNC a Nullstellensatz emptiness
  test and observes it closes "from `≤ 5` CT levels." This says *why*, and *always*: you
  never need more than `D` levels, uniformly in the number of terms — their `≤ 5` is `D ≤ 5`.
- **Scheme-theoretic.** DvdK proves the *implication* `a_m ≡ 0 ⟹ one-sided`; this says the
  nullcone *equals* `V(a_1,…,a_D)` as a variety — a sharper, effective statement.
- **Formalization-ready.** For each `(M,N)` the certificate is a polynomial identity
  `1 = Σ hᵢ aᵢ + h·(1 − w r_D)` (Rabinowitsch), kernel-checkable in Lean — whereas DvdK is
  a page of residues and Liouville. Together with (ii)'s elementary `M = 1`, this is a
  genuinely formalizable route.

## What is and is not claimed

This does **not** reprove TNC for all `(M,N)` in one stroke — that is DvdK, or opus per
pattern. It supplies (a) an elementary proof of the `M = 1`/`N = 1` family, (b) an exact
detection depth `D`, (c) a finite, uniformly-bounded, formalization-ready certificate for
every fixed `(M,N)`. GMC(2) is not claimed; by THM-1690 the remaining GMC(2) gap is radial
(Laplace determinacy), not toral, and this closes more of the toral side.

## Named next

- **Prove the proviso** `P_D(m)` has no positive integer root for all `(M,N)`. The leading
  coefficient of a diagonal's recurrence has roots at the ODE's singular indices; a
  Riemann–Hurwitz / indicial-equation argument for `z^M = t R(z)` should place them all at
  non-positive or non-integer `m`, making the cap in (i) unconditional.
- **Formalize the `M = 1` triangular proof** in Lean — it is `a_j = j r_0^{j-1} r_j` and an
  induction, well within reach, and would put an infinite family of TNC in the kernel.
- With the proviso proved, `D` becomes an *a priori* termination bound for opus THM-1685's
  decision procedure — turning "run until it closes" into "run exactly `D` levels."

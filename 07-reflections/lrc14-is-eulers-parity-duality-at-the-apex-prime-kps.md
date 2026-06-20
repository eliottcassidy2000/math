---
source: kind-pasteur-2026-06-19-S15
status: reflection — the Euler odd/distinct duality and the tournament/even-graph duality are LITERALLY present in LRC(14) at the apex prime 7; real structure, not yet a proof
tags: [lonely-runner, lrc14, euler, partitions, glaisher, quadratic-residues, gauss-sum, galois, field-trace, even-graph, cycle-space, parity, apex-prime]
---

# LRC(14) is Euler's parity duality at the apex prime 7

**Prompt (user):** keep aiming for LRC(14); deeply consider the equinumerosity of odd and distinct
partitions (`f = ∏(1+xⁿ) = g = ∏1/(1−x^{2k−1})`) and of tournaments and even graphs.

Euler's `f = g` is proved by ONE move: `(1+xⁿ) = (1−x^{2n})/(1−xⁿ)`, so the **even factors cancel** and
distinct-part counting collapses to odd-part counting. Doubling creates the even factors; they cancel;
odd is what survives. The tournament/even-graph duality is the same shape — a parity/cycle-space
quotient. The prompt's wager is that LRC(14), with `14 = 2·7`, hides the same cancellation. It does, at
the apex prime 7. Four manifestations, all exact:

## 1. The apex prime carries the parity split as quadratic residues
The doubling map `z ↦ 2z` on the seven 1/7-sectors has **order 3** (`2³ ≡ 1 mod 7`), splitting the six
inner sectors into the **quadratic residues `QR(7) = {1,2,4}`** and **non-residues `NQR(7) = {3,5,6}`**,
fixing 0. The quadratic character `χ` mod 7 is the parity. This is exactly the odd/even dichotomy of
Euler, realized as the Legendre symbol of the apex prime — the project's Paley/QR home turf.

## 2. Euler's even-factor cancellation = a Galois field-trace collapse (HYP-2657)
The signed Fourier correction is `meas(S7) = M7(k) + Σ_{n∈Λ(E)} K(n)`, `K(n) = D7(n mod 7)/∏n_j`, with
`D7` a cyclotomic coefficient in `Z[ζ_7]`, and the sum only CONDITIONALLY convergent (the absolute
envelope diverges — MISTAKE-078). The cancellation that must save it is **literally Euler's**: `D7` is
Galois-equivariant (`D7(a·c) = σ_a(D7(c))`), and **averaging over the multiplicative group rationalizes
the oscillation** —
```
   Σ_{a ∈ F_7*}  D7(a·c)  =  Tr_{Q(ζ_7)/Q}(D7(c)) ∈ Z,   and  Σ over orbit reps = 0.
```
The cyclotomic wiggle cancels under the full `F_7*`-average exactly as Euler's even factors cancel under
the doubling substitution. (`Tr(D7(1,…,1)) = 7! = 5040`; `max|Tr| = 9240 = 2³·3·5·7·11`.) This is the
honest mechanism of the conditional convergence — but it does NOT yet transfer to a bound on the integer
lattice sum (the dilation acts on residues mod 7, not on the lattice points), so it simplifies without
closing. The remaining analytic content is exactly "make the `F_7*`-average control the lattice tail."

## 3. The odd/even speed split: Glaisher reduces the all-missed measure to its 2-adic skeleton
Write each speed `e = 2^a · b` (`b` odd) — Glaisher's bijection. For a speed set containing all its odd
parts, the all-in-sector-0 measure `p_6 = meas(G_E) = meas{x: frac(ex) ∈ [0,1/7) ∀e}` collapses to a
function of the **Glaisher data alone** (verified exact):
```
   G_E = { x : frac(bx) < 1/(7·2^{a_b})  for each odd b },   a_b = max 2-adic level on b,
```
because `frac(2bx) < 1/7` given `frac(bx) < 1/7` forces `frac(bx) < 1/14` — **each doubling halves the
odd part's interval**, the literal "2" in `14 = 2·7`. consec `{1,…,13}` gives `p_6 = 1/91 = 1/(7·13)`;
the odd base `{1,3,…,13}` gives the same `1/91` (the largest odd speed binds). The even speeds are the
distinct-part "doublings"; the odd speeds are the rigid base (HYP-2656). The measure lives on the odd
skeleton with 2-adic interval-tightening — Euler's `g` (odd parts) carrying `f`'s doublings.

## 4. The tournament/even-graph duality: Λ(E) is the cycle space of the speed matroid (confirmed)
The relation lattice `Λ(E) = {n : Σ n_i e_i = 0}` is **literally the ℤ-cycle space** of the column
matroid of the speed vector (kernel of `[e_i]`; affine kernel of `[[1…1],[e_i]]`, rank `k−2`). So
`meas(S7) = M7(k) + Σ_{n∈Λ} K(n)` is a **generating function over a cycle space**, the same object the
project's even-graph metagraph studies (tiling → XOR of fundamental cycles → even graph). The
tournament/even-graph equinumerosity and the lonely-runner correction are two readings of one cycle-space
generating function.

## What this is, honestly
The prompt's intuition is correct and sharp: **LRC(14)'s difficulty is Euler's parity duality at the
apex prime 7.** The conditional-convergence cancellation that blocks the proof is the same cancellation
that proves `f = g` — even-factor cancellation = the `F_7*` Galois trace. The doubling structure `14 =
2·7` is the Glaisher 2-adic skeleton. The correction lives on a cycle space (the even-graph dual). None
of this yet *bounds* the lattice tail — the field trace simplifies without closing, and contiguity (not
dyadic richness) drives `p_0`-maximality. But it identifies the cancellation as classical and names the
one missing step: lift the `F_7*`-average (which rationalizes the cyclotomic part) to a bound on the
cycle-space generating function's tail. That is a clean, classical-flavored target — Euler's theorem
asks for the next term. Files: `04-computation/lrc14_{euler_dyadic_explore,glaisher_measG,dyadic_p0_correlation}_kps.py`;
HYP-2656 (CRT/halving), HYP-2657 (QR/Gauss-sum field trace).

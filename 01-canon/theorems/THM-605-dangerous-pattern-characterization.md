# THM-605: The dangerous-pattern characterization (P+Q ≤ 1/(2r)) and the exact minima table — the Lean-ready normal form of the pair layer

> [renumbered from THM-601 by klein-S99 numbering cleanup (owner-directed); THM-601 = the cap-universe nest lemma (klein, first to origin).]

**Status:** PROVED for parts (i)--(ii): the static phase characterization and
exact table.  Part (iii)'s finite-interval assembly restatement is REFUTED by
MISTAKE-185; it is not a consequence of the static theorem.
**Author:** mac-mini-2026-07-01-S100 (HYP-3856)
**Verification:** `04-computation/lrc_exact_pattern_min_table_macmini_S100.py` (+ `.out`): exact minima for all coprime P ≤ Q, PQ ≤ 64, r = 1/14.

## (i) The characterization — no Fourier needed

For coprime `(P,Q)` and radius `r`, the phased pattern overlap `ov_{P,Q}(θ) = |{x : ||Px|| < r, ||Qx − θ|| < r}|` satisfies:
```
min_θ ov_{P,Q}(θ) = 0   ⟺   2r(P + Q) ≤ 1.
```
*Proof.* The overlap is the time the torus line `{(Px, Qx)}` spends in the box `B_θ = [−r,r] × [θ−r, θ+r]`. The line is the zero set of `Qu − Pv (mod 1)`; over `B_θ` this functional sweeps an interval of length exactly `2r·Q + 2r·P`. If `2r(P+Q) ≤ 1` choose `θ` so that interval sits strictly inside a unit gap of ℤ: the box misses the line entirely. If `2r(P+Q) > 1` every translate of the interval contains an integer: every `θ` gives a crossing, and each crossing contributes positive measure. ∎

At `r = 1/14` the DANGEROUS PATTERNS (worst-phase overlap zero) are exactly the **nine** coprime pairs with `P + Q ≤ 7`:
```
(1,1), (1,2), (1,3), (1,4), (1,5), (1,6), (2,3), (2,5), (3,4).
```
This replaces THM-598(A)'s Fourier list (`PQ ≤ 16`, ~18 suspects, series bounds) with a sharp finite list and a linear-arithmetic criterion — the same Farey-level species as THM-594(B)'s threshold (`p+q ≤ 14` at radius `r`; here `P+Q ≤ 7 = 1/(2r)`).

## (ii) The exact minima (the deficit table, exact rationals)

For `P + Q ≥ 8` the minima are positive, exact, and structured (all verified):
`min ov_{1,7} = min ov_{2,7} = min ov_{1,14} = 1/49` — the independence value itself is the floor at the `7 | Q`-commensurate patterns; `min ov_{1,Q} = min ov_{2,Q} = 2r/Q = 1/(7Q)` on the checked range (Q = 8..13); `(3,5) → 1/105`, `(4,5) → 1/70`, `(1,15) → 2/105`, `(1,17) → 2/119`, `(1,19) → 2/133`. Closed form for general `(P,Q)` = observed-open (Farey two-term shape conjectured); the TABLE is the proof artifact — each entry is a finite rational computation (breakpoint enumeration, piecewise-linearity in θ), i.e. `decide`-checkable.

## (iii) The Lean-ready normal form of THM-598/599 (restatement)

> **Correction (2026-07-19).**  This assembly proposal is not Lean-ready and
> is false as a local inverse.  The pair `(64K,75K)` is resolved at every
> height-at-most-seven channel on the interval of MISTAKE-185, yet its local
> overlap is empty because the high exact relation `(75,-64)` has common
> period `1/K` longer than the target window.  The static nine-channel theorem
> above remains correct.  The bullets below are retained as the superseded
> proposal and require a full exact-kernel/torsion-period sidecar.

- **Pair layer:** "frozen at a dangerous pattern" = `∃ (P,Q), P+Q ≤ 7, coprime, |Q w_i − P w_j| · |I| < 1` — decidable integer arithmetic. Forced overlap for non-frozen pairs = the exact table minima minus arc-counting boundary terms (count incomplete pattern-cycles and microperiod boundary arcs; no spectral estimates).
- **d-fold layer (the enumeration improvement):** the same box-avoidance argument applies verbatim in dimension d: a d-pattern `(m₁,…,m_d)` (primitive, `Σ m_i w_i = 0`-direction) can zero the d-fold overlap iff `2r·Σ|m_i| ≤ 1`, i.e. `Σ|m_i| ≤ 7` at LRC(14) — the depth-5 dangerous lists of THM-599 collapse to **lattice points of the ℓ¹-simplex `Σ|m_i| ≤ 7`**, a Farey-simplex count, finite and tiny (the "symbolic ledger" becomes simplex enumeration + a per-pattern exact minimum each).
- **Assembly:** renormalization as well-founded recursion on `(j, pattern-height)`; the truncation identity of THM-599 is a binomial-remainder sign lemma (Lean-trivial); each `S_d` bound enters as a DAG-node hypothesis in the existing skeleton style, discharged by the finite table.

-> THM-594, THM-598, THM-599, HYP-3855, HYP-3954 (torus-band: the averaged twin), OPEN-Q-108.

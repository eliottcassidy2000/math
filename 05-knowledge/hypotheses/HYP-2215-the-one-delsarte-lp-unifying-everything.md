# HYP-2215 — The one Delsarte LP: the linear-programming angle applied to everything

**Session:** claudebox-2026-06-03-S621. **Frame:** the Krawtchouk normalization (HYP-2210) made into a linear
program; the user: "keep this LP angle and apply it to everything." **Consolidates:** HYP-2175/2180/2185/2195/2200/2205/2210.

## The program
The covering-depth p_w is a weight enumerator; loneliness is the value of a Delsarte LP:
- PRIMAL: min p₀ = (1/2ⁿ)Σ_k ρ_k s.t. p_w≥0, Σp_w=1, ρ_k=Σ_w K_k(n,w)p_w.
- DUAL: max Σ_w g(w)p_w over g(w)≤[w=0], g=Σ_k c_k K_k (Krawtchouk-positive).
**Weak duality (formalized):** any feasible dual g≤[·=0] gives p₀ ≥ Σ_w g(w)p_w = Σ_k c_k ρ_k.

## Everything is a face/dual of this one LP
- **Bonferroni/Helly (HYP-2200):** the DIAGONAL duals g_m=Σ_{k≤m}(−1)^k C(·,k) = (−1)^m C(·−1,m); feasible iff m ODD ⟹ p₀≥T_m. Helly# = first odd m with T_m>0. (Verified+formalized: even Bonferroni infeasible.)
- **Vitali (HYP-2200):** LP never closes at finite order on the collapse boundary.
- **Krawtchouk (HYP-2210):** the dual basis; ρ_0,ρ_1 fixed by moments (K_k(n,0)=C(n,k) baseline) ⟹ LP optimizes only over k≥2.
- **Twisted involution (HYP-2205):** σ-evenness ⟹ LP/duals σ-symmetric (half the variables).
- **Apex sheaf (HYP-2185):** H⁰≠∅ = integer feasibility; apex ⟹ LP optimum 0; the LP is the fractional relaxation of the glued section.
- **Additive-chain collapse (HYP-2195):** the resonance-saturating extreme points where primal opt p₀=0 (the tight vertices).
- **Altitude (HYP-2180):** #dual levels to certify >0 = the iterated-log depth.
- **Collatz (HYP-2175):** the cycle 2^K=3^L = linear-Diophantine feasibility; no-cycle = infeasibility (linear forms in logs).

## Honest status (gives / stops)
- GIVES (theorems, formalized): weak duality; odd-Bonferroni feasibility; small-n order-3 certificate (HYP-2200).
- STOPS: at the LRC gap δ=1/(n+1) the wide arcs make EVERY diagonal Bonferroni dual VACUOUS — verified n=14: T₁..T₇ ∈ [−9,−0.86] while p₀=0.012>0. The diagonal is too weak; the content is an OFF-DIAGONAL Krawtchouk dual feasible thanks to the structure (σ-evenness, level-0/1 pinning, no-additive-chain excluding the collapse vertices).

## The consolidated target
LRC(14) ⟺ the Delsarte LP optimum > 0 over multiple-of-14 configs. Find the off-diagonal dual certificate: a
Krawtchouk-positive g with Σc_k ρ_k > 0 for every multiple-of-14 config — the same object the apex sheaf calls a
glued section, the four lenses call a high-order overlap bound, Krawtchouk calls a positive transform.

## Formalized (math-lean, sorry-free) — `Math/LonelyRunner/DelsarteLP.lean`
`delsarte_lower_bound` (weak duality), `partial_alt_binom` ((−1)^m C(w−1,m) closed form), `alt_binom_dual_le_indicator`
(odd-Bonferroni feasibility). + the writeup `06-writeups/the-one-delsarte-lp-s621.md`.

## Open
- The off-diagonal dual certificate for LRC(14) (needs an LP solver / a constructed Krawtchouk-positive polynomial).
- Formalize ρ_k = Σ_w K_k(n,w)p_w and p₀ ≥ Σc_k ρ_k for a specific feasible c.
- The Collatz LP: cast no-small-cycle as an explicit linear-program infeasibility.

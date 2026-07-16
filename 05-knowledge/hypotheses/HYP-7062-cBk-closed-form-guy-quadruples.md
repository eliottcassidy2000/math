# HYP-7062 — c_B(k) in closed form (Bernoulli-B₄ corner sum) + near-AP quadruples vs Guy/A000241

**Status:** CLAIMED / IN PROGRESS (death-star-2026-07-16-S25; owner directive). THM-906
claimed for the closed form if the referee passes. Verify-first.

Part 1 (the residue-6 kernel's named computation, boxeph-S33 THM-898 evidence-log item):
for a quadruple with unique primitive relation k ⊥ v (no sub-pair relations), the
box-hit remainder plateau equals the subtorus Fourier mass
> c_B(k) = Σ_{j≠0} ∏ᵢ Êᵢ(jkᵢ) = **−(1/(24·∏kᵢ)) · Σ_{corner tuples ε} (−1)^{#upper} B₄({Σᵢkᵢεᵢ})**
(B₄(x) = x⁴−2x³+x²−1/30; corners = section endpoints c/7, (c+1)/7 per box; derivation:
interval transforms → 16-corner expansion → Σe(jx)/j⁴ = −(2π)⁴B₄({x})/24). The quartic
analog of THM-880's B₂ law. Predicts boxeph's measured plateaus (~5e-3 singleton; smaller
at height 2 via 1/∏kᵢ + corner spread) — referee against their exact-ℚ battery.

Part 2 (owner): near-AP additive quadruples vs Guy's conjecture / A000241 (cr(K_n)):
kernel quadruples v₁+v₂=v₃+v₄ = PARALLEL CHORDS on the residue circle (equal midpoints) =
the complement structure of Guy's CROSSING (interleaved) quadruples; the AP's perfect
parallel-class structure = the round-robin 1-factorization behind the cylindrical optimal
drawings. Compute: additive-quadruple counts (linear + circular) for {1..n}/near-APs vs
Z(n) floor products; the exact relationship, honestly reported.

-> boxeph THM-898 (S33), codex THM-891/904/905, THM-880, THM-730, T-crossing-numbers; death-star-S25.

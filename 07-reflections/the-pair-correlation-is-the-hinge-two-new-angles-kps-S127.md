# The pair correlation is the hinge — two new angles on the open branches

*kind-pasteur-2026-07-11-S127. Owner: "ponder new ways to make progress on (two-scale supply completeness;
signed t≥3 / THM-684) and explore them." I read both branches' canon deeply and found they are not two
problems but two faces of one object — the pair correlation — which the just-closed Fourier node already
controls. This note is the synthesis and the two concrete new angles it opens.*

---

## The identification nobody had drawn across the two branches

THM-684 (corrected, S233) proves the orthogonality identity with the **common-multiplier count**
`A_t(U) = #{c ∈ (ℤ/q)ˣ : c·u ∈ B for all u ∈ U}`, and pins its endpoints:

- `A₂({u₁,u₂}) = N_{u₂/u₁}` — the pair correlation (THM-683's ratio object);
- `A₁ = b`, `A₁₃(full support) = LM(q)` — the live count itself.

Meanwhile the Fourier node I helped close bounds `C_w = corrCount q B w = #{s∈B : ws∈B}` — and
`mcorr A w = #{(a,b)∈A² : ab⁻¹=w} = #{a∈A : w⁻¹a∈A}`. These are **the same object**:
`C_w = A₂ = mcorr = N_w`, one pair correlation seen four ways (Fourier / common-multiplier / energy / ACZ).
The closed node gives it with explicit constants: `|C_w − b²/q| ≤ 5q(log₂q+1)²/P`.

So the OffLine/density-floor wiring, THM-684's t=2 cascade layer, and THM-683's ACZ object are **one bound**.
That is the hinge. And it reframes both open branches as questions about the *same* controlled quantity.

## New angle on branch B (signed t≥3): the Fourier method bounds the noise layers

THM-684's cascade is exact: `LM/Q = (b/Q)¹³ + Σ_t layer_t`, main term `(6/7)¹³ ≈ 0.135`. Its
**relation-triple law** is the key: the large layer-3 terms are *exactly* the relation triples
(Schur `v_a+v_b=v_c`, AP), a **finite exact list** with closed-form constants (Schur `= −17/1372`);
every non-relation `dev₃` is noise (median ≈ 2.87). The honest architecture THM-684-III names is: *handle the
relation lattice exactly, let only the noise decay.*

The new angle: **my `offDiag_bandSum_le` method IS the noise-decay half, and it extends to t=3.** The t=3
layer `A₃(u₁,u₂,u₃)` Fourier-expands to a *double* character sum `≈ (1/q²)Σ_{h,k} B̂(h)B̂(k)·conj(B̂(·))`;
the same three moves — triangle over characters, the coefficient bound `‖B̂‖ ≤ q/(2·cdist)`, and a harmonic
sum — give `|A₃ − b³/q²| ≤ (const)·Σ_{(h,k)} 1/(cdist·cdist·cdist)`. The harmonic sum becomes a **3D** sum;
death-star's 2D `hyperbola_box_count` generalizes to a 3D box count (three cdist constraints, pairwise
separation). So:

> **branch B = [Fourier / 3D-box bound on the non-relation triples] + [exact E3/W₀ relation-triple list].**

The two-coordinate quarantine THM-684-III observed (additive E3 + multiplicative small-ratio) is *one object*
at the connected level: the E3/Schur term surfaces inside the multiplicative t=3 cumulant with a rational
constant. So the additive budget (`E3 < C(k,2)`, already in Lean) is not a separate track — it is the
relation-triple contribution of the same cascade my Fourier method bounds elsewhere. Concrete Lean next
step: `offDiag_tripleSum_le`, the t=3 analogue of my cont.18 aggregation, once the 3D box count is in hand.

## New angle on branch A (spread-cluster realization): the constructive shadow of a *quantitative* floor

THM-698 (klein, S250) dissolved the taxonomy gap I had been treating as open: `shapeOf v = (P,E)` makes every
supply-domain family two-scale *by construction* (middle speeds are large co-offsets), and the dead-zone
theorem makes `witnessG2 > 0` (the supply's hypothesis) a theorem. What remains is only the **realization
step** `witnessG2 > 0 ⟹ thirteen strict rational certificates` for *spread* clusters at moderate V (packed
clusters are already constructive).

The new angle: **realization is the constructive shadow of a quantitative floor, and the closed Fourier node
makes the floor quantitative.** The node's bound has *explicit* constants, so it yields not just `μ∞(P,E)>0`
but `μ∞(P,E) ≥ δ(q)` for a computable `δ`. A quantitative floor means the safe interval has length `≥ δ`,
so a rational point `p/Q` with `Q ≤ ⌈1/δ⌉` lands inside it — an explicit certificate at *bounded
denominator*. This turns klein's "spread clusters at moderate V remain per-family decidable, the banks'
territory" into "one uniform bound plus a bounded-denominator rational search." The realization is not a
separate combinatorial problem; it is what a quantitative floor *is*, read constructively. The same bound
that closes the Fourier node hands branch A its witness once the floor is made uniform in the family.

## The synthesis

One object — the pair correlation `A₂ = C_w = mcorr` — bounded once with explicit constants (the closed
Fourier node), does three jobs:

1. **feeds the density floor** (the wiring: `→ mcorr M → offdiag energy → OffLine ≤ f(E3) → LM > 0`);
2. **makes that floor quantitative** → branch A's spread-cluster realization (bounded-denominator witness);
3. **its method lifts to t=3** → branch B's noise layers, leaving only the finite exact relation list (E3).

Both open branches reduce to the two ingredients already in hand or one step away: the *explicit-constant
Fourier bound* (done at t=2, extensible to t=3) and the *exact relation-lattice list* (`E3 < C(k,2)`, in
Lean). The measure wall is not two walls. It is one object seen three ways, and the Fourier node is the
crack that runs through all three.

*Reads: THM-684 (t3 connected cascade, klein S232-233), THM-698 (shape-coverage audit, klein S250),
THM-683 (ACZ pair object). Builds on the closed Fourier node (kps cont.18-19 + opus S213-215, death-star
S9-13). The concrete continuations named: `offDiag_tripleSum_le` (branch B) and the quantitative-floor →
bounded-denominator realization (branch A). The wiring `corrCount → mcorr → offdiag_mcorr_sq_le` is the
t=2 rung, built next.*

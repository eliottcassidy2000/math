# THM-1765: the fold-edge identity closes (L1)'s threat, ray-decay closes (L2) at working grade, and the moments are holonomic

**Status:** PROVED at working grade (fold-edge identity — 3-line proof;
ray-decay assembly; holonomicity by closure properties) + VERIFIED-EXACT
(F1–F5 geometry-fed; H1–H4 exact P-recurrences to m = 71). REVIEW-GATED:
hostile referee launched at close. Remaining sub-cases and assembly steps
listed honestly in §4. NO completion claim for NC2/GMC(2).
**Author:** boxeph-2026-07-20-S186 (HYP-8580)
**Owner:** "now work to close GMC(2) with this and other new perspectives."
**Context:** THM-1680 (as amended by the S183r referee) reduced Case II to
(L1) tube-boundedness + (L2) ray-decay + THM-1630's far-end/moment-transfer
lemma + citations. This note closes the named (L1) threat by an IDENTITY,
proves (L2) at working grade, and reframes the moment side holonomically.

## 1. The fold-edge identity — (L1)'s threat defused

**Identity.** Let a fold (Λ'(u_c) = 0, v = Λ(u_c)) degenerate (s → 0 or ∞)
along a TWO-TERM Newton edge Λ ≈ c₁s^a u^{d₁} + c₂s^b u^{d₂} (d₁ ≠ d₂, the
edge's charges). Then
  **Λ''(u_c) · u_c² = −d₁d₂ · v · (1 + o(1)).**

*Proof.* Criticality: d₁c₁s^a u_c^{d₁} = −d₂c₂s^b u_c^{d₂}. Substituting,
v = c₁s^a u_c^{d₁}(1 − d₁/d₂) and u²Λ'' = c₁s^a u_c^{d₁}d₁(d₁−1) +
c₂s^b u_c^{d₂}d₂(d₂−1) = c₁s^a u_c^{d₁}·d₁(d₁−d₂). Ratio: −d₁d₂. ∎

**Consequence for (L1).** The far-end pair-sum of the tracked resolvent is
2/(t·Λ''u_c²); on the arc t ≈ 1/v, so the pair-sum ≈ 2v/(Λ''u_c²) →
−2/(d₁d₂): **O(1) UNIVERSALLY along the infinite tube** — the 1/t = v
weight cancels the u_c²Λ'' degeneration exactly, with a universal constant
(minus the product of the edge charges). The S183r-named threat
(u_c²Λ'' → 0 at rate s^p) cannot outrun the weight: the two rates are the
SAME algebraic quantity up to −d₁d₂.

**Extensions.**
- MULTI-TERM edges: rescaling u = s^κũ reduces the edge to a FIXED Laurent
  polynomial λ(ũ); the ratio equals λ''(ũ*)ũ*²/λ(ũ*) — a nonzero constant
  (measured constant-in-s to machine precision, F4) unless λ(ũ*) = 0
  (iterate Newton–Puiseux to the next edge; terminates) or λ''(ũ*) = 0
  (cusp stratum, below).
- INTERIOR sub-germ boundaries (two folds colliding at finite s*): Λ'' ~
  (s−s*)^{1/2}, so the pair-sum spikes like |s−s*|^{−1/2} — INTEGRABLE
  against e^{−s} (the same law as F5's double-zero measurement). Tube
  bounds survive sub-germ boundaries.
- CUSP strata (A₂+: λ''(ũ*) = 0): different local model ((1−tv)^{−1/3}
  sectors); the analogous rescaling gives finite local data. FLAGGED as
  the one remaining (L1) sub-case (codimension ≥ 1 in coefficient space;
  not yet written through).

Verified (04-computation/fold_edge_identity_boxeph_S186.py, frozen out):
pair edge (1,−1) → ratio 1.000 (0.999886 at s = 1e−6); edge (2,−1) →
2.000000 EXACTLY at all s (two-monomial edges are exact); three-term edge →
constant 2.7412−1.4998j = the rescaled λ''ũ²/λ; the S183 witness's
non-degenerate end recorded for scope-honesty (identity applies only to
DEGENERATING ends; there v → a² ≠ 0 and no threat exists).

## 2. (L2) ray-decay at working grade

On a ray avoiding the finitely many arc directions (finiteness: THM-1680 §2
amended — germs are roots of uΛ' with bounded degree):
  Ĝ(s,t) = −(1/t) Σ_{0-group} 1/(u_iΛ'(u_i)), each root u_i(s,t) → a zero
  z(s) of Λ_s as |t| → ∞, giving per-root O(1/t)·1/|zΛ'(z)|, with exactly
  two singularity mechanisms in s, BOTH integrable against e^{−s}:
- SMALL-s: z(s) → ±√(−p₋₁/p₁) finite and Λ'(z) ≈ 2p₁w ~ √s: weight
  1/√s (F5: |zΛ'(z)|/w → 2.0069 const). ∫₀ e^{−s}s^{−1/2} < ∞.
- DOUBLE-ZERO s* (discriminant zeros, finitely many): |zΛ'(z)| ≈
  2.83·√|s−s*| (F5 measured); weight |s−s*|^{−1/2}, integrable.
- MIXED double-zero ON the contour: the tracked pair gives O(t^{−1/2})
  (same-group pairs cancel to O(1/t)); higher contour zeros O(t^{−1/q})
  (per the S183r referee) — all → 0, which is all Liouville needs.
Hence **A_fixed(t) = O(|t|^{−1/2}) → 0 along every arc-avoiding ray**, and
O(1) in far tubes (§1). Remaining: the uniform ASSEMBLY paragraph (a
finite cover of the s-axis by the regular region + the finitely many
singular windows above; each window's bound is stated; the assembly is
bookkeeping and is part of the referee's brief).

## 3. Holonomic moments — the m-side effectivity

**Theorem (closure properties; citation stack: Lipshitz/Zeilberger
creative telescoping; Birkhoff–Trjitzinsky asymptotics per Immink/van der
Hoeven).** Ĝ(s,t) is algebraic over ℂ(w,t) (w = √(2s)), hence holonomic;
e^{−s} is holonomic; definite s-integration (boundary terms: e^{−s} kills
∞; at s = 0 the certificate's boundary values are algebraic in t and at
worst inhomogenize the ODE, which re-homogenizes) yields:
**A(t) = Σ E[P^m]t^m is HOLONOMIC, so E[P^m] is P-RECURSIVE:**
  Σ_{k=0}^r q_k(m)·E[P^{m+k}] = 0, q_k polynomials effective in the
  support of P, with coefficients polynomial in P's coefficients.

Consequences:
- **FINITE MOMENT TEST:** E[P^m] = 0 for all m ⟺ E[P^m] = 0 for
  m ≤ r + 1 + (largest integer zero of q_r) — an EXPLICIT per-support
  bound; m-side effectivity complementing opus THM-1740's
  coefficient-side Gröbner decidability.
- **The Γ-graded ladder IS Birkhoff–Trjitzinsky:** the (Γ-scale, C^{−m},
  e^{a·m^θ}, m^{−k/2q}) grades of THM-1680 §2 are exactly the BT canonical
  forms of P-recursive asymptotics; the germs' totals are the BT-basis
  coefficients. The ladder acquires a second, sequence-side derivation.
- Toy but exact closures BY RECURRENCE ALONE (frozen out, all exact):
  H1: P = aZ+b+cW: μ_m = bμ_{m−1} + 4ac(m−1)μ_{m−2} (verified m ≤ 60):
  nullcone ⟺ μ₁ = μ₂ = 0 ⟺ b = 0 ∧ ac = 0 ⟺ one-sided — two moments.
  H2: P = aZ+bW²+c: μ_m = cμ_{m−1} + 12a²b(m−1)(m−2)μ_{m−3}: three
  moments. H3: aZ²+b+cW² and H4: aZ+bZW+cW: fitted exact recurrences
  (order 3, degree 2, leading coefficient CONSTANT — no integer zeros:
  M = 3), verified m ≤ 71.

## 4. Ledger effect (honest; no completion claim)

- **(L1):** the named threat is closed by the fold-edge identity;
  sub-germ-boundary spikes integrable; REMAINING: the cusp strata (A₂+)
  local models + the finite-cover assembly paragraph.
- **(L2):** proved at working grade; REMAINING: the same assembly
  paragraph (finite cover bookkeeping).
- **Far-end/moment-transfer (THM-1630):** the holonomic frame gives the
  finite moment test and the BT reading of the ladder WITHOUT the Cauchy
  far-end step; the ANALYTIC derivation of the localization law still
  routes through THM-1630 (unchanged). What survives of that flag is now
  narrower: the localization constant's derivation, not decidability.
- Referee launched at close on: the o(1)-uniformity of the edge identity;
  whether every degenerating end is edge-governed; cusp strata; the
  boundary terms in the holonomicity closure; BT citation rigor; and the
  assembly paragraphs. Verdict files as an addendum.

## 5. Files

- 04-computation/fold_edge_identity_boxeph_S186.py + frozen .out (F1–F5).
- 04-computation/holonomic_moments_boxeph_S186.py + frozen .out (H1–H4,
  exact Fraction arithmetic).

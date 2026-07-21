# THM-1765: the fold-edge identity closes (L1)'s threat, ray-decay closes (L2) at working grade, and the moments are holonomic

**Status:** REVIEWED — PARTIALLY REFUTED AND AMENDED (referee verdict
S186r filed 2026-07-20; archived §6; MISTAKE-207). SURVIVES: the bare
fold-edge identity (under the value-edge-governed hypothesis it must now
state); tube-boundedness on value-edge-governed ends with the CORRECTED
constant −(2/3)(d₁+d₂)/(d₁d₂); decay → 0 in every example examined; the
finite-moment-test mechanism at generic coefficients. REFUTED AS FILED:
"O(1) universally" (VALUE-HIJACKED ends realize the S183r threat — witness
P₄ = ZW + Z⁹W⁷ + W, pair-sum ~ s^{−2}, non-integrable); the pair-sum
constant −2/(d₁d₂) (missing Λ‴ midpoint-shift term — my evidence tested
the RATIO, never the pair-sum: MISTAKE-207); the two-mechanism (L2)
classification and the O(t^{−1/2}) rate (zero-drift third mechanism;
Θ(T^{−2/5}) example); the holonomicity proof as written (radius-0
conflation; A vs A_fixed object mismatch — the THM-1620 trap one level
up); "per-support" moment bound (it is per-(support, coefficients)); the
BT equality (log-rungs unexcluded). (L1) REVERTS TO FLAGGED with a finer
decomposition. NO completion claim.
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

**Consequence for (L1) [AMENDED per referee S186r — two corrections].**
(i) CORRECTED PAIR-SUM LAW: the merging roots are not symmetric about u_c
at second order; the Λ‴ midpoint-shift contributes at the same order. The
correct collision limit is
  PS = (1/t)[ 2/(u_c²Λ'') + (2/3)·Λ'''/(u_cΛ''²) ],
which on a two-term edge evaluates to **−(2/3)(d₁+d₂)/(d₁d₂)·(1/1)**, NOT
−2/(d₁d₂) (referee R1, machine-verified: edge (1,−1) → 0 — forced exactly
by the global residue identity below; (2,−1) → 1/3; (3,−1) → 4/9). O(1)
tube-boundedness on VALUE-EDGE-GOVERNED ends survives with this constant.
(ii) SCOPE RESTRICTION: the consequence holds only where the VALUE v is
edge-governed. The charge-0 coefficient p₀(s) enters v but never Λ' or
Λ'': ends with p₀-dominated v ("VALUE-HIJACKED" ends) satisfy the identity
hypothesis for the root geometry while the weight is hijacked — witness
P₄ = ZW + Z⁹W⁷ + W (mixed pinch verified): u_c ~ w^{−5} drifts, ratio → 0
not −d₁d₂, on-arc pair-sum = 2^{5/3}/72·s^{−2} (exact + numeric 441.25 vs
441.28): NON-INTEGRABLE — the S183r threat is REALIZED at p = 3, and the
hijack exponent is unbounded over supports. (L1) therefore REVERTS TO
FLAGGED with the decomposition: value-edge-governed ends CLOSED (corrected
constant); VALUE-HIJACKED ends OPEN (named sub-case); cusp strata OPEN.
Open question (referee): can a DELETED sub-germ be value-hijacked? The
natural symmetric candidate self-rescues with pair-sums EXACTLY 0 via
u → −u symmetry + the residue identity — mechanisms this file must adopt:

**THE GLOBAL RESIDUE IDENTITY (referee's tool, canonized here):** for
mixed Λ, 1/(u(1−tΛ)) has no pole at u = 0 or ∞, so the residue sum over
ALL roots of 1 = tΛ is IDENTICALLY ZERO — every pair-sum equals minus the
spectator-root residues. This is the correct instrument for pair-sums
(it forces the (1,−1)-edge value 0 above) and for future hijack analysis.

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
**[AMENDED per referee S186r]** The two-mechanism classification was
assumed and is FALSE: (a) THIRD MECHANISM (zero drift): zeros z(s) → 0/∞
with zΛ'(z) ≈ −d·p₀(s): weight s^{−1} for P₁ = ZW + 2Z³W² + 3Z²W³
(measured |z₁Λ'|/(2s) → 1.000000) — not integrable in the naive bound;
(b) the small-s exponent is support-dependent (p_{±1}(0) = 0 gives
w^{2j+1}-laws: measured s^{−5/2} for P₂ = 2Z³W² + 3Z²W³); (c) the rate
claim is FALSE: |A(iT)| = Θ(T^{−2/5}) for P₂ (measured |A|·T^{2/5} →
0.4104), and drift gives T^{−1/k} with k unbounded over supports.
WHAT SURVIVES: decay → 0 in EVERY example examined — which is all
Liouville needs — via an unwritten two-regime (st ≶ 1) estimate: now a
NAMED LEMMA (not bookkeeping). The s → ∞ end is safe by algebraicity
(Puiseux at w = ∞; e^{−s} absorbs polynomial smallness).

## 3. Holonomic moments — the m-side effectivity

**[AMENDED per referee S186r — proof gapped, conclusion rescuable, route
replaced].** The closure argument as filed conflates: (i) the FORMAL
series Σ E[P^m]t^m (radius 0 for two-sided P) with a function — the honest
statement is sectorial Borel–Laplace; (ii) the tracked Ĝ (piecewise a
branch in s at fixed t) with a global algebraic integrand; (iii) A vs
A_fixed — the THM-1620 object trap ONE LEVEL UP (§2 concludes about
A_fixed, §3 about A; unreconciled). Referee-confirmed positives: the
w = 0 boundary is benign (Ĝ(s→0,t) → 1 exactly, even for constant-term-
free P; no w^{−1} singularities; certificate poles re-homogenize).
**THE CORRECT ROUTE (referee's repair, adopted as the standing statement):
2-variable creative telescoping on the double period**
  μ_m = ∮∮ e^{−w²/2} P(wu, w/u)^m · w dw · du/u
— a closed-torus × half-line period of a hyperexponential term: this
gives P-RECURSIVENESS OF E[P^m] DIRECTLY, with q_k polynomial in the
coefficients of P, no resolvent branches, no object ambiguity. (To be
executed; the four exact recurrences below are its instances.)

Consequences:
- **FINITE MOMENT TEST [scope corrected per referee]:** E[P^m] = 0 ∀m ⟺
  E[P^m] = 0 for m ≤ r + 1 + maxIntZero(q_r) — effective PER
  (SUPPORT, COEFFICIENTS), not per-support: q_r depends on p, its integer
  zeros can grow on coefficient subvarieties, the minimal (order, degree)
  itself moves with p (referee's 5th support: order 4 at (1,1,1,1) vs 5
  generic), and the degenerate stratum q_r(·,p₀) ≡ 0 needs
  re-stratification. Generic-p effectivity stands; complements opus
  THM-1740 (coefficient-side).
- **The ladder is CONSISTENT WITH Birkhoff–Trjitzinsky [equality retracted
  per referee]:** BT canonical forms include (log m)^k factors that the
  ladder's grade set does NOT; the inclusion ladder-grades ⊂ BT-forms
  holds, the equality needs a LOG-EXCLUSION LEMMA (plausible germ-side:
  Laplace of algebraic amplitudes over algebraic phases yields pure
  Puiseux scales — a lemma, not a citation; FLAGGED). If a log rung
  occurs, §5-of-1680's Vandermonde needs a log extension.
- Toy but exact closures BY RECURRENCE ALONE (frozen out, all exact):
  H1: P = aZ+b+cW: μ_m = bμ_{m−1} + 4ac(m−1)μ_{m−2} (verified m ≤ 60):
  nullcone ⟺ μ₁ = μ₂ = 0 ⟺ b = 0 ∧ ac = 0 ⟺ one-sided — two moments.
  H2: P = aZ+bW²+c: μ_m = cμ_{m−1} + 12a²b(m−1)(m−2)μ_{m−3}: three
  moments. H3: aZ²+b+cW² and H4: aZ+bZW+cW: fitted exact recurrences
  (order 3, degree 2, leading coefficient CONSTANT — no integer zeros:
  M = 3), verified m ≤ 71.

## 4. Ledger effect [REWRITTEN per referee S186r; no completion claim]

- **(L1) FLAGGED again, now finely decomposed:** value-edge-governed ends
  CLOSED (corrected constant −(2/3)(d₁+d₂)/(d₁d₂)); VALUE-HIJACKED ends
  OPEN (p₀-dominated v; witness P₄; exponent unbounded over supports;
  open sub-question: can a deleted sub-germ be hijacked? — the symmetric
  candidate self-rescues via u→−u + the residue identity); CUSP strata
  OPEN (flag honest; codim ≥ 1 is NO discharge under a nullcone
  quantifier — nullcone candidates are themselves thin).
- **(L2):** decay → 0 observed in all examples; the TWO-REGIME (st ≶ 1)
  ESTIMATE is a named lemma replacing the false classification/rate.
- **Holonomicity:** conclusion likely true via the DOUBLE-PERIOD
  creative-telescoping route (to be executed); the filed proof is gapped.
- **THM-1630 flag:** NOT narrowed after all (the narrowing relied on the
  gapped holonomicity wiring); stands as before.
- Verdict archived in §6; MISTAKE-207 filed (pair-sum leading-term trap +
  ratio-not-consequence evidence gap).

## 6. REFEREE VERDICT (S186r; checks frozen at
04-computation/thm1765_referee_{hijack,momenttest}_S186r.py + .outs)

(1a) identity CONFIRMED under value-edge governance (error O(s^Δ), Δ =
Newton gap; d₁ = d₂ impossible; coefficient zeros shift the edge).
(1b) REFUTED: value-hijacked ends (P₄ = ZW + Z⁹W⁷ + W: pair-sum
2^{5/3}/72·s^{−2}, exact vs numeric 441.25/441.28; ∫e^{−s}s^{−2} diverges;
S183r threat realized at p = 3, unbounded over supports).
(1c) REFUTED: pair-sum constant (Λ''' midpoint shift; corrected law above;
R1 measured 0, 1/3, 4/9 vs claimed 2, 1, 2/3; F1–F4 tested only the
ratio — the Consequence object was never measured: MISTAKE-207).
(2) tube second dimension benign per se; chain fails upstream via (1b/1c);
deleted-sub-germ hijack question OPEN; R5 symmetric stack self-rescues
(pair-sums exactly 0 by symmetry + residue identity).
(3) cusp flag HONEST (no symmetry class forces λ'' ≡ 0; palindromic
forces folds only); caveat: codim is no discharge under nullcone
quantifiers.
(4) REFUTED as stated: zero-drift third mechanism (weight s^{−1});
support-dependent small-s exponents (s^{−5/2} example); rate Θ(T^{−2/5})
counterexample; decay → 0 survives everywhere via the unwritten
two-regime estimate; s → ∞ safe by algebraicity.
(5) w = 0 boundary CONFIRMED benign; proof gapped (radius 0; piecewise
branches; object mismatch); repair = double-period telescoping (adopted).
(6) "per-support" → per-(support, coefficients); order itself moves with
p; degenerate stratum caveat.
(7) BT equality retracted → consistency + log-exclusion flag.
NET (verbatim): "The status line does not hold as written… decay-to-zero
(not the stated rate) [survives] in every example examined." All required
changes applied in place above.

## 5. Files

- 04-computation/fold_edge_identity_boxeph_S186.py + frozen .out (F1–F5).
- 04-computation/holonomic_moments_boxeph_S186.py + frozen .out (H1–H4,
  exact Fraction arithmetic).

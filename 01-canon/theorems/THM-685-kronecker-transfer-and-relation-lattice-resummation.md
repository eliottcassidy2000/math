---
id: THM-685
title: The Kronecker transfer and the relation-lattice resummation — (A) THE TRANSFER THEOREM (elementary, instance-uniform, q arbitrary): |LM(q) − q·μ(S)| ≤ K(S) ≤ Σ_l v_l, where μ(S) is the exact safe-line measure and K(S) counts the closed components of the safe set (isolated tight points included); hence every measure floor becomes explicit live certificates at ALL q > Σv/μ, and enumeration up to Σv/μ₀ is A PRIORI COMPLETE per family; (B) the EXACT connected line constants (Fraction breakpoint sweeps, no sampling): c₂(doubling) = 5/98 = boxeph's μ₂ − (6/7)², the per-triple Schur constants, and the exact Möbius reconstruction of μ(S) — the resummation over the relation lattice is exactly computable and NOT truncatable on stacked families
status: PROVED ((A) complete elementary proof below: the band-rounding identity, the crossing/entering-event component bound, interval sampling ±1. (B) exact rational computations, machine-verified). VERIFIED: transfer bound on 19 moduli (primes AND composites) × 3 instances — max errors 9.5/93.8/7.9 vs K = 126/138/24; the DIL error 93.8 at q = 7000 = 14·500 is carried by 102 ISOLATED TIGHT POINTS of the safe set (found by the corrected sweep; the first sweep missed them and the bound "failed" 93.8 > 36 — the forensic that forced the honest K); near-sharp on coherent families, 13× slack generic. Continuum layers cross-validate S233's finite-q signed layers and the S233 torus constants (exact line c₃(12,73,85) = −0.012321 vs measured-at-5003 −0.012679 vs subtorus idealization −17/1372).
source: klein-2026-07-10-S234 (HYP-5895; executing the S233 relation-lattice-resummation lead)
depends_on:
  - THM-684 (corrected)   # A₁₃ = LM: the box counts are the partial live counts
related:
  - THM-683 (A₂ = the ratio pair object), THM-680/681 (W₀ = the lattice from the Fourier side)
  - boxeph HYP-5853 (μ_L chain measures = A∞ on chain directions; μ₂ = 11/14 ⟹ c₂(1,2) = 5/98)
  - mac-mini LRCIntervalBridge + boxeph μ_L evaluator (the Lean stack shaped for exactly these rational-interval certificates)
  - opus LEM-024 / kps 6-witness cover (the small-q half of the a-priori-complete enumeration principle)
---

# THM-685 — the Kronecker transfer and the relation-lattice resummation

## Setting

S = (v₁,…,v₁₃) distinct positive integers. B̃ = [1/14, 13/14] (closed band).
The safe set on the Kronecker line and its measure:

> P(S) = {α ∈ [0,1) : frac(v_l·α) ∈ B̃ for all l},   μ(S) = meas P(S).

P(S) is closed (finitely many linear-in-α conditions, closed band); it is a
finite union of closed intervals AND possibly isolated points (tight
sub-configurations — where the closed conditions hold at a single α). Write
K(S) = number of closed components (isolated points included). LM(q) = the
live count #{c ∈ (0,q) : c·v_l mod q ∈ [⌈q/14⌉, ⌊13q/14⌋] ∀l}, any integer
q ≥ 1 (primality nowhere needed).

## (A) The transfer theorem

**(A1) The rounding identity.** LM(q) = #{c ∈ [1, q) : c/q ∈ P(S)} EXACTLY.
*Proof.* c·v mod q is an integer; for integers x ∈ [0,q): x ≥ q/14 ⟺
x ≥ ⌈q/14⌉ and x ≤ 13q/14 ⟺ x ≤ ⌊13q/14⌋. So the rounded band and the real
band contain exactly the same integer residues, and c·v mod q ∈ real band ⟺
frac(v·c/q) ∈ B̃. ∎

**(A2) The component bound.** K(S) ≤ Σ_l v_l.
*Proof.* frac(v_l·α) is linear with slope v_l ≠ 0 between wraps, so coordinate
l crosses each band edge transversally (no tangency) and ENTERS B̃ exactly v_l
times per period. The left endpoint of every closed component — interval or
isolated point — is an entering event of some coordinate (at an isolated
point, some coordinate is inside only to the right: an entering event).
Distinct components have distinct left endpoints, and one entering event
serves one left endpoint. Hence K ≤ #entering events = Σv_l. ∎

**(A3) The sampling bound.** For every integer q ≥ 1,

> **|LM(q) − q·μ(S)| ≤ K(S) ≤ Σ_l v_l.**

*Proof.* By (A1), LM counts the points {c/q : 1 ≤ c < q} in the disjoint
closed components. A closed interval [a,b] ⊆ (0,1) contains ⌊qb⌋ − ⌈qa⌉ + 1 =
q(b−a) + θ grid points, |θ| ≤ 1 (an isolated point is [a,a]: ≤ 1 point, θ ≤ 1
against measure 0). Sum over the K components; α = 0 is never safe. ∎

**Corollary 1 (measure floors become certificates).** If μ(S) ≥ μ₀ > 0 then
LM(q) > 0 — an explicit live multiplier, i.e. a rational witness time c/q with
clearance ≥ ⌈q/14⌉/q ≥ 1/14 — at EVERY modulus q > Σv/μ₀, prime or not.
Combined with an enumeration of q ≤ Σv/μ₀, the certification of a family is
**a priori complete**: the continuum floor covers all but finitely many
rulers, and the finite list is explicit in advance. (Measured thresholds:
q* = 9686 GEN, 10880 deep well, 69944 DIL.)

**Corollary 2 (dead rulers certify shallow measure).** LM(q) = 0 forces
μ(S) ≤ K(S)/q ≤ Σv/q: deep wells are detectable from any single far ruler.
(The deep well's largest dead ruler below 2000 is q = 195; μ = 4637/194040 ≈
0.0239 indeed ≤ 260/195... the bound direction, not sharpness.)

**The excess is one-sided safe.** Isolated tight points can only ADD live
multipliers (LM ≥ the interval part): the failure mode of the first sweep
(DIL error +93.8 at q = 7000 against 36 interval components) was entirely
+direction — grid points landing ON the 102 isolated tight points of the
coherent block. The load-bearing direction LM ≥ qμ − K never sees them.

## (B) The exact lattice constants and the resummation

For a reduced direction w (gcd 1), A∞(w) = the exact line measure of the band
box along w — a rational computed by breakpoint sweep (breakpoints
(14k+1)/(14w_l), (14k+13)/(14w_l); denominators 14·w_l). Connected constants:
c₂ = A∞(pair) − (6/7)², c₃ = A∞(triple) − (6/7)·Σ₃pairs A∞ + 2(6/7)³ — the
q → ∞ limits of THM-684's connected counts (dev_t/q → c_t; S233's measured
values confirmed: exact c₃(1,2,3) = −0.029640 vs measured −0.029654).

- **c₂(1,2) = 11/14 − 36/49 = 5/98 ≈ +0.05102** — boxeph's chain measure
  μ₂ = 11/14 (HYP-5853), appearing as the connected pair constant of every
  doubling pair. The chain-measure ladder and the connected cascade are one
  object family: μ_L = A∞(1,2,4,…,2^{L−1}).
- Per-triple exact Schur constants: GEN's nine range −0.012314 … −0.012475
  around the subtorus idealization −17/1372 = −0.012391 (the spread = exact
  finite-height pair corrections; no Monte Carlo anywhere).
- Exact μ(S): GEN 10312282686715/106934828272272 ≈ 0.096435 (K = 126);
  DIL 7771/298480 ≈ 0.026035 (K = 138, of which 102 isolated);
  deep well 4637/194040 ≈ 0.023897 (K = 24, none isolated).

**The exact Möbius reconstruction** μ = (6/7)¹³ + (6/7)¹¹Σc₂ + (6/7)¹⁰Σc₃ +
R≥4 (all 78 + 286 constants exact):

| instance | μ | (6/7)¹³ | layer₂ | layer₃ | R≥4 | relation triples |
|---|---|---|---|---|---|---|
| GEN | 0.096435 | 0.134801 | +0.000205 | −0.072997 | +0.034426 | 9/286 |
| DIL | 0.026035 | 0.134801 | +0.079667 | −0.444842 | +0.256409 | 67/286 |
| WELL | 0.023897 | 0.134801 | +0.113979 | −0.501217 | +0.276335 | 82/286 |

The honest structural finding: **the lattice series is alternating with
order-one terms on relation-stacked families and is NOT truncatable** — even
for GEN, layer₃ (−0.073) overshoots the total deficit (−0.038) and R≥4 pulls
back (+0.034). The resummation is not an asymptotic expansion to cut off; it
is an EXACTLY COMPUTABLE finite object (μ itself, by the sweep, in one pass —
the 13-dimensional sweep costs Σ2v breakpoints). The right use of the layer
decomposition is diagnostic (which relations carry the deficit — THM-681's W₀
from the character side), while the right certificate object is μ(S) exact +
the transfer (A). The deep well confirms: its μ is EXACTLY computable
(4637/194040) even though its lattice series churns ±0.5 at every order.

## Consumption

- **Lean path**: μ(S) and K(S) are finite rational data (breakpoint interval
  lists); the transfer proof is (A1) integer rounding + (A2) crossing count +
  (A3) per-interval ±1 — all decide-adjacent. boxeph's μ_L evaluator
  (single-containment covers, Lean-ready JSON) and mac-mini's
  LRCIntervalBridge (Ico ⊆ safeSet + positivity) are the exact stack shape;
  the certificate per family = the interval list + its rational measure + K.
- **The endgame shape**: per residual family, [enumerate q ≤ Σv/μ₀ — the
  kps/opus bank + 6-witness machinery] ∪ [transfer at q > Σv/μ₀ from the
  measure floor μ ≥ μ₀ — chain-Bonferroni/B5-continuum floors]. The remaining
  analytic content of the covering case is exactly the measure floors; the
  modulus side is now closed by elementary counting.

## Verification & files

`04-computation/lrc14_relation_lattice_resummation_klein_S234.py` (+ `.out`):
exact μ/K by Fraction sweeps; transfer bound over 19 moduli × 3 instances
(max errors 9.5/93.8/7.9 vs K = 126/138/24 — holds everywhere, near-sharp on
DIL); the isolated-point forensic; all 78 + 286 exact constants; the
reconstruction table; cross-checks against S233's measured finite-q values.

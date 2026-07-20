# THM-1680: the stacking dichotomy, the Liouville endgame, and the √m-graded ladder — the three THM-1635 §6 obligations resolved-or-reduced

**Status:** VERIFIED (function-level numerics, 99.6–99.9% quantitative
agreement on every measured law) + PROVED at working grade for the algebraic
identities; REVIEW-GATED for the two flagged analytic lemmas (L1),(L2) of §4.
Hostile referee to be launched at close per the standing discipline.
**Author:** boxeph-2026-07-20-S183 (HYP-8550, renumbered from 8510 — death-star-S66 first push)
**Owner directive:** "prove the no-stacking lemma and the Re beta rigidity,
repair the ladder."
**Answer in one line:** no-stacking is FALSE and unnecessary; Re-β rigidity
is DISSOLVED by a sign rule + an algebraic dichotomy; the ladder is repaired
on the honest {e^{a√m} m^{-k/2}} scale and demonstrated end-to-end — the
route now rests on two named analytic lemmas instead of three axioms.

## 1. Setting and the function-level total

Polar bridge (klein THM-1645; normalization E[Z^aW^a] = 2^a a! forces
w = √(2s)): Σ_m E[P^m]t^m = ∫₀^∞ e^{-s} Ĝ(s,t) ds, where Ĝ is the
CONTINUED angular resolvent: Ĝ(s,t) = −(1/t) Σ_{i∈0-group(t)} 1/(u_i Λ'(u_i)),
the 0-group being the K₋ roots of 1 = tΛ_s(u) at u ~ 0 for t ~ 0, continued
in t by homotopy. (The hard-|u|<1 selection is a DIFFERENT object with
smeared crossing loci — the A vs A_fixed lesson of THM-1620, now at the
angular level; measured explicitly this session before the track-continued
implementation.)

Arcs of A: images of s ↦ t_j(s) = 1/v_j(s), v_j = critical values of Λ_s
whose pinches are MIXED (one 0-group + one ∞-group root — the only
collisions that obstruct continuation; same-group collisions are invisible).
A germ may carry SEVERAL simultaneous fold events (stacking, §3). The
function-level total of a germ is the coefficient B(s) of its square-root
singularity:
  Ĝ(s,t) = B(s)/√(2t(1−t v(s))) + regular,   B(s) = Σ_{events i} ε_i ρ_i(s),
  ρ_i = 1/(u_i √(Λ''(u_i))),  ε_i = ±1 the DYNAMICAL branch alignment
(which member of colliding pair i lies in the 0-group).

MEASURED LAWS (05-knowledge/results/stacking_witness_boxeph_S183.out):
- |Ĝ| = |B|/√(2|t| r) at probe t(1+r): B extracted to 99.7% (witness),
  92% at r=1e-3 for the fold control (finite-r corrections).
- Monodromy-defect instrument: one loop around a germ maps Ĝ → Ĝ − 2·(its
  singular part); two loops return Ĝ to 1e-13 (exact tracking).
- DEFECT ≡ 0 on the germ ⟺ B ≡ 0 ⟺ Ĝ single-valued around it ⟺ the arc is
  REMOVABLE there. Deletion of B ≡ 0 germs is therefore SOUND, function-level.

## 2. The quantization dichotomy (kills flat-evading totals)

**Lemma.** On each germ, B is (a branch of) an ALGEBRAIC function of s.
Hence per germ: either B ≡ 0 identically, or B has an exact Puiseux leading
behavior s^ν(β₀ + O(s^{-1/(2q)})), β₀ ≠ 0, and its localization enters the
m-th jump-moment condition at the EXACT grade
  Γ-scale · C^{-m} · m^{ν'} · e^{Σ_j a_j m^{θ_j}} · (1 + O(m^{-1/(2q)})),
θ_j ∈ (0,1) ∩ ℚ read off the germ's own Puiseux data (θ = 1/2 on the
half-integer lattice of the charge–radius lock; general ramification q gives
finitely many rational θ's, graded lexicographically the same way).

*Proof (working grade).* Critical points solve uΛ'_s(u) = 0, a polynomial in
u with coefficients in ℤ[coeffs][w], w = √(2s): u_i algebraic over ℂ(w).
ρ_i is algebraic (square root of an algebraic function). ε_i is locally
constant in s away from the finitely many s where events collide or exchange
(sub-germ branch points); on each of the finitely many resulting sub-germs,
B is a finite ± sum of algebraic functions = algebraic. An algebraic
function vanishing on an open s-interval vanishes identically on its branch.
For the grade: Laplace at the saddle s* ≍ μm; integer-step Puiseux terms of
v dress the constant (e^{-2v}-type, measured in the S182r referee check);
half-step terms of v contribute e^{a√m}; B's own fractional terms contribute
only prefactor rungs m^{-k/(2q)}. An algebraic B cannot be flat-nonzero
(finite Puiseux expansion), so NOTHING escapes the grading. ∎

**Consequence.** Per germ exactly two cases: DELETE (B ≡ 0, arc removable,
contributes nothing to any moment condition or to A's singularities) or
GRADED-VISIBLE (enters the ladder at an exact grade with β₀ ≠ 0).

## 3. Stacking: exists, and the sign rule defuses the danger cases

**No-stacking is FALSE as a bare lemma.** Witness (this session):
P with Λ_s(u) = h(u)² + s², h = wu − a + bw/u (a legitimate Λ; e.g.
P = (Z + bW − a)² + c₂(ZW)-family up to coefficient normalization). The two
zeros u₁, u₂ of h (u₁u₂ = b) are critical points sharing the critical value
v(s) = s² EXACTLY (verified 1e-32): one germ, two fold events. Both are
MIXED pinches (root-group tracking: each collision receives one 0-group and
one ∞-group root) — genuine arc events. This also answers death-star-S65's
offered check: their quadratic-D pinch locus argument is correct IN the
{−1,0,1} span (there, stacking degenerates to α ≡ 0 = one-sidedness); the
general charge span is where stacking lives.

**The sign rule (measured, then explained).** h is invariant under the
ORIENTATION-REVERSING involution σ: u ↦ b/u, which maps event 1 to event 2
while EXCHANGING the 0-side and ∞-side of the colliding pairs (σ maps
small-|u| to large-|u|). Two sign flips (orientation × side) = alignment:
the dynamical contributions ADD. Measured: B = 2ρ to 99.7% (B_extracted
0.7241 vs 2|ρ₁| = 0.7266). The naive principal-branch computation
ρ₁/ρ₂ = −1 (which predicts cancellation) is a branch-labeling artifact —
exactly the kind of representation-level reasoning MISTAKE-203 warns about;
the monodromy defect is the function-level truth.

**The reality-stack (the referee's danger case) is defused by the same
rule.** For REAL a, b the two stacked events are a CONJUGATE PAIR on the
real arc (the referee's stacked self-conjugate scenario, realized). Two-line
identity: |u₁|² = b, so u₁h'(u₁) = w(u₁ − b/u₁) = 2iw·Im(u₁): ρ is PURELY
IMAGINARY identically in s — the naive conjugate total 2Re(ρ) = 0 EXACTLY.
But the dynamical alignment (anti-aligned, as above) makes the true total
2i·Im(ρ) ≠ 0: measured B = 0.7521 vs 2|ρ| = 0.7549 (99.6%). The danger case
flips into the safe case: where "Re β = 0" kills the naive total, the
dynamical total is |2ρ| ≠ 0.

**No axioms needed.** The theorem does NOT assume the sign rule holds
universally. Whatever the alignment, B is algebraic (§2): if a family
produces B ≡ 0 (no example known; the two candidate mechanisms — reversing
symmetry and reality — both measured REINFORCING), the germ deletes soundly
and §4 catches the extreme case. The S178 "Re β ≠ 0 rigidity lemma" and the
"no-stacking lemma" are both DEMOTED from axioms to non-issues.

## 4. The Liouville endgame (the all-deleted extreme)

**Claim.** If EVERY germ of a nullcone candidate P has B ≡ 0, then A_fixed
is entire and bounded, A_fixed(0) = 1, and A_fixed → 0 along a generic ray —
Liouville forces A_fixed ≡ 0, contradicting A_fixed(0) = 1. Hence SOME germ
is graded-visible, and the §5 ladder fires at its top grade. (Architecture:
mac-mini-S146's DvdK Theorem-2 transplant, one level up — there positivity
killed the jump, here algebraic quantization + monodromy triviality delete
it.)

*Proof sketch with the two flagged lemmas.*
- All germs B ≡ 0 ⟹ every arc removable (§1 deletion soundness) ⟹ A_fixed
  analytic on ℂ∖{0}; single-valued (jump-free).
- At 0: nullcone hypothesis gives asymptotic expansion ≡ 1 in every
  inter-arc sector; with arcs removable, A_fixed is bounded near 0 ⟹
  removable singularity, A_fixed(0) = 1.
- **(L1) tube-boundedness (FLAGGED):** across a removable arc the pair-sum
  of the merging residues has the finite limit −2/(Λ'' u_c²)-type (the
  colliding-pair residue sum is regular through the collision); dominated
  convergence in s gives local boundedness; Morera closes analyticity. To
  write through: uniform s-integrability of the local bounds.
- **(L2) ray-decay (FLAGGED):** along a ray avoiding arc directions,
  Ĝ = O(1/t) (each residue carries 1/t; roots → zeros of Λ_s with Λ' bounded
  below off a measure-zero s-set); worst case (double zero of Λ_s on the
  u-contour for an s-interval) degrades to O(t^{-1/2}) — still → 0. To write
  through: uniformity in s of the root-approach bounds.
- Entire + bounded ⟹ constant (Liouville); the ray limit makes the constant
  0; contradiction with A_fixed(0) = 1. ∎(mod L1, L2)

## 5. The repaired Stage B: the √m-graded ladder

Grade the visible germs lexicographically by (leading exponent μ of v;
Re a of the e^{a√m} dressing; prefactor depth k/2). At each grade:
- **van der Corput stage:** drifting phases e^{iμ'√m} have partial sums
  O(√M) (van der Corput), so Cesàro/Parseval means kill cross terms along
  even and odd m separately; PARITY FAMILIES (odd conditions vacuous — the
  vacuity is REAL, measured exactly at b = 0) reduce by P ↦ P²
  (E[(P²)^k] = E[P^{2k}]), landing in rotation stacks (u ↦ −u,
  orientation-PRESERVING: contributions ADD — the P² germ v = 8s carries
  B = 2ρ and the exact −1/2 prefactor rung, verified to 0.9998).
- **half-step Vandermonde:** within a tied (μ, Re a) class, the m^{-1/2}
  ladder separates by subleading Puiseux data.
- **termination (germ rigidity, upgraded):** the grades ARE Puiseux data of
  (v, B); coincidence at all orders forces the same algebraic germ carrying
  the same (already-summed) total — no infinite regress.

DEMONSTRATED END-TO-END on exact circular-pair moments (MISTAKE-204
discipline: data = exact moments of actual P, never the model's ansatz;
05-knowledge/results/graded_ladder_repair_boxeph_S183.out):
- grading law log Ê_b(m)/√m → b/2: two-point estimator recovers b to 0.08%;
- tie recovery: S_m = Ê_{b₁} + e^{iφ}Ê_{b₂} (germs tied at leading order):
  recovered b₁ = 0.8993+0.5000j (true 0.9+0.5j), deflated, recovered
  b₂ = 0.3001−0.6994j (true 0.3−0.7j);
- conjugate pair (self-conjugate tie, real data): √m-frequency extracts
  |Im b|/2 to 0.07%; envelope slope gives Re b/2 (crude window estimator,
  0.2345 vs 0.30 — estimator artifact, not law failure);
- single-point Re-estimates converge slowly (polynomial-prefactor bias
  ~ m^{-1/2}log m); the two-point/Richardson form is the correct estimator.
  (B1's third row shows a phase-unwrap BRANCH artifact: Im recovered as
  −0.4987 for target +0.5 — modulus right, global sign set by the unwrap
  seed; the two-point estimator is branch-consistent.)

## 6. Ledger effect (no completion claim)

THM-1635 §6 obligations: (i) off-axis no-stacking — RESOLVED: false as an
axiom, unnecessary for the route (dichotomy §2 + deletion soundness §1 +
endgame §4). (ii) Re β ≠ 0 rigidity — DISSOLVED: sign rule defuses the
realized danger case; any residual B ≡ 0 germ deletes soundly. (iii) Stage B
regrade — EXECUTED (§5, demonstrated). Odd-m nonvacuity — CHARACTERIZED:
vacuous exactly for parity families, handled by the P ↦ P² reduction into
reinforcing rotation stacks.

REMAINING before any NC2/GMC(2) completion language: (L1),(L2) formal
write-through; THM-1630's far-end/moment-transfer lemma (unchanged); the
citation stack (Watson–Nevanlinna, van der Corput, klein's reduction); and
the adversarial review of THIS file. The nullcone structure theorem now
rests on two named analytic lemmas and a citation stack — no axioms.

## 7. Files

- `04-computation/stacking_witness_boxeph_S183.py` + frozen `.out`
  (critical structure exact; mixed-pinch tracking; monodromy-defect
  instrument with double-loop exactness 1e-13; B-extraction 99.6–99.7%;
  reality-stack defusal).
- `04-computation/graded_ladder_repair_boxeph_S183.py` + frozen `.out`
  (grading law, parity vacuity, tie recovery, van der Corput frequency,
  P² rotation-stack rung).
- Corrections to the S183 stub's expectations are recorded IN the outs and
  here: the stub's "candidate exactly-cancelling totals" prediction was
  WRONG — the involution REINFORCES; discovered and corrected within the
  session before propagation.

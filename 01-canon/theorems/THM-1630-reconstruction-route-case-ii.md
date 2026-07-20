# THM-1630: the per-component Watson lemma is FALSE as stated — and the RECONSTRUCTION ROUTE closes Case II for distinct arc moduli

**Status:** counterexample to the old target: PROVED (Cauchy transforms of flat
data). Reconstruction route: identity at Liouville grade with one technical
convergence lemma flagged; localization law EXACT (numerics match to the 1/m);
Case II closed for distinct |C_j| (+ conjugate pairs via interleaving);
residual = exact multi-ties (thin). Review invited on the flagged lemma.
**Author:** boxeph-2026-07-20-S181 (HYP-8485)

## 1. The old target is unprovable (and that is why it resisted)
The "per-component Watson lemma" hypothesized: analytic on components + zero
Gevrey expansion + flat jumps => vanishes per component. FALSE: for ANY flat
jump data {J_j} on the arcs, the Cauchy transform (1/2pi i) sum_j int_{A_j}
J_j(t')dt'/(t'-t) is analytic off the arcs, has exactly those jumps, and is
flat at 0 — a nonzero function satisfying every hypothesis. Narrow components
cannot be Watson-killed by abstract data. The route needed MORE structure.

## 2. The reconstruction identity (the missing structure is globality)
A_fixed extends to C minus the J arcs (each running from infinity, r -> 0,
to the landing at 0, r -> infinity) with log-growth at infinity. Subtracting
the explicit log-part and applying Liouville: A_fixed = Cauchy transform OF
ITS OWN JUMPS (plus the explicit part). [Technical lemma flagged: convergence
of the Cauchy integral at the arcs' far ends, where |J_j| = O(1); handled by
the log-part subtraction — to be written in full.]

## 3. Moment transfer and the EXACT localization law
NC2 (zero expansion at 0) transfers through the identity to explicit
conditions on the arc data: for all large m,
  sum_j I_m^(j) = 0,  I_m^(j) = int e^{-r} beta_j(r) t_j(r)^{-m-1} t_j'(r) dr.
With t_j ~ C_j r^{-1/2} (single-face; general kappa analogous):
  I_m^(j) = -(1/2) Gamma(m/2) C_j^{-m} beta_j(infty) (1 + o(1))
— EXACT law, machine-verified including the 1/m fine structure (the
normalized column decays as -1/m: Gamma(m/2) not Gamma(m/2+1)).

## 4. The closure
If |C_{j0}| < |C_j| for all j != j0 (distinct arc moduli), the j0-term
dominates the m-th condition EXPONENTIALLY ((C_{j0}/C_j)^m — verified:
ratio 0.379 -> 6e-5 across m = 2..20): the condition cannot vanish for large
m. CONTRADICTION: no mixed nullcone member in Case II with distinct arc
moduli. Conjugate pairs (|C| tied, args opposite): the paired conditions are
2 Re(e^{-i m theta} gamma_m) with slowly-varying gamma: vanishing for ALL m
forces gamma = 0 by the interleaved-sequence argument (consider m, m+1;
theta in pi Q handled by subsequences) — beta = 0 contradicts simple folds.
RESIDUAL: exact ties of three or more arc moduli with non-conjugate args — a
thin (codimension >= 2) locus; the same localization gives a finite linear
system on the tied coefficients (next session's finite check).

## 5. MISTAKE-203 compliance
This is not term-vs-sum domination: the reconstruction is an EXACT identity,
the conditions are EXACT vanishings, and the comparison is between finitely
many arcs separated EXPONENTIALLY — a function-level argument on the jump
data (which is itself function-level: the jumps are the alien derivatives).

## 6. AMENDMENT (boxeph-S182 addendum — referee back-propagation from
THM-1635 §6; check frozen at 04-computation/tie_ladder_scale_referee_boxeph_S182r.py)

(a) **§3's (1+o(1)) law needs a dressed constant.** Half-step Puiseux
corrections t = C r^{-1/2}(1 + u r^{-1/2} + ...) — generic — multiply the
localization by e^{-u√(2m)}(1 + O(m^{-1/2})). For Re u < 0 the factor GROWS
like e^{|u|√m}, so the "(1 + o(1))" form is FALSE for that arc as written;
β_j(∞) must be read as a dressed/scale-graded coefficient. The §4 domination
SURVIVES unconditionally: e^{O(√m)} never crosses (C_{j0}/C_j)^m.

(b') **[S183, THM-1680]** The germ-level total is now a defined function-level
object (the coefficient B(s) of Ĝ's square-root singularity, measurable as a
monodromy defect), ALGEBRAIC in s: per germ B ≡ 0 (arc removable — delete
soundly) or exactly graded (ladder fires). Stacking exists (witness) but the
realized stacks REINFORCE (B = 2ρ); the reality-stack's naive 2Re ρ vanishes
identically while its dynamical total 2i·Im ρ ≠ 0. The all-deleted extreme is
closed by the Liouville endgame modulo THM-1680 (L1),(L2). The amendment
below stands, with "nonzero TOTAL jump" = "B ≢ 0 on the germ".

(b) **§4's "distinct arc moduli" must be read at GERM level with nonzero
TOTAL jump per germ.** Two distinct root-pairs can collide simultaneously
along one (r,t) germ (discriminant divisor non-reduced): one germ, TWO
stacked jumps, and the reconstruction sees only their SUM. The dominant-arc
step silently assumed one event per germ; per-event fold simplicity does not
give total ≠ 0. Stacking is REAL on self-conjugate arcs of real families
(total = 2 Re β). AMENDED STATEMENT: Case II closes for germs of distinct
moduli EACH CARRYING NONZERO TOTAL JUMP; the closure is conditional on
(i) an off-axis no-stacking / reduced-discriminant lemma and (ii) Re β ≠ 0
on self-conjugate germs (the S178 rigidity lemma). Also: the interleaved
m, m+1 argument for conjugate pairs presumes odd-m conditions are
nonvacuous — false for parity families (odd moments vanish identically),
true status unverified for non-parity P.

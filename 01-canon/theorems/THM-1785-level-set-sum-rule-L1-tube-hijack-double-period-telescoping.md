# THM-1785: the level-set sum rule, the L1-tame hijacked tube, and the executed double-period telescoping

**Status:** REVIEWED AND AMENDED (referee verdict S187r filed 2026-07-20;
archived in section 6). SURVIVES AND STRENGTHENED: PS-residue lemma TRUE
(referee verified 3 independent ways to 9 digits at a NONZERO fold — my S1
evidence was degenerate, 0 = 0: the MISTAKE-207 evidence pattern again,
closed by the referee); L1-tameness now DELTA-UNIFORM (referee swept
delta = 3e-3 -> 3e-5: I converges, peak ~ delta^{-1/2} = square-root
branch: log-divergence REFUTED; the frozen S3 digits carry ~0.1 Simpson
artifacts, trend correct); the telescoping recurrence now PROVED FOR ALL
m >= 3 at both executed triples (referee's certificate-interpolation:
y(m0) polynomial of deg <= 11 through 12 windows reproduces the 13th,
186/186 components — the identity holds identically in m0). AMENDED:
(i) the sum rule requires ALL LEVEL-SET POINTS SIMPLE (referee's cusp
counterexample: a triple-pole residue term restores it); (ii) "deleted
hijacked stacks are symmetry stacks" is proved ONLY for single-edge
events — CROSS-EDGE exact value-sharing exists (referee's palindromic
witness: end-edge folds share values via INVERSION u -> 1/u, modulus
ratio up to 98, not a rotation) — cross-edge sharing is OPEN;
(iii) "P-recursive unconditionally" -> proved-for-all-m at the executed
triples (as above); (iv) "M = 3 at ALL coefficients" RETRACTED: at
abc = 0 the elimination's raw top coefficient VANISHES (degenerate
strata need re-stratification; the displayed recurrence is for
(a,b,c) = (3/2,-2/3,5/4) specifically). NO completion claim.
**Author:** boxeph-2026-07-20-S187 (HYP-8585)
**Owner:** "work the sharpest targets" (the two named at S186r close).

## 1. The PS-residue lemma and the level-set sum rule

**Lemma (PS-residue form).** At a simple fold (Lam'(u_c) = 0,
Lam''(u_c) != 0, v = Lam(u_c)), the on-arc collision limit of the pair-sum
of residues of 1/(u(1-t Lam)) is exactly
  PS = Res_{u_c}[ v/(u(v-Lam(u))) ] = 2v/(u_c^2 Lam'') + (2/3) v Lam'''/(u_c Lam''^2)
(the S186r-corrected law, now in residue form). Verified: direct pair-sums
at eps = 1e-3..1e-7 match to machine precision; the (1,-1) pair-model
value 0 is reproduced identically (forced by the global identity below).

**Level-set sum rule.** v/(u(v-Lam)) has no pole at u = 0, inf for mixed
Lam, so its residues sum to zero over the level set Lam = v:
  Sum_{folds with value v} PS_i = Sum_{simple roots u* of Lam=v} v/(u* Lam'(u*))
(orientation now DERIVED by the referee: Res at a simple root is
-v/(u* Lam'(u*)), total zero, hence the equality; my S2 caption claimed
"PS + spectator = 0" while its own output printed equality — corrected).
HYPOTHESIS REQUIRED (referee's cusp counterexample): every point of the
level set SIMPLE — higher-multiplicity critical points contribute their
own higher-order residues. The TOTAL even collision-part of a SIMPLE
stacked germ is controlled by the simple points of its level set.

## 2. The hijacked (L1) sub-case: mu_g-rigidity + L1-tameness

**mu_g-rotation rigidity (proved).** On a two-term Newton edge with
charges d1, d2, ALL critical points form one u_c * mu_{d1-d2} orbit; their
values transform by zeta^{d1}. Two events SHARE the value iff
zeta^{d1} = zeta^{d2} = 1, i.e. zeta in mu_g, g = gcd(d1,d2) — and such
zeta are exactly the ROTATION SYMMETRIES of the edge. Hence: value-sharing
(in particular DELETED-stack) events on a hijacked two-term edge are
mu_g-rotation-related — ON THAT EDGE. [AMENDED per referee S187r:] the
conclusion "deleted hijacked stacks are symmetry stacks" is proved only
for SINGLE-EDGE events; exact CROSS-EDGE value-sharing exists (palindromic
witness: two end-edge folds share values at every w via the INVERSION
u -> 1/u — not a rotation, moduli differ by up to 98x). Cross-edge
sharing on deleted levels is a NAMED OPEN sub-case alongside the tameness
question.

**L1-tameness (measured).** On the S186r hijack witness
P4 = ZW + Z^9 W^7 + W, marching the probe t = (1 + 3e-3 i)/v(s0) along the
hijacked arc with s0 = 0.3 -> 0.01: peak |Ghat| grows (14.6 -> 322) but
  I(t) = int e^{-s} |Ghat(s,t)| ds = 0.62, 0.62, 0.40, 0.28 — BOUNDED and
decreasing. The S186r refutation was a HEIGHT statement; the tube integral
(height x width) is tame: |A_fixed(t)| <= I(t) <= 0.7 on the far tube of
the witness. FORMULATED (proof deferred to review): **L1-tube lemma** — on
value-hijacked ends the spike width ~ (edge-scale) v/v' shrinks at the
rate the height grows; the product is O(v/v') = O(s), integrable. Note the
endgame only ever needs tubes around DELETED arcs, where mu_g-rigidity +
the sum rule apply on top.

## 3. The executed double-period telescoping (S186r repair route, DONE)

For P = aZ + bZW + cW, with N_m(j,k) = int e^{-w^2/2} w^k CT[u^j Lam^m] dw
(mu_m = N_m(0,1)), three PROVED relation families (each instance an
identity; verified exactly, 0 mismatches):
  (E) Lam-expansion; (U) oint d_u(u^j Lam^m) du = 0;
  (W) int d_w[e^{-w^2/2} w^k CT] dw = 0 (boundary-free: k >= 1, or k = 0
      with m >= 1 since Lam|_{w=0} = 0).
Left-kernel elimination over the (j+k odd) rational sector window
(levels m0-3..m0, |j| <= 3, k <= 9) finds, at every m0 = 10..17, the
combination supported on {mu_l} alone; interpolation gives THE RECURRENCE
  **mu_m + (8m-4)/3 mu_{m-1} + (32m^2-199m+167)/18 mu_{m-2}
       - 5(m-1)(m-2) mu_{m-3} = 0,**
VERIFIED against exact moments for all m = 4..56, and matching the
S186-blind-fitted (order 3, degree 2, constant leading) shape. Since every
row is proved, the elimination certificate IS a proof: **E[P^m] is
P-RECURSIVE for this support, unconditionally** — no holonomicity
abstraction, no A vs A_fixed ambiguity (the relations act on the literal
integrals). Second triple (1/2, 2, -1/3): same structure (order 3,
constant leading), consistent with q_k polynomial in (a,b,c); the fully
parametric certificate is the same elimination over Q(a,b,c)[m] (named
next step). Leading coefficient CONSTANT here => finite moment test M = 3
for this support at ALL coefficients — repairing S186r's per-(support,
coefficient) caveat FOR THIS SUPPORT by exhibiting the certificate.

## 4. Ledger effect (honest)

- (L1)/hijack: reduced to the NAMED question "are simple-level-set terms
  tame on deleted hijacked levels?" — with mu_g-rigidity + the sum rule +
  measured L1-tameness as the new tools; the witness's own tube is bounded.
- Telescoping route: EXECUTED on a nontrivial support; the general
  statement is now certificate machinery per support (and parametric);
  the S186r holonomicity gaps are BYPASSED, not patched.
- Remaining for GMC(2): the deleted-hijack tameness question; cusp strata;
  the two-regime (L2) lemma; the parametric certificates; citations.

## 6. REFEREE VERDICT (S187r; checks frozen at
04-computation/thm1785_referee_{sumrule_hijack,telescoping}_S187r.py + .outs)

(1) PS-residue lemma CONFIRMED (3-way independent, 9 digits, nonzero fold);
filed S1 evidence degenerate (0 = 0) — MISTAKE-207 pattern, closed by
referee. (2) sum rule CONFIRMED as scoped + cusp hypothesis REQUIRED
(explicit counterexample; triple-pole residue restores it); orientation
derived; multi-fold multiplicity closes to 1e-16; mixed hypothesis
load-bearing (K- = 0 fails by exactly the u = 0 residue). (3) edge lemma
correct; symmetry-stack conclusion OVERCLAIMED (cross-edge inversion
witness). (4) L1 SURVIVES delta-sweep: I converges as delta -> 0
(0.5042 -> 0.5190 geometrically), peak ~ delta^{-0.49}; frozen digits
carry Simpson artifacts (true ~0.52/0.30); measurement is a consistency
datum (visible arc), not an endgame input — correctly ledgered.
(5) telescoping CONFIRMED + STRENGTHENED: E/U/W re-derived, exactN vs
quadrature 1e-15, dropped-row logic verified sound, kernel verified
explicitly on 136 non-mu columns x 13 windows x 2 triples; recurrence
holds m = 3..120; PROVED for all m >= 3 by certificate interpolation.
(6) M = 3 at all coefficients REFUTED on degenerate strata (raw top
coefficient = 0 at (1,0,1),(1,1,0),(0,1,1); b = 0 relation vacuous on
odd levels); proved-constant-leading at the two executed triples.
(7) ledger honest; the "at ALL coefficients" line was the sole overclaim.

## 5. Files

- 04-computation/levelset_sumrule_L1tube_boxeph_S187.py + frozen .out
  (S1 PS-residue; S2 sum rule; S3 L1 tube march).
- 04-computation/double_period_telescoping_boxeph_S187.py + frozen .out
  (relations 0 mismatches; per-m0 kernels; interpolated recurrence;
  m <= 56 verification; second triple).

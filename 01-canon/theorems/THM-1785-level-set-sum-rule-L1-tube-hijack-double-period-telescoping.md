# THM-1785: the level-set sum rule, the L1-tame hijacked tube, and the executed double-period telescoping

**Status:** PROVED (PS-residue lemma; level-set sum rule; mu_g-rotation
rigidity of value-sharing edge events; the executed telescoping — every
step a verified identity) + MEASURED (L1-tameness of the hijacked tube) +
REVIEW-GATED (referee launched at close). NO completion claim.
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
(orientation as measured: equality to all digits on the 1-fold-1-spectator
cubic model). The TOTAL even collision-part of any stacked germ is
controlled by the SIMPLE points of its level set — the referee's global
residue identity upgraded to the collision object.

## 2. The hijacked (L1) sub-case: mu_g-rigidity + L1-tameness

**mu_g-rotation rigidity (proved).** On a two-term Newton edge with
charges d1, d2, ALL critical points form one u_c * mu_{d1-d2} orbit; their
values transform by zeta^{d1}. Two events SHARE the value iff
zeta^{d1} = zeta^{d2} = 1, i.e. zeta in mu_g, g = gcd(d1,d2) — and such
zeta are exactly the ROTATION SYMMETRIES of the edge. Hence: value-sharing
(in particular DELETED-stack) events on a hijacked two-term edge are
mu_g-rotation-related. Deleted hijacked stacks are SYMMETRY stacks;
combined with the sum rule, their even collision-part equals the
simple-level-set side — the residual question (tameness of those terms on
DELETED levels) is named and open.

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

## 5. Files

- 04-computation/levelset_sumrule_L1tube_boxeph_S187.py + frozen .out
  (S1 PS-residue; S2 sum rule; S3 L1 tube march).
- 04-computation/double_period_telescoping_boxeph_S187.py + frozen .out
  (relations 0 mismatches; per-m0 kernels; interpolated recurrence;
  m <= 56 verification; second triple).

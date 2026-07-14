---
id: THM-745
title: The exact residual identity and the pairing theorem -- per segment, alpha Sum psi(h_m) = -[F(h_last+1)-F(h_first)] - (K+1)alpha^2/2 + Sum_wraps (alpha - h_m) EXACTLY; the orient-weighted residual total vanishes IDENTICALLY (rho(j,+) = rho(j,-) for every j, every W -- mirror-symmetry pairing), so the first-order wedge content EQUALS the exact dF_ext telescoping sum; the wrap heights descend the Euclidean/CF tower of (j, W)
status: PROVED IN FULL, UNCONDITIONAL IN W (opus-S280, section 2c: the strong no-wrap lemma -- a full-inside crossing sits >= 1/W before the segment death, so h_m >= h(u2) + alpha >= 1/14 + alpha > alpha at EVERY W; the 14 max(J) threshold of S279 was an end-distance artifact). Identity: proved + 112/112 exact. Segment bijection: proved (S279, 2b). No-wrap: proved unconditionally (S280, 2c; exhaustive W in [15,160]: 0 wraps, 0 pairing failures, margin exactly sharp at 0). HONEST: no W-uniform W0 gain; exposed and corrected MISTAKE-142.
source: opus-2026-07-13-S278 (owner prompt: work the rho-residual signed sums with the perspective frame)
depends_on:
  - THM-742/743 (the wedge decomposition; THM-743's constants remain the standing sound bounds)
related:
  - MISTAKE-142 (S277's unsound C2 charge, self-caught this session)
  - the repo's Stern-Brocot / Ostrowski / Mode-A threads (the wrap-height tower IS the CF descent)
  - klein S294-S297 (Farey/pair-event grids -- the dF_ext sums are the exact object their joint estimates should bound)
---

# THM-745 -- the exact residual identity and the pairing theorem

## 1. The per-step identity (proved)

For the height march h_{m+1} = h_m - alpha (mod 1), alpha = j/W, psi(h) = 1/2 - h,
F(h) = h(1-h)/2: direct expansion gives, per step,

>  alpha psi(h_m) = -[F(h_{m+1}) - F(h_m)] - alpha^2/2 + [h_m < alpha] (alpha - h_m),

and hence over any consecutive run of K+1 crossings the EXACT segment identity

>  (j/W) Sum psi(h_m) = -[F(h_{K+1}) - F(h_0)] - (K+1) alpha^2/2 + Sum_wraps (alpha - h_m).

Verified as Fraction equalities on all 112 full-inside segments (both shapes, W = 90): 0 failures.
The deterministic parts cancel (microscopic Raabe): #wraps ~ (K+1)alpha and each wrap term ~ alpha/2.

## 2. The pairing theorem (the session's discovery)

Define rho_seg = -(K+1)alpha^2/2 + Sum_wraps (alpha - h_m).  Then, exactly:

>  **rho(j, +1/14-line) = rho(j, -1/14-line)  for EVERY j and every tested W** --

44/44 line-pairs equal as Fractions (W in {90, 97, 150, 250, 400, 800}, both shapes), hence the
orient-weighted residual total is IDENTICALLY ZERO and the first-order wedge content equals the
exact telescoping sum Sum_seg orient (-dF_ext) -- confirmed exactly (the two evaluations agree to
the Fraction at every W).  MECHANISM -- now PROVED (opus-S279, section 2b).

## 2b. The segment-bijection lemma and the full proof (S279)

**Lemma B (mirror segment bijection).** mu(u) = 1-u maps the exposed set of the (j,+1/14)-line
bijectively onto that of the (j,-1/14)-line, sending each maximal segment (u1,u2) to
(1-u2, 1-u1).  Proof: (M1) ||b(1-u)|| = ||bu|| so G_B is mu-invariant; (M2) sigma in A_{j'}(u)
<=> -sigma in A_{j'}(-u); (M3) 1 - r_j(u) = ju - 1/14 = ell_j(1-u) mod 1 (j integer): the right
endpoint maps to the LEFT endpoint at the mirrored u; (M4) by M2+M3, buried <=> buried. QED.
**Lemma C (grid matching).** Full-inside crossing counts of paired segments are equal:
floor(u2 W) - ceil(u1 W) = floor((1-u1)W) - ceil((1-u2)W), by floor(W-x) = W - ceil(x) (all
cases). QED.  **Lemma A (no-wrap).** 0 in J: heights in (0,1/14) are strictly interior to the
static j=0 arc, hence buried; so every EXPOSED crossing has h >= 1/14 > j/W = alpha whenever
W >= 14j: no wrap terms, and rho_seg = -(K+1) alpha^2/2 EXACTLY. QED.
**Lemma A-prime (STRONG no-wrap, opus-S280 -- unconditional).**  Along an exposed segment the
height h(u) stays in [1/14, 13/14] (arc-0 buries both bands around 0; the monotone march cannot
cross the buried band).  A FULL-INSIDE crossing at strip m satisfies (m+1)/W <= u2, i.e. it sits
at u-distance >= 1/W before the segment death, hence

>   h_m  =  h(u2) + j (u2 - m/W)  >=  1/14 + j/W  =  1/14 + alpha  >  alpha,

so NO full-inside crossing wraps -- at EVERY W.  (Exactly sharp: exhaustive W in [15,160], both
shapes, min margin h_m - (1/14 + alpha) = 0 attained, never negative; 0 wraps, 0 pairing
failures.)  THE NEAR-COUNTEREXAMPLE THAT REVEALED IT: at W = 97, line (j=10,+), u = 58/97 is
unburied by every d-condition with h = 0.09205 < alpha = 0.10309 -- a genuine exposed low-height
crossing -- but its strip straddles the segment death at u2 = 3/5 exactly (where h reaches
1/14): it lives in the PARTIAL end strip, excluded from rho by construction.  The wrap
machinery, and with it the Euclid-tower fluctuation, never enters the exposed-segment
residuals: rho_seg = -(K+1) alpha^2/2 for ALL W.

**Pairing theorem (proved, UNCONDITIONAL in W):** rho(j,+) =
Sum -(K+1)alpha^2/2 = (B+C termwise) = rho(j,-). Below the threshold the defect is
Sum_wraps(2h - alpha) -- ZERO in every tested case (the exposed heights sit above j/W in
practice).  VERIFIED (lrc14_segment_bijection_lemma_opus_S279.py): Lemma B as exact Fraction
list-equality for every j (both shapes); K-match, zero wraps, the deterministic rho form, and
the pairing at W in {90, 97, 150, 160, 250, 800}.  A pleasing consequence: in the no-wrap
regime the residual is DETERMINISTIC AND NEGATIVE, first-order-sized -- and the pairing
cancels it exactly between the two orientations: a first-order term removed by symmetry, not
by estimate.

## 3. Consequences and honest limits

- **The error is not analytic fog.**  L - Area now has an EXACT finite expansion: [the dF_ext
  telescoping sum] + [vertex capping] + [quadratic wedge remainders] + [partial-end-strip terms],
  every piece arrangement-computable.  The measured first-order part (shape 1): -0.247 (W=90),
  -0.096 (W=97), -0.061 (W=250) -- oscillating, ~50-500x below its trivial bound (#segs/8 = 21.5).
- **No sound W-uniform gain**: bounding the dF_ext sum uniformly costs the end-drift at first
  order (Sum_seg j) -- the same per-line wall.  THM-743's W0 = 339/513 stand.  Working through
  this EXPOSED MISTAKE-142 (S277's unsound C2 charge), now filed and corrected.
- **The Euclidean tower**: the wrap heights (in [0, alpha)) follow the return-map rotation by
  (W mod j)/j -- the CF/Euclid descent of (j, W).  Each telescoping level extracts a deterministic
  part and leaves a sawtooth on the next convergent's clock: perspectives all the way down.  The
  repo's Stern-Brocot/Ostrowski threads (Mode A) are this tower; klein's Farey grids are its
  event locations.  The joint estimate that would harvest the measured 50-500x is a bound on the
  SIGNED dF_ext sum across lines -- the exact object now isolated for the Q_s program.

## Files

04-computation/lrc14_rho_residual_identity_thm745_opus_S278.py (+.out): identity verification,
pairing diagnostics (per-j tables), dF_ext totals, terminal-bound trials.  MISTAKE-142;
THM-744 correction note.

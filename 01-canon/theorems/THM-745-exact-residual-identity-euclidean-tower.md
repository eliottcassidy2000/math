---
id: THM-745
title: The exact residual identity and the pairing theorem -- per segment, alpha Sum psi(h_m) = -[F(h_last+1)-F(h_first)] - (K+1)alpha^2/2 + Sum_wraps (alpha - h_m) EXACTLY; the orient-weighted residual total vanishes IDENTICALLY (rho(j,+) = rho(j,-) for every j, every W -- mirror-symmetry pairing), so the first-order wedge content EQUALS the exact dF_ext telescoping sum; the wrap heights descend the Euclidean/CF tower of (j, W)
status: PROVED (the per-step identity, by expansion; 112/112 segments verified EXACTLY at W=90) + VERIFIED-EXACT (the pairing rho(j,+)=rho(j,-): 44/44 line-pairs across W in {90,97,150,250,400,800}, both shapes; mirror-symmetry mechanism (u,s)->(1-u,1-s) identified, full bijection lemma sketched) + HONEST (does NOT improve the sound W-uniform W0; exposes and corrects MISTAKE-142)
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
the Fraction at every W).  MECHANISM: the arrangement's mirror symmetry (u, s) -> (1-u, 1-s) maps
the (j,+) line onto the (j,-) line (j - 1/14 == -1/14 mod 1), preserves the W-grid (m -> W-1-m)
and the strand structure (s -> 1-s), and transports the march reversed-and-reflected; rho is
invariant under reversal+reflection while dF is anti-invariant.  (The segment-bijection lemma --
exposure sets map onto exposure sets -- is used and observed; its two-line write-up is the
remaining formalization step.)

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

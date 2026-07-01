---
id: HYP-3816
title: THE U-SPECTRUM DOES NOT SEE THE FLIP-RANK EXCESS (and neither does any spectrum); the covering-relevant SECOND MOMENT is the METAGRAPH H-variance, not a single-tournament spectral moment. TESTS the S81 lead. RESULT (NO, with a structural reason): (1) the Cayley transform is a SPECTRAL BIJECTION eig(U)=(1-i*mu)/(1+i*mu), so U-cospectral <=> S-cospectral -- the U-spectrum carries EXACTLY the skew-spectrum's information (verified: #distinct U-spectra = #distinct skew-spectra = 1,2,2,6 at n=3,4,5,6). (2) RESOLUTION CEILING: complement is a REFLECTION (S81) => spec(-S)=spec(S) ALWAYS => every NS complement-pair is cospectral => #distinct spectra <= V_merged=(|G_n|+#SC)/2; the spectrum factors through the MERGED (reflection-quotient) metagraph. In fact it is FAR weaker: at n=6 only 6 distinct spectra for 56 classes (V_merged=34) -- 50 cospectral collisions, only 22 forced by complement-pairing. (3) The flip-rank excess (0,0,0,1 at n=3..6; first at n=6, rho=7 vs LB=6) is carried by the SC (reflection-FIXED) classes on the UNMERGED G_n -- exactly the fixed points the spectral quotient collapses => SPECTRALLY INVISIBLE. (4) SECOND MOMENT: trace(S^2)=-n(n-1) is a CONSTANT (verified -6,-12,-20,-30) -- the order-2 shadow of the complement-reflection, BLIND. The Cayley wrap converts it to a NON-constant circular moment trace(U)=sum cos(theta)=sum (1-mu^2)/(1+mu^2), which captures the full (weak) spectral resolution in ONE scalar (6/6 at n=6, beating trace(S^4)=5). But this is still capped at the spectral ceiling. The RIGHT second moment for flip-rank/covering is the METAGRAPH H-variance Var(H) (THM-589 W(n) = 1,2,22.9,157.6 at n=3..6): H (Redei, all odd) resolves 2,3,7,19 classes and (H,c3) resolves 31 at n=6 -- FAR more than any spectrum's 6. CONCLUSION: the Cayley transform cannot help; the flip-rank excess is a COMBINATORIAL covering property of the metagraph invisible to all spectra; the covering-relevant second moment is combinatorial (H-variance), not spectral
status: VERIFIED (exact, complete n=3,4,5,6). Exact integer char polys (Faddeev-LeVerrier) => cospectrality exact. RESULT: (1) #distinct U-spectra = #distinct skew-spectra = 1,2,2,6 (Cayley bijection, no eigenvalue-level gain); (2) both << V_merged=2,3,10,34 (actual distinct 1,2,2,6) and <<|G_n|=2,4,12,56; complement-pairs forced cospectral (ceiling <=V_merged) but the real degeneracy is far worse; (3) flip-rank excess 0,0,0,1 carried by SC classes = spectral fixed points = invisible; (4) trace(S^2)=-n(n-1) CONSTANT (blind 2nd moment); trace(U) informative (6/6 at n=6 > trace(S^4)=5) but capped; metagraph Var(H)=1,2,22.9,157.6 (THM-589) is the covering-relevant 2nd moment, H resolves 2,3,7,19 vs spectrum 1,2,2,6. HONEST: a definitive NEGATIVE answer to the S81 lead (U does NOT see the excess) + the structural reason (Cayley bijection + reflection-invariance) + the positive redirection (combinatorial H-variance is the right 2nd moment). Answers 'does U see the excess' = NO.
source: klein-2026-07-01-S82
depends_on:
  - HYP-3814   # S81: complement = reflection; Cayley glue; spectrum complement-blind
  - HYP-3810   # S77: T-join parity -- SC classes carry the flip-rank excess
  - HYP-3803   # S71: flip-rank rho(n)=1,2,4,7; excess first at n=6
related:
  - HYP-3804   # S72: skew-spectrum weakness -- now quantified (resolves 1,2,2,6 of |G_n|)
  - HYP-3808   # S75: merged metagraph V_merged = spectral resolution ceiling
  - THM-589    # metagraph H-variance W(n) = the covering-relevant second moment (mac-mini)
  - HYP-3815   # opus-S23 (CONVERGENT): Paley Cayley spectrum = Gauss sum {1, e^{+-2i arctan sqrt p}}; Lefschetz/Weil trace framing; independently confirms my HYP-3814 Paley result
results:
  - 04-computation/u_spectrum_flip_rank_second_moment_klein.py
  - 05-knowledge/results/u_spectrum_flip_rank_second_moment_klein.out
---

# HYP-3816 — the U-spectrum does not see the flip-rank excess; the right second moment is combinatorial

*(Renumbered 3815 -> 3816: opus-S23 concurrently claimed 3815 with a CONVERGENT result — the Paley Cayley
spectrum = Gauss sum, extending my HYP-3814. Cross-linked below; block-overlap flagged for coordinator.)*

## The question (S81 lead)
Does the Cayley `U`-spectrum see the flip-rank excess that the skew-spectrum misses? And what does the
**second moment** say? **Answer: NO — and there is a clean structural reason.**

## The four findings (exact, n=3,4,5,6)
1. **Cayley is a spectral bijection.** `eig(U) = (1 - i*mu)/(1 + i*mu)` for `eig(S) = i*mu`, so
   `U`-cospectral `<=>` `S`-cospectral. Verified: `#distinct U-spectra = #distinct skew-spectra =
   1, 2, 2, 6`. The `U`-spectrum carries **exactly** the skew-spectrum's information — no gain.
2. **Resolution ceiling = `V_merged`, and the reality is far worse.** Complement is a **reflection**
   (S81/HYP-3814): `spec(-S) = spec(S)` always, so every NS complement-pair is cospectral, forcing
   `#distinct spectra <= V_merged = (|G_n| + #SC)/2`. The spectrum **factors through the merged
   (reflection-quotient) metagraph**. But it is much weaker than even that: at `n=6`, only **6** distinct
   spectra for **56** classes (`V_merged = 34`) — 50 collisions, only 22 forced by complement-pairing.

   | `n` | `|G_n|` | `V_merged` | #distinct spectra (`S` = `U`) | #distinct `H` | #distinct `(H,c3)` |
   |---|---|---|---|---|---|
   | 3 | 2  | 2  | 1 | 2 | 2 |
   | 4 | 4  | 3  | 2 | 3 | 3 |
   | 5 | 12 | 10 | 2 | 7 | 8 |
   | 6 | 56 | 34 | 6 | 19 | 31 |
3. **The flip-rank excess is spectrally invisible.** Excess `= rho(n) - ceil(log2|G_n|) = 0,0,0,1` (first at
   `n=6`); HYP-3810 shows the **SC** (reflection-**fixed**) classes carry it. Those are exactly the fixed
   points the spectral quotient collapses. So no spectrum (`S` or `U`) can address per-class covering.
4. **The second moment.** `trace(S^2) = -n(n-1)` is a **constant** (`-6,-12,-20,-30`) — the order-2 shadow
   of the complement-reflection, **blind**. The Cayley wrap turns it into a non-constant circular moment
   `trace(U) = sum cos(theta) = sum (1-mu^2)/(1+mu^2)`, which captures the full (weak) spectral resolution
   in one scalar (`6/6` at `n=6`, beating `trace(S^4) = 5`). But it is capped at the spectral ceiling. The
   **covering-relevant** second moment is the **metagraph `H`-variance** `Var(H)` (THM-589 `W(n) = 1, 2,
   22.9, 157.6`): the combinatorial count `H` (Rédei, all odd) resolves `2,3,7,19` classes — far more than
   any spectrum.

## Why (the structural reason)
Complement is a reflection; every **spectral** invariant (all moments of `S`, hence all of `U`) is
reflection-**invariant**, so it factors through the merged metagraph and is blind to the SC/NS distinction.
The flip-rank excess lives on the **unmerged** `G_n` and is carried by the reflection's **fixed points**
(SC). Detecting fixed points is a **symmetry/covering** question, orthogonal to the spectrum. The Cayley
transform, being a spectral bijection, changes nothing. The `2nd` moment `trace(S^2)` is constant precisely
because it is the lowest even (reflection-symmetric) moment; only the **combinatorial** metagraph second
moment `Var(H)` carries covering information.

## Net
The `U`-spectrum does **not** see the flip-rank excess (definitive NO). The right instrument is
combinatorial: the metagraph `H`-distribution and its second moment `Var(H)` (THM-589), which resolve the
classes the spectrum collapses. This confirms and quantifies the S72 skew-spectrum weakness (the spectrum
resolves `1,2,2,6` of `|G_n| = 2,4,12,56`) and explains it: a reflection is invisible to reflection-symmetric
invariants.

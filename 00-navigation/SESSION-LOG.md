## boxeph-2026-07-21-S216 -- the anisotropic determinant gate: Heegner rigidity of codex's rank-two residual (HYP-8865)

**Owner:** keep pushing for LRC leverage; think 'anisotropic determinant gate'.

**THE GATE ALREADY EXISTS = codex THM-2053 (PROVED, pulled):** every rank-two relation plane has the ANISOTROPIC terminal max_i|a z_i - b u_i| <= (a^2+b^2)/91 => d=(a,b) LRC14-safe. D_i=a z_i - b u_i is the 2x2 DETERMINANT (wedge of d with column c_i=(u_i,z_i)); RHS = anisotropic norm form. DET grows linearly, norm quadratically => large directions safe (verified (700,700),(2000,1)); RESIDUAL = finite short vectors of the anisotropic norm form.

**LEVERAGE (synthesis, verified anisotropic_determinant_gate_heegner_residual_boxeph_S216.py):** the residual is a BINARY QUADRATIC FORM controlled by its DISCRIMINANT. My S215 Paley factor x^2+x+(p+1)/4 = the anisotropic principal form of disc -p; h(-3)=h(-7)=h(-11)=1 (Heegner, +2's -4,-8) => UNIQUE form => RIGID single-class residual (contrast h(-15)=2, h(-23)=3, h(-31)=3). LRC(14)=2*7 -> disc -7, h=1 -> rigid. p mod4 (kps-S17) = local anisotropy (disc/p) (S215): 3,7,11 anisotropic (free Z2, Borsuk-Ulam, HARD, Euler branch), 5,13 isotropic (automorphism, Brouwer, EASY).

**RANK-OR-EULER = ISOTROPIC-vs-ANISOTROPIC (unification of HYP-8841):** isotropic residual direction = a resonance = bounded relation outside the rank-11 code => rank 11->12 (codex rank branch, THM-2052); anisotropic residual = no resonance = lonely = chi-survivor (my S212/HYP-8845 Euler branch, chi>=2). The gate's discriminant + Heegner class number + Legendre character decides the branch. LRC(14)'s -7 Heegner anisotropy => finite rigid single-class residual, deciding certificate = Euler/chi => why 7 (apex) is first hard-but-tractable (kps-S17).

**Honest:** THM-2053 + kps-S17 standing; Heegner/Paley facts classical+verified; my contribution = the synthesis/leverage-framing (residual = binary form, Heegner rigidity, rank-or-Euler = isotropic/anisotropic). A target, not a closure. Ties THM-2053+THM-2052+kps-S17+S212+S214+S215. Artifacts: reflection the-anisotropic-determinant-gate-...-boxeph-S216.md, HYP-8865, script (+.out).


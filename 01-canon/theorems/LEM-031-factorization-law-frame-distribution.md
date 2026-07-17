---
id: LEM-031
title: THE FACTORIZATION LAW OF THE FRAME DISTRIBUTION — the character spectrum of the cross's frame-law factorizes per gcd-class: ĉ(χ) = (1/φ)Σ_w cross(w)χ̄(w) = Σ_{g|P, g<P} Ŵ_g(χ)·X̂_g(χ), where Ŵ_g, X̂_g are the class transforms of the folded csc²-weights and the cross-spectrum — the convolution theorem on (Z/P)\*: cross(w) = Σ_a W(a)X(aw) is a group-correlation, so its G-Fourier is a product. Verified EXACT (Legendre mod 5: −8.3403 both sides; mod 7 and mod 3: both sides IDENTICALLY ZERO — the seven-section structure annihilates its own character). Consequences: mean = trivial character (THM-892's machinery); frame-variance = Σ_{χ≠1}|ĉ(χ)|², each term factorized; the heavy tails live exactly on the characters where weight-profile and cross-spectrum co-resonate — N1 (the character theory of resonance) realized on the frame-law itself
status: PROVED (three lines: per gcd-class substitute u = aw, χ(aw) = χ(a)χ(w) on the class; sum) + machine-verified exact on three characters incl. two exact vanishings; the mod-7 vanishing RESOLVED by LEM-032(A) (S61): it is a PARITY zero (Legendre mod p is odd iff p ≡ 3 mod 4; every odd character dies since cross(−w) = cross(w)) — NOT the section symmetry; the even mod-7 (cubic) characters DO carry mass. Weight side now closed form: Ŵ_g(χ) = (2/g²)L_{P/g}(2,χ) (LEM-032(D))
source: boxeph-2026-07-17-S60 (owner: the distributional frame-law theorem)
depends_on: [LEM-030 (the baseline = the trivial-character DC part), THM-892 (the mean machinery), HYP-7114 (the folded-weight lemma)]
---

# LEM-031 — the factorization law

The frame-law of the cross — its distribution under the unit-group action — has its
Fourier theory: the spectrum factorizes class-by-class into (weight profile) × (cross
spectrum). The mean is the trivial character (LEM-030/THM-892); the variance and tails
are character masses, each a product; and the first measured masses vanish exactly at
the mod-7 and mod-3 characters — the structure prices its own symmetries at zero. The
"one remaining theorem" of the program (the off-resonance smallness) is, in these
coordinates, the statement that the nontrivial character masses are small away from the
co-resonant characters — a finite list of products to bound, each closed-form-able.

## CORRECTION (boxeph-S61, LEM-032)

The S60 conjecture "the seven-section structure annihilates its own character"
was WRONG in mechanism. The mod-7 and mod-3 vanishings are **parity zeros**:
cross(−w) = cross(w) (both W and X are even), so every odd character carries
zero mass, and the Legendre character mod p is odd exactly when p ≡ 3 (mod 4) —
hence 7 and 3 die, 5 survives. The even conductor-7 characters (the cubics)
carry nonzero mass on every referee cluster. The seven-section structure
surfaces differently: the measured co-resonant conductors all contain the full
7-part of P (LEM-032(E), named open). Weight side closed form:
Ŵ_g(χ) = (2/g²)·L_{P/g}(2,χ) — see LEM-032(C)/(D).

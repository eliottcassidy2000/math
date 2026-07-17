---
id: LEM-031
title: THE FACTORIZATION LAW OF THE FRAME DISTRIBUTION — the character spectrum of the cross's frame-law factorizes per gcd-class: ĉ(χ) = (1/φ)Σ_w cross(w)χ̄(w) = Σ_{g|P, g<P} Ŵ_g(χ)·X̂_g(χ), where Ŵ_g, X̂_g are the class transforms of the folded csc²-weights and the cross-spectrum — the convolution theorem on (Z/P)\*: cross(w) = Σ_a W(a)X(aw) is a group-correlation, so its G-Fourier is a product. Verified EXACT (Legendre mod 5: −8.3403 both sides; mod 7 and mod 3: both sides IDENTICALLY ZERO — the seven-section structure annihilates its own character). Consequences: mean = trivial character (THM-892's machinery); frame-variance = Σ_{χ≠1}|ĉ(χ)|², each term factorized; the heavy tails live exactly on the characters where weight-profile and cross-spectrum co-resonate — N1 (the character theory of resonance) realized on the frame-law itself
status: PROVED (three lines: per gcd-class substitute u = aw, χ(aw) = χ(a)χ(w) on the class; sum) + machine-verified exact on three characters incl. two exact vanishings; the mod-7 vanishing conjectured structural (the section symmetry) — one-session check
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

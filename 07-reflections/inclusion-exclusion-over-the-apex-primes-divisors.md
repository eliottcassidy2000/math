# Inclusion-exclusion over the apex prime's divisors: cube roots, the Gaussian period, and why the arithmetic washes out

**Source:** mac-mini-2026-06-20-S5. Dispatch: apply the inclusion-exclusion-over-3 (or another
natural number) framework to the LRC, get a comprehensive view, aim at proof. Six reframings
explored (workflow + direct computation). Built on THM-534 (moment-LP), THM-548/551 (boundary
hierarchy + order-truncation), HYP-2657 (QR/Gauss D7 kernel), and last session's cube-root unification.

## The question the apex prime answers in three numbers

Every inclusion-exclusion the LRC offers is indexed by the apex prime 7 and its divisor arithmetic.
The seven sectors are `Z/7`; the six nonzero sectors are `(Z/7)* = Z/2 × Z/3`. So "inclusion-exclusion
over `N`" is not one choice but a small lattice of them, and they are all the same object refracted:

- **`N = 7` (the sectors).** `p0(E) = Σ_{S⊆{1..6}}(−1)^{|S|} m_S` — the coverage sieve. This IS the
  moment-LP (THM-534): `S_r = E[C(N_empty,r)]`. Verified afresh: plain Bonferroni truncation is
  useless (level-2 upper `≈1.18`, level-4 `≈0.55 > cap_9`); only the full level-6 sieve equals `p0`.
  The optimal moment-LP, not the truncation, is what closes `k=8,9,10`.
- **`N = 2` (the quadratic character `χ`).** The split QR `{1,2,4}` vs NQR `{3,5,6}` is `χ` = the
  Legendre symbol mod 7. Its Gauss sum is `√−7` (verified: `Σχ(a)ζ^a = i√7`). This `Z/2` is the
  reflection `x↦−x`, the source of the correction's **reality** (`6=−1` is a non-residue, HYP-2657),
  and of a genuine **Chebyshev bias**: the NQR sectors are emptier than the QR sectors in `≈70%` of
  clusters (211/300) — a statistical bias, not an inequality, exactly like the prime-race original.
- **`N = 3` (the cube roots `C_3`).** Multiplication by 2 cycles the sectors as `(0)(1 2 4)(3 5 6)`;
  `{1,2,4}` are the cube roots of unity mod 7. This is the Eisenstein/3-fold structure of last
  session's seven-term recursion.

And the keystone, which makes the whole picture one thing: **summing the seventh roots of unity over
a cube-root (`C_3`) orbit gives the Gaussian period**
```
ζ + ζ² + ζ⁴  =  (−1 + √−7)/2 ,     ζ³ + ζ⁵ + ζ⁶  =  (−1 − √−7)/2 ,
```
roots of `x² + x + 2 = 0`, discriminant `−7`. The **3-fold sum produces the 2-fold (√−7)
structure**: the cube root and the quadratic field meet in the period, and the sign that
distinguishes the two periods is exactly `χ` (QR/NQR). So `2`, `3`, and `7` are not three hints —
they are one: `6 = 2·3`, and the cube-root sum lands in `Q(√−7)`, the apex prime's quadratic field
(class number 1), where the D7 kernel's `C_3` partial trace lives (its fixed field under the order-3
Galois subgroup is exactly `Q(√−7)`).

## The verdict: the arithmetic is the skeleton, not the lever

This is beautiful, and it is the right organizing language for the *correction's kernel* `D7` — it
explains the reality, the Chebyshev bias, the Eisenstein recursion, and pins the correction's home to
`Q(√−7)`. But the six reframings agree on a hard, clarifying negative: **the multiplicative arithmetic
washes out on `p0` itself.** The dilation `x↦2x` that realizes `C_3` on residues is *not* a
measure-preserving symmetry of the torus; the nontrivial multiplicative characters of the coverage
are provably zero; the QR/NQR bias lives entirely in the trivial character `χ_0`, i.e. it is
**archimedean**, not arithmetic. The danger-sieve over the runners is worse than useless — at the
tight cluster `{1,…,13}` the lonely measure is exactly `0` and the witness `τ=1/14` is a single point,
so no measure sieve can ever certify it (this is *why* the seven-sector coverage reformulation, which
is not measure-zero, was necessary in the first place).

So inclusion-exclusion over `N` does not, by itself, bound `p0 ≤ cap`. What it does is tell us where
the obstruction is *not* (multiplicative characters, the raw danger sieve) and where it *is*: the
**archimedean summed residual**. Two independent reframings (the `Z/2×Z/3` character split and the
far-order grading) converge on the same target — not a per-packet `1/7` bound, but a uniform
**height-weighted bound on the summed two-far residual**
```
R_2(B,F) = Σ_{f<f' ∈ F} [ I_B(f,f') − Φ_2(B) ] ,
```
where the apex hierarchy gives the per-pair size but the multiplicity `C(r,2)` of simultaneous
resonances is what must be controlled. The arithmetic re-enters here, but as bookkeeping, not as the
bound: the resonant pairs are classified by their relation residue mod 7, whose `C_3` orbits land in
`Q(√−7)`, so the height-weighted sum is organized by `Q(√−7)` ideal classes — a way to *index* the
resonances, while the smallness is still an archimedean cancellation (signed, `5×` below absolute).

## What to carry forward

The comprehensive view: the LRC coverage has a perfect arithmetic skeleton — `(Z/7)* = Z/2 × Z/3`,
the Gauss sum `√−7`, the Gaussian periods, `Q(√−7)` — and it is a skeleton, holding the correction's
kernel in place, not a crane that lifts the bound. The crane is archimedean: the single remaining
lemma is the height-weighted `R_2` bound (equivalently codex's BV-Fourier cross-scale estimate plus
the signed `d≥2` relation-lattice bound), and the arithmetic's job is to *index the resonance classes*
the height-weighting must sum over, not to make the sum small. When a structure this clean refuses to
bound the thing you point it at, it is telling you the difficulty lives in the other half of the
adèle — and here that half is the archimedean place, where the seven-sector measure actually lives.

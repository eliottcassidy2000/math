# The QR(-1) gate: the quadratic character of Z/q governs the LRC decorrelation constant

*kind-pasteur-2026-06-22-S34. A reflection on THREAD 2 of the q-uniform witness route.*

## The observation

For the doubled-prime family n = 2q, the LRC witness route's decorrelation is
carried by a single spectral object, the **Fourier-correlation constant**

> c^F_q = Σ_{m≠0} |chat(m)| |ghat(m)|,

where `chat(m) = a(m)/k` is the cluster's normalized autocorrelation and
`|ghat(m)| = |sin(πm/q)|/(π|m|)` is the Fourier coefficient of the width-1/q
sector kernel — the canonical "1/m decay modulated by the Z/q grid" (it vanishes
on q|m, peaks at m ≡ ±1 mod q).

Decompose c^F_q over the residue m mod q into a **quadratic-residue shell** and a
**non-residue shell**. The split is

> **exactly 50/50  ⟺  −1 is a non-residue mod q  ⟺  q ≡ 3 (mod 4).**

Verified, zero mismatches, on primes q = 3,5,7,11,13,17,19,23,29,31. The LRC(14)
prime q = 7 lies in the balanced class (7 ≡ 3 mod 4), and its residue set
QR(7) = {1,2,4} is exactly the Fano line / the Hamming(7,4) cyclic shift set.

## Why it is forced (and why it is not a coincidence)

The proof is a two-line involution argument. `|ghat|` is even in m; for a cluster
symmetric about its centroid the autocorrelation satisfies `a(m) = a(−m)`. So the
term `t(m) = |chat(m)||ghat(m)|` is even, `t(m) = t(−m)`. Group by class
`r = m mod q ∈ F_q^×`. Negation `r ↦ −r = q − r` is a fixed-point-free involution
on F_q^× that carries equal t-mass. The only question is whether it **swaps**
QR ↔ NQR or **preserves** each. It swaps iff −1 ∉ QR, i.e. iff
`(−1|q) = (−1)^{(q−1)/2} = −1`, i.e. iff q ≡ 3 mod 4. When it swaps, the two
shells are images of each other and carry identical mass — exactly 50/50.

This is the same Legendre-symbol fact (`−1` is a square mod p iff p ≡ 1 mod 4)
that decides whether p is a sum of two squares, whether Gauss sums are real, and
the sign in quadratic reciprocity. Here it surfaces as the balance condition of a
**measure-theoretic** decorrelation constant on the circle. The arithmetic of
Z/q is not a metaphor for the equidistribution — it is the mechanism.

## What it answers

The owner's THREAD-2 question (4) asked whether the decorrelation goes **through**
the Z/q sector structure (Gauss sums, QR for q=7) or **around** it. The committed
answer (HYP-2854/2856) said "through the Farey grid", generically. The QR(-1) gate
is the sharp form: it goes through the **quadratic character** of Z/q, and the
character's value at −1 is the exact switch. "Through", precisely.

It also dissolves the Cayley-Dickson red herring. One expects a property-loss
discontinuity at q = 7 (n = 14, the octonion level: 14 = 2·7 → 18 = 2·9 → 24 = 8·3
mirrors ℝ→ℂ→ℍ→𝕆). But every analytic constant of the family — c_q, c^F_q, R'_q,
the widening margin φ_q/m_P — moves **smoothly** through q = 7. The only fork is
2-periodic in q mod 4, governed by `(−1|q)`. The relevant tower is not
Cayley-Dickson; it is the splitting of −1, which alternates with period 4 and is
**orthogonal** to the doubling tower. q = 7 is special not because octonions lose
associativity, but because 7 ≡ 3 mod 4 puts −1 among the non-residues, balancing
the Fano line against its complement.

## The neighbouring number

There is a second, softer surprise next to the gate. The quasi-independence ratio

> R'_q = meas(coverSet^c ∩ G_P) / (meas(G_P)·(1−p0))

runs 3/4, 0.886, 1.005, 1.045 for q = 3,5,7,9: it **approaches 1 from below and
crosses 1 at q = 7**. Below the LRC prime the G_P holes sit preferentially inside
the covered region (anti-correlated, R' < 1); at q = 7 the lonely set and the
small-speed obstruction are almost exactly independent (R' = 1.005); above, faintly
positively correlated. The first open case is, to three decimals, the
independence threshold of its own proof. The decorrelation the route needs is not
merely *present* at n = 14 — it is *centred* there.

## The transferable principle

When a clean cancellation, a forced balance, or an exact 50/50 appears in an
analytic estimate over Z/n, look for the quadratic character before reaching for
the heavier machinery (Weyl sums, additive energy, Vitali coverings). A symmetric
weight plus the involution m ↦ −m reduces a measure identity to the single bit
`(−1|q)`. The same move recurs across the project — path-reversal, complement,
relabeling involutions in tournament proofs all collapse a sum to a parity. Here
the involution is negation on F_q^× and the parity is `(q−1)/2 mod 2`.

Files: `04-computation/lrc_witness_2q_qr_dichotomy_kpswf12.py`,
`04-computation/lrc_witness_2q_fourier_corr_kpswf12.py`, and their `.out`s in
`05-knowledge/results/`. Logged as a refinement of HYP-2854.

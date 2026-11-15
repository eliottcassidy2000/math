# Generic single-point lines: an exact formula, the measured law, and the precise missing estimate

Context: in a carrier row G(u) = c u + k p(u) + φ(u mod 2), every collision has a two-ended core (x, y) of span a,
a common suffix of length j and a prefix w. It satisfies the lattice-line equation (1) of
[TWISTED_LATTICES.md](../chanpak_cone_collision_placement_20261115/TWISTED_LATTICES.md):

    d w + δ r + [w odd] − ε [r even] = ((c/k) 2^j + 2^j − 1) δ + p(x) − p(y),

where δ = x − y, d = rev_a x − rev_a y, ε = (x mod 2) − (y mod 2), and r = rev_j(suffix).
[THM-5232](../../../01-canon/theorems/THM-5232-reverse-complement-collisions-every-odd-core-with-its-reverse-complement-collides-on-an-anti-diagonal-so-the-row-collision-bound-is-sharp-and-one-explicit-family-carries-linear-energy.md)
removed the explicit lattice-ray families. What remains of HYP-9752 is the **generic single-point lines**: cores
whose line meets the (w, r) box in at most one point per parity class.

## 1. Fixed signature: p(x) − p(y) is affine in the shared bits (PROVED)

Write σ_i = x_i − y_i ∈ {−1, 0, 1} for the signed difference and Z = {k : σ_k = 0} for the shared positions, with
shared bits z_k. Then

    p(x) − p(y) = c₀(σ) + Σ_(k∈Z) c_k(σ) z_k,
    c_k(σ) = − Σ_(j>k) σ_j 2^(j−k−1) − Σ_(i<k) σ_i 2^(k−i−1) − σ_(k+1) − σ_(k−1).

*Proof.* The recursion unrolls to

    p(u) = Σ_(i<j) (1 − u_i) u_j 2^(j−i−1) + #{i : u_i = 0, u_(i+1) = 1}.

In p(x) − p(y):
- a pair (i, j) with both positions shared cancels;
- a pair with i shared and j not contributes (1 − z_i) σ_j 2^(j−i−1);
- a pair with j shared and i not contributes −σ_i z_j 2^(j−i−1);
- the adjacency term is handled the same way.

Collecting the coefficient of each z_k gives the formula. ∎

- *Check.* Exhaustive for a ≤ 7, sampled for a ≤ 10: 38,999 evaluations, 0 violations (`affine_lemma.py`).
- *Reading.* The first sum is the part of δ above k, divided by 2^(k+1). The second is the top part of d, scaled.
  So the coefficients are the partial forward and reversed differences of the core: two binary sequences running in
  opposite directions.

## 2. The exact count is a sum of Riesz products (PROVED)

Fix a, j and a generic signature σ, so δ, d, ε and the c_k are fixed. Then the cores with signature σ are exactly
the z ∈ {0,1}^Z, and the collision count with w ∈ [W, 2W) and r ∈ [0, 2^j) is

    T_σ = Σ_(e_w, e_r ∈ {0,1}) ∫₀¹ D_W^(e_w)(θd) · D_(2^j)^(e_r)(θδ) · Π_(k∈Z) (1 + e(−θ c_k)) · e(−θ M_(σ,e)) dθ.

- D are the Dirichlet kernels of the two progressions (w and r in a fixed parity class).
- M_(σ,e) = ((c/k)2^j + 2^j − 1)δ + c₀ − e_w + ε[e_r = 0].

The z-sum factors into the **Riesz product** R_σ(θ) = Π_(k∈Z)(1 + e(−θ c_k)), because the shared bits enter
linearly. For j = 0 there is no r, and the count is

    T_σ = #{z : d | A_σ + Σ_k c_k z_k, quotient in the box}
        = (1/|d|) Σ_(h mod d) e(h A_σ/d) R_σ(h/d) · (box condition).

## 3. The measured law: generic cores obey random divisibility (NUMERICAL)

Take the row 2U + p, where every collision class occurs. Compare the exact generic collision count (all cores with
d ≠ δ and |d| > 8) with the prediction Σ_(core, r, e_w) 1[w_real ∈ [W, 2W)]/|d|.

Note: d is always odd and the right side always even, so w's parity matches automatically. An earlier version that
charged a factor 1/2 for parity was off by exactly 2.

| span a | suffix j | W | actual | predicted | ratio |
|---:|---:|---:|---:|---:|---:|
| 12 | 0 | 2³ | 152 | 145.5 | 1.045 |
| 13 | 0 | 2³ | 335 | 291.9 | 1.147 |
| 13 | 0 | 2⁵ | 86 | 77.9 | 1.104 |
| 11 | 2 | 2³ | 879 | 843.7 | 1.042 |
| 12 | 2 | 2⁵ | 712 | 623.4 | 1.142 |
| 10 | 4 | 2³ | 264 | 264.3 | 0.999 |
| 10 | 4 | 2⁵ | 1,906 | 1,744.0 | 1.093 |

- The largest-W cells have only 13–27 hits and ratios near 1.4 (noise).
- Divisibility itself is unbiased. In the strip, d | rhs at rate 0.00301 against mean(1/d) = 0.00276. The share of
  pairs with gcd(d, rhs) > 1 is 0.195, against 1 − 8/π² ≈ 0.189 for a random odd d and a random even rhs.
- So the generic part of HYP-9752 behaves exactly like the random model, within about 10%.

## 4. Why the proof must cancel across signatures, and the estimate that would finish it

Weight signatures by their number of cores, 2^|Z|. A typical signature then has |Z| ≈ a/2, while |d| ≈ 2^a.
- *Per-signature counts are below 1.* The main term of a single σ is about 2^|Z|/|d| ≈ 2^(−a/2), so a core-by-core
  asymptotic is impossible: each T_σ is 0 or 1.
- *Bounds on single σ are too weak.* Summing |R_σ| over σ gives at best Σ_σ (4/π)^|Z| = 3.27^a against a main term
  of order 2^a. Cauchy–Schwarz in θ gives √(W2^j)·3.414^a. That reproduces the concurrent 0.3928·L span threshold
  and no more.

What remains is a large-sieve statement over the moduli d_σ with digital weights.

**Hypothesis (generic-line large sieve; HYP-9782).** For the j = 0 sector, uniformly in W,

    Σ_(σ generic, span a)  (1/|d_σ|) Σ_(h ≢ 0 mod d_σ) e(h A_σ/d_σ) · R_σ(h/d_σ) · B_σ(h)  =  o( Σ_σ 2^|Z(σ)| β_σ/|d_σ| ),

where B_σ is the box condition written in Fourier form and β_σ its density. The analogue for j ≥ 1 adds the r-kernel.
- Together with THM-5232's families and the concurrent short-span bounds, it gives HYP-9752 in the form
  N(l, w) ≪ 2^l/(1 + w).
- The data of §3 say it holds with constant about 1.1.

**Why it is plausible.** Off the resonant frequencies the factors |1 + e(φ)| of a Riesz product have geometric mean
1, since ∫₀¹ log|2cos πφ| dφ = 0 (the Mahler measure of 1 + z is 1). So a typical |R_σ(h/d)| is e^(o(a)), not
2^|Z|. Large values occur only when h c_k/d is close to an integer for most k. That forces the forward and reversed
partial differences to resonate, which is the lattice-ray structure THM-5232 already removed.

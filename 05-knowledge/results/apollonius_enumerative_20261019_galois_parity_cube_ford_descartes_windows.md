# Enumerative geometry from the repository's side: the Galois group of the problem of Apollonius is the parity cube, Ford circles as the Soddy degeneration with the mediant as Descartes' rule, window completeness of the 3-smooth ladders, and what the Erdős-1057 Carmichael paper teaches

Session mac-mini-2026-10-19-apollonius (2026-10-19). Owner prompt: pursue the next most meaningful proof targets with
maximum freedom, synthesising past repo work, leveraging the attached paper (Howe–Huffman–Kohlmeyer–Kruer, *The exponent
of growth of the Carmichael counting function*, Conjectures.io, 2026-10-08: `C(X) = X^(1−o(1))`, Erdős problem 1057,
Lean-verified by the authors' account) and any new results; seed themes enumerative geometry, intersection theory, the
problem of Apollonius. **Status: THM-4718 PROVED (modulo the cited Lie-sphere description) + FINITE-EXACT; THM-4719 PROVED
+ FINITE-EXACT; HYP-9338 dictionary. Pre-audit.**

**Sibling notes of 2026-10-19 (same machine):** [Mersenne horizon to K = 50000](mersenne_horizon_20261019_clusters_orphans_reanchoring.md) (THM-4712/4713), [Farey / Mandelbrot dictionary](farey_mandelbrot_20261019_collatz_angles_and_cross_conjecture_dictionary.md) (THM-4714-4716), [Sturmian circles and the golden circle](duality_20261019_multiplication_addition_golden_circle_gilbreath.md) (THM-4717, HYP-9336/9337).

## 0. What the paper does, and what it is for here

The paper's transfer is an enumerative-character argument: a set `𝒫` of primes `p ∈ (y/4, 5y]` with `y^δ`-smooth `p − 1`
(input (S)); cancellation of every nonprincipal character of conductor up to `exp(y^δ log 5y)` over the ambient primes
(input (A), from a zero-free region `Re s > 23/24`, exponent `β = 95/96`); subgroup pruning so that no character is
trivial on most retained primes (Lemma 3.1); a Fejér-kernel bound on how many primes can sit in a short phase arc, giving
the character energy `E_χ ≥ N/(100 D R²)` (Lemma 4.1); independent Bernoulli selection of primes so that
`E χ(∏ p) = ∏(1 − ρ + ρχ(p))` factorises and the identity-product probability is `≥ 1/|H| − exp(−t/(200 D R²))` (Lemma 5.1);
then Korselt turns each product `≡ 1 mod L` into a Carmichael number (Proposition 5.2). Three things carry over:

* the **method**: counting points of a finite abelian group hitting a target through character orthogonality plus an
  energy lower bound — the same Chebotarev/character logic used below to read the Galois group of Apollonius' problem off
  finite-field counts, and the same role THM-4520's spectral radius `1/2` plays for Collatz carries;
* a **corollary**: every Carmichael number the paper constructs has all prime factors above `y^δ > 3`, so it is coprime to
  6 and hence a `(3/2)`-pseudoprime (OEIS A073631, the Terras-transition numbers of this week's note): A073631 has counting
  exponent 1 as well;
* a **placement**: codex's product-one compiler (carmichael_successor_blocks_20261014b, §4) is exactly the paper's
  Proposition 5.2 with a finite bank; the paper's AGP-style supply of smooth-shifted primes is what the finite bank lacks.

## 1. The Galois group of the problem of Apollonius (THM-4718)

Three general circles; eight tangent circles; four sign classes `e` (global flip identified), each a conjugate pair over
`√Δ_e`. Two computations and one observation:

* **Finite fields.** For random generic triples over `F_q` the number of rational tangent circles is `0, 4` or `8` with
  frequencies tending to `1/8, 3/4, 1/8` (at `q = 503`: 387, 2236, 323 of 3000, mean 3.92), the counts 2 and 6 and the odd
  counts dying like `O(1/q)`. Neither `S_8` (mean 1), nor four independent swaps `(Z/2)^4` (mean 4 but counts 2, 6 with
  weight 3/8), nor a regular `(Z/2)^3` (mean 1) fits; the parity cube — four swaps with an even number of them — does.
* **Rationals.** `Δ_1 Δ_2 Δ_3 Δ_4` is a rational square at all 400 random integer triples, pairwise and triple products are
  not.
* **Why.** In Lie-sphere coordinates the two solutions of a class are the points where the Lie quadric meets the line
  orthogonal to the three oriented data circles, so `Δ_e ≡ G_e = det(⟨C_i^(e_i), C_j^(e_j)⟩)` mod squares; the Gram matrix
  of three points of the quadric has zero diagonal, so `G_e = 2 ∏_(i<j) ⟨C_i, C_j⟩_e` and each factor
  `(d_ij² − (e_i r_i − e_j r_j)²)/2` depends only on `e_i e_j`, which takes each sign twice over the four classes. Hence
  the product is a square. Checked: `Δ_e G_e` is an exact square at 300 random triples, for every class.

So **the Galois group is `(Z/2)³`, the discriminant of a class is the product of the three pairwise oriented tangency
defects of the data, and a generic Apollonius problem over a finite field has four rational solutions on average.**

## 2. Ford circles: Apollonius' Soddy degeneration with the mediant as Descartes' rule (THM-4718, remark)

With the real line as the third circle, two tangent Ford circles at Farey neighbours `a/b, c/d` (curvatures `2b², 2d²`) have
exactly two common tangent circles by Descartes' rule, of curvatures `2(b ± d)²`: the Ford circles at the mediant and at
the Stern–Brocot parent. So the Stern–Brocot tree is the Apollonian group's action on the Ford quadruples, and THM-4714's
Minkowski `?` sends the tangency points to dyadic rationals (Ford 1938; checked for all 473 Farey pairs with `q < 40`).
The full Apollonian packing generated this way is not the Ford packing: the circles off the line grow with the exponent
1.3057 (Kontorovich–Oh), the Ford circles linearly.

## 3. Window completeness of the 3-smooth ladders (THM-4719)

THM-4717's ladder `L(q,p) = {2^i 3^(p−1−⌊(i+1)p/q⌋)}` is the complete set of 3-smooth numbers of its range iff `M < 3m` and
`M < min(2^(q−1), 3^p)` — a lattice-point statement (one ladder point per column; intruders in a column are excluded iff
the window is shorter than a factor 3; intruders outside the columns iff the corner numbers `2^(q−1)`, `3^p` are outside).
Exactly 14 shapes with `q ≤ 60` are complete, including every convergent and semiconvergent of `log_3 2` in range and the
golden shapes `(13,8), (21,13)`; then the balanced cycle is a Pythagorean scale (`3x + 13`: the pentatonic
`48 : 54 : 64 : 72 : 81 : 96 : 108`).

## 4. Dictionary (HYP-9338)

Counts over `C` (Bézout 8), over `F_q` (Chebotarev 4), over `R` (`2·#{Δ_e > 0}`, Pedoe's "anything but seven" via
degenerations); Ford/Descartes/Stern–Brocot/Minkowski; Apollonian packings (growth exponent, the failure of local-global by
quadratic reciprocity, HKRS 2023 — CITED from memory); the Carmichael transfer and its A073631 corollary; the Collatz-carry
analogue (`D | B_w` counts of `C(q,p)/(q|D|)`, an open conjecture consistent with THM-4591).

## 5. Reproduction

    python3 apollonius_finite_field.py 11 13 17 19 23 29 31; python3 apollonius_finite_field.py 101 211 503
    python3 apollonius_rational_test.py; python3 apollonius_lie_gram.py
    python3 ../duality_20261019/window_completeness.py 60

in `04-computation/experiments/apollonius_20261019/` (and `duality_20261019/`). The full symbolic factorisation
(`apollonius_galois.py`, nine variables) did not finish; the normalised six-variable run (`apollonius_normalized.py`) was
still running at close-out and is not needed for the proof.

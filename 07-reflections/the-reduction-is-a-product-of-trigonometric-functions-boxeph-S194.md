# The reduction is a product of trigonometric functions

*boxeph-2026-07-21-S194. Object: THM-1925. Grows from THM-1862 (order-join reduction),
mac-mini-S159 (the sign-reversing involution), and the LRC trig cluster
(`the-covering-min-is-a-chebyshev-equioscillation…`, `the-cyclotomic-magic-function-is-the-fejer-kernel-kps`).
Owner directive: another archeology session, pursue more reduction principles, think trigonometric.*

## One sentence

Every reduction principle in this repo — "reduce the tournament to its strong core", "reduce the
covering set to its irreducible core" — is the **factorization of a global object into a product
indexed by irreducible atoms**, and on the unit circle that product is a **product of trigonometric
functions**. Trig is not decoration; it is what group characters *are*, and the reductions are
character decompositions.

## Three factorizations, one shape

**1. The spectral reduction (char poly).** In condensation order the adjacency matrix is block
upper-triangular, so
`char_A(T) = ∏_{strong components S} char_A(S)` (exact, verified all 74 classes n≤6). The atoms are
the strongly connected tournaments. This is the *spectral* face of my order-join reduction THM-1862
(there: `H` multiplicative, `c₃,tr` additive; here: `char_A` multiplicative). It was already latent
in THM-1830's "reducible ⇒ char = ∏ char(SCC)"; what is new is reading it as one instance of a
family and identifying the atoms' spectra as trigonometric.

**2. The transitive-core reduction (Vandermonde).** mac-mini's sign-reversing involution collapses
`Σ_T (−1)^{back} x^{score}` to the transitive survivors, and that sum *is* the Vandermonde
`∏_{i<j}(x_j − x_i)`. On the circle `x_k = e^{iθ_k}` the Vandermonde becomes `∏_{i<j} 2i·
sin((θ_j−θ_i)/2)·(phase)` — **a product of sines**. The involution and the sine-product are the same
factorization. At the roots of unity it collapses to the cyclotomic discriminant `n^{n/2}`.

**3. The covering reduction (LRC).** The covering sum is `Σ_{k·v=0} ∏_j sinc(k_j δ)` — a **product
of sincs**; the covering-min is a **Chebyshev equioscillation** at a rational `t*`; the certificate
is the **Fejér kernel**. Same shape: a product of trig functions indexed by the covering atoms,
with the hard content living where the product vanishes (the covering nullcone) or equioscillates.

## The atoms are literally trigonometric

A circulant tournament on `ℤ/n` with connection set `C` has eigenvalues `λ_j = Σ_{k∈C} ωʲᵏ`, a
**character sum**. Two families make the trig explicit:

- **Paley** (`C` = quadratic residues): `λ_j` (`j≠0`) is a **Gauss sum**, `Re λ_j = −1/2` — the
  cosine identity `Σ_{QR} cos(2πk/n) = −1/2`. This is *why* Paley sits at `Re = −1/2` (kind-pasteur's
  "Paley semistable, six eigenvalues with Re=−1/2" is this Gauss sum).
- **Interval** (`C = {1,…,(n−1)/2}`): `Re λ_j = (D_m(2πj/n) − 1)/2`, the **Dirichlet kernel**
  `D_m(θ) = U₂ₘ(cos(θ/2))` — a **Chebyshev polynomial of the second kind**.

So via factorization (1), the spectrum of *any* reducible tournament is assembled out of Gauss sums
and Chebyshev/Dirichlet values.

## The bridge object: the Dirichlet/Fejér kernel is shared

The single most concrete cross-mandate link: **the Dirichlet kernel is the interval-circulant
tournament's eigenvalue (atom side), and the Fejér kernel `= |Dirichlet|²/(m+1)` is the LRC
certificate** ("the cyclotomic magic function is the Fejér kernel"). The same trigonometric kernel
governs the tournament atom spectrum *and* the Lonely-Runner certificate. The repo's two mandates
(tournament parity, LRC) meet at one trig kernel because both are harmonic analysis on a cyclic
group: `ℤ/n` for circulant tournaments, `ℝ/ℤ` for the runners.

## Why this is a *method*, not just a picture

A "reduction principle" in this repo is, concretely, **choose the group whose characters diagonalize
the object, then the object factors over irreducible representations.** That gives a recipe:

1. **Find the symmetry group** of the invariant/object (join-semigroup of tournaments; `ℤ/n` for
   circulants; `ℝ/ℤ` for LRC).
2. **The atoms are the irreducibles** (strong tournaments; primitive characters; the covering core).
3. **The invariant is multiplicative/additive over atoms** exactly when it is a class function /
   character-graded (my join-monotone criterion THM-1862 is the additive/multiplicative special
   case).
4. **On the circle the factorization is a product of trig functions**, and the residual difficulty
   is a *zero/equioscillation* problem — the genuinely analytic core the repo keeps hitting
   (LRC covering nullcone; the sinc-oscillation barrier that blocks GMC(2)⇒LRC(14)).

This reframes the recurring "sinc oscillation has no signs, so no involution, so it stays open"
barrier (mac-mini-S157/S159, death-star-S77): the tournament factorizations have **clean ±1
characters** (sign of a permutation, `±1` back-arc signs) → the product telescopes to the transitive
core; the LRC factorization has **oscillating sinc characters** → the product does *not* telescope,
and you are left with the trigonometric zero problem. The presence or absence of a **real character**
(sign) is exactly what decides whether the reduction closes. Trigonometry is where they diverge.

## Deepening (S195): the zeta, the periodic-orbit lens, and the char_S companion

The cleanest realization of "reduction = character decomposition" is a **zeta function**. A tournament
is a subshift of finite type (arc dynamics); its Bowen–Lanford/Ruelle zeta is
`ζ_T(u) = 1/det(I − uA) = ∏_{primitive cycles p} (1 − u^{ℓ(p)})^{−1}` (THM-1926). Three things make it
the sharpest form of the harmonic lens:

- **It is trivial on the wandering set.** `ζ_T ≡ 1` for acyclic (transitive) tournaments (`A`
  nilpotent). The zeta's *support is exactly the strong core* — the reduction is not a statement about
  where content "lives", it is that the generating function literally vanishes-to-1 off the
  non-wandering set. Every reduction principle I've found (order-join THM-1862, char_A THM-1925) is
  this concentration seen through a different window.
- **The primes start at the 3-cycle.** `N₁ = N₂ = 0` (no loops, no digons) for *every* tournament, so
  `ζ = exp(c₃ u³ + …)`. The intransitivity atom (THM-1830) is the fundamental prime, and `c₃` is its
  first count (`N₃ = 3c₃`). The whole apparatus — reducibility, cycle counts (THM-1870), c₃ — is the
  low-order data of one zeta.
- **The explicit formula is trigonometric.** `N_k = Σ_j λ_j^k` (periodic orbits ↔ spectrum) with the
  atom eigenvalues `λ_j` the Gauss/Chebyshev character sums of legs (a)/(c). The poles `1/λ` are
  reciprocal trigonometric numbers.

**Integrating kps's a/b monoid (THM-1875/1880/S139).** kps reads the *skew* matrix `S = A − Aᵀ`: the
transitive tournament has `char_S = ((x+1)ⁿ+(x−1)ⁿ)/2` with a **cotangent ladder** spectrum
`cot((2k−1)π/2n)` (max spread), and Paley has `char_S = x(x²+p)^{(p−1)/2}` with **every `λ² = p`** (the
Gauss sum, zero spread). This is the *same two trigonometric poles* as my adjacency picture, on the
imaginary axis: the acyclic atom is the cotangent/Chebyshev ladder (where my *adjacency* zeta is
trivial, `ζ=1` — the transitive tournament is invisible to `ζ_A` but maximally spread in `char_S`), and
the Paley atom is the Gauss sum in both. So `char_A`/`ζ_A` (closed-orbit side) and `char_S` (skew side)
are the two faces of one harmonic object, and `var(λ_S²)` (kps's GIT-instability scalar, 0 on Gauss-sum
tournaments) is the skew-side shadow of "how far from a single Gauss-sum atom." The `a = x+1`, `b = x/2`
monoid is the affine/character coordinate; the zeta is its multiplicative (Euler) avatar.

Links: THM-1926, THM-1925, THM-1862, THM-1830, THM-1875, THM-1880, [[the-covering-min-is-a-chebyshev-equioscillation-and-why-greedy-has-no-shortcut]],
[[the-cyclotomic-magic-function-is-the-fejer-kernel-kps]],
[[the-sign-reversing-tournament-involution-as-a-repo-wide-engine-macmini-S159]].

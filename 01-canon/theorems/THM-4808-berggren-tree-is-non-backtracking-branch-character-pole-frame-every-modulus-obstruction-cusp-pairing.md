---
id: THM-4808
title: "The Berggren tree is a non-backtracking walk, and its branch character at the integer level. (G) On P1 of any field of characteristic not 2 (and over Q), the three Berggren branches U(t) = 2 - 1/t, A(t) = 2 + 1/t, D(t) = t + 2 on Euclid slopes t = m/n form exactly the Hashimoto non-backtracking operator of the graph whose vertices are the Klein orbits {t, -t, 1/t, -1/t} and whose dart reversal is t -> 2 - t. Over Q this graph is a 4-regular tree with one orbifold point: its vertices are the primitive Pythagorean triangles up to leg swap and sign, plus the degenerate triangle (1, 0, 1), and its edges are the Berggren parent-child pairs. The Berggren ternary tree is its non-backtracking walk from the dart (1,0,1) -> (3,4,5). The branch character chi = (-1)^(#A) = det of the word matrix satisfies: (a) chi = alpha n - gamma m, where (alpha, gamma) is the image of the pole (1, 0); (b) #A is computed by an automaton on the regular continued fraction of m/n (it counts the even heads); (c) for every modulus M the triple of (16M+2, 4M+1) is congruent to (3, 4, 5) mod M and has chi = -1, so no function of residues modulo a fixed M, and in particular no Dirichlet or Hecke character of finite conductor, determines chi; (d) sign(a - b) sees only the trailing A-run, and on the middle ray chi = -(a - b), the Pell norm sign; (e) the depth-k level identity: sum of chi x^#D y^#U over words of length k equals (x + y - 1)^k; (f) all non-root parabolic rays pair off with opposite chi and interleaved heights, so sum_{c <= X} chi = 1 + N_U(X) + N_D(X) + E(X), with the root-ray counts N_U + N_D = (1/2 + 1/sqrt 2) sqrt X + O(1) and E(X) a signed sum of bounded local terms."
status: >
  PROVED, elementary: (G), (a)-(f). The tree enumeration of primitive triples is classical (Berggren 1934, Barning
  1963, Hall 1970) and CITED; its theta-group form (U = D S, D = T^2) is classical as well. No priority is claimed for
  the non-backtracking reading. A literature check found no source for it, but the check was not exhaustive.
  FINITE-EXACT:
  (G) the operator identity B = NB(G_q) for every prime 5 <= q < 120. The Ihara-Bass correspondence lambda^2 - mu
  lambda + 3 = 0 was checked numerically; exact factorisation in PARI/GP shows that every eigenvalue other than 3 and
  +-1 has modulus exactly sqrt 3 for all primes q <= 211. G_q is Ramanujan (max nontrivial |mu| <= 2 sqrt 3) for every
  prime 5 <= q < 1000 except q = 331, 569, 719, 809, 887, where the excess is at most 1.33% (q = 719). The 331
  exception is confirmed exactly: a degree-82 factor has a root off the circle.
  (a), (b), (d) on all 456,137 primitive triples with m <= 1500; (c) for M = 1..2000; (f) the exact identity at
  X = 10^3, 10^5, 10^6, 10^7, and S(X) to 10^10.
  CITED: uniform expansion of the quotients G_q follows from Selberg's 3/16 theorem for congruence subgroups (Lubotzky,
  Discrete Groups, Expanding Graphs and Invariant Measures, 1994), since <U, A, D> generates the extended theta group,
  of index 3 in PGL2(Z). The converse bit of Paley 7 is twin_pythagorean_124_20261026/finite_field section 5.
  The 16-state cover of actual Paley tournaments, the determinant as a monoid evaluation at (1, -1, 1), the
  non-descent of the bit through the radical, and the mod-15 forward observer are CITED from the concurrent
  growing_sieve_path_20261026/field (audited PASS). That packet keeps the bit as an extra state and does not show
  non-descent through residues; (c) supplies that, at every modulus. Pre-audit at filing.
session: mac-mini-2026-10-25-pyth (2026-10-26)
scripts:
  - 04-computation/experiments/berggren_sieve_parity_20261026/branch_parity.py (+ .out): (a) pole frame, (b) automaton, (c) every-modulus family, (d) sign law
  - 04-computation/experiments/berggren_sieve_parity_20261026/berggren_hashimoto2.py (+ .out): (G) multiset out-law, B == NB(G_q), Ihara roots, Ramanujan check q < 120
  - 04-computation/experiments/berggren_sieve_parity_20261026/berggren_ihara.gp (+ .out), berggren_exact_331.gp (+ .out), berggren_charpoly_exact.py (+ .out): exact spectra
  - 04-computation/experiments/berggren_sieve_parity_20261026/berggren_ramanujan_scan.py (+ .out): primes q < 1000
  - 04-computation/experiments/berggren_sieve_parity_20261026/ray_pairing.py (+ .out): (f) the exact identity
  - 04-computation/experiments/berggren_sieve_parity_20261026/parity_height_sum.c, parity_height_remainder.c, parity_height_big.c (+ .out): S(X), E(X) to 10^10
  - 04-computation/experiments/berggren_sieve_parity_20261026/parity_candidates.py, parity_xor_search.py, companion.py (+ .out): no character formula (XOR search over 18 arithmetic features), companion box
related:
  - 05-knowledge/results/twin_pythagorean_124_20261026/finite_field/README.md     # the F8 marking, F1 (PSL/PGL swap), the mod-7 endpoint hostile
  - 05-knowledge/results/twin_pythagorean_124_20261026/residue_obstruction/README.md  # the CRT twin obstruction (fixed moduli)
  - 05-knowledge/results/twin_pythagorean_124_20261026/tree/RESULT.md             # L(H) = J on the U and D rays
  - 01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md
  - 05-knowledge/results/growing_sieve_path_20261026/field/README.md            # 16-state Paley cover, monoid evaluation, radical non-descent
  - 05-knowledge/results/growing_sieve_path_20261026/sieve/README.md            # W1-W3: the exact growing-wheel twin count
  - 05-knowledge/hypotheses/HYP-9418-legal-path-sub-observers-berggren-sieve-levels-and-the-branch-character.md
---

# THM-4808 — The Berggren tree is a non-backtracking walk; the branch character at the integer level

## Conventions

- **Euclid pairs.** A pair `(m, n)` has `m > n ≥ 1`, `m` and `n` coprime and of opposite parity. Its triple is
  `Φ(m,n) = (a, b, c) = (m² − n², 2mn, m² + n²)`, and its slope is `t = m/n`.
- **Branches.** They act on columns as
  `U = [2 −1; 1 0]`, `A = [2 1; 1 0]`, `D = [1 2; 0 1]`, and on slopes as `U(t) = 2 − 1/t`, `A(t) = 2 + 1/t`,
  `D(t) = t + 2`.
- **Words.** Every primitive triple is `Φ(M_w (2,1)ᵀ)` for a unique word `w` (applied root first). This is the
  classical tree (Berggren, Barning, Hall).
- **Branch character.** It is `χ(w) = (−1)^#A(w) = det M_w`. Over `F_7` it is the swap bit of the two Steiner systems,
  and with the pole fixed it is the converse bit of Paley 7 (finite_field §§3–5).
- **Klein group and reversal.** `V = {t, −t, 1/t, −1/t}`, generated by `ι: t ↦ −t` and `S: t ↦ −1/t`. The reversal
  is `τ = D∘ι: t ↦ 2 − t`. It is an involution because `ιDι = D⁻¹`.
- **Theta-group form.** The branches factor through it as `U = D∘S`, `A = D∘S∘ι` and `D = D`.

## (G) The non-backtracking structure

For every point `t` of `P¹(K)` (`K` a field of characteristic ≠ 2, or `Q`), let `M(t)` be the multiset
`[Ut, At, Dt, τt]`. It is the multiset `D(V·t)` (each element counted `4/|V·t|` times). It depends only on the orbit
`V·t`, and it contains `τt`.

Make a graph `G`:
- its vertices are the `V`-orbits `O`;
- vertex `O` has outgoing darts `D(O)`;
- `tail(f)` is the orbit `O` with `f ∈ D(O)`, and `head(t) = V·t`.

Then:
- `τ(t) ∈ D(V·t)`;
- `tail(τt) = head(t)` and `head(τt) = tail(t)`, because `V·Dιt = V·D⁻¹t`.

Hence

    B = P_U + P_A + P_D  =  the Hashimoto non-backtracking operator of G:
    B[e, f] = (multiplicity of f in M(e)) − [f = τ(e)].

*Proof.* `Ut = D(St)`, `At = D(Sιt)`, `Dt = D(t)`, `τt = D(ιt)`, and `{t, ιt, St, Sιt} = V·t` with multiplicity. ∎

**Over Q.** A dart `t` (a slope) runs from the triangle with slope-orbit `V·D⁻¹t` to the triangle with slope-orbit
`V·t`. At `(3,4,5)`, whose slope orbit is `{±2, ±1/2}`, the four outgoing darts are:
- `4` → `(15,8,17)`;
- `5/2` → `(21,20,29)`;
- `3/2` → `(5,12,13)`;
- `0` → the degenerate triangle `(1,0,1)`.

So the Berggren tree is the non-backtracking walk, from the dart `(1,0,1) → (3,4,5)`, on a 4-regular tree. That tree's
vertices are the unsigned primitive triangles plus `(1,0,1)`, which is an orbifold vertex carrying a doubled edge and a
doubled loop. The **pole** `t = ∞` of the F8 marking is that degenerate vertex.

**Mod q (q an odd prime).** `G_q` is the reduction: about `(q+1)/4` vertices of degree 4, with small orbits counted
with multiplicity.
- **Ihara–Bass.** Every eigenvalue `λ ≠ ±1` of `B` satisfies `λ² − μλ + 3 = 0` for an eigenvalue `μ` of the
  adjacency of `G_q`.
- **Exact spectra, q ≤ 211.** Every eigenvalue other than `3` and `±1` has modulus exactly `√3`. Equivalently,
  `|μ| ≤ 2√3` for `μ ≠ 4`, which is the Ramanujan property.
- **Primes below 1000.** `G_q` is Ramanujan except at `q = 331, 569, 719, 809, 887`, with a largest ratio of 1.0133.
- **Uniform expansion in `q`** is CITED: Selberg's 3/16 theorem plus Lubotzky's transfer, for a finite generating set
  of a finite-index subgroup of `PGL₂(Z)`.

**Why this matters for sieves.** A growing modulus sees the tree through `G_q`. Expansion makes the depth-`k` residues
equidistribute in every residue class of slopes. Hence fixed moduli see only congruence-forced predicates (part (c)
and the CRT twin obstruction). A level growing with height is what an affine sieve on the tree uses. See HYP-9418 (4)
and the note.

## (a) The pole frame

Let `(α, γ) = M_w (1, 0)ᵀ` be the image of the pole. Then

    χ(w) = det[(α, γ) | (m, n)] = αn − γm.

*Proof.* `det[M(1,0)ᵀ | M(2,1)ᵀ] = det M · det[(1,0)ᵀ | (2,1)ᵀ] = det M`. ∎

Both candidates `α ≡ ±n⁻¹ (mod m)` lie in `(0, m)`. The bit is which side of `m/n` the pole image `α/γ` lies on: the
penultimate convergent of the theta-group expansion. FINITE-EXACT, `m ≤ 1500`: `0 ≤ α ≤ m` and `0 ≤ γ ≤ n` always.

## (b) The continued-fraction automaton

Write the regular continued fraction `m/n = [c₀; c₁, …, c_L]`, with `c_L ≥ 2` when `L ≥ 1`. Run the digit automaton:

- **Even head `h` with digits remaining.** Do `(h−2)/2` steps of `D⁻¹`, then one `A⁻¹`. Count one, and the next head is
  `c₁`.
- **Odd head with at least two further digits.** Do `(h−1)/2` steps of `D⁻¹`, then `c₁` steps of `U⁻¹`. The next head
  is `c₂ + 1`.
- **Stop** at an even integer slope, or at `[odd; c₁]`.

`#A` equals the number of counted even heads.

*Proof.* The three inverse branches act on digits as follows:
- `D⁻¹: [c₀; …] ↦ [c₀−2; …]` when `t > 3`;
- `A⁻¹: [2; c₁, …] ↦ [c₁; …]`;
- `U⁻¹: [1; c₁, c₂, …] ↦ [1; c₁−1, c₂, …]` when `c₁ ≥ 2`, and `[1; 1, c₂, …] ↦ [c₂+1; …]`. ∎

Exact on all 456,137 triples with `m ≤ 1500`.

**Not a character.** No XOR combination of up to five of 18 arithmetic features matches the parity on `m ≤ 120`. The
features were Jacobi symbols of `(m, n, m ± n, c)`, classes of `a, b, c` modulo 16, sign of `a − b`, and inverses
modulo `m` or `n`. The best two-feature agreement is 0.61. Part (c) turns the residue case into a theorem.

## (c) Every fixed modulus loses χ

For every `M ≥ 1`:

    (m, n) = (16M + 2, 4M + 1) ≡ (2, 1)  (mod M),       word  A · D^(M−1) · U · D,       χ = −1,

and `Φ(m, n) ≡ (3, 4, 5)` componentwise mod `M`.

*Proof.* Descend:
- `t = (16M+2)/(4M+1) > 3`, so `D⁻¹` gives `(8M, 4M+1)`;
- `1 < t < 2`, so `U⁻¹` gives `(4M+1, 2)`;
- `M − 1` steps of `D⁻¹` give `(5, 2)`;
- `A⁻¹` gives `(2, 1)`. ∎

The root has `χ = +1`. So `χ` is not a function of `(a, b, c)`, of `(m, n)`, or of `t` modulo any fixed `M`. In
particular it is not a finite product of Dirichlet or Hecke characters of finite conductor evaluated at `m + ni`, the
sides, or the slope.

This is the every-modulus form of the packet's mod-7 depth-two hostile (`UU` against `A`). The packet's hostile
compares projective conic points; here the Euclid vectors agree exactly.

Group-theoretically, `R = D U D⁻¹ A = [−1 4; 0 1]` fixes `(2, 1)` and has determinant `−1`. The word above is
`R mod M`.

## (d) What the sign of a − b sees

Put `Q(m,n) = a − b = (m−n)² − 2n² = N((m−n) + n√2)` and `P(m,n) = a + b > 0`. Then

    Q(A·) = −Q,     Q(D·) = P > 0,     Q(U·) = −P < 0,     Q(2,1) = −1.

Write `w = w′ X A^r`, with `X` the last non-`A` letter. Then `sign(a − b) = (−1)^r · (+1 if X = D, −1 if X = U or w′X`
is empty`)`. On the middle ray `A^k`, `|a − b| = 1` and `χ = −(a − b)`, the sign of the Pell norm. Away from that ray
the sign records only the trailing `A`-run, while `χ` records every `A`.

## (e) Level identity

    Σ_{|w| = k} χ(w) x^#D(w) y^#U(w) = (x + y − 1)^k.

In particular, at every depth exactly one more triple has even than odd branch parity.

## (f) Cusp pairing and the √X bias

For a node `w = (m, n)` with character `s`, compare the rays started at its children.

**D-rays.** `wU·D^j = (2m − n + 2jm, m)` and `wA·D^j = (2m + n + 2jm, m)` share the second coordinate. Their heights
interleave:

    c(wU D^j) < c(wA D^j) < c(wU D^(j+1)),

because `2n < 2m`.

**U-rays.** `wD·U^j = (m+2n + j(m+n), n + j(m+n))` and `wA·U^j = (2m+n + j(m+n), m + j(m+n))` lie on the same line
`x − y = m + n`, and they interleave because `n < m < m + 2n`.

The paired rays carry opposite characters. These rays are exactly the rays based at nodes ending in `U`, `D` or `A`.
So every parabolic ray except the root's two is paired. Every node is the root, an `A`-child, or an interior point of
exactly one ray. Hence

    S(X) := Σ_{c ≤ X} χ = 1 + N_U(X) + N_D(X) + E(X),
    N_U(X) = #{k ≥ 1 : (k+2)² + (k+1)² ≤ X},   N_D(X) = #{k ≥ 1 : (2k+2)² + 1 ≤ X},
    E(X)   = Σ_w χ(w) (δ_D(w,X) + δ_U(w,X) − [c(wA) ≤ X]),     δ ∈ {0, 1},

and `N_U + N_D = (1/√2 + 1/2)√X + O(1)`. ∎

The identity was checked exactly at `X = 10³, 10⁵, 10⁶, 10⁷`. Numerically:
- `S(X)/√X = 1.2158, 1.1532, 1.2421, 1.2194, 1.2941` at `X = 10⁸, 10⁹, 2·10⁹, 5·10⁹, 10¹⁰`, against
  `1/2 + 1/√2 = 1.2071`;
- `E(X)/√X` changes sign repeatedly (10 sign changes over 22 sampled heights from 100 to `10⁹`), with
  `|E| ≤ 0.21√X` at every sampled height from `10³` to `10¹⁰` (0.30 at `X = 100`).

The two pre-registered predictions at `2·10⁹` and `5·10⁹` were confirmed.

The two root rays are the D ray, which carries the twin images `T_p` and the centre identity `1 + p(p+2) = (p+1)²`,
and the U ray, which carries the prime-leg vertices `P(p)` and the classical abc family `P(3^(2^j))`. The parity bias
of the whole tree is located on these two rays and on no others.

## Exact non-connections (typed)

- **7 and 21.** Mod 7, with the pole fixed, `χ` is the converse bit of Paley 7, whose automorphism group has order 21
  (finite_field §5). `H(Tᵒᵖ) = H(T)` for every tournament. So no statement about Hamilton-path counts, in particular
  THM-1370's exclusion of 7 and 21, can detect `χ`. The numbers 7 and 21 meet `χ` only as the point count and group
  order of that converse. **NUMEROLOGY beyond this.**
- **abc.** The abc carriers on the tree are the twin centre identity and the classical `9^(2^j)` family. Both lie on
  the root rays, so `χ = +1` throughout, and abc data cannot see `χ` there. Off the rays, (c) rules out every
  residue-level arithmetic proxy.
- **F1.** On the free branch-word monoid, `χ` is the evaluation of the commutative word monomial at `(1, −1, 1)`, and
  augmentation loses it. This, and the non-descent of `χ` through radical support, are growing_sieve_path field §5.
  (c) adds the residue version: `χ` is not a function of any fixed-modulus data. So the F1-style multiplicative
  evaluations of the free prime monoid cannot see it at any finite level. Berggren words and prime factorisations are
  two free monoids with incompatible gradings. This is an exact loss, not a bridge.

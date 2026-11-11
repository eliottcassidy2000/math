# Stein modules of Keller pairs: Gaussian integration by parts between GMC(2) and JC(2) (kind-pasteur-2026-11-11-S24)

Owner directive (2026-11-11): merge prior work on the planar Jacobian conjecture with the proof of the
two-dimensional Gaussian moments conjecture, and define conjectures between GMC(2) and JC(2) that could form a
chain excluding high-degree counterexamples that rest on an undiscovered topological anomaly.

Canon: THM-5138, HYP-9668. Pre-registration: `PREREGISTRATION.md` (commits d701fad1cd and 90e0a6855d; predictions
S2-S6 all CONFIRMED). **JC(2) remains OPEN.**

Audits. First audit (`audit/AUDIT.md`, of the checkpoint 86ab0707d6 and of a draft addendum): PASS WITH CORRECTIONS,
no number wrong, six major corrections; all are applied in this version and logged in MISTAKE-764. Sections 5.2-5.4
and 6.3 were written after the first audit. **Second audit pending at this commit** (it will be
`audit/AUDIT-2.md`): until then Proposition 2.4, Lemma 5.2 to Corollary 5.7 and Proposition 6.3 are unaudited.

## 0. What comes out

GMC(2) and JC(2) live in the same ring `C[x, y]` and are statements about the same pair of operators.

For a Keller pair `(f, g)` (`{f, g} = 1`) take the dual fields `d_u = {., g}`, `d_v = {f, .}` (THM-4620), a centre
`z = (a, b)`, `U = f - a`, `V = g - b`, and the two operators `d_u - V`, `d_v - U`: Gaussian integration by parts in
the coordinates `(f, g)`. The **Stein kernel** `M_z` is the sum of their images and the **Stein space** is
`S_z = C[x, y] / M_z`.

1. For the standard pair `(x, y)` the Stein space is one-dimensional, its functional is the Gaussian expectation, and
   **GMC(2) (THM-2022) is the statement that the Stein kernel is a Mathieu-Zhao subspace** (this form of GMC is
   Zhao's). The same holds for every automorphism, by transport.
2. **Dimension count (Theorem 3.1; PROVED from a CITED comparison theorem).** For a Keller pair of generic degree
   `d`: `dim S_z = c + sum_i delta_i (deg_u S_i + deg_v S_i) = c + ((deg f + deg g)/gcd) sum_i delta_i rho_i >= 2d - 1`,
   where the `S_i` are the components of the non-properness curve. So **JC(2) says: the Stein functional of every
   Keller pair is unique.** A counterexample of generic degree `d` would carry at least `2d - 1` independent
   functionals that all satisfy the Gaussian integration by parts rules in `(f, g)`.
3. **Finite chaos (Theorem 4.1; sharp form due to the first audit).** If every polynomial has a finite Hermite
   expansion in `(f, g)` modulo the Stein kernel at one centre, then `S_z = 0` or the map is an isomorphism. For a
   Keller pair `S_z != 0`, so finite chaos at any single centre forces an automorphism.
4. **The Stein transform is injective (Theorem 5.4; PROVED, elementary).** The spaces `S_z` are the fibres of the
   single module `C[x, y]` over the polynomial ring `C[f - d_v, g - d_u]`, and this module is torsion-free. Hence a
   nonzero polynomial lies in `M_z` only for centres in a countable union of curves. Consequences:
   - **detection (Corollary 5.5):** at one centre that is algebraically independent over the field of coefficients,
     a Keller pair is an automorphism as soon as some iterate `d_u^r x` and some iterate `d_u^s y` lie in `M_z`
     (one direction, two test polynomials, first hits). So K-GMC, CV and SN each detect at such a centre;
   - **JC(2) <=> the Stein kernel of every Keller pair is a Mathieu-Zhao subspace (Corollary 5.7).** This is
     Zhao's Image Conjecture for the potentials `(f - a)(g - b)`; GMC(2) is the pair `(x, y)`.
5. **The ladder is not intermediate.** `automorphism => K-GMC => CV => SN => first hits` holds at every centre, and
   at a centre in general position every rung is equivalent to the automorphism property, because by item 4 it is
   then equivalent to the classical criterion "the dual field is locally nilpotent". Required of all Keller pairs,
   each rung is JC(2) again. What the Gaussian side adds to the classical criterion is confined to special centres.
   A two-direction version of detection with a moving centre, conditional on item 2, was pushed first by Session269
   (CONCURRENT, section 1).
6. **Examples (FINITE-EXACT, sections 6-7).** On the two planar non-proper etale maps available (a double cover on
   `C^* x C`, a slice of the `C^3` Keller map of THM-1300 on `C^2` minus a hyperbola) the count is confirmed
   (dimensions 2 and 11; linear potentials 2, 4, 6), finite chaos fails, and SN fails at every centre in general
   position.
7. **Open (HYP-9668).** Detection at a prescribed centre (it is false on surfaces); the Mathieu property for
   rank-one potentials with nonzero class of 1. One part of the hypothesis as first filed was false (section 8).

## 1. What exists

- **GMC(2)**: two real Gaussians, i.e. one pair `(Z, W = conj Z)` with `L(Z^a W^b) = a! [a = b]`. THM-2022 (PROVED;
  Lean endpoint `GMC2Main.gmc2`) through NC2 and one-variable DvdK (THM-2067, THM-2101, THM-2111). GMC(n) is false
  for `n >= 3` (THM-1480, THM-1490, THM-1500, THM-3290). No arrow from GMC(2) to JC(2) is in canon (MISTAKE-237;
  THM-2801: GMC(2) reaches the special image conjecture in dimension one only).
- **Images of first-order operators (the home of the conjectures used here).** W. Zhao, *Images of commuting
  differential operators of order one with constant leading coefficients*, J. Algebra 324 (2010) 231-247
  (arXiv:0902.0210): Conjecture 1.3 (Image Conjecture: the image of a commuting family is a Mathieu subspace),
  Theorem 2.5 (reduction to `{d_i - d_i(q)}`), Corollary 2.11 (the quotient is finite dimensional). With
  `{f, g} = 1`, `M_z = Im(d/dx - W_x) + Im(d/dy - W_y)` for `W = (f - a)(g - b)`, so **K-GMC of section 5 is the
  instance `q = (f - a)(g - b)` of the Image Conjecture** for two variables over `C`, and "`S_z` is finite
  dimensional" is Zhao's Corollary 2.11. GMC as a statement about such images: Zhao (J. Pure Appl. Algebra 214
  (2010)); Derksen, van den Essen, Zhao (Israel J. Math. 219 (2017): `GMC(2n) => SIC(n)`). The nearest published
  equivalence of JC(2) with a Mathieu statement about a Keller pair: van den Essen, Wright, Zhao, *Images of locally
  finite derivations of polynomial algebras in two variables*, J. Pure Appl. Algebra 215 (2011) 2130-2134
  (arXiv:1004.0521), Theorem 4.3 (JC(2) iff `Im D` is a Mathieu subspace for every divergence-free derivation `D`
  with 1 in its image). In the repo: the vacuum-Mathieu programme HYP-8350 (boxeph-2026-07-20-S166,
  `INDEX-HISTORICAL-THROUGH-2026-07-21.md:251`; "GMC(2) is a vacuum-Mathieu statement for A_1") and
  `THM-1830-conditional-map-what-assuming-the-2d-nullcone-conjecture-buys.md` (a superseded programme).
- **Dual fields** of a Keller map, the Weyl-algebra module they define (Bavula's twisted module), the kernel ladder,
  "cokernel of a dual field = relative de Rham cohomology", "`(f, g)` is an automorphism iff `D_f` is locally
  nilpotent iff `D_f` is surjective iff the generic fibre of `f` is a line" (KNOWN: Gurjar-Masuda-Miyanishi 2013,
  Cerveau 1997, Y. Stein 1989): THM-4620, THM-4626. The untwisted statements of this note (potential `W = f`) are
  there.
- **The `C^3` Keller map** with a cubic fibre: THM-1300; its Jelonek quartic and master identities: THM-1310,
  THM-1335; its coordinate fibrations: THM-4621. **Non-properness sets**: THM-2633 (the inertia groups normally
  generate the monodromy group; `1 <= k_D <= d - 1`), THM-4122 (statement (4), Nguyen Van Chau 2004: every component
  has normalization `A^1` and pole orders `rho_i (deg f, deg g)/gcd(deg f, deg g)`).
- **Session267** (same owner prompt, 2026-11-11 12:02-12:20; CONCURRENT;
  `05-knowledge/results/session267_gaussian_planar_20261111/`): consequences of NC2 for families, the row repair
  operator of Session254 as Gaussian rotation ("integration by parts gives `E[R P] = 0`", a one-dimensional
  rotational quotient), a positive cone theorem (a homogeneous sector of the Mathieu property), a shear-word
  compiler, a wreath law for the fixed-sheet gate, a conditional implication G5. These concern the standard
  Gaussian, not a Keller pair; there is no Stein functional of a Keller pair, no twisted cohomology, no count.
- **Session269** (2026-11-11 13:03-13:08; CONCURRENT, built on the pre-audit checkpoint of this note;
  `05-knowledge/results/session269_centres_ports_switches_20261111/jacobian/proof.md`): Theorem 1 (first hits in
  both pure directions, for every monomial, on an open set of centres, give mixed finite chaos at very general
  centres; proof by constructible membership, a Baire argument and a compiler of rational certificates through the
  connection `d/da + V`, `d/db + U`), Theorem 2a (one centre algebraically independent over the field of
  coefficients), Conditional Corollaries 3-4 (with the count and the finite-chaos criterion of this note: detection
  with a moving centre and the collapse of the universally quantified ladder to JC(2)). These were on `origin/main`
  before the corresponding statements here: **first-pusher credit for moving-centre detection in two directions and
  for the collapse of the ladder is Session269's.** The constructibility step of Theorem 5.4 below is theirs.
- **CITED.** A. Dimca, M. Saito, *On the cohomology of a general fiber of a polynomial map*, Compositio Math. 85
  (1993), Theorem 0.3: for a polynomial `W` on `C^n`, `H^k(Omega^., d - dW ^)` is the reduced cohomology of the
  generic fibre in degree `k - 1`. Already cited in the repo (`logic_lenses_audit_A_20261012.md:32`, as the theorem
  behind Cerveau 1997 and Gurjar-Masuda-Miyanishi 2013, Lemma 3.7, which is the case `W = f` of Theorem 3.1). For a
  smooth affine variety `X`: C. Sabbah, Nagoya Math. J. 141 (1996), Theorem 1.1 (`H^k(X, d - dW ^)` is
  `H^k(X, W^(-1)(t))`, generic `t`) and Corollary 1.2 (Euler characteristic form); C. Hien, C. Roucairol, Bull. SMF
  136 (2008) (period pairing). I read no proof; the first audit confirmed the statements from the sources, and they
  are checked on the examples below.

Searches. Repo: canon and notes for "twisted de Rham", "bifurcation set" (nothing), "Stein identity" (only the
Student-Stein rows of THM-4315 / THM-4328, another object), "Dimca" (the audit line above), "Mathieu" (HYP-8350).
Literature: one sweep by the first audit (sources in `audit/AUDIT.md`, section 6) found the items above and did not
find the product potential `(f - a)(g - b)`, the count, the bound `2d - 1`, or the finite-chaos criterion. No
priority is claimed for anything here.

Names. "Stein" refers to Charles Stein's lemma `E[g'(Z)] = E[Z g(Z)]`; it is neither the Stein kernel of Stein's
method nor Yosef Stein (1989) of the Jacobian literature.

## 2. The Stein module

Setting. `X` a connected smooth affine complex surface, `B = O(X)`, `F = (u, v) : X -> C^2` etale. For `X = C^2`
this is a Keller pair (scaled so that `{u, v} = 1`). The dual fields `d_u, d_v` are the derivations of `B` with
`d_u u = d_v v = 1`, `d_u v = d_v u = 0`; they commute. Centre `z = (a, b)`; `U = u - a`, `V = v - b`, `W = UV`.

    T_u := d_u - V,   T_v := d_v - U,   M_z := T_u B + T_v B,   S_z := B / M_z.

Classes in `S_z` are written `[h]_z`, or `[h]` when the centre is fixed.

**Proposition 2.1 (twisted de Rham; PROVED).** For `omega = alpha du + beta dv`,
`(d - dW ^) omega = (T_u beta - T_v alpha) du ^ dv`. Hence `S_z = H^2(Omega^.(X), d - dW ^)`.

*Proof.* `d omega = (d_u beta - d_v alpha) du ^ dv` and `dW ^ omega = (V du + U dv) ^ omega = (V beta - U alpha) du ^ dv`. QED

**Proposition 2.2 (identities; PROVED).** For all `k, h` in `B` and all `i, j, m, n >= 0`:

- (a) `[V k] = [d_u k]`, `[U k] = [d_v k]`;
- (b) `[U^i V^j] = i! [i = j] [1]`;
- (c) (iterated integration by parts) `[h V^m] = [d_u^m h]`, `[h U^m] = [d_v^m h]`;
- (d) (Hermite) `[d_u^m d_v^n h] = [h H_(n,m)(U, V)]`,
  `H_(n,m)(U, V) = sum_j (-1)^j j! C(n,j) C(m,j) U^(n-j) V^(m-j)` (Ito's polynomials);
- (e) (rotation) `(U d_u - V d_v) B` lies in `M_z`.

*Proof.* (a) is the definition. (b): `[U^i V^j] = [d_v (U^(i-1) V^j)] = j [U^(i-1) V^(j-1)]`, and `[V^j] = [d_u V^(j-1)] = 0`,
`[U^i] = 0` for `i, j >= 1`. (c): `[h V^m] = [d_u (h V^(m-1))] = [(d_u h) V^(m-1)]`. (d): put
`Psi_h(s, t) = [h e^(sV + tU)]` (formal). By (a), `d/ds Psi_h = Psi_(d_u h) + t Psi_h` and
`d/dt Psi_h = Psi_(d_v h) + s Psi_h`, so `Phi_h = e^(-st) Psi_h` satisfies `d/ds Phi_h = Phi_(d_u h)`,
`d/dt Phi_h = Phi_(d_v h)`, `Phi_h(0, 0) = [h]`, hence `Phi_h = sum [d_u^m d_v^n h] s^m t^n / (m! n!)`; and
`e^(sV + tU - st)` is the generating function of the `H_(n,m)`. (e): `[U d_u h] = [d_v d_u h] = [d_u d_v h] = [V d_v h]`. QED

Check: (b) for `i + j <= 4` and (d) for `m + n <= 3` on the slice of section 6 (15 and 21 cases, 0 failures; the
audit: 21 and 27).

**Proposition 2.3 (automorphisms; PROVED from THM-2022).** If `F` is an isomorphism onto `C^2`, then at every centre
`S_z = C [1]`, the Stein functional is `U^i V^j -> i! [i = j]`, i.e. `h -> L(h o F^(-1))` centred at `z`, and `M_z`
is a Mathieu-Zhao subspace of `B`.

*Proof.* `B = C[U, V]`; by 2.2 (b) `S_z` is spanned by `[1]`; the functional `U^i V^j -> i! [i = j]` kills
`T_u (U^i V^j) = i U^(i-1) V^j - U^i V^(j+1)` and likewise `T_v`, and is not zero. So `M_z` is its kernel, and
GMC(2), applied to the variables `X_0 = U`, `X_1 = V`, says that this kernel is a Mathieu-Zhao subspace. QED

**Proposition 2.4 (the Stein spaces are the fibres of one module; PROVED).** The operators

    E_a := u - d_v,      E_b := v - d_u

on `B` commute, and `M_z = (E_a - a) B + (E_b - b) B`. So `S_z` is the fibre at `z` of `B` regarded as a module over
the polynomial ring `C[E_a, E_b]`, and `d_u, d_v` act on this module as `d/da, d/db`
(`[d_u, E_a] = 1`, `[d_v, E_b] = 1`, `[d_u, E_b] = [d_v, E_a] = 0`).

*Proof.* `[u - d_v, v - d_u] = d_u(u) - d_v(v) = 0`, and `(E_a - a) h = -T_v h`, `(E_b - b) h = -T_u h`. QED

For the standard pair the heat operator `e^(d_x d_y)` carries `E_a, E_b` to the multiplications by `x, y`, so `B` is
free of rank one over `C[E_a, E_b]` and `[h]_z = (e^(d_x d_y) h)(a, b) [1]_z`: the class of `h` at `z` is the
Gaussian expectation of `h` centred at `z`. (The first audit found this description of the family of Stein spaces:
its family module `G` over the plane of centres is `B` with `a = E_a`, `b = E_b`; `G` is Bavula's twisted module
pulled back by the automorphism `(U, V, d_U, d_V) -> (U - d_V, V - d_U, d_U, d_V)` of the Weyl algebra, hence
holonomic.)

## 3. The dimension count

**Theorem 3.1 (PROVED from the CITED comparison theorem).** Let `F = (f, g)` be a Keller pair on `C^2` of generic
degree `d = [C(x,y) : C(f,g)]`, with non-properness set `S` (empty or a curve with components `S_i`), and let
`delta_i = d - #F^(-1)(p)` for a generic point `p` of `S_i`. For every centre `z = (a, b)`:

    dim S_z = b_1(Y) = c + sum_i delta_i (deg_u S_i + deg_v S_i) = c + ((deg f + deg g)/G) sum_i delta_i rho_i,

where `Y` is a generic fibre of `W = (f - a)(g - b)`, `c >= 1` the number of its connected components (`c = 1` for
generic `z`), `deg_u S_i`, `deg_v S_i` the degrees in `u` and `v` of an equation of `S_i`, `G = gcd(deg f, deg g)`
and `rho_i >= 1` the integers of THM-4122 (4). Moreover

    sum_i delta_i deg_v S_i = d - 1 + b_1(generic fibre of f),   sum_i delta_i deg_u S_i = d - 1 + b_1(generic fibre of g),

so `dim S_z >= 2d - 1`, and `dim S_z = 1` iff `F` is an automorphism. In particular `S_z != 0` at every centre.

*Proof.* By Dimca-Saito, `S_z = H^2(C^2, d - dW ^)` is `H^1(Y)` (the case `W = f` is Gurjar-Masuda-Miyanishi 2013).
`Y = F^(-1)(H_t)`, `H_t = {(u-a)(v-b) = t}`, a copy of `C^*`. Over `H_t \ S` the map `F : Y -> H_t` is finite etale
of degree `d`; over a point `p` of `H_t cap S` it has `d - delta_p` preimages. So
`chi(Y) = d chi(H_t \ S) + sum_p (d - delta_p) = - sum_p delta_p`, and `b_1(Y) = c - chi(Y)`. By THM-4122 (4) each
`S_i` is the image of a polynomial map `t -> (phi_i(t), psi_i(t))` with `deg phi_i = rho_i deg f / G`,
`deg psi_i = rho_i deg g / G`; so the closure of `S_i` in `P^1 x P^1` meets the boundary only at
`(infinity, infinity)`, never at the base points `(a, infinity)`, `(infinity, b)` of the pencil of curves `H_t` of
bidegree (1, 1), and for generic `t` the curve `H_t` meets `S_i` transversally in
`deg_u S_i + deg_v S_i = rho_i (deg g + deg f)/G` generic points. The same count for a vertical line `u = u_0` gives
`chi(f^(-1)(u_0)) = d - sum_i delta_i deg_v S_i`; the generic fibre of `f` is connected (THM-4620 (5): the ring of
constants of `D_f` is `C[f]`), so `b_1 = 1 - chi >= 0`. If `d = 1` then `F` is birational and etale, hence an
automorphism, `S` is empty and `Y` is `C^*`. If `d >= 2` then `dim S_z >= 1 + 2(d - 1) >= 3`. For generic `z` and
`t` the curve `Y` is a general member of the linear system spanned by `1, f, g, fg`, which has no base point and
maps `C^2` onto a surface; by Bertini it is irreducible, so `c = 1`. QED

**Corollary 3.2.** JC(2) holds iff every Keller pair has a one-dimensional Stein space at one centre. A
counterexample of generic degree `d` has `dim S_z >= 2d - 1` at every centre; with the cited bound `d >= 5` of
`CORE-PAPERS-PLANAR-JC.md` (Orevkov, Domrina) this is at least 9. (The audit notes that Zoladek's bound `d >= 6`
would give 11; that source is not in `CORE-PAPERS-PLANAR-JC.md` and is not used.)

**Remark 3.3 (surfaces).** For a smooth affine surface `X` the cited statement (Sabbah, Theorem 1.1) and the exact
sequence `H^1(X) -> H^1(Y) -> H^2(X, Y) -> H^2(X) -> H^2(Y) = 0` give
`dim S_z = b_1(Y) - rank(H^1(X) -> H^1(Y)) + b_2(X)`, with the same Euler characteristic count for `b_1(Y)`. Here
the position of the centre matters: it must avoid the lines `u = a`, `v = b` through boundary points
`(a, infinity)`, `(infinity, b)` of the closure of `S` (on the slice of section 6, `b = 0` is excluded).

**Remark 3.4 (degrees).** The count depends on the degrees of `f` and `g` through the reduced degree pair
`(deg f, deg g)/G` only, and the two fibre counts give
`(d - 1 + b_1(f-fibre))/deg f = (d - 1 + b_1(g-fibre))/deg g = (sum_i delta_i rho_i)/G`.

## 4. Finite chaos

Say that `F` has **finite chaos** at `z` if for every `h` in `B` there is `N` with `[d_u^m d_v^n h]_z = 0` whenever
`m + n >= N` (by 2.2 (d): the Hermite expansion of `h` in `(U, V)` is finite modulo the Stein kernel).

**Theorem 4.1 (PROVED; the sharp form and the trace argument are the first audit's).** Let `F : X -> C^2` be etale
of generic degree `d`, `X` connected. If `F` has finite chaos at `z`, then `S_z = C [1]`; and if `[1]_z != 0`, then
`F` is an isomorphism onto `C^2`. So finite chaos at `z` holds iff `S_z = 0` or `F` is an isomorphism.

*Proof.* Define `Theta : B -> S_z (x) C[s, t]`, `Theta(h) = sum [d_u^m d_v^n h] s^m t^n / (m! n!)`. Then
`Theta(d_u h) = d/ds Theta(h)`, `Theta(d_v h) = d/dt Theta(h)`, and by 2.2 (a) and the Leibniz rule
`[d_u^m d_v^n (V h)] = [d_u^(m+1) d_v^n h] + n [d_u^m d_v^(n-1) h]`, so

    Theta(V h) = (t + d/ds) Theta(h),      Theta(U h) = (s + d/dt) Theta(h).

The operators `s + d/dt`, `t + d/ds`, `d/ds`, `d/dt` on `C[s, t]` satisfy the relations of `U, V, d_u, d_v`; call
this module `P`. It is generated by 1 and `(s + d/dt)^i (t + d/ds)^j 1 = s^i t^j +` lower terms, so `P` is free of
rank one over `C[U, V]`, with `d_u, d_v` acting as `d/dU, d/dV`: `P = C[U, V]`.

*`Theta` is onto.* Let `N` be a subspace of `S_z (x) P` stable under the four operators and
`S' = N cap (S_z (x) 1)`. Let `n` in `N` and let `sigma s^I t^J` be a monomial of `n` of largest total degree. No
other monomial of `n` is divisible by `s^I t^J`, so `d^I/ds^I d^J/dt^J n = I! J! sigma` and `sigma` is in `S'`. The
element `(s + d/dt)^I (t + d/ds)^J sigma = sigma s^I t^J +` terms of smaller degree lies in `N`; subtracting it
lowers the number of monomials of top degree. By induction `N = S' (x) P`. For `N = Theta(B)` the constant term of
`Theta(h)` is `[h]`, so `S' = S_z`.

*The trace.* Let `K = Frac(B)`, `k = C(u, v)`. For a linear form `l` on `S_z`, `tau = (l (x) 1) Theta` is a
`C[U, V]`-linear map `B -> C[U, V]` that commutes with `d_u, d_v`; extend it to `K -> k`. The trace form of `K/k`
is nondegenerate, so `tau(h) = Tr(kappa h)` for one `kappa` in `K`; the trace commutes with the derivations, so
`Tr((d_u kappa) h) = d_U tau(h) - tau(d_u h) = 0` for all `h`, hence `d_u kappa = d_v kappa = 0` and `kappa` is a
constant. With `h = 1`: `tau = l([1]) Tr / d`. So `Theta(h) = [1] (x) Tr(h)/d`. If `[1] = 0` then `Theta = 0`, and
`Theta` is onto, so `S_z = 0`. If `[1] != 0` then `S_z = C [1]` and `Tr(h)` lies in `C[U, V]` for every `h` in `B`;
the power sums `Tr(h^j)` are polynomials, so `h` is integral over `C[u, v]`; `F` is finite and etale, a connected
covering of `C^2`, an isomorphism. Conversely, for an isomorphism `d_u, d_v` are locally nilpotent. QED

**Corollary 4.2.** (i) (PROVED from the CITED count, used only through `S_z != 0`.) A Keller pair has finite chaos
at one centre, any centre, iff it is an automorphism. (ii) (PROVED; no cited theorem.) The same at one centre that
is algebraically independent over the field of coefficients, where `[1]_z != 0` by Theorem 5.4.

Consistency (audit): the double cover of section 6.1 has dimension 2 and no finite chaos at centres with `a != 0`;
at the centres `(0, b)` it has `S_z = 0` and finite chaos holds trivially; the open immersion `C^* x C -> C^2` has
dimension 1 at (2, 3) and 0 at (0, 3).

## 5. The ladder and detection

### 5.1 One centre

For an etale `F` and a centre `z`:

- **(K-GMC)** `M_z` is a Mathieu-Zhao subspace of `B`: if `[P^m] = 0` for all `m >= 1` then for every `h`,
  `[h P^m] = 0` for `m` large. This is the instance `q = (u - a)(v - b)` of Zhao's Image Conjecture (section 1).
- **(CV), cone vanishing**: for every `h` there is `N` with `[h U^k V^(m+k)] = 0` and `[h U^(m+k) V^k] = 0` for
  all `m >= N`, `0 <= k <= m`.
- **(SN), one-sided Stein nilpotence**: for every `h`, `d_u^m h` and `d_v^m h` lie in `M_z` for `m` large
  (`SN_u`, `SN_v` for one direction).
- **(FH), first hit** (Session269): for the `h` under test, some `d_u^r h` with `r >= 1` lies in `M_z`.

**Proposition 5.1 (PROVED).** At every centre: isomorphism `=>` K-GMC `=>` CV `=>` SN `=>` FH, and finite chaos
`=>` K-GMC. For isomorphisms K-GMC is GMC(2).

*Proof.* Isomorphism `=>` K-GMC is Proposition 2.3. K-GMC `=>` CV: for `c` in `C` put `P_c = V (1 + c UV)`. By
2.2 (b), `[P_c^m] = sum_k C(m,k) c^k [U^k V^(m+k)] = 0` for `m >= 1`. So for every `c` there is `m_0(c)` with
`[h P_c^m] = 0` for `m >= m_0(c)`. The sets `Z_N = {c : [h P_c^m] = 0 for all m >= N}` are Zariski closed, increase
with `N` and cover `C`; a countable union of finite sets does not, so some `Z_N` is `C`, and then every coefficient
`[h U^k V^(m+k)]` of the polynomial `c -> [h P_c^m]` vanishes for `m >= N`. The other half uses `U (1 + c UV)`.
CV `=>` SN is `k = 0` and 2.2 (c). Finite chaos `=>` K-GMC by Theorem 4.1 (`M_z = B`, or an isomorphism). QED

What K-GMC does not give in this way is charge zero: the classes `[h (UV)^k]`, or the diagonal Hermite classes
`[h H_(n,n)(U, V)]`, for large `k, n`. The element `UV` is not null (`[(UV)^k] = k! [1]`), and the families
`V (1 + c (UV)^j)` reach large radial exponents only at nonzero charge. Finite chaos asks for charge zero as well.
This is the split of the proof of GMC(2), where the angular layer is DvdK and the gap was radial (THM-1645).

### 5.2 All centres at once

`K = Frac(B)` is a finite extension of `k = C(u, v)`, and `d_u, d_v` extend to `K`.

**Lemma 5.2 (a valuation at infinity of the target plane; PROVED).** Let `nu = -deg` be the valuation of `C(u, v)`
given by the total degree in `(u, v)`. Every extension `w` of `nu` to `K` satisfies

    w(d_u phi) >= w(phi) + 1,      w(d_v phi) >= w(phi) + 1      for all phi in K.

*Proof.* Put `tau = 1/u`, `s = v/u`. Then `C(u, v) = C(s)(tau)`, `nu` is the `tau`-adic valuation (a polynomial of
degree `n` with top form `p_n` is `tau^(-n) (p_n(1, s) + tau (...))`), and on `C(u, v)`

    d/dv = tau D_2,      d/du = - tau D_1 - s tau D_2,

where `D_1 = tau d/dtau` and `D_2 = d/ds` are the derivations with `D_1 s = 0`, `D_1 tau = tau`, `D_2 s = 1`,
`D_2 tau = 0` (check on `u = 1/tau`, `v = s/tau`). The completion of `K` at `w` is a complete discretely valued
field of residue characteristic zero containing `C(s)((tau))`; its residue field `k'` is a finite extension of
`C(s)` and lifts uniquely to a subfield containing `C(s)` (Hensel); there is a uniformizer `pi` with
`pi^e = c tau`, `c` in `k'`, `e` the ramification index (extract an `e`-th root of a 1-unit); so the completion is
`k'((pi))`. On `k'((pi))` define

    D~_1 (sum c_j pi^j) = sum (j/e) c_j pi^j,      D~_2 (sum c_j pi^j) = sum (c_j' + (j/e)(c'/c) c_j) pi^j,

where `'` is the extension of `d/ds` to `k'`. These are derivations, they do not lower the valuation, and they
restrict to `D_1, D_2` on `C(s)((tau))` (`D~_1 tau = tau`, `D~_2 tau = 0`, by `tau = pi^e / c`). A derivation of a
field extends uniquely to a finite separable extension, so `D~_1, D~_2` are the extensions of `D_1, D_2`, and on `K`
the extensions of `d/du, d/dv` are `d_u = -tau D~_1 - s tau D~_2`, `d_v = tau D~_2`. As `w(tau) = 1`, `w(s) = 0`,
both raise `w` by at least one. QED

**Lemma 5.3 (torsion-freeness; PROVED).** For a nonzero polynomial `P(a, b)` of degree `n` and `phi != 0` in `K`:
`w(P(E_a, E_b) phi) = w(phi) - n`. In particular `P(E_a, E_b) h != 0` for `h != 0` in `B`: the module `B` over
`C[E_a, E_b]` is torsion-free.

*Proof.* Multiplication by `u` or by `v` lowers `w` by exactly one, and `d_u, d_v` raise it by at least one. Expand
`E_a^i E_b^j = (u - d_v)^i (v - d_u)^j` into words: the word without derivative is `u^i v^j`, every other word
shifts `w` by at least `-(i + j) + 2`. So `P(E_a, E_b) phi = P_n(u, v) phi +` terms of valuation at least
`w(phi) - n + 1`, where `P_n` is the top form of `P`; and `w(P_n(u, v)) = -n` because `P_n(1, s) != 0`. QED

**Theorem 5.4 (the Stein transform is injective; PROVED).** Let `X`, `F` and `h != 0` in `B` be defined over a
subfield `k_0` of `C` (for a Keller pair: `f, g, h` in `k_0[x, y]`). If `a, b` are algebraically independent over
`k_0`, then `h` is not in `M_(a,b)`. Consequently, for every `h != 0` in `B` the set of centres `z` with
`[h]_z = 0` is contained in a countable union of proper algebraic subsets of `C^2`.

*Proof.* Let `B_0` be the `k_0`-form of `B`, with a `k_0`-basis `e_1, e_2, ...`, and `B_L` the span of
`e_1, ..., e_L`; the operators `E_a, E_b` preserve `B_0`. Suppose `h = (E_a - a) p + (E_b - b) q` with `p, q` in
`B`; then `p, q` lie in `C (x) B_L` for some `L`. In coordinates this is a linear system `N(a, b) (p, q) = h` whose
matrix has entries in `k_0[a, b]` and whose right side has entries in `k_0`. Since `a, b` are algebraically
independent over `k_0`, the field `k_0(a, b)` is a field of rational functions in two indeterminates, and the ranks
of `N` and of the augmented matrix do not change under the extension from `k_0(a, b)` to `C`. So the system has a
solution over `k_0(a, b)` with `a, b` indeterminates (this step is Session269's, section 3a there). Clearing
denominators, there are a nonzero `P` in `k_0[a, b]` and `p', q'` in `k_0[a, b] (x) B_0` with

    P(a, b) h = (E_a - a) p' + (E_b - b) q'      in C[a, b] (x) B.

Apply the linear map `Phi : C[a, b] (x) B -> B`, `a^i b^j (x) k -> E_a^i E_b^j k`. It kills
`(E_a - a)(a^i b^j (x) k) = a^i b^j (x) E_a k - a^(i+1) b^j (x) k` and likewise `(E_b - b)(...)`, because `E_a` and
`E_b` commute. Hence `P(E_a, E_b) h = 0`, against Lemma 5.3. For the last assertion take for `k_0` a countable
field of definition of `X`, `F`, `h`: the centres that are not algebraically independent over `k_0` form the union
of the countably many curves `{P = 0}`, `P` in `k_0[a, b]`. QED

Checks (`runs/inject.txt`; residues modulo the prime play the role of a centre in general position): for the
standard pair, the double cover and the slice, every monomial tested (21, 36, 31) has a nonzero class at three
pseudo-random centres and at (2, 3); at the special centres the classes do vanish: 18 of 21 at (0, 0) for the
standard pair (`x^i y^j`, `i != j`), all at (0, 3) for the double cover (`S_z = 0`), `[d_u x] = 0` at (0, 0) on the
slice.

### 5.3 Detection

**Corollary 5.5 (detection from first hits in one direction; PROVED).** Let `(f, g)` be a Keller pair with
coefficients in `k_0`, and let `z` be a centre whose coordinates are algebraically independent over `k_0`. The
following are equivalent:

1. `(f, g)` is an automorphism;
2. K-GMC at `z`; 3. CV at `z`; 4. `SN_u` at `z`, for `h = x` and `h = y` only;
5. FH at `z` in the direction `d_u`, for `h = x` and `h = y`: some `d_u^r x` and some `d_u^s y` (`r, s >= 1`) lie
   in `M_z`;
6. `d_u` is locally nilpotent.

The same holds with `d_v` in place of `d_u`, and with "at `z`" replaced by "at every centre of a set that is not
contained in a countable union of proper algebraic subsets" (a nonempty open set, a set of positive measure).

*Proof.* 1 => 2 => 3 => 4 => 5 is Proposition 5.1. 5 => 6: `d_u^r x` is a polynomial with coefficients in `k_0`
that lies in `M_z`; by Theorem 5.4 it is zero. Likewise `d_u^s y = 0`. A derivation that is nilpotent on generators
is locally nilpotent. 6 => 1 is THM-4620 (5) (KNOWN: the kernel of `d_u` is `C[g]`, `f` is a slice, so
`C[x, y] = C[g][f]`). For the last sentence: such a set contains a centre that is algebraically independent over a
countable field of coefficients. QED

**Corollary 5.6 (surfaces; PROVED).** Let `F : X -> C^2` be etale, defined over `k_0`, and `z` algebraically
independent over `k_0`. Then `SN_u` holds at `z` iff `d_u` is locally nilpotent on `B`, and SN holds at `z` in both
directions iff `F` is an isomorphism onto `C^2`; so K-GMC at `z` holds iff `F` is an isomorphism. More precisely,
if `d_u^m h != 0` for all `m` then `[d_u^m h]_z != 0` for all `m`.

*Proof.* The first and the last statement are Theorem 5.4 applied to the elements `d_u^m h`, `h` in `B_0`. If
`d_u` and `d_v` are both locally nilpotent, the slice theorem gives `B = Ker(d_u)[u]`, then
`Ker(d_u) = (Ker d_u cap Ker d_v)[v]`, and `Ker d_u cap Ker d_v = C` on a connected surface; so `B = C[u, v]`. QED

**Corollary 5.7 (JC(2) as an image conjecture; PROVED).** The following are equivalent:

1. JC(2);
2. for every Keller pair `(f, g)`, `Im(d/dx - (fg)_x) + Im(d/dy - (fg)_y)` is a Mathieu-Zhao subspace of `C[x, y]`
   (K-GMC at the centre 0 for all Keller pairs; equivalently at all centres, since `(f - a, g - b)` is a Keller pair);
3. CV, or SN, or `SN_u`, for all Keller pairs and all centres;
4. for every Keller pair, `C[x, y]` has rank one over `C[f - d_v, g - d_u]` (this item: PROVED from the CITED
   count and the CITED fact that a holonomic module is a vector bundle over a dense open set).

GMC(2) is statement 2 for the pair `(x, y)`. Statement 2 is the Image Conjecture of Zhao for the potentials `fg`;
the Image Conjecture itself fails for some other potentials (section 8). The equivalence of 1 with 2-3 uses GMC(2)
(THM-2022), Proposition 5.1, Theorem 5.4 and the classical criterion; it does not use the comparison theorem.
Session269's Conditional Corollary 4 is the two-direction form of this collapse, conditional on Theorem 3.1.

The rungs by generic degree, "`K-GMC_D`: every Keller pair of generic degree at most `D` satisfies K-GMC at every
centre", are restatements of "no Keller pair of generic degree `2 <= d <= D`"; `K-GMC_4` is the cited theorem of
Orevkov and Domrina and carries no further information.

| property of a Keller pair at one centre `z` | forces an automorphism |
|---|---|
| finite chaos | at every `z` (Corollary 4.2) |
| K-GMC, CV, SN, `SN_u`, first hits of `x` and `y` in one direction | at every `z` algebraically independent over the coefficients (Corollary 5.5) |
| the same | at a prescribed `z`: OPEN on `C^2`; false on surfaces (section 6.3: SN and CV hold at a special centre with `[1]_z != 0`; double cover at `a = 0`: `S_z = 0`) |

### 5.4 Periods (a second route; the draft audited as "Theorem 5.2")

`audit/AUDIT.md`, section 7, audited a draft (kept as `SECOND-PROOF-PERIODS.md`) that reaches two-direction
detection on an open set of centres differently: `B` over `C[E_a, E_b]` with `D_a = d_u`, `D_b = d_v` is a holonomic
module, hence a vector bundle with flat connection over a dense open set `Omega`; for a flat local frame
`lambda^(i)` of the dual bundle the periods `pi_h^(i)(z) = lambda_z^(i)([h]_z)` satisfy
`d/da pi_h = pi_(d_u h)`, `d/db pi_h = pi_(d_v h)`; SN on an open set makes them polynomials (Taylor series, an
uncountability argument), which is finite chaos, and Theorem 4.1 concludes. It is superseded by Corollary 5.5 as a
route to detection, and it is CONCURRENT with Session269's Theorem 1, which was pushed first and is elementary.
What it adds is the reading of the Stein functionals as flat sections: on the double cover
`d/da lambda_+-(x^k y^l) = ((k + l)/2) lambda_+-(x^(k-2) y^l) = lambda_+-(d_u x^k y^l)` and
`lambda_+-(x^(-2)) = 1/(2a)` is not a polynomial in `a`.

## 6. The planar non-proper etale maps in canon (FINITE-EXACT)

`scripts/stein.py` computes `O(X) / M` by exact linear algebra modulo the prime 2147483629, in the frame `dx ^ dy`
(`M' = Im(d/dx - W_x) + Im(d/dy - W_y)`; the class of `h` in `S_z` is the class of `h J`, `J` the Jacobian of
`(u, v)`), with elimination that puts pivots in the monomials of largest weight; `c(N, N')` is the codimension in
weight at most `N` of the span of the generators of weight at most `N'`. It is an upper bound for the rank of the
weight-`N` part in `S_z`, and a nonzero normal form certifies non-membership relative to the truncation `N'` only;
the first audit re-ran everything with an independent engine (another basis, two other primes, three filtrations).

**Controls** (dimension of `H^2(C^2, d - dW ^)`): `xy`: 1; `x^3 + y^3`: 4; `x + x^2 y` (no critical point; fibre
`C^*`): 1; `x^2 y`: 1; `fg` for two automorphisms: 1; audit: `x^2 y^2`: 2; `x^3 y^3`: 3; `x^3 y^2`: 1. All as the
cited theorem says.

**6.1 A double cover on `C^* x C`.** `(u, v) = (x^2, y/x)`, Jacobian 2; etale of degree 2 onto `C^* x C`, not
proper over `C^2` (the line `u = 0` is missed). Here `d_v = x d/dy` is locally nilpotent and
`d_u = (1/2x) d/dx + (y/2x^2) d/dy` is not.

**Proposition 6.1 (PROVED by hand).** For `W = (x^2 - a)(y/x - b)`, `a != 0`, the Stein space is two-dimensional,
with basis of functionals

    lambda_+-(x^k y^l) = l! Res_(x = +-sqrt(a)) [ x^(k+l+1) e^(b (x^2 - a)) / (x^2 - a)^(l+1) ],

and `lambda_+-(x^(-2) V^m) = (-1)^m m! / (2 a^(m+1))` for all `m >= 0`, while `[V^m] = 0` for `m >= 1`. So for
`a != 0`, `M_z` is not a Mathieu-Zhao subspace and `SN_u` fails; `SN_v` holds trivially. At the centres `(0, b)`,
`S_z = 0`.

*Proof.* On a generator `(d/dy - W_y)(x^k y^l)` the two terms cancel for `l >= 1` and the term for `l = 0` is regular
at `x = +-sqrt(a)`. On `(d/dx - W_x)(x^k y^l)` the expression under the residue is the `x`-derivative of
`l! x^(k+l+1) e^(b(x^2-a)) / (x^2-a)^(l+1)` (direct computation with `W_x = y + a y/x^2 - 2bx`,
`W_y = (x^2 - a)/x`). The two functionals are exchanged by `(x, y) -> (-x, -y)` and are independent. With
`s = x^2`: `lambda_+-(x^(-2-l) y^l) = (l!/2) Res_(s=a) [s^(-1) e^(b(s-a)) / (s-a)^(l+1)]`, and summing the
binomial expansion of `V^m = (y/x - b)^m` leaves `(1/2) (-1)^m m! / a^(m+1)`. Upper bound (audit): with
`D_x = d/dx - W_x`, `D_y = d/dy - W_y`,
`D_y(x^(i-1) y^(j+1)) - D_x(x^i y^j) = 2a x^(i-2) y^(j+1) - (i - j - 1) x^(i-1) y^j - 2b x^(i+1) y^j` lowers the
degree in `y` modulo `M'`, and `D_y(x^i) = a x^(i-1) - x^(i+1)` reduces the powers of `x` to `1, x`. At `a = 0`:
`x^(i+1) y^j = j x^i y^(j-1)` modulo `M'`, so everything is zero. QED

The two functionals are one half of the local Gaussian functionals at the two preimages of the centre (audit:
`lambda_+- = (1/2) L_(q_+-)`, `L_q(h) = sum_i i! h_ii` for the Taylor expansion `h = sum h_ij U^i V^j` at `q`;
the sum is finite because `B` is polynomial in `V`). So here there is no extra functional: `SN_u` is broken by the
local Gaussians themselves, because `1/u` is singular on the missed line
(`L_q(x^(-2) V^m) = (d/du)^m (1/u)` at `u = a`).

Machine (`runs/s2_toy.txt`, centre (2, 3)): `c(N, N') = 2` for `N < N'`; classes of `x^(-2) V^m` and `x^(-1) V^m`
nonzero for `1 <= m <= 6`; classes of `x U^m`, `x^(-2) U^m` zero. Topology: fibre `P^1` minus four points,
`b_1 = 3`, `b_1(X) = 1` (injective into `H^1(Y)`: `dx/x` has a residue at a puncture), `b_2(X) = 0`:
`3 - 1 + 0 = 2`.

**6.2 A slice of the `C^3` Keller map.** For `F = (F_1, F_2, F_3)` of THM-1300, the surface `F_1 = e` (`e != 0`) is
`X_e = C^2 \ {1 + xy = 0}` and `(F_2, F_3)` restricts to

    u = y + 3e x w,   v = x w + x w^2 - e x^3 w^3,   w = 1/(1 + xy),   Jacobian -2 w^3,

an etale map of generic degree 3 onto its image in `C^2`. Its non-properness curve is the slice of the Jelonek
quartic of THM-1310, `27 e^2 v^2 - 18 e u v + 16 e + u^3 v - u^2 = 0`, of bidegree (3, 2) in `(u, v)`, with two
sheets missing over a generic point (`scripts/slice_pred.py`; fibre sizes 3 off the curve and 1 on it: audit).

**Proposition 6.2 (FINITE-EXACT; `e = 1`).**

| potential | points of the curve on the generic fibre image | predicted dimension | computed |
|---|---|---|---|
| `(u - 2)(v - 3)` | 5 | `1 + 2*5 - 1 + 1 = 11` | 11 (`N <= 36`, `N' = 40, 46`; audit: `4 <= N <= 60`, `N' <= 70`) |
| `u` | 2 | `1 - (3 - 4) = 2` | 2 |
| `v` | 3 | `1 - (3 - 6) = 4` | 4 |
| `u + v` | 4 | `1 - (3 - 8) = 6` | 6 |

Other centres: 11 at (5, -7), (1, 1), (-3, 4), (7, 2), (-1, -2); 9 at (0, 0) and (audit) at (2, 0), (4, 0), where
`b = 0` is the horizontal asymptote `v ~ 1/u` of the curve and one intersection point is lost at infinity
(Remark 3.3). At the centre (2, 3): `[U^i V^j] = i! [i = j] [1]` for `i + j <= 4`; the normal form of 1 is not
zero; the normal forms of `w V^m`, `w U^m` (`m <= 5`) and `x V^m`, `x U^m` (`m <= 4`) are not zero (audit:
`m <= 9` in the direction `V`, `m <= 12` in the direction `U`, modulo the generators of weight at most 70, identical
at `N' = 46, 58, 70`); every Hermite class `[h H_(n,m)(U, V)]` computed for `h = w` and `h = x` (`n, m <= 4`) has a
nonzero normal form, while for `h = 1` only `(n, m) = (0, 0)` has. At the six centres above with dimension 11 the
normal forms of `w V^m`, `w U^m`, `x V^m`, `x U^m` (`m <= 3`) are not zero (`runs/sn_centres.txt`).

What is proved on the slice. (i) Finite chaos fails at every centre where `S_z != 0` (Theorem 4.1), in particular
where the cited count gives 11. (ii) `d_u` and `d_v` are not locally nilpotent (they do not kill the unit `w`), so
by Corollary 5.6, at every centre algebraically independent over `Q`, SN fails in both directions, and so do CV and
K-GMC: for at least one of the generators `h = x, y, w` all the classes `[d_u^m h]_z` are nonzero, and likewise for
`d_v`. (iii) At the rational centre (2, 3) the failure of SN, CV, K-GMC is not proved; no vanishing is observed.

Reading (HEURISTIC). `11 = 3 + 8`, where 3 is the number of preimages of the centre (the critical points of `W`):
three local Gaussians, and eight functionals attached to the missing sheets. The computations do not identify the
functionals that break SN.

Observation (audit, FINITE-EXACT). At the excluded centre (0, 0) the potential `uv` is invariant under
`(x, y) -> (-x, -y)` and some classes do vanish in the direction `d_u`: `[w V^m] = 0` for `m = 4, 6, 8`,
`[x V^m] = 0` for `m = 1, 3, 5, 7, 9`, `[y V^m] = 0` for `m = 3, 5, 7, 9`, while the classes `h U^m` stay nonzero.
Special centres behave differently.

(`runs/s3_s5_slice.txt`, `runs/s6_linear.txt`, `runs/centres.txt`, `runs/hermite_slice.txt`, `runs/inject.txt`;
`audit/a7_slice_stability.out`, `audit/a11_slice_centre00.out`.)

**6.3 A special centre where SN and CV hold without an isomorphism.** `X = C^2 \ {1 + xy = 0}`, `F` the inclusion
`(u, v) = (x, y)`, centre `z = (0, 0)`, `W = xy`. The map is not an isomorphism and `d_u = d/dx` is not locally
nilpotent on `B = C[x, y, w]`, `w = 1/(1 + xy)`.

**Proposition 6.3 (PROVED by hand; machine check `runs/special_centre.txt`).** At this centre `S_z` has dimension
2, `[1]_z != 0`, and for every `h` the classes `[d_u^m h]_z`, `[d_v^m h]_z` and `[h U^k V^(m+k)]_z`,
`[h U^(m+k) V^k]_z` vanish for `m` large: SN and CV hold. Finite chaos fails.

*Proof.* Grade `B` by the charge `deg x - deg y` (`w` has charge 0): `B_k = x^k R` for `k >= 0`, `B_k = y^(-k) R`
for `k < 0`, `R = C[t, 1/(1 + t)]`, `t = xy`. The operators `d/dx - y` and `d/dy - x` shift the charge by -1 and
+1, so `M_z` is graded. For `k >= 1`: `(d/dy - x)(x^(k-1) phi) = x^k (phi' - phi)` and
`(d/dx - y)(x^(k+1) phi) = x^k (k phi + (d/dt - 1)(t phi))`, so `M_z` contains `x^k R`; likewise for `k <= -1`. In
charge 0, `M_z cap R = (d/dt - 1)(t R)`, of codimension `dim R/(d/dt - 1) R + dim R / t R = 1 + 1` (`d/dt - 1` is
injective on `R`), and it does not contain 1 (`psi' - psi = 1` has no solution in `t R`). Now `d_u^m h` has charge
`k - m` and `h U^k V^(m+k)` has charge (charge of `h`) `- m`: nonzero for `m` large. Finite chaos fails by
Theorem 4.1. QED

Machine: `c(N) = 2`; the classes of `d_u^m h`, `d_v^m h` (`h = w, w^2, xw, yw, x^2 w`, `m <= 6`) are zero exactly
when the charge is not zero; the classes of `(d_u d_v)^m w` (`m <= 3`) are not zero. At the centre (2, 3) of the
same map the dimension is 3 (as Remark 3.3 predicts: `1 + 2 - 1 + 1`) and none of these classes vanishes. The
centre (0, 0) is special because the missed curve `uv = -1` is itself a member of the pencil of hyperbolas
`uv = t`: no hyperbola of the pencil meets it. For a Keller pair on `C^2` this cannot happen (every component of
the non-properness curve leaves the plane at `(infinity, infinity)` only, section 3).

## 7. Where the extra functionals sit (HEURISTIC; not used in any proof)

For a Keller pair restrict to a generic horizontal line `v = v_0`. The curve `C = g^(-1)(v_0)` is etale over the
`u`-line with punctures over the points of `S` on that line, `d_u` is `d/du` on `C`, and the fibrewise twisted
cohomology `H^1(C, d - (v_0 - b) du)` has dimension `sum_i delta_i deg_u S_i` (the missing sheets; a puncture of
ramification `e` carries `e` of them). Integrating fibrewise residues along a cycle `gamma` in the `v_0`-line gives
Stein functionals of the form

    lambda(h) = int_gamma r_h(v_0) e^(-(s(v_0) - a)(v_0 - b)) dv_0,

where `u = s(v_0)` is a branch of `S` and `r_h` is a finite expression in the Laurent coefficients of `h` at the
puncture. For such a functional `lambda(h V^m)` is a one-variable moment sequence along the non-properness curve; it
should fail to vanish for large `m` because of the branch points of `S` over the `v`-line. In the language of
D-modules: `B / T_u B` is a holonomic module on the `v`-line of generic rank `sum_i delta_i deg_u S_i >= d - 1`, and
`SN_u` says that all of its de Rham cohomology is carried by a quotient supported at `v = b`; the exponential
factors `e^(-(s(v) - a)(v - b))` should leave no room for that. Section 5 proves the conclusion at centres in
general position by another argument; the picture is only a guess at what the `2d - 2` further functionals of a
counterexample are (one-variable exponential periods on `S`, the family of THM-3031), and it may matter for
prescribed centres.

## 8. Remarks

- **Rank-one potentials, and a failure of the general image statement.** `x + x^2 y` and `x^2 y` have a
  one-dimensional `H^2(C^2, d - dW ^)` without being Gaussian. The first version of HYP-9668 asked whether the kernel
  of every such functional is Mathieu-Zhao. **That is false** (first audit): for `W = x^2 y`,

      1 = 2 (d/dy - W_y)(y) - (d/dx - W_x)(x),

  so 1 lies in `M = Im(d/dx - W_x) + Im(d/dy - W_y)`, while `M` is the kernel of the nonzero functional
  `x^k y^l -> l! [k = 2l + 1]` (value 1 on `x`); a Mathieu-Zhao subspace that contains 1 is the whole ring. The same
  happens for every `x^p y^q` with `p, q >= 1`, `p != q`: for weights `(alpha, beta)` with `p alpha + q beta = 1`,
  `(d/dx - W_x)(alpha x) + (d/dy - W_y)(beta y) = alpha + beta - W`, so `W = alpha + beta` modulo `M` along a whole
  line of weights, and `M` is proper because `H^2` has dimension `gcd(p, q)`. The two operators commute and have
  order one with constant leading coefficients, so this is also an instance where the statement of Zhao's
  Conjecture 1.3, read for an arbitrary potential, fails in characteristic zero. I do not claim this is new: the
  one search made (first audit) found no record of it, and van den Essen's survey (arXiv:1006.5801) reports that no
  counterexample in characteristic zero was known in 2010; it needs its own sourced note before anyone relies on
  it. The repaired question (HYP-9668 (E')) asks for the Mathieu property when `H^2` is one-dimensional **and the
  class of 1 is not zero** (`x + x^2 y`: the functional takes the value -1 on 1).
- **Other dimensions.** The definitions make sense for Keller maps of `C^n` with a nondegenerate quadratic form
  `Q`: `M = sum_i Im(d_(u_i) - Q_i(F))`, `dim S` is the reduced Betti number `b_(n-1)` of `F^(-1)(Q = t)`. But GMC(n)
  is false for `n >= 3`: the bottom rung exists only in the plane. For the map of THM-1300 the full
  three-dimensional Stein space was not computed.
- **Composition (ANALOGY).** For etale `G : X' -> X` and `F : X -> C^2`, `chi(Y_(F o G)) = deg(G) chi(Y_F) -`
  (sheets of `G` missing over `Y_F`): the pair (degree, Euler characteristic) is transported by an affine rule
  through a composition, in the same way as the state (offset, slope) of THM-5130 through a closer. Nothing is
  proved with this.

## 9. What this does not do

It does not prove JC(2), and it excludes no Keller pair: Corollaries 5.5 and 5.7 are reformulations, and at centres
in general position the Gaussian rungs say no more than the classical criterion "the dual field is locally
nilpotent". There is no statement here that lies strictly between GMC(2) and JC(2): every rung, asked of all Keller
pairs, is JC(2), and asked of the pair `(x, y)` it is GMC(2). Detection at a prescribed centre is not proved on
`C^2` and is false on surfaces. Theorem 3.1 rests on a cited comparison theorem that I did not re-prove. The
examples of section 6 live on surfaces other than `C^2`; they show what the count and the failure of the rungs look
like when sheets are missing, not that this must happen on `C^2`.

## 10. Reproduce

    cd scripts
    python stein.py                                   # controls
    python slice_pred.py                              # Jelonek quartic, slice map, intersection count
    python stein_runs.py toy 12 18 24 30              # S2          -> runs/s2_toy.txt
    python stein_runs.py slice 40 46                  # S3, S4, S5  -> runs/s3_s5_slice.txt   (about 90 s)
    python stein_more.py s6                           # S6          -> runs/s6_linear.txt
    python stein_more.py centres ; python stein_more.py auto ; python stein_more.py hermite 40 ; python stein_more.py identity
    python stein_centres.py                           # SN at six centres on the slice -> runs/sn_centres.txt
    python stein_inject.py                            # Theorem 5.4 sanity checks      -> runs/inject.txt
    python stein_special_centre.py                    # Proposition 6.3                -> runs/special_centre.txt
    python c3check.py                                 # THM-1300: det J = -2 and the three-point fibre

Audit programs: `audit/` (first audit), `audit/audit2/` (second audit).

# Stein modules of Keller pairs: Gaussian integration by parts between GMC(2) and JC(2) (kind-pasteur-2026-11-11-S24)

Owner directive (2026-11-11): merge prior work on the planar Jacobian conjecture with the proof of the
two-dimensional Gaussian moments conjecture, and define conjectures between GMC(2) and JC(2) that could form a
chain excluding high-degree counterexamples that rest on an undiscovered topological anomaly.

Canon: THM-5138, HYP-9668. Pre-registration: `PREREGISTRATION.md` (pushed before the runs, commits d701fad1cd and
90e0a6855d; all five predictions S2-S6 CONFIRMED). **PRE-AUDIT DRAFT.** JC(2) remains OPEN. Nothing here proves or
refutes a rung of the ladder that was not already known; what is new is the ladder itself, one criterion that is
equivalent to JC(2), a dimension count, and exact computations on the only planar non-proper etale maps in canon.

## 0. What comes out

GMC(2) and JC(2) live in the same ring `C[x, y]` and are statements about the same pair of operators.

For a Keller pair `(f, g)` (`{f, g} = 1`) take the dual fields `d_u = {., g}`, `d_v = {f, .}` (THM-4620), a centre
`(a, b)`, `U = f - a`, `V = g - b`, and the two *Stein operators* `d_u - V`, `d_v - U`: they are what Gaussian
integration by parts looks like in the coordinates `(f, g)`. The **Stein kernel** `M_F` is the sum of their images
and the **Stein space** is `S_F = C[x, y] / M_F`; its dual consists of the linear functionals that behave like the
Gaussian expectation with respect to `(f, g)`.

1. For the standard pair `(x, y)` the Stein space is one-dimensional, the functional is the Gaussian expectation,
   and **GMC(2) (THM-2022) is the statement that the Stein kernel is a Mathieu-Zhao subspace.** The same holds for
   every automorphism, by transport.
2. **Dimension count (Theorem 3.1).** For a Keller pair of generic degree `d` and a generic centre,
   `dim S_F = c + sum_i delta_i (deg_u S_i + deg_v S_i) >= 2d - 1`: one functional per connected component of the
   generic fibre of `UV`, plus one for every sheet that is missing over a point where a generic hyperbola meets the
   non-properness curve `S`. So **JC(2) says: the Stein functional of every Keller pair is unique.** A counterexample
   of generic degree `d` would carry at least `2d - 1` Gaussian-like functionals on the same polynomial ring.
3. **Finite chaos (Theorem 4.1).** A Keller pair is an automorphism iff every polynomial has a finite Hermite
   expansion in `(f, g)` modulo the Stein kernel: the classes of `d_u^m d_v^n h` vanish for `m + n` large. This is
   weaker in appearance than local nilpotence of the dual fields (the derivatives need only land in a subspace of
   finite codimension) and it is equivalent to JC(2).
4. **The ladder (section 5).** For all Keller pairs:
   `JC(2) <=> finite chaos <=> Stein uniqueness  =>  K-GMC (the Stein kernel is Mathieu-Zhao)  =>  cone vanishing
   =>  one-sided Stein nilpotence`, and GMC(2) is K-GMC for automorphisms. The three weaker rungs are new
   conjectures that lie between the theorem GMC(2) and the open problem JC(2); whether any of them implies JC(2)
   (the *detection* question) is open.
5. **What an anomaly would have to look like (sections 6, 7).** On the two planar non-proper etale maps available
   (a double cover on `C^* x C`, and a slice of the `C^3` Keller map of THM-1300 on `C^2` minus a hyperbola) the
   count is confirmed exactly (dimensions 2 and 11; linear potentials 2, 4, 6) and the lowest rung already fails:
   Gaussian moments see the missing sheets. The extra functionals of a counterexample are not Gaussian at all; they
   sit on the non-properness curve, and the obstruction they create does not depend on the degree of the map.

## 1. What exists

- **GMC(2)**: two real Gaussians, i.e. one pair `(Z, W = conj Z)` with `L(Z^a W^b) = a! [a = b]`. THM-2022 (PROVED;
  Lean endpoint `GMC2Main.gmc2`) through NC2 (a polynomial all of whose positive moments vanish has charges
  `a - b` of one strict sign) and one-variable DvdK (THM-2067, THM-2101, THM-2111). GMC(n) is false for `n >= 3`
  (THM-1480, THM-1490, THM-1500, THM-3290). No arrow from GMC(2) to JC(2) is in canon (MISTAKE-237; THM-1830 lists
  one as a superseded programme; THM-2801: GMC(2) reaches the special image conjecture in dimension one only).
- **Dual fields** of a Keller map, the Weyl-algebra module they define on `K[x]` (Bavula's twisted module), the
  kernel ladder, "cokernel of a dual field = relative de Rham cohomology", and "`(f, g)` is an automorphism iff
  `D_f` is surjective iff the generic fibre of `f` is a line" (KNOWN: Gurjar-Masuda-Miyanishi, Cerveau, Stein):
  THM-4620, THM-4626. The untwisted statements of this note (potential `W = f`) are there.
- **The `C^3` Keller map** with a cubic fibre: THM-1300; its Jelonek quartic and master identities: THM-1310,
  THM-1335; its coordinate fibrations: THM-4621; a plane family: THM-1345. **Non-properness sets**: THM-2633
  (inertia generates the monodromy; `1 <= k_D <= d - 1` finite branches), THM-4122 (components are polynomial
  curves).
- **Session267** (same owner prompt, 2026-11-11 12:02-12:20; CONCURRENT): consequences of NC2 for families
  (`L(Q P^m) = 0` once `m > deg Q`), the row repair operator of Session254 as Gaussian rotation, a rational-cone
  criterion for rigid debts, a shear-word compiler for inverse coordinates, a wreath law for the fixed-sheet gate,
  open bridges A, B, C. It has no Stein functional of a Keller pair, no twisted cohomology and no dimension count.
  The Student-Stein rows of THM-4315 / THM-4328 are Stein operators of a different object (source-normal rows).
- **CITED** (not in the repo before): A. Dimca, M. Saito, *On the cohomology of a general fiber of a polynomial
  map*, Compositio Math. 85 (1993): for a polynomial `W` on `C^n`, `H^k(Omega^., d - dW ^)` is the reduced
  cohomology of the generic fibre in degree `k - 1`. For a smooth affine variety `X`: `H^k(X, d - dW ^)` is
  `H^k(X, W^(-1)(t))` for generic `t` (C. Sabbah, Nagoya Math. J. 141 (1996); M. Hien, Invent. Math. 178 (2009),
  for the period pairing with rapid-decay homology). I read neither proof for this session; the statements are
  checked on ten examples below, including non-tame potentials and surfaces with `b_1, b_2 != 0`.

A grep of canon and notes for "twisted de Rham", "Dimca", "bifurcation set" returned nothing; "Stein identity"
returns only the Student-Stein rows.

## 2. The Stein module

Setting. `X` a smooth affine complex surface, `B = O(X)`, `F = (u, v) : X -> C^2` etale. For `X = C^2` this is a
Keller pair (scaled so that `{u, v} = 1`). The dual fields `d_u, d_v` are the derivations of `B` with
`d_u u = d_v v = 1`, `d_u v = d_v u = 0`; they commute. Centre `(a, b)`; `U = u - a`, `V = v - b`, `W = UV`.

    T_u := d_u - V,   T_v := d_v - U,   M_F := T_u B + T_v B,   S_F := B / M_F.

Classes in `S_F` are written `[h]`.

**Proposition 2.1 (twisted de Rham; PROVED).** For `omega = alpha du + beta dv`,
`(d - dW ^) omega = (T_u beta - T_v alpha) du ^ dv`. Hence `S_F = H^2(Omega^.(X), d - dW ^)`.

*Proof.* `d omega = (d_u beta - d_v alpha) du ^ dv` and `dW ^ omega = (V du + U dv) ^ omega = (V beta - U alpha) du ^ dv`. QED

**Proposition 2.2 (identities; PROVED).** For all `k, h` in `B` and all `i, j, m, n >= 0`:

- (a) `[V k] = [d_u k]`, `[U k] = [d_v k]`;
- (b) `[U^i V^j] = i! [i = j] [1]`;
- (c) (Cameron-Martin) `[h V^m] = [d_u^m h]`, `[h U^m] = [d_v^m h]`;
- (d) (Hermite) `[d_u^m d_v^n h] = [h H_(n,m)(U, V)]`,
  `H_(n,m)(U, V) = sum_j (-1)^j j! C(n,j) C(m,j) U^(n-j) V^(m-j)`;
- (e) (rotation) `(U d_u - V d_v) B` lies in `M_F`.

*Proof.* (a) is the definition. (b): `[U^i V^j] = [d_v (U^(i-1) V^j)] = j [U^(i-1) V^(j-1)]`, and `[V^j] = [d_u V^(j-1)] = 0`,
`[U^i] = 0` for `i, j >= 1`. (c): `[h V^m] = [d_u (h V^(m-1))] = [(d_u h) V^(m-1)]`. (d): put
`Psi_h(s, t) = [h e^(sV + tU)]` (formal). By (a), `d/ds Psi_h = Psi_(d_u h) + t Psi_h` and
`d/dt Psi_h = Psi_(d_v h) + s Psi_h`, so `Phi_h = e^(-st) Psi_h` satisfies `d/ds Phi_h = Phi_(d_u h)`,
`d/dt Phi_h = Phi_(d_v h)`, `Phi_h(0, 0) = [h]`, hence `Phi_h = sum [d_u^m d_v^n h] s^m t^n / (m! n!)`; and
`e^(sV + tU - st)` is the generating function of the `H_(n,m)`. (e): `[U d_u h] = [d_v d_u h] = [d_u d_v h] = [V d_v h]`. QED

Check: (b) for `i + j <= 4` and (d) for `m + n <= 3` on the slice of section 6 (15 and 21 cases, 0 failures).

**Proposition 2.3 (automorphisms; PROVED from THM-2022).** If `F` is an isomorphism onto `C^2`, then
`S_F = C [1]`, the Stein functional is `U^i V^j -> i! [i = j]`, i.e. `h -> L(h o F^(-1))` centred at `(a, b)`, and
`M_F` is a Mathieu-Zhao subspace of `B`.

*Proof.* `B = C[U, V]`; by 2.2 (b) `S_F` is spanned by `[1]`; the functional `U^i V^j -> i! [i = j]` kills
`T_u (U^i V^j) = i U^(i-1) V^j - U^i V^(j+1)` and likewise `T_v`, and is not zero. So `M_F` is its kernel, and
GMC(2) is the statement that this kernel is a Mathieu-Zhao subspace. QED

## 3. The dimension count

**Theorem 3.1 (PROVED from the cited comparison theorem).** Let `F = (f, g)` be a Keller pair on `C^2` of generic
degree `d = [C(x,y) : C(f,g)]`, with non-properness set `S` (empty or a curve with components `S_i`), and let
`delta_i = d - #F^(-1)(p)` for a generic point `p` of `S_i`. For a centre `(a, b)` such that the closure of `S` in
`P^1 x P^1` contains neither `(a, infinity)` nor `(infinity, b)` (all centres outside finitely many lines):

    dim S_F = b_1(Y) = c + sum_i delta_i (deg_u S_i + deg_v S_i),

where `Y` is a generic fibre of `W = (f - a)(g - b)`, `c` the number of its connected components, and
`deg_u S_i`, `deg_v S_i` the degrees in `u` and `v` of an equation of `S_i`. Moreover

    sum_i delta_i deg_v S_i = d - 1 + b_1(generic fibre of f),   sum_i delta_i deg_u S_i = d - 1 + b_1(generic fibre of g),

so `dim S_F >= 2d - 1`, and `dim S_F = 1` iff `F` is an automorphism.

*Proof.* By Dimca-Saito, `S_F = H^2(C^2, d - dW ^)` is `H^1(Y)`. `Y = F^(-1)(H_t)`, `H_t = {(u-a)(v-b) = t}`, a
copy of `C^*`. Over `H_t \ S` the map `F : Y -> H_t` is finite etale of degree `d`; over a point `p` of `H_t cap S`
it has `d - delta_p` preimages. So `chi(Y) = d chi(H_t \ S) + sum_p (d - delta_p) = - sum_p delta_p`, and
`b_1(Y) = c - chi(Y)`. In `P^1 x P^1` the curve `H_t` has bidegree (1, 1) and its only points at infinity are
`(a, infinity)`, `(infinity, b)`, the base points of the pencil; for the centres allowed and generic `t` it meets
each `S_i` transversally in `deg_u S_i + deg_v S_i` generic points. The same count for a vertical line `u = u_0`
gives `chi(f^(-1)(u_0)) = d - sum_i delta_i deg_v S_i`; the generic fibre of `f` is connected (THM-4620 (5): the
ring of constants of `D_f` is `C[f]`), so `b_1 = 1 - chi >= 0`. If `d = 1` then `F` is birational and etale, hence
an automorphism, `S` is empty and `Y` is `C^*`. If `d >= 2` then `dim S_F >= 1 + 2(d - 1) >= 3`. QED

**Corollary 3.2.** JC(2) holds iff every Keller pair has a one-dimensional Stein space at one such centre. A
counterexample of generic degree `d` has `dim S_F >= 2d - 1`; with the cited bound `d >= 5` of
`CORE-PAPERS-PLANAR-JC.md` (Orevkov, Domrina) this is at least 9.

**Remark 3.3 (surfaces).** For a smooth affine surface `X` the cited statement gives
`dim S_F = b_1(Y) - rank(H^1(X) -> H^1(Y)) + b_2(X)`, with the same Euler characteristic count for `b_1(Y)`.

**Remark 3.4 (degree).** The count involves only `c`, the `delta_i` and the bidegrees of the `S_i`. It does not see
the degrees of `f` and `g`: replacing `(f, g)` by `(f + g^N, g)` changes the potential, not the cover.

## 4. Finite chaos

Say that `F` has **finite chaos** at the centre `(a, b)` if for every `h` in `B` there is `N` with
`[d_u^m d_v^n h] = 0` whenever `m + n >= N` (by 2.2 (d): the Hermite expansion of `h` in `(U, V)` is finite modulo
the Stein kernel).

**Theorem 4.1 (PROVED; uses Theorem 3.1).** Let `F : X -> C^2` be etale of generic degree `d` with `S_F` finite
dimensional. If `F` has finite chaos at `(a, b)` then `dim S_F <= d`. Consequently a Keller pair on `C^2` has finite
chaos at one centre allowed in Theorem 3.1 iff it is an automorphism.

*Proof.* Define `Theta : B -> S_F (x) C[s, t]`, `Theta(h) = sum [d_u^m d_v^n h] s^m t^n / (m! n!)`. Then
`Theta(d_u h) = d/ds Theta(h)`, `Theta(d_v h) = d/dt Theta(h)`, and by 2.2 (a) and the Leibniz rule
`[d_u^m d_v^n (V h)] = [d_u^(m+1) d_v^n h] + n [d_u^m d_v^(n-1) h]`, so

    Theta(V h) = (t + d/ds) Theta(h),      Theta(U h) = (s + d/dt) Theta(h).

The operators `s + d/dt`, `t + d/ds`, `d/ds`, `d/dt` on `C[s, t]` satisfy the relations of `U, V, d_u, d_v`; call
this module `P`. It is generated by 1, and `(s + d/dt)^i (t + d/ds)^j 1 = s^i t^j +` lower terms, so `P` is free
of rank one over `C[U, V]`.

*Submodules.* Let `N` be a subspace of `S_F (x) P` stable under the four operators and `S' = N cap (S_F (x) 1)`.
Let `n` in `N` and let `sigma s^I t^J` be a monomial of `n` of largest total degree. No other monomial of `n` is
divisible by `s^I t^J`, so `d^I/ds^I d^J/dt^J n = I! J! sigma` and `sigma` is in `S'`. The element
`(s + d/dt)^I (t + d/ds)^J sigma = sigma s^I t^J +` terms of smaller degree lies in `N`; subtracting it from `n`
lowers the number of monomials of top degree. By induction, `n` is in `S' (x) P`; so `N = S' (x) P`.

Apply this to `N = Theta(B)`: the constant term of `Theta(h)` is `[h]`, and constant terms of elements of
`S' (x) P` lie in `S'`, so `S' = S_F` and `Theta` is onto. Hence `B` has a quotient that is free of rank
`dim S_F` over `A = C[U, V]`. Tensoring with the fraction field of `A`, a space of dimension `d` maps onto one of
dimension `dim S_F`. For a Keller pair and an allowed centre, `dim S_F >= 2d - 1 > d` unless `d = 1`. Conversely
for an automorphism `B = C[U, V]` and `d_u, d_v` are locally nilpotent. QED

**Corollary 4.2.** JC(2) <=> every Keller pair has finite chaos <=> every Keller pair has a unique Stein functional.

What the criterion adds to "the dual fields are locally nilpotent": the high derivatives are only required to lie in
`M_F`, which has finite codimension and contains every `U^i V^j` with `i != j`, every `U^i V^i - i!`, and the image
of the rotation field.

## 5. The ladder

For a Keller pair `F` and a centre:

- **(K-GMC)** `M_F` is a Mathieu-Zhao subspace of `C[x, y]`: if `[P^m] = 0` for all `m >= 1` then for every `h`,
  `[h P^m] = 0` for `m` large.
- **(CV), cone vanishing**: for every `h` there is `N` with `[h U^k V^(m+k)] = 0` and `[h U^(m+k) V^k] = 0` for
  all `m >= N`, `0 <= k <= m`.
- **(SN), one-sided Stein nilpotence**: for every `h`, `d_u^m h` and `d_v^m h` lie in `M_F` for `m` large.

**Proposition 5.1 (PROVED).** automorphism <=> finite chaos => K-GMC => CV => SN. For automorphisms K-GMC is
GMC(2).

*Proof.* The first equivalence is Theorem 4.1 and the first implication is Proposition 2.3. K-GMC => CV: for
`c` in `C` put `P_c = V (1 + c UV)`. By 2.2 (b), `[P_c^m] = sum_k C(m,k) c^k [U^k V^(m+k)] = 0` for `m >= 1`. So
for every `c` there is `m_0(c)` with `[h P_c^m] = 0` for `m >= m_0(c)`. The sets
`Z_N = {c : [h P_c^m] = 0 for all m >= N}` are Zariski closed, increase with `N` and cover `C`; a countable union of
finite sets does not, so some `Z_N` is `C`, and then every coefficient `[h U^k V^(m+k)]` of the polynomial
`c -> [h P_c^m]` vanishes for `m >= N`. The other half uses `U (1 + c UV)`. CV => SN is `k = 0` and 2.2 (c). QED

What K-GMC does not give in this way is the radial sector: the classes `[h (UV)^k]`, or the diagonal Hermite
classes `[h H_(n,n)(U, V)]`, for `k, n` large. Finite chaos asks for them; the Mathieu property is blind to `UV`,
which is not null (`[(UV)^k] = k! [1]`). This is the same split as in the proof of GMC(2), where the angular layer
is DvdK and the gap was radial (THM-1645).

**The conjectures between GMC(2) and JC(2)** (HYP-9668): K-GMC, CV, SN for all Keller pairs. Each follows from
JC(2), and each is a theorem for automorphisms. **Detection** asks whether SN (the weakest) already forces an
automorphism. If it does, all four statements are equivalent to JC(2), and JC(2) would follow from any proof of the
Mathieu property that works in the coordinates of an arbitrary Keller pair.

## 6. The planar non-proper etale maps in canon (FINITE-EXACT)

`scripts/stein.py` computes `O(X) / M` by exact linear algebra modulo the prime 2147483629, in the frame `dx ^ dy`
(`M' = Im(d/dx - W_x) + Im(d/dy - W_y)`; the class of `h` in `S_F` is the class of `h J`, `J` the Jacobian of
`(u, v)`), with elimination that puts pivots in the monomials of largest weight; `c(N, N')` is the codimension in
weight at most `N` of the span of the generators of weight at most `N'`.

**Controls** (dimension of `H^2(C^2, d - dW ^)`): `xy`: 1; `x^3 + y^3`: 4; `x + x^2 y` (no critical point; fibre
`C^*`): 1; `x^2 y` (critical line; fibre `C^*`): 1; `fg` for two automorphisms: 1. All as the cited theorem says.

**6.1 A double cover on `C^* x C`.** `(u, v) = (x^2, y/x)`, Jacobian 2; etale of degree 2 onto `C^* x C`, not
proper over `C^2` (the line `u = 0` is missed). Here `d_v = x d/dy` is locally nilpotent and
`d_u = (1/2x) d/dx + (y/2x^2) d/dy` is not.

**Proposition 6.1 (PROVED by hand; the dimension is FINITE-EXACT and agrees with Remark 3.3).** For
`W = (x^2 - a)(y/x - b)`, `a != 0`, the Stein space is two-dimensional, with basis of functionals

    lambda_+-(x^k y^l) = l! Res_(x = +-sqrt(a)) [ x^(k+l+1) e^(b (x^2 - a)) / (x^2 - a)^(l+1) ],

and `lambda_+-(x^(-2) V^m) = (-1)^m m! / (2 a^(m+1))` for all `m >= 0`, while `[V^m] = 0` for `m >= 1`. So `M_F` is
not a Mathieu-Zhao subspace and SN fails in the direction `d_u`; in the direction `d_v` it holds trivially.

*Proof.* On a generator `(d/dy - W_y)(x^k y^l)` the two terms cancel for `l >= 1` and the term for `l = 0` is regular
at `x = +-sqrt(a)`. On `(d/dx - W_x)(x^k y^l)` the expression under the residue is the `x`-derivative of
`l! x^(k+l+1) e^(b(x^2-a)) / (x^2-a)^(l+1)` (direct computation with `W_x = y + a y/x^2 - 2bx`,
`W_y = (x^2 - a)/x`). The two functionals are exchanged by `(x, y) -> (-x, -y)` and are independent. With
`s = x^2`: `lambda_+-(x^(-2-l) y^l) = (l!/2) Res_(s=a) [s^(-1) e^(b(s-a)) / (s-a)^(l+1)]`, and summing the
binomial expansion of `V^m = (y/x - b)^m` leaves `(1/2) (-1)^m m! / a^(m+1)`. QED

Machine (`runs/s2_toy.txt`, centre (2, 3)): `c = 2` at every truncation; classes of `x^(-2) V^m` and `x^(-1) V^m`
nonzero for `1 <= m <= 6`; classes of `x U^m`, `x^(-2) U^m` zero. Topology: fibre `P^1` minus four points,
`b_1 = 3`, `b_1(X) = 1`, `b_2(X) = 0`: `3 - 1 + 0 = 2`.

**6.2 A slice of the `C^3` Keller map.** For `F = (F_1, F_2, F_3)` of THM-1300, the surface `F_1 = c` (`c != 0`) is
`X_c = C^2 \ {1 + xy = 0}` and `(F_2, F_3)` restricts to

    u = y + 3c x w,   v = x w + x w^2 - c x^3 w^3,   w = 1/(1 + xy),   Jacobian -2 w^3,

an etale map of generic degree 3 onto its image in `C^2`. Its non-properness curve is the slice of the Jelonek
quartic, `S_c : 27 c^2 v^2 - 18 c u v + 16 c + u^3 v - u^2 = 0`, of bidegree (3, 2) in `(u, v)`, with two sheets
missing over a generic point (`scripts/slice_pred.py`).

**Proposition 6.2 (FINITE-EXACT; `c = 1`).**

| potential | points of `S_1` on the generic fibre image | predicted `b_1(Y) - 1 + 1` | computed |
|---|---|---|---|
| `(u - 2)(v - 3)` | 5 | `1 + 2*5 = 11` | 11 (stable for `N <= 36`, `N' = 40, 46`) |
| `u` | 2 | `1 - (3 - 4) = 2` | 2 |
| `v` | 3 | `1 - (3 - 6) = 4` | 4 |
| `u + v` | 4 | `1 - (3 - 8) = 6` | 6 |

Other centres: 11 at (5, -7), (1, 1), (-3, 4); 9 at (0, 0), where `b = 0` is the horizontal asymptote `v ~ 1/u` of
`S_1` and one intersection point is lost at infinity (as Theorem 3.1 predicts for an excluded centre).
At the centre (2, 3): `[U^i V^j] = i! [i = j] [1]` for `i + j <= 4`; `[1] != 0`; the classes of `w V^m`, `w U^m`
(`m <= 5`) and `x V^m`, `x U^m` (`m <= 4`) are all nonzero; every Hermite class `[h H_(n,m)(U, V)]` computed for
`h = w` and `h = x` (`n, m <= 4`, weight at most 30: 22 classes each) is nonzero, while for `h = 1` only
`(n, m) = (0, 0)` is. (`runs/s3_s5_slice.txt`, `runs/s6_linear.txt`, `runs/centres.txt`, `runs/hermite_slice.txt`.)

So on this surface the Stein space is `3 + 8`: three is the number of preimages of the centre, and eight further
functionals come with the missing sheets; SN, CV, K-GMC and finite chaos all fail, in both directions.

## 7. What an anomaly would be (HEURISTIC; not used in any proof)

For a Keller pair restrict to a generic horizontal line `v = v_0`. The curve `C = g^(-1)(v_0)` is etale over the
`u`-line with punctures `p` over the points of `S` on that line, `d_u` is `d/du` on `C`, and the fibrewise twisted
cohomology `H^1(C, d - (v_0 - b) du)` has dimension `sum_i delta_i deg_u S_i`: one residue functional
`h -> Res_p (h e^(-(u - a)(v_0 - b)) du)` for each puncture. Integrating these residues along a cycle `gamma` in the
`v_0`-line gives Stein functionals of the form

    lambda(h) = int_gamma r_h(v_0) e^(-(s(v_0) - a)(v_0 - b)) dv_0,

where `u = s(v_0)` is a branch of `S` and `r_h` is a finite expression in the Laurent coefficients of `h` at the
puncture. For such a functional `lambda(h V^m) = int_gamma (v_0 - b)^m r_h e^(...) dv_0`: a one-variable moment
sequence along the non-properness curve. It vanishes for all large `m` only if `(v_0 - b)^N r_h e^(...)` extends
holomorphically across `gamma`, which fails at the branch points of `S` over the `v`-line when `b` is generic. In
the language of D-modules: `B / T_u B` is a holonomic module on the `v`-line of generic rank
`sum_i delta_i deg_u S_i >= d - 1`, SN in the direction `d_u` says that all of its de Rham cohomology is carried by a
quotient supported at `v = b`, and the irregularity at infinity of its exponential factors
`e^(-(s(v) - a)(v - b))` leaves no room for that when `a` is generic.

This is the content of the detection conjecture (HYP-9668 (D)). It would say that the extra functionals of a
counterexample are one-variable exponential periods on `S` (the family of THM-3031 and the factorial-conjecture
lane), and that the obstruction is the same at every degree: nothing in it refers to `deg f` or `deg g`.

## 8. Remarks

- **Other dimensions.** The definitions make sense for Keller maps of `C^n` with a nondegenerate quadratic form
  `Q`: `M_F = sum_i Im(d_(u_i) - Q_i(F))`, `dim S_F` is the reduced Betti number `b_(n-1)` of `F^(-1)(Q = t)`, and
  finite chaos gives `dim S_F <= d` by the same proof. But GMC(n) is false for `n >= 3`: the bottom rung of the
  ladder exists only in the plane. For the map of THM-1300 the full three-dimensional Stein space was not computed.
- **Rank-one potentials.** `x + x^2 y` and `x^2 y` have a one-dimensional `H^2` without being Gaussian. Whether the
  kernel of such a functional is Mathieu-Zhao is a natural extension of GMC(2) (HYP-9668 (E)); for monomial
  potentials it is a weighted-charge statement of DvdK type.
- **Composition (ANALOGY).** For etale `G : X' -> X` and `F : X -> C^2`, `chi(Y_(F o G)) = deg(G) chi(Y_F) -`
  (sheets of `G` missing over `Y_F`): the pair (degree, Euler characteristic) is transported by an affine rule
  through a composition, in the same way as the state (offset, slope) of THM-5130 through a closer. Nothing is
  proved with this.

## 9. What this does not do

It does not prove JC(2), K-GMC, CV or SN for a single Keller pair that is not an automorphism (none is known), and
it does not prove detection. Theorem 3.1 rests on a cited comparison theorem that I did not re-prove. The examples
of section 6 live on surfaces other than `C^2`; they show what the count and the failure of the rungs look like when
sheets are missing, not that this must happen on `C^2`.

## 10. Reproduce

    cd scripts
    python stein.py                                   # controls
    python slice_pred.py                              # Jelonek quartic, slice map, intersection count
    python stein_runs.py toy 12 18 24 30              # S2          -> runs/s2_toy.txt
    python stein_runs.py slice 40 46                  # S3, S4, S5  -> runs/s3_s5_slice.txt   (about 90 s)
    python stein_more.py s6                           # S6          -> runs/s6_linear.txt
    python stein_more.py centres ; python stein_more.py auto ; python stein_more.py hermite 40 ; python stein_more.py identity
    python c3check.py                                 # THM-1300: det J = -2 and the three-point fibre

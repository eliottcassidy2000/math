---
id: THM-4700
title: "The unit-path ring of the six-colour field F = Q(zeta15, sqrt-11, w): clock_duality's B2 holds in F, so the component of 0 in the unit-distance graph is the intersection of the localizations at the conjugation-fixed primes, which lie exactly over the primes p = 14 mod 15 with (p/11) = -1 (8 or 24 of them per p); 1/n is unit-reachable iff n has no such prime factor (least excluded prime 29); short reciprocal walks, with the two-step reciprocals given by four Pell equations; the Heron criterion for triangles in a CM field; and the census of coefficient triangles of the cubic recipe (only equilateral and golden), so its b- and c-fibres admit no translation matching while single reflection links always exist"
status: >
  PROVED: statements 1-3 and 5. Statement 4 is PROVED for the six realizable triangles (explicit constructions) and
  FINITE-EXACT for the fourteen excluded ones: each has an explicit quadratic non-residue certificate at a degree-one
  place.
  CITED (classical, with proofs retained in the cited repo notes):
  - the Dedekind overring lemma (proved in section 2 of cyclotomic_unit_component_20261011.md, clock_duality_oct15);
  - integral bases of composita with coprime discriminants (classical; the cubic order O_(Q(w)) = Z[1, w, w^2/2] was
    computed in cubic_adelic_20261016e.md section 3).
  KNOWN: statement 1 is clock_duality_oct15's B2 (cubic_adelic_20261016e.md section 3) transported from
  F' = Q(zeta5, sqrt-11, w) to F = F'(sqrt-3). The proof needs only zeta3, so the content new here is:
  - the fixed-prime census of F (correcting the brief-2 count "8 per p", MISTAKE-635);
  - the rational and reciprocal statements of 2;
  - the Heron criterion in this form, and its triangle census;
  - the consequences for the cubic recipe in 5.
  The exact identities were checked in the degree-48 field. Self-audited only (PRE-AUDIT).
source: mac-mini-2026-10-18-unitpath
depends_on: [clock_duality_oct15 B1-B3 (cubic_adelic_20261016e.md), cyclotomic_unit_component_20261011.md section 2 (overring lemma), golden_contact_circuit_20261017i.md section 4 (cubic recipe u^2 = a + b u + c Y), THM-4690 statement 5 (19-adic separation), THM-4558 (Lemma R, residue colourings), merge_18_19_11_22_collatz_hn_20261018.md section 2.1 (fixed places of F)]
verification:
  - 04-computation/experiments/unit_path_ring_20261018/unit_path_ring_six_colour.py (+ .out): generators of O_F in Gamma, every walk of statement 2, the reflection link, recipe denominators (17 exact checks)
  - 04-computation/experiments/unit_path_ring_20261018/heron_census.py (+ .out, .json): the twenty side multisets of statement 4
  - 04-computation/experiments/unit_path_ring_20261018/matching_criterion.py (+ .out, .json): 1 - 4N for a, b, c; spindle pair denominators
  - 04-computation/experiments/unit_path_ring_20261018/two_step_reciprocals.out: the n <= 3000 with L(1/n) = 2
  - 04-computation/experiments/unit_path_ring_20261018/explicit_inverses.py (+ .out): the first explicit walks (session start)
  - 04-computation/experiments/merge_18_19_11_22_20261018/fixed_places_golden_field.py (+ .out): the fixed places for p <= 2099 (PARI-confirmed to 1300 by audit F3)
---

# THM-4700 — the unit-path ring of the six-colour field

## Setting

- `j = sqrt(-11)` and `f(w) = w^3 + 2w^2 - 4w - 4`.
- `F = Q(zeta15, j, w)` has degree 48. clock_duality's cubic recipe keeps its coordinate factors here
  (golden_contact_circuit_20261017i.md section 4).
- **F is a CM field.** Complex conjugation `c` sends `zeta -> zeta^-1`, `j -> -j` and `w -> w`. Its fixed field `F+`
  is totally real of degree 24, and `F = F+(sqrt-3)`. Since `c` commutes with every complex embedding, `z conj(z) = 1`
  says `|z| = 1` in every embedding at once.
- `U = {z in F : z conj(z) = 1}` is the set of unit vectors.
- `Gamma = sum Z U` is the set of points joined to 0 by a finite unit walk in the unit-distance graph `G(F)`.
  It is a ring, because products of unit vectors are unit vectors.
- `L(x)` is the least number of unit vectors with sum `x`.
- `N(x) = x conj(x)`.
- `S` is the set of primes `p` with `p = 14 (mod 15)` and `(p/11) = -1`:
  `S = {29, 149, 239, 359, 479, 569, 659, 809, 1019, ...}`.
- The cubic circuit is `u^2 = a + b u + c Y`, with `a = conj(1 + zeta5)` and `b`, `c` the coefficients of denominator 19.
  The factors on axes 0, 1, 4 are `A + aM`, `C + bM`, `C + cM`, where `M` is the Moser spindle
  `{0, 1, eta, 1 + eta, rho, rho eta, rho(1 + eta)}`, `eta = 1 + zeta3` and `rho = (5 + j)/6`.
- Write `d_k = 1/|k|` for `k = a, b, c`. Two coefficients `m, m'` interact:
  - through the u-square family when `|m - m'| = 1`;
  - through a translation matching of the axis-`k` fibres when `|m - m'| = d_k`.

## Statement

**1. The unit-path ring (B2 in F).** `O_F ⊆ Gamma` and

    Gamma = ⋂_{c(P) = P} (O_F)_P.

The `c`-fixed primes lie exactly over the primes of `S`, and none of them is ramified. Over `p` in `S` there are:
- **8** fixed primes if `f` has at most one root mod `p`. The residue field is `F_(p^2)`, or `F_(p^6)` when `f` is
  irreducible mod `p`.
- **24** fixed primes if `f` splits completely mod `p`. The first such prime is `p = 2549`, with roots 548, 2146, 2402.

The components of `G(F)` are the cosets `x + Gamma`.

**2. Reciprocals.** `Gamma ∩ Q = Z[1/p : p not in S]`.
- So `1/n` is unit-reachable iff no prime factor of `n` lies in `S`. Every prime `p < 29` has `1/p` in `Gamma`, and
  `1/29` is not in `Gamma`.
- Explicit walks, all verified exactly, with `delta_g = conj(g)/g`:

| Point | Walk | Length |
|---|---|---|
| `1/2` | `v + conj v`, `v = (1 + sqrt-15)/4` | 2 |
| `1/5` | `v + conj v`, `v = (1 + 3j)/10` | 2 |
| `1/13` | `v + conj v`, `v = (1 + 15 sqrt-3)/26` | 2 |
| `1/3` | `delta + delta^-1 + 2`, `g = (1 + j)/2` | 4 |
| `1/4` | `delta + delta^-1 + 2`, `g = (1 + sqrt-15)/2` | 4 |
| `1/7` | `delta + delta^-1 + 2`, `g = (1 + 3 sqrt-3)/2` | 4 |
| `1/19` | `delta + delta^-1 + 2`, `g = (1 + 5 sqrt-3)/2` | 4 |
| `1/23` | `2 delta + 2 delta^-1 - 3`, `g = (9 + j)/2` | 7 |
| `1/17` | `18 delta + 18 delta^-1 - 28`, `g = (11 + sqrt-15)/2` | 64 |
| `1/11` | clock_duality's connector `11 sum beta_r - 38` | 82 |

- **Two-step reciprocals.** For `n >= 2`, `L(1/n) = 2` iff `2n = x` for a solution of `x^2 - 3k y^2 = 1` with
  `k ∈ {1, 5, 33}`. The class `k = 165` would also be allowed, but `x^2 - 495 y^2 = 1` has no solution with `x` even.
  The values are `n = 2, 5, 13, 122, 181, 1985, 2521, 7562, ...`, from `k = 5, 33, 1, 5, 1, 33, 1, 5`.

**3. The Heron criterion.** Let `R_1, R_2, R_3` in `F+` be totally positive, and suppose some `z_1` in `F` has
`N(z_1) = R_1`. Then there is `z_2` in `F` with `N(z_2) = R_2`, `N(z_1 - z_2) = R_3` and `z_2/z_1` not in `F+` iff

    H = 2(R_1R_2 + R_2R_3 + R_3R_1) - R_1^2 - R_2^2 - R_3^2   (= 16 · area^2)   lies in 3·(F+)^{*2}.

In particular, for `x` not in `{0} ∪ 2U`: `L(x) = 2` iff `N(x)(4 - N(x))/3` lies in `(F+)^{*2}`.

**4. Coefficient triangles of the cubic recipe.** Take the twenty multisets of three sides from `{1, d_a, d_b, d_c}`.
Exactly six are triangles with vertices in `F`:
- the four equilateral ones;
- the two golden ones, `(1, 1, d_a)` and `(1, d_a, d_a)`.

The other fourteen are not realizable. This includes every triangle with a side `d_b` or `d_c` and two other sides
that are not both equal to it.

**5. Consequences for the cubic recipe.**
- **(a) One component.** Every point of `A + aM`, `C + bM` and `C + cM` lies in `Gamma`: the power-basis denominators
  involve only 2, 3 and 19. So the fourteen translated coefficient copies of `C` lie in the zero component of `G(F)`.
  THM-4690(5)'s 19-adic separation is therefore a statement about step alphabets that are integral at the `h = 7` primes
  above 19.
  - Those primes are `c`-fixed in `F' = Q(zeta5, j, w)`, where `1/19` is not in `Gamma_(F')` (clock_duality's B3).
  - In `F`, the factor `sqrt-3` splits each of them into a conjugate pair. The unit vector
    `delta_19 = (1 - 5 sqrt-3)/(1 + 5 sqrt-3)` has a pole at one prime of the pair and a zero at the other.
- **(b) No translation matching on the b- and c-axes.** For units `q ≠ q'`, the copies `C + bq` and `C + bq'` are never
  joined by the matching `x + bq ~ x + bq'`, since that would be a triangle `(1, 1, d_b)`. The same holds for `c`.
  For `a` the matching exists exactly when `q' = zeta10^(±1) q`.
- **(c) Reflection links.** For every unit difference `s` of `C`, every unit `q` and `k = b, c`, put

      q' = (conj(k)/k) · conj(q) · (kq - s)/(conj(kq) - conj(s)).

  Then `q'` is a unit and `|s + k(q' - q)| = 1`. So a single unit edge joining copy `q` to a new copy `q'` is always
  available.
- **(d) Connector length.** No walk of one or two unit steps joins two copies of `C + bM` or of `C + cM`. The one-step
  case is clock_duality's observation; the two-step case is THM-4690's connector census, with certified non-squares.
  The least connector length is therefore at least 3, and finite by (a).

## Proofs

**1.** *Generators.* The following identities are checked exactly (`unit_path_ring_six_colour.py`):
- `zeta15` is in `U`;
- `1/2 = v + conj(v)` with `v = (1 + sqrt-15)/4`;
- `1/11 = 11(beta_1 + ... + beta_4) - 38`, where `beta_r = (zeta5^r - 3)/(zeta5^-r - 3)` is in `U`; this is
  clock_duality's connector;
- `j = 6 rho - 5`;
- `w = (19(u + conj u) + 6)/11`, where `u = (w + 3 - j)/(w + 3 + j)` is in `U`; this is clock_duality's trace
  identity B1(3).

Since `Gamma` is a ring, it contains `(1 + j)/2`, `w` and `w^2/2`.

*Integral basis.* The discriminants of `Q(zeta15)`, `Q(j)` and `Q(w)` have pairwise disjoint supports `{3, 5}`, `{11}`
and `{2, 37}`. So `O_F = Z[zeta15] ⊗ Z[(1 + j)/2] ⊗ Z[1, w, w^2/2]`, and hence `O_F ⊆ Gamma`.

*Overring step.* The overring lemma gives `Gamma = ⋂ (O_F)_P` over the primes `P` with `Gamma ⊆ (O_F)_P`.
- If `c(P) = P`, every `z` in `U` has `2 v_P(z) = v_P(z conj z) = 0`.
- If `c(P) ≠ P`, pick `t` in `c(P)` but not in `P`. Then `t/conj(t)` is in `U` and has a pole at `P`.

This is clock_duality's argument verbatim.

*The fixed primes.* Put `K = Q(zeta15, j)`, with `Gal(K/Q) = (Z/15)^× × {±1}` and `c = (-1, -1)`.
- **Unramified `p`.** The decomposition group is generated by `Frob_p = (p mod 15, (-11/p))`, and `(-11/p) = (p/11)`.
  It contains `c` iff some odd power of `p` is `-1` mod 15 and `(p/11) = -1`. In `(Z/15)^×` only `p = 14` generates a
  subgroup containing `-1`.
- **Ramified primes 2, 3, 5, 11, 37.** None is fixed in `K`:
  - 2 and 37: `⟨2⟩ = {1, 2, 4, 8}` and `⟨7⟩ = {1, 4, 7, 13}` miss 14;
  - 3 and 5: they split in `Q(j)` into a pair swapped by `c`, while the decomposition group acts trivially on `j`;
  - 11: the decomposition group fixes `zeta5`.

  A fixed prime of `F` restricts to a fixed prime of `K`, so none of these is fixed in `F` either.

Now take `p` in `S`. Then `c = Frob_p`, all 8 primes of `K` above `p` are fixed, and `c` acts on each residue field
`F_(p^2)` as `Frob_p`. Since `f` is separable mod `p` (`p ∤ 148`), Hensel's lemma makes the primes of `F` above such a
`P` correspond to the irreducible factors of `f` over `F_(p^2)`. These are the orbits of `Frob_p^2` on the roots of `f`,
and `c` permutes them through `Frob_p`. So such a prime is `c`-fixed iff its orbit is `Frob_p`-stable:
- a root in `F_p` gives a fixed prime;
- the two roots of a quadratic that is irreducible over `F_p` are exchanged;
- an irreducible cubic stays irreducible over `F_(p^2)` and gives one fixed prime.

Hence there are `8 × (number of roots of f in F_p)` fixed primes, or 8 if `f` is irreducible mod `p`. A separable cubic
never has exactly two roots in `F_p`. The search over `p <= 20000` gives 2549 as the first prime of `S` with three roots.

**2.** A fixed prime `P` over `p` in `S` is unramified, so `v_P(1/p) = -1`; by 1, `1/p` is not in `Gamma`. If `p` is not
in `S`, no fixed prime lies over `p`, and `1/p` is in every `(O_F)_P` occurring in 1. This gives `Gamma ∩ Q`.
- The walks are exact identities.
- **The `delta` trick.** If `g` is in `F^*`, put `t = g + conj g` and `n = g conj g`. Then `delta = conj(g)/g` is in `U`
  and `delta + delta^-1 = t^2/n - 2`. With `t = 1` this gives `1/n = delta + delta^-1 + 2`. For `1/23` and `1/17`,
  Bezout combines `t^2/n` with integers.
- **Two-step reciprocals.** Statement 3 with `N = 1/n^2` asks for `(4n^2 - 1)/3` to lie in `(F+)^{*2}`. A positive
  rational is a square in `F+` iff its square class lies in `{1, 5, 33, 165}`:
  - the real quadratic subfields of `F+` are those of `K+`, because `[F+ : K+] = 3` is odd;
  - and those are `Q(sqrt5)`, `Q(sqrt33)`, `Q(sqrt165)`.

  If `3 | n`, the square class contains 3, so this cannot happen. Otherwise `4n^2 - 1 = 3k y^2`, which gives four
  Pell equations `x^2 - 3k y^2 = 1` with `x = 2n` even. Their fundamental solutions are `(2, 1)`, `(4, 1)`, `(10, 1)`
  and `(89, 4)`. In the first three, `x` is even exactly at the odd powers. For `k = 165`, the recurrence
  `x' = 89x + 1980y` keeps `x` odd, so there is no solution. The list `n <= 12000` was cross-checked by brute force
  (`two_step_reciprocals.out`).

**3.** Put `y = z_2/z_1`. Then `y + conj y = (R_1 + R_2 - R_3)/R_1 =: T` and `y conj y = R_2/R_1`, so `y` is a root of
`Y^2 - T Y + R_2/R_1`. Its discriminant is `T^2 - 4R_2/R_1 = -H/R_1^2`. The root `y` lies in `F` iff `-H` is a square
in `F`.

If `z^2` lies in `F+` with `z` in `F`, then `z` lies in `F+` or in `sqrt-3·F+`. This is because
`F = F+ ⊕ sqrt-3 F+`, and `(x + sqrt-3 y)^2` lies in `F+` iff `xy = 0`. So the following are equivalent:
- `y` lies in `F` but not in `F+`;
- `sqrt(-H)` lies in `sqrt-3·F+`;
- `H` lies in `3(F+)^{*2}`.

Conversely, such a `y` has norm `R_2/R_1`, and `z_2 = y z_1` works.

For the corollary, `L(x) = 2` iff `x = q_1 + q_2` with `q_1`, `q_2` in `U`. That is the triangle `(0, q_1, x)` with
squared sides `(1, 1, N(x))`, and then `H = N(4 - N)`.

**4.** *The six triangles.*
- Equilateral: `(z, zeta6 z)` for any side.
- `(1, 1, d_a)`: `(1, zeta10)`. Here `N(1 - zeta10) = 2 + zeta5^2 + zeta5^3 = 1/N(a)` exactly in `F`, because
  `(2 + zeta5 + zeta5^4)(2 + zeta5^2 + zeta5^3) = 1`.
- `(1, d_a, d_a)`: `(1/a, -zeta5/a)`.
- `(d_a, d_a, d_a)`: equilateral.

*The fourteen exclusions.* In each, `H/3` is a quadratic non-residue at a degree-one place of `F` where it is a unit,
so it is not a square even in `F`. The primes used are 1021, 1171, 1741 and 2011 (`heron_census.out`). Ten of the
fourteen also fail total positivity of `H` in some embedding, which is a second, independent obstruction.

**5.**
- (a) is the denominator computation in `unit_path_ring_six_colour.py`, together with 1 and 2. Pole and zero of
  `delta_19`: `(1 + 5 sqrt-3)/2` generates one of the two primes of `Q(sqrt-3)` above 19, and `c` exchanges them.
- (b) is 4.
- (c) **Reflection.** The two circles `|z| = |k|` and `|z - (kq - s)| = 1` meet at `kq` (since `|s| = 1`). The second
  meeting point is the reflection of `kq` in the line through 0 and `kq - s`. That point is
  `((kq - s)/conj(kq - s)) · conj(kq)`, which lies in `F`. Divide by `k`. The script checks it for `b`, `c` and every
  unit of `M`.
- (d) is cited.

## Remarks

**Walk-length geometry.**
- `L(x) <= 1` is a norm condition: `U` is the conic `P^1(F+)`, through `u_lambda = (lambda + sqrt-3)/(lambda - sqrt-3)`.
- `L(x) <= 2` is one square-class condition, by statement 3.
- `L(x) <= 3` asks for some `lambda` in `F+ ∪ {∞}` with `L(x - u_lambda) <= 2`. Clearing denominators, this is a point
  with `y ≠ 0` on the quartic

      3 y^2 = A_x(lambda) · (4 lambda^2 + 12 - A_x(lambda)),
      A_x(lambda) = N(x - 1) lambda^2 + 2 sqrt-3 (x - conj x) lambda + 3 N(x + 1),

  a curve of genus one when the quartic is squarefree.

So three-step connectivity between coefficient copies is a question about rational points on explicit genus-one curves
over the degree-24 field `F+`. This is why a finite direction screen (THM-4690's 450 named directions) cannot decide it.

**What does not follow.**
- Statements 1 and 5(a) say nothing about colourings. `chi(G(F)) = chi(G(Gamma))`, and the bound `chi(G(F)) <= κ(29)
  <= 7` with `κ(29) >= 6` is unchanged (brief-2 note, THM-4558).
- Connecting copies is cheap (5(c)), but rigidity is not. Translation matchings and mixed triangles are excluded on the
  b- and c-axes (5(b), statement 4).
- So any coefficient context that couples the fibre shifts must use cycles of length at least 4 that mix the distances
  `1` and `d_b` (or `d_c`).
- This sharpens the design target recorded in golden_contact_circuit_20261017i.md section 4 and in HYP-9289
  ("difference-forcing layers"). It does not refute either.

**Credits.**
- clock_duality_oct15: B1-B3, the overring proof, the `t/conj(t)` pole argument and the 82-step connector.
- The cubic recipe and the spindle `M` are clock_duality's.
- The six-colour field's fixed places are the brief-2 note's (corrected here).
- The Heron criterion is the CM-field form of the classical fact that lattice triangles of the Eisenstein plane have
  area in `sqrt3·Q`.

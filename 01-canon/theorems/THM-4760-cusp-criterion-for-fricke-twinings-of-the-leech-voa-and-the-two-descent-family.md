---
id: THM-4760
title: "Cusp criterion for Fricke twinings of the Leech VOA, and the 2-descent family. For v in Lambda (x) Q of order e with e v^2/2 in Z, Z_v - 24 = Theta_Lambda^[v]/eta^24 - 24 is the normalized Hauptmodul of Gamma_0(e)+e iff (0') Z_v has rational coefficients, (i) Lambda + v has a unique vector of norm 2 - 2/e, (ii) exactly 24 of norm 2, and (iii) at every intermediate cusp the twined shallow vectors of Lambda + cv cancel (Phi_c divides their residue polynomial). (0') is certified through q^1 by the vector-valued Sturm lemma. On the cusps adjacent to the 23 deep-hole cusps the criterion holds exactly for the 23 (h+1)-points (shell 0) and for seven shell-1 points, carried by A1 in A1^24, D4, D6, D8, E7 in D10E7^2, D12 and E8 in D16E8, of orders 4, 8, 12, 16, 20, 24, 32, which twine V_Lambda by 24 + T_4A, T_8A, T_12H, T_16C, T_20F, T_24I, T_32A; shell 2 gives none. Doubling realises Borcherds' 2-descent: for the nine Niemeier lattices with -1 in W the orbifold of 2v is the source N_1 (h = 2h_1 + 2) and 2v has a cancelling pair of minimal vectors, while for the other even h the point 2v is the (h/2+1)-point of a lattice of Coxeter number h/2. The other shell-1 points of the nine (carriers D10 in D10E7^2, D16 in D16E8, E8 in E8^3, D24) satisfy (i)-(iii) but are irrational (fields containing Q(sqrt 5), Q(sqrt 2), Q(sqrt 2), Q(sqrt 3)). The selection is Galois acting on Kac classes: for D_n at order 2n the shell-1 class is stable iff n in {2, 4, 6, 8, 12}; for E8 at order 32 the shell-1 class is stable and the principal one is not; for E7 at order 20 Galois swaps the two, which the swap of the two E7's of D10E7^2 undoes"
status: >
  PROVED (the criterion), plus FINITE-EXACT certificates and CITED class names.
  - The criterion (1) is PROVED: classical theta transformations of lattice cosets, the Weil representation of
    M = {x : x.v in Z}, and the principal-part argument on X_0(e) (the argument of THM-4759 (7), with rationality
    made explicit and shallow multiples allowed).
  - The rationality certificate (0) is PROVED: an isometry of the discriminant form commutes with the Weil
    representation, and the vector-valued Sturm lemma (weight 12) reduces invariance to vectors of norm <= 2.
  - The verdicts (3)-(6) are FINITE-EXACT: exact integer enumeration of all vectors of norm <= 2 in the e cosets
    Lambda + bv (fricke_cusp_certificate.py, fricke_cusp_certificate.out).
    - Shell 0 is complete.
    - Shell 1 is complete: its classes are the component orbits, every orbit occurs in the sample, and the counts
      are proportional to the component numbers.
    - Shell 2 is complete for shells with at most 200 vectors and sampled (200 vectors) otherwise.
    - Shells 2-8 of A15D9, A17E7, D10E7^2, D12^2, A24, D16E8, E8^3 and D24 (orders 18-39 and 47-55) are complete up
      to 150 vectors per shell (fricke_cusp_certificate_deep.out).
  - The Monster names in (4) are PROVED given the OEIS eta quotients for T_X (CITED; fricke_monster_classes.py):
    Ligozat's cusp orders and W_e-symmetry make each quotient the Hauptmodul of Gamma_0(e)+e, and uniqueness gives
    Z_v - 24 = T_X.
  - The doubling table (5) is VERIFIED: exact root enumeration of the Kneser neighbours (fricke_cusp_orbifolds.py).
  - The structural explanation of (5), V_Lambda^{orb(g^2)} = V_N^{orb(g_N^{e/2})} with g_N^{e/2} conjugate to the
    principal involution of N, is a SKETCH (standard cyclic-orbifold duality plus the vector-representation
    count below). It is not used in the proofs of (1)-(4).
  Prior art (third S53 agent search, 2026-10-22).
  - Lattice side of the shell-1 cusps: KNOWN. Chenevier-Lannes, Ergebnisse 69 (2019), Prop. III.4.3.4, shows that
    Lambda is a (2h + 1 - rho.alpha)-neighbour of every Niemeier lattice with roots. Taking alpha the highest root
    gives the (h+2)-neighbour, i.e. the orbifold of v is N, for all 23 types.
  - Borcherds' pairs h_2 = 2h_1 + 2: KNOWN (Borcherds' thesis Prop. 4.5.3; Chenevier-Lannes §III.4.4).
  - Method analogues of the criterion, CITED: Tuite, CMP 166 (1995) §4.1; Paquette-Persson-Volpato (2017) §6; van
    Ekeren-Moller-Scheithauer, arXiv 1704.00478, Lemmas 4.5-4.8.
  - The principal-part method is standard (Conway-Norton; THM-4759 (7)).
  - NOT FOUND: the criterion for lattice-VOA translations; any twining of V_Lambda by a translation equal to 24 +
    T_X for X = 4A, 8A, 12H, 16C, 20F, 24I, 32A; the cancellation mechanism; the level-2 / Galois reading.
  Session opus-2026-10-22-S53.
session: opus-2026-10-22-S53
source: 05-knowledge/results/leech_deep_hole_swap_theorem_20261022.md
scripts:
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_cusp_certificate.py (+ .out: shells 0-2, all 23 types; _deep.out: shells 0-8 of the large-h types; _small.out: shells 0-6 of the 14 types with h <= 14)
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_monster_classes.py (+ .out) (Ligozat certificates for 4A, 8A, 12H, 16C, 20F, 24I, 32A)
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_cusp_orbifolds.py (+ .out) (orbifolds of the shell-1 points and their multiples)
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_cusp_carriers.py (+ .out) (which component carries each shell-1 class)
  - 04-computation/experiments/deep_holes_theta_20261022/root_multisets.py (+ .out) (Galois action on the root spectra, (4'))
  - 04-computation/experiments/deep_holes_theta_20261022/kac_galois_classes.py (+ .out) (exact Galois action on the level-2 Kac classes P, C, (4'))
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_cusp_distance.py (+ _A2.out, _A3.out) (cusps at distance 2 and 3 from the deep-hole cusps)
  - 04-computation/experiments/deep_holes_theta_20261022/fricke_shell_test.py (+ .out), fricke_shell_all.py (+ .out) (q-series checks through q^2; q^4 for A1^24, D4^6)
related:
  - THM-4759 (deep-hole twining; (7), (10): this theorem supplies the rationality input of their proofs)
  - THM-4758 (swap theorem; the same vector-valued Sturm lemma)
  - HYP-9368 (genus-zero inner automorphisms of V_Lambda; part (B))
---

# THM-4760 — cusp criterion for Fricke twinings, and the 2-descent family

## Notation

- Λ is the Leech lattice. A point v ∈ Λ ⊗ Q has order e if ev ∈ Λ and no smaller multiple is.
- v has *type 0* if e·v²/2 ∈ Z. Then every norm of Λ + v lies in (2/e)Z.
- The twining is Z_v = Θ_Λ^{[v]}/η²⁴, where Θ_Λ^{[v]} = Σ_{λ∈Λ} e^{2πiλ·v} q^{λ²/2}. It is the graded trace of
  e^{2πiv(0)} on V_Λ, and its constant term is 24.
- For u, w ∈ Λ ⊗ Q put F(u, w) = Σ_{x∈Λ+u} e^{2πiw·x} q^{x²/2}. So Θ_Λ^{[v]} = F(0, v) and θ_{Λ+v} = F(v, 0).
- The *cusps adjacent to the deep-hole cusp* z of a Niemeier lattice N are the null vectors w = y + z′ + (y²/2)z,
  for y ∈ N, in II_{25,1} = N ⊕ U. They are again of type N, and −w·ρ = (h/2)|y − ρ_N/h|² (THM-4759 (11)).
- *Shell s* consists of the w with −w·ρ = h + 1 + s.
  - Shell 0 is the (h+1)-points (y = 0 and its glue translates).
  - Shell 1 is y = θ_j, the highest root of a component j, up to glue translation.
- For each w, v = c_w = proj_Λ(w)/m is the Leech point of order m = −w·ρ.

## Statement

**(1) The cusp criterion (PROVED).** Let v have order e and type 0. Then Z_v − 24 is the normalized Hauptmodul of
Γ_0(e) + e (in particular X_0(e)/W_e has genus 0) if and only if all four of the following hold:
- **(0′) Rationality.** Θ_Λ^{[v]} has rational coefficients. Equivalently, Θ^{[dv]} = Θ^{[v]} for every unit d mod e.
- **(i) Simple pole at 0.** Λ + v has a unique vector of norm < 2, and its norm is 2 − 2/e.
- **(ii) Fricke constant.** Λ + v has exactly 24 vectors of norm 2.
- **(iii) No other poles.** Take any c | e with 1 < c < e and any norm t < 2 occurring in Λ + cv.
  - On that norm, e(v·x) runs through one coset r_0 + (e/c)Z of Z/e.
  - Put P(T) = Σ_x T^{s(x)}, with s(x) = (e(v·x) − r_0)/(e/c) mod c.
  - Then Φ_c divides P.

**(2) The rationality certificate (PROVED).** Put M = {x ∈ Λ : x·v ∈ Z}. Then:
- M*/M ≅ (Z/e)² with q(a, b) = ab/e. Here Λ/M is the line b = 0, Λ + bv is the column b, and
  x ∈ Λ + bv has label a(x) = e(v·x) − a_0·b, where a_0 = e v²/2.
- For a unit d, the map α_d(a, b) = (da, b/d) is an isometry of M*/M, so it commutes with the Weil representation.
- **Sufficient condition for (0′), call it (0).** The theta vector of M is α_d-invariant for every unit d. By the
  vector-valued Sturm lemma (weight 12, a form vanishing to order > 1 at ∞ is 0) this holds iff, for every b and
  every t ≤ 2,

        #{x ∈ Λ + bv : x² = t, a(x) = a}  =  #{x ∈ Λ + (b/d)v : x² = t, a(x) = da}.

  Then θ_{(a,0)} = θ_{(da,0)}, and so Θ^{[dv]} = Θ^{[v]}.
- **Necessary condition for (0′).** If Θ^{[dv]} = Θ^{[v]}, then applying S gives θ_{Λ+dv} = θ_{Λ+v}. So a unit d
  with a different norm profile of Λ + dv below norm 2 *proves* that Z_v is irrational.

**(3) Shell 0 (FINITE-EXACT).** All 23 (h+1)-points satisfy (0), (i), (ii) and (iii). This gives a second proof of
the twining 24 + T_{(h+1)+(h+1)} of THM-4759 (10), and it supplies the rationality input that the proof given there
leaves implicit.

**(3′) Deeper shells and farther cusps (FINITE-EXACT, sampled where large).**
- No cusp in shell 2 of any type satisfies the criterion. Nor does any in shells 2-6 of the 14 lattices with h ≤ 14
  (fricke_cusp_certificate_small.out), or in shells 2-8 of the eight with h ≥ 16. For h ≤ 7 this includes every root
  cusp y = α, at height 2h + 1 − ht α (Chenevier-Lannes' neighbours).
  The orders reached include every missing genus-zero Fricke level between 18 and 50: 18, 21, 25, 27, 29, 35, 36, 39,
  49 and 50.
- Every cusp w has w·z = −A for some A ≥ 1, w = Az′ + y + (y²/2A)z, of height (h/2A)|y − Aρ_N/h|². At A = 2 (14 types
  with h ≥ 10) and A = 3 (D12², A24, D16E8, E8³, D24), the cusps with |y − Aρ/h|² ≤ 4.5 have heights 4–53, sampled.
  Every class that passes is one of the 30 of (3) and (4), recognised by its orbifold and series; no new order
  appears.
- The same Leech point is reached from several deep holes. For example, the D4⁶ shell-1 point (8A) and the A6⁴
  (h+1)-point (also 8A) both appear at A = 2 from D6⁴; they are inequivalent (2v shallow versus deep).

**(4) Shell 1 (FINITE-EXACT; class names CITED).** Exactly seven shell-1 classes satisfy the criterion, one in each
of the lattices below. The column "carrier" is the component j with y = θ_j (identified by the excess
(h/2)|y_j − ρ_j/h|² − (h/2)|ρ_j/h|² = 1).

| order e = h + 2 | N | carrier | Z_v − 24 = q⁻¹ + … | class |
|---|---|---|---|---|
| 4 | A1²⁴ | A1 | 276q + 2048q² + 11202q³ + 49152q⁴ | 4A |
| 8 | D4⁶ | D4 | 36q + 128q² + 386q³ + 1024q⁴ | 8A |
| 12 | D6⁴ | D6 | 14q + 36q² | 12H (A058486) |
| 16 | D8³ | D8 | 8q + 16q² | 16C (A058516) |
| 20 | D10E7² | E7 | 5q + 10q² | 20F (A058555) |
| 24 | D12² | D12 | 4q + 6q² | 24I (A058579) |
| 32 | D16E8 | E8 | 2q + 4q² | 32A (A058629) |

In all seven, the only shallow proper multiple is 2v. Its coset Λ + 2v has exactly two vectors of norm < 2, namely
±μ with μ² = 2 − 4/e, and P = 1 + T. The eigenvalues of g on the two lowest states of the g²-twisted sector are
opposite, and that cancels the pole at the cusp 1/2.

Every other shell-1 class fails, and each failure is a proof of non-Hauptmodularity:
- **Odd h** (A2¹², A4⁶, A6⁴, A8³, A12², A24): n_2 = 23, and the point is irrational.
- **Even h outside the nine lattices of (5)**: 2v has a single vector of norm < 2, so P = 1 and there is a pole at
  1/2.
- **The other carriers of the nine lattices of (5)**: D10 in D10E7², D16 in D16E8, the E8's of E8³, and D24. These
  pass (i), (ii) and (iii) but are IRRATIONAL. Only the units d ≡ ±1 mod 10, ±1 mod 8, ±1 mod 8 and ±1 mod 12
  respectively preserve θ(Λ + dv) below norm 2. So the coefficient fields contain Q(√5), Q(√2), Q(√2) and Q(√3).
  For E8³ at d = 3, Λ + 3v has minimum 15/8, below the Fricke norm 31/16.

**(4′) Why these seven: Galois permutes Kac classes (VERIFIED on root spectra).** On a component R of Coxeter number
h, the N-side element of a shell-1 point has one of two Kac classes of order h + 2:
- P = ρ/(h + 2), Kac type (3; 1, …, 1), on every non-carrier;
- C = (ρ + θ)/(h + 2), Kac type 2 at the node(s) next to α_0 and 1 at every other node, on the carrier.

The Galois power x ↦ dx must permute the multiset of component classes up to Aut(N). On the root eigenvalue
spectra (root_multisets.py, an independent necessary check), this works out as follows.

| component, order | P stable | C stable | behaviour |
|---|---|---|---|
| D_n, 2n | yes (the degrees 2, 4, …, 2n − 2, n are stable) | iff (n − 2)d ≡ ±(n − 2) or ±2 mod 2n for all units d, i.e. n ∈ {2, 4, 6, 8, 12} among the Niemeier D_n | n = 10, 16, 24 are unstable |
| E8, 32 | no | yes | |
| E7, 20 | — | — | d ≡ ±3 mod 10 swaps P and C |

**Exact Kac classes (kac_galois_classes.py).** Here dx is reduced into the fundamental alcove by W ⋉ Q in exact
arithmetic, and its Kac coordinates are compared with P and C. Conjugacy holds either exactly or up to an
extended-diagram symmetry; the latter is realised by glue.
- **E8, order 32.** C is fixed by every unit. P goes to a third class when d ≡ ±3 mod 8. So the E8-carrier point of
  D16E8 is conjugate to all its Galois powers inside W ⋉ Q, with no glue, and is rational for a structural reason.
  E8³ is not.
- **E7, order 20.** σ_3 swaps P and C exactly, and the swap of the two E7's of D10E7² undoes it. The units ±7 act up
  to the E7 diagram symmetry.
- **D_n, order 2n.** P is fixed up to a diagram symmetry for every n. C is fixed up to a diagram symmetry for n = 4,
  6, 8, 12, and goes to another class for n = 10, 16, 24.
- **A1, order 4.** σ_3 swaps P and C, and the all-ones Golay word (the glue (1/2)Σα_i) undoes it on all 24
  components at once.
- **E6, order 14.** The classes scatter.

The glue compatibility in the "up to a diagram symmetry" cases is not checked; rationality itself is proved by the
certificate (0).

This reproduces all nine cases:
- A1²⁴, D4⁶, D6⁴, D8³ and D12² are stable.
- D10E7² with an E7 carrier is stable up to the swap of its two E7's (its umbral Z/2).
- D16E8 with the E8 carrier is stable, being an E8 C next to a D16 P.
- The D10 carrier (an unstable C), the D16 carrier, E8³ (two unstable E8 P's) and D24 (an unstable C) are not.

**(5) Doubling realises Borcherds' 2-descent (VERIFIED, all 23 types).** For every shell-1 point v of a Niemeier
lattice N with even h, the orbifold of V_Λ at v is V_N, and 2v has order (h + 2)/2. Two cases occur.
- **The nine lattices with −1 ∈ W(N)**: A1²⁴, D4⁶, D6⁴, D8³, D10E7², D12², D16E8, E8³, D24. Here 2v has exactly two
  minimal vectors ±μ at its own Fricke norm 2 − 4/(h+2). Its orbifold is Borcherds' 2-descent source N_1, with
  h = 2h_1 + 2:

  | N | Λ-orbifold of 2v |
  |---|---|
  | A1²⁴ | Λ |
  | D4⁶ | A1²⁴ |
  | D6⁴ | A3⁸ |
  | D8³ | D4⁶ |
  | D10E7² | A7²D5² |
  | D12² | D6⁴ |
  | D16E8, E8³ | D8³ |
  | D24 | D12² |

- **The other eight even-h lattices**: A3⁸, A5⁴D4, A7²D5², A9²D6, A11D7E6, E6⁴, A15D9, A17E7. Here 2v has a
  *unique* minimal vector at norm 2 − 4/(h+2), and deep proper multiples, so it is a strict Fricke point. Its
  orbifold is the lattice of Coxeter number h/2: A1²⁴, A2¹², A3⁸, A4⁶, A5⁴D4, A5⁴D4, A7²D5², A8³ respectively. By
  HYP-9368 (A) it is that lattice's (h/2 + 1)-point. (A) is verified or sampled at orders 3, 4 and 6, sketched at
  the primes, and open at the orders 9 and 10 that occur here.

In both cases every multiple cv with c ≥ 3 is a deep hole of a lattice whose Coxeter number equals its order (an
A-type lattice: A1²⁴, A2¹², A3⁸, A4⁶, A5⁴D4, A7²D5², A8³, A11D7E6, A15D9).

**(6) The scale of the 2-descent (PROVED, elementary).** Put K = h + 2, the order of a shell-1 point. Borcherds'
relation h = 2h_1 + 2 is K(N) = 2K(N_1), and (5) is the statement that doubling halves the order. More generally,
the Fricke descent of THM-4759 (8), H(N) = p·H(X_p) + 1 with H = h + 1, becomes multiplication by p on the scale
L_p = (p − 1)H + 1. The genuine towers are:
- p = 2, L_2 = K: 2 → 4 → 8 → 16 → 32 (Λ → A1²⁴ → D4⁶ → D8³ → D16E8, E8³), 6 → 12 → 24 → 48 (A3⁸ → D6⁴ → D12²
  → D24), and 10 → 20 (A7²D5² → D10E7²);
- p = 3, L_3 = 2h + 3: 3 → 9 → 27 (Λ → A2¹² → E6⁴), and 21 → 63 (A8³ → E8³);
- p = 5, L_5 = 4h + 5: 5 → 25 → 125 (Λ → A4⁶ → E8³).

The shell-1 Fricke family of (4) is the rational part of the 2-towers, on the K-scale. Collatz's x ↦ 3x + 1 is
conjugate to x ↦ 3x by x ↦ 2x + 1 in the same way. That conjugacy, and the towers of powers of 3 and 5, are
NUMEROLOGY beyond the exact relation H = pH′ + 1. The E6⁴ entry L_3 = 27 = 3³ is the "27" of the S53 prompt.

## Proof

**(1), the theta calculus.** For x ∈ Λ + u, e^{πix²} = e^{2πix·u}e^{−πiu²}. Hence

    F(u, w)(τ + 1) = e^{−πiu²} F(u, u + w).

Poisson summation for the unimodular Λ of rank 24 gives

    F(u, w)(−1/τ) = τ¹² e^{2πiu·w} F(w, −u).

So the pair (u, w) transforms as a row vector, and for γ = (a b; c d) ∈ SL_2(Z)

    F(0, v)|_{12}γ = ε(γ) F(cv, dv),   |ε(γ)| = 1.

**(1), Γ_0(e)-invariance.** By type 0, the lattice M has level e, and det M = e² is a square. For γ ∈ Γ_0(e) the
Weil representation therefore acts on the cosets M + δ, δ ∈ Λ/M, by θ_{M+δ}|γ = θ_{M+aδ}. There is no character,
since (−1)^{12}e² is a square, and no phase, since q(δ) ∈ Z. Hence

    Θ^{[v]}|γ = Θ^{[dv]} = σ_d(Θ^{[v]}),

and under (0′) this is Θ^{[v]}.

**(1), poles at the cusps.** Every cusp of Γ_0(e) is equivalent to some a/c with c | e. There Z_v∘γ is a root of
unity times F(cv, dv)/η²⁴. It has a principal part iff some norm t < 2 of Λ + cv has a nonzero twined sum
Σ_x e^{2πid v·x}.

On a fixed norm, c(v·x) ≡ x²/2 + c²v²/2 mod 1 is constant. So e(v·x) lies in one coset r_0 + (e/c)Z, and the sum
is ζ_e^{dr_0} P(ζ_c^d). As a runs over the units mod c, so does d. Since P has integer coefficients, P(ζ_c^d) = 0
for all such d iff Φ_c | P. This is condition (iii).

**(1), the cusps ∞ and 0.**
- At c = e, F(0, v) = 1 + …, so there is a simple pole q⁻¹.
- At c = 1, Z_v(−1/τ) = θ_{Λ+v}(τ)/η(τ)²⁴. By type 0 the norms below 2 lie in {2 − 2/e, 2 − 4/e, …}. So (i) holds
  iff this has the single negative term q^{−1/e}, a simple pole in the width-e parameter, with coefficient 1.

**(1), "if".** Put G = Z_v∘W_e = θ_{Λ+v}(eτ)/η(eτ)²⁴, which is q⁻¹ + n_2 + O(q) by (i). G is Γ_0(e)-invariant,
because W_e normalises Γ_0(e). The difference G − Z_v is holomorphic at:
- ∞ and 0, where the leading coefficients 1 cancel;
- the other cusps, by (iii), since W_e permutes them.

So G − Z_v is constant, equal to n_2 − 24, which is 0 by (ii). Then Z_v is W_e-invariant with simple poles exactly
at the two cusps ∞ and 0, which W_e exchanges. So it has degree 1 on X_0(e)/W_e, and Z_v − 24 is the normalized
Hauptmodul.

**(1), "only if".**
- X_0(e)/W_e is defined over Q and has the rational cusp ∞. If it has genus 0 it is P¹_Q, so it has a uniformizer
  defined over Q. By the q-expansion principle that uniformizer has rational coefficients at ∞, and hence so does
  the normalized Hauptmodul. This gives (0′).
- A Hauptmodul has no poles at the intermediate cusps, which gives (iii).
- The simple pole at 0 with leading coefficient 1 gives (i).
- G = Z_v gives (ii).

**(2).**
- α_d preserves q(a, b) = ab/e, and O(M*/M) commutes with the Weil representation ρ. So α_d^*Θ_M lies in M_12(ρ).
- Vector-valued Sturm lemma. Every component of ρ(γ)F has order > 1 at ∞ whenever every component of F does. So
  each component F_i, a scalar form of weight 12 on Γ(N), vanishes at every cusp to order > N·(12/12) in the
  width-N parameter. The valence formula on Γ(N) then forces F_i = 0.
- The label a(x) = e(v·x) − a_0 b follows from v ≡ a_0 μ + ν mod M, where μ, ν generate the two isotropic lines with
  μ·ν = 1/e. Indeed q(v) = a_0/e, and B((a_0, 1), (a, b)) = (a + a_0 b)/e.
- The necessary condition is the S-transform F(0, w)(−1/τ) = τ¹²F(w, 0).

**(3)-(5).** These are exact enumerations, described under Scripts. For each cusp the certificate:
- builds v from w;
- enumerates every vector of norm ≤ 2 in the e − 1 cosets Λ + bv, b ≠ 0 (shell_vecs, exact integer norms);
- tests (0) for every unit d, together with the necessary condition, and then (i), (ii) and (iii), using PARI's
  polcyclo for Φ_c | P.

The doubling table computes the Kneser neighbour Λ_u + Zw of each multiple u = cv. Its roots are the norm-2 vectors
x ∈ Λ + ju with (x − jw)·u ∈ Z, split into irreducible components.

**(4), class names.** fricke_monster_classes.py takes the OEIS eta quotients E_X:
- 12H: (η(3τ)η(4τ)/(η(τ)η(12τ)))⁴;
- 16C: η(4τ)¹⁰/(η(τ)²η(2τ)³η(8τ)³η(16τ)²);
- 20F: (η(4τ)η(5τ)/(η(τ)η(20τ)))²;
- 24I: η(4τ)⁴η(6τ)⁴/(η(τ)η(2τ)²η(3τ)η(8τ)η(12τ)²η(24τ));
- 32A: (η(2τ)η(16τ))³/(η(τ)²η(4τ)η(8τ)η(32τ)²);
- 4A, 8A: the W_e-symmetric quotients (η(2τ)²/(η(τ)η(4τ)))²⁴ and (η(2τ)η(4τ)/(η(τ)η(8τ)))⁸.

For each, it checks:
- Ligozat's conditions (a Γ_0(e)-function);
- the cusp orders, which are −1 at ∞ and at 0 and ≥ 0 elsewhere;
- r_δ = r_{e/δ}, so that E_X∘W_e = E_X exactly.

So E_X − c_X is the Hauptmodul of Γ_0(e) + e, and it equals Z_v − 24 by uniqueness. Its coefficients also match
both the twinings and the OEIS terms.

**(5), why doubling lands on N_1 (SKETCH).**
- The dual of ⟨g⟩ ⊃ ⟨g²⟩ under cyclic-orbifold duality gives V_Λ^{orb(g²)} = V_N^{orb(g_N^{e/2})}.
- For the shell-1 class, g_N = e^{2πix(0)} with x = (ρ_N + θ_j)/(h + 2), so g_N^{e/2} = e^{πi(ρ_N + θ_j)(0)}.
- On a D_n component, ρ + θ = (n, n−1, n−3, …, 1, 0) has as many odd entries as ρ. So this involution is conjugate
  to the principal involution e^{πiρ(0)}, whose fixed roots are those of even height: Borcherds' X_2.
- The cancelling pair ±μ is then the lowest-weight space of the g²-twisted sector, on which g acts by ±i·const.

## Remarks

- **Two kinds of Fricke points.** The 23 (h+1)-points are strict: all proper multiples are deep. The seven
  2-descent points are not: their double is a "doubled" Fricke point with a ± pair. Both kinds now satisfy the
  single criterion (1).
- **Levels 0, 1, 2 (the THM-4759 (11) analogy extended; VERIFIED on root spectra).** Write the N-side element of a
  cusp as (ρ̂ + λ̂)/(h + k) for a level-k weight λ̂ of each component.
  - The deep hole is level 0 (λ̂ = 0, Kostant's principal element).
  - The (h+1)-point is the level-1 vacuum (λ̂ = Λ_0, Kac (2; 1, …, 1)).
  - A shell-1 point is level 2, with the adjoint θ̂ (the level-2 weight with finite part θ) on the carrier and the
    vacuum 2Λ_0 on every other component. This is the Kac type of (4′).

  Rationality is then Galois stability of this level-2 configuration under the Galois action on affine primaries,
  σ_d(λ̂) = w(d(λ̂ + ρ̂)) − ρ̂ (Coste-Gannon 1994; Fuchs-Gato-Rivera-Schellekens-Schweigert 1994; CITED as the
  framework). The computed instances (root_multisets.py) are:
  - D_n at level 2 always fixes the vacuum, and fixes the adjoint iff n ∈ {2, 4, 6, 8, 12};
  - E8 at level 2 is Ising-like: σ_3 swaps 2Λ_0 ↔ Λ_1 and fixes the adjoint Λ_8;
  - E7 at level 2: σ_3 swaps the vacuum class {2Λ_0, 2Λ_7} with the adjoint class {Λ_1, Λ_6} (classes up to the
    affine diagram symmetry).

  The level-1 vacuum is always Galois-fixed, which is why all 23 (h+1)-points are rational. The swap theorem holds at
  levels 0 and 1 and fails at level 2, where the Fricke property survives only for the Galois-stable configurations.
- **Two classes, one twining.** Several orders carry two inequivalent Fricke points with the same Hauptmodul:
  4 (A2¹² (h+1)-point and A1²⁴ shell-1), 7 (A5⁴D4, D4⁶), 8 (A6⁴ (h+1)-point and D4⁶ shell-1), 11, 13, 19 and 31.
  Equal twining does not determine the point; the orbifolds differ.
- **Prediction (registered in HYP-9368).** The realised Fricke orders are the 23 values 3–17, 19, 20, 23, 24, 26,
  31, 32 and 47. No Fricke Hauptmodul twining of V_Λ by a translation should exist at the other genus-zero Fricke
  levels: 18, 21, 25, 27, 29, 35, 36, 39, 41, 49, 50, 59, 71. Prime orders other than h + 1 are excluded by the
  dimension argument of HYP-9368 (A). The adjacent cusps of (3′) already exclude 18 through 50 on the shells tested.
- **Rationality is the real filter.** At shell 1, the poles cancel for all nine lattices with −1 ∈ W. What selects
  seven of the 13 classes is whether Galois permutes the component Kac classes, (4′). The answer is not "D-type":
  E7 and E8 carriers occur, and the D10, D16 and D24 carriers are irrational.

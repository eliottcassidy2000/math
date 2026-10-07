---
id: THM-4602
title: "Tournaments as chirotopes: for a tournament with skew sign matrix B, every 4-Pfaffian b_ij b_kl - b_ik b_jl + b_il b_jk is +-1 or +-3, and |Pf| = 3 exactly on 4-sets containing exactly one 3-cycle; hence T is the sign pattern of the Pluecker coordinates of n vectors in R^2 (a rank-2 chirotope, a point of Gr(2,n) up to signs) iff all 4-Pfaffians are +-1 iff T is locally transitive, and the transitive tournament is the totally positive part (where Conway-Coxeter friezes live); over F_p, p = 3 mod 4, reading 'positive' as 'quadratic residue', the Legendre pattern of all p+1 points of P^1(F_p) is the Paley tournament plus a sink up to switching, every SL2 frieze strip over F_p has as Legendre pattern a switched induced subtournament of it, and Singer frieze strips realise the whole class"
status: >
  PROVED (elementary; the 4-vertex statement is an exhaustive check of all 64 labelled 4-tournaments). FINITE-EXACT checks: 800 random
  tournaments and 800 planar vector configurations for n = 5..8; Singer frieze strips for 13 primes p <= 83 (isomorphism to Paley + sink found by
  backtracking); 200 random frieze strips. KNOWN ingredients: rank-2 oriented matroids are realizable, and locally transitive = circular tournaments
  (classical); the Paley skew-Hadamard construction from P^1(F_p) (Paley 1933); SL2-friezes over finite fields (Morier-Genoud 2021). The packaging
  as a Legendre frieze dictionary is new to the repo; no priority claim for the ingredients.
session: mac-mini-2026-10-07-twoanchor
source: 05-knowledge/results/twoanchor_reset2_friezes_20261007.md (section 8.3)
scripts:
  - 04-computation/experiments/twoanchor_20261007/tournament_chirotope.py (ALL CHECKS PASSED)
related:
  - THM-438 (Paley cluster integrals are Catalan: plane-tree Euler tours, the same Catalan numbers as Conway-Coxeter friezes)
  - THM-4552/4553 (PSL(2,7) > Borel > torus; the Paley tournament QR7)
  - 05-knowledge/results/frieze_collatz_reset2_20261007.md (positive marked minors of Collatz carries)
---

# THM-4602 — tournaments as chirotopes; Legendre frieze strips

## Setting

* A tournament T on [n] is encoded by the skew matrix `B` with `b_ij = +1` if i → j and `−1` otherwise (`b_ji = −b_ij`).
* For vectors `v_1, …, v_n` in a plane over a field with a "positivity" character χ, the χ-*pattern* orients i → j iff `χ(det(v_i, v_j)) = +1`. Over R, χ is the sign; over F_p, the Legendre symbol.
* This is a tournament whenever χ(−1) = −1, i.e. over R, or over F_p with p ≡ 3 mod 4.

## Statements

1. **Four-point criterion (PROVED).**
   * For every 4-set {i<j<k<l}, `Pf(B_ijkl) = b_ij b_kl − b_ik b_jl + b_il b_jk ∈ {±1, ±3}`.
   * `|Pf| = 3` iff the 4-set spans exactly one 3-cycle, i.e. a 3-cycle with a vertex dominating it or dominated by it.
2. **Rank two = local transitivity (PROVED).** The following are equivalent:
   * (a) T is a sign pattern `sign det(v_i, v_j)` of vectors in R^2: a realisable rank-2 chirotope, i.e. the sign vector of a point of Gr(2,n) with nonzero Plücker coordinates;
   * (b) every 4-Pfaffian is ±1, i.e. every 3-term Grassmann–Plücker relation has terms of both signs;
   * (c) T is locally transitive: every out- and in-neighbourhood is transitive.

   The totally positive part (all `p_ij > 0` for i < j) is exactly the transitive tournament in its order. This is the locus of positive SL2 frieze patterns (Conway–Coxeter) and of positive cluster coordinates of type A.
3. **Legendre chirotopes (PROVED).** Let p ≡ 3 mod 4.
   * (a) For representatives `(1, x)`, `x ∈ F_p`, and `(0, 1)`, the Legendre pattern of P^1(F_p) is the Paley tournament on F_p with ∞ a sink.
   * (b) Rescaling a representative by λ switches its vertex iff `χ(λ) = −1`. So the switching class is an invariant of the point set.
   * (c) Every SL2 frieze strip over F_p (distinct points with consecutive determinants 1) therefore has as Legendre pattern a switched induced subtournament of Paley + sink.
   * (d) The PGL2 Singer cycle orders all p+1 points. Normalised to consecutive determinant 1, it gives a frieze strip whose pattern is the full Paley + sink class (verified for p = 3, …, 83).
   * (e) Constant-quiddity (elliptic Chebyshev) friezes cover only (p+1)/2 points: the elliptic torus of SL2(F_p) has image of order (p+1)/2 in PSL2.

## Proof sketch

* (1) Exhaustive over the 64 labelled 4-tournaments, or directly: the three products are ±1, so `|Pf| = 3` iff they are equal, iff the GP relation is sign-violated. Checking the four isomorphism types (TT4, strong, 3-cycle + source, 3-cycle + sink) gives the stated characterisation.
* (2) (a) ⇒ (b) is the GP relation `p_ij p_kl − p_ik p_jl + p_il p_jk = 0`. (b) ⇒ (a) is the realisability of rank-2 uniform chirotopes (points on a circle in angular order, classical). (b) ⟺ (c) follows from (1): a 4-set with one 3-cycle is precisely a vertex whose out- or in-neighbourhood contains a 3-cycle.
* (3) `det((1,x),(1,y)) = y − x` and `det((1,x),(0,1)) = 1`; `det(λv, w) = λ det(v, w)`. ∎

## Reading

* **Over R.** Positivity of all frieze entries gives the transitive tournament. Rank two in general gives the locally transitive (circular) tournaments.
* **Over F_p.** "Positivity = quadratic residue" gives the Paley tournament, which our tournament work (THM-438, R(p) → e; QR7) treats as the extremal object. In this sense the Paley tournament is the finite-field totally positive tournament. A χ-positive frieze over F_p (all entries residues) is a transitive subtournament of Paley + sink, so its width is bounded by the largest transitive subtournament, of logarithmic size.
* **Catalan dictionary (KNOWN bijections).** THM-438's leading patterns are Euler tours of plane trees, counted by C_k. Via triangulations these are Conway–Coxeter friezes.

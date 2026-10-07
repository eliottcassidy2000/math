---
id: THM-4602
title: "Tournaments as chirotopes: every 4-Pfaffian b_ij b_kl - b_ik b_jl + b_il b_jk of a tournament's skew sign matrix is +-1 or +-3, with |Pf| = 3 exactly on 4-sets containing one 3-cycle; so T is the sign pattern of the Pluecker coordinates of n vectors in R^2 (a rank-2 chirotope) iff all 4-Pfaffians are +-1 iff T is locally transitive (KNOWN: Babai-Cameron 2000, Gunderson-Semeraro 2017); positive real SL2 friezes are the slice p_(i,i+1) = p_(1n) = 1 of the totally positive (transitive) part and Conway-Coxeter friezes are its integer points; over F_p (p = 3 mod 4) the Legendre pattern of P^1(F_p) is Paley + sink up to switching (KNOWN, Paley 1933), every SL2 frieze strip over F_p gives a switched induced subtournament of it, the PGL2 Singer ordering closes as a frieze of width p-2 with 2-periodic quiddity, constant quiddities cover at most (p+1)/2 (elliptic), (p-1)/2 (hyperbolic) or p (parabolic -2) points, and chi-positive configurations have at most tt(QR_p) + 1 points"
status: >
  PROVED (elementary). KNOWN as statements: (1) is the local-order / vortex / 4-graph criterion (Babai-Cameron 2000 section 3; Knuth 1992;
  Gunderson-Semeraro JCTB 126 (2017), Fact 18, Lemma 21); (2) is Babai-Cameron 2000 Lemma 3.3 (with Brouwer 1980, Lachlan 1984) in
  chirotope language; (3)(a)-(b) are Gunderson-Semeraro 2017 Def. 17 / Thm 19 and arXiv:2204.10775 Lemma 3.4, Prop. 3.5, going back to
  Paley 1933. New to the repo: the Pfaffian / frieze-strip phrasing, the slice description of positive friezes inside the transitive part,
  the closure of Singer strips (3)(d), the constant-quiddity census (3)(e), and the chi-positive bound tt(QR_p) + 1 (3)(f). FINITE-EXACT:
  64 labelled 4-tournaments; LT = rank-2 = Pfaffian-+-1 counts (n-1)! 2^(n-1) for n = 4..9 (audit B, exhaustive with realisation
  certificates); Singer strips for the 17 primes p = 3 mod 4 up to 131; tt(QR_p) + 1 attained for p <= 43. Independently audited 2026-10-07
  (audit B); corrections applied (slice, Singer closure, constant quiddities, Paley analogy, Catalan coincidence, KNOWN typing); MISTAKE-586.
session: mac-mini-2026-10-07-twoanchor
source: 05-knowledge/results/twoanchor_reset2_friezes_20261007.md (section 8.3)
scripts:
  - 04-computation/experiments/twoanchor_20261007/tournament_chirotope.py (ALL CHECKS PASSED)
  - 04-computation/experiments/twoanchor_20261007/audit_B/ (chirotope_audit.py, lt_enum.c, legendre_audit.py + outputs, REPORT.md)
related:
  - THM-438 (Paley cluster integrals have leading coefficient C_k, a signed Moebius / free-cumulant sum, not plane-tree tours (MISTAKE-060/061); the same numbers count Conway-Coxeter friezes)
  - THM-4552/4553 (the PSL(2,7) > Borel > torus ladder; by (3)(b) and Babai-Cameron Lemma 3.1 an odd-order group such as the Borel fixes a tournament in the Legendre switching class, while PSL(2,p) of even order fixes none)
  - 05-knowledge/results/frieze_collatz_reset2_20261007.md (positive marked minors of Collatz carries)
---

# THM-4602 — tournaments as chirotopes; Legendre frieze strips

## Setting

* A tournament T on [n] has the skew matrix `B`, with `b_ij = +1` if i → j and `−1` otherwise.
* For vectors `v_1, …, v_n` in a plane over a field with a "positivity" character χ, the *χ-pattern* orients i → j iff `χ(det(v_i, v_j)) = +1`. Over R, χ is the sign; over F_p it is the Legendre symbol.
* The χ-pattern is a tournament when χ(−1) = −1: over R, or over F_p with p ≡ 3 mod 4.

## Statements

1. **Four-point criterion.**
   * For every 4-set {i<j<k<l}, `Pf(B_ijkl) = b_ij b_kl − b_ik b_jl + b_il b_jk ∈ {±1, ±3}`.
   * `|Pf| = 3` iff the 4-set spans exactly one 3-cycle, i.e. a 3-cycle together with a vertex dominating it or dominated by it.

2. **Rank two = local transitivity.** The following are equivalent:
   * (a) T is a sign pattern `sign det(v_i, v_j)` of vectors in R^2, i.e. a realisable uniform rank-2 chirotope (the sign vector of a point of Gr(2,n) with nonzero Plücker coordinates);
   * (b) every 4-Pfaffian is ±1, i.e. every 3-term Grassmann–Plücker relation has terms of both signs;
   * (c) T is locally transitive.

   **Positivity.** The totally positive part (all `p_ij > 0` for i < j) consists exactly of the configurations whose pattern is the transitive tournament in its order.
   * Positive real SL2 friezes with n points (width n−3) are the points of this part normalised by `p_(i,i+1) = 1` and `p_(1n) = 1`. This is an (n−3)-dimensional slice, not the whole part.
   * For odd n the slice meets every torus orbit once, so positive friezes ≅ Gr⁺(2,n)/T.
   * For even n it meets only the orbits with `p_12 p_34 ⋯ p_(n−1,n) = p_23 ⋯ p_(n−2,n−1) p_(1n)`, each in a one-parameter family.
   * Conway–Coxeter friezes are the positive *integer* friezes, C_(n−2) of them.
   * The Plücker coordinates are the cluster variables of type A_(n−3).

3. **Legendre chirotopes, p ≡ 3 mod 4.**
   * (a) With representatives `(1, x)` (x ∈ F_p) and `(0, 1)`, the Legendre pattern of P^1(F_p) is the Paley tournament on F_p with ∞ a sink.
   * (b) Rescaling a representative by λ switches its vertex iff `χ(λ) = −1`. So the switching class is an invariant of the point set.
   * (c) Every SL2 frieze strip over F_p (distinct points, consecutive determinants 1) has as Legendre pattern a switched induced subtournament of Paley + sink. Any ordering of all p+1 points gives the whole class.
   * (d) **Singer closure.** Let g be a PGL2 Singer cycle: projective order p+1, with `det g` a non-residue. The ordering it induces closes: `v_(p+1) = χ(det g)·v_0 = −v_0`. So it is a closed SL2 frieze of width p−2 over F_p, with no zero entry and 2-periodic quiddity (a, b), `ab = tr(g)²/det(g)`. Random orderings close in about 1/(p−1) of cases.
   * (e) **Constant quiddities.** A constant quiddity covers as many points as the order of its matrix in PSL2(F_p):
     * elliptic: at most (p+1)/2;
     * hyperbolic: at most (p−1)/2;
     * parabolic (−2): a closed frieze through p points, all but its fixed point, whose pattern is the switching class of QR_p.

     No constant quiddity covers all p+1 points.
   * (f) **χ-positive configurations.**
     * A configuration with all Plücker coordinates residues has transitive Legendre pattern.
     * After an SL2(F_p) change of coordinates it is a transitive subtournament of Paley + sink, with ∞ as the sink.
     * So it has at most tt(QR_p) + 1 points, and this is attained (brute force, p ≤ 43).
     * The proved bound is tt(QR_p) ≤ 2√p + 1. Logarithmic growth is observed (tt(QR_199) = 11) but unproved.

## Proof sketch

* **(1)** The three products are ±1, so `|Pf| = 3` iff they are equal, iff the Grassmann–Plücker relation is sign-violated. Checking the four isomorphism types (TT4, strong, 3-cycle + source, 3-cycle + sink) gives the characterisation.
* **(2)**
  * (a) ⇒ (b) is the Grassmann–Plücker relation `p_ij p_kl − p_ik p_jl + p_il p_jk = 0`.
  * (b) ⇒ (a) is the realisability of rank-2 uniform chirotopes (angular order of ±v_i).
  * (b) ⟺ (c) follows from (1).
  * The slice statement follows from the frieze normalisation and the torus action `p_ij ↦ t_i t_j p_ij`.
* **(3)**
  * (a), (b), (c): `det((1,x),(1,y)) = y − x`, `det((1,x),(0,1)) = 1`, and `det(λv, w) = λ det(v, w)`.
  * (d): write `v_k = C_k g^k v_0`. Consecutive determinant 1 gives `C_(k+1)/C_(k−1) = 1/det g`. Together with `g^(p+1) = det(g)·I` this gives `v_(p+1) = χ(det g) v_0`.
  * (e): the orbit of a point under the quiddity matrix in PSL2.
  * (f): send the last vector to (0, 1) by SL2(F_p), then apply (a). ∎

## Reading

* **Over R.** The standard chart (1, x), (0, 1) has the transitive pattern: the real line in its order, with ∞ a sink. Rank-2 configurations in general give the locally transitive (circular) tournaments.
* **Over F_p.** The same chart has the Paley + sink pattern. In this sense (ANALOGY) the Paley tournament, the extremal object of our tournament work (THM-438, QR7), is the finite-field counterpart of the order of the real line. The counterpart of total positivity is different: χ-positive configurations are small transitive subtournaments, as in (3)(f).
* **Catalan coincidence.** THM-438's leading coefficient is C_k, but as a signed Möbius sum over even-series patterns, not a count of plane-tree tours (MISTAKE-060/061). Conway–Coxeter friezes with k+2 points are also counted by C_k. No bijection is claimed.
* **Mutation.** Quiver mutation at a source or sink reverses the arrows at that vertex. That is switching, which (3)(b) shows is the natural symmetry of Legendre chirotopes.

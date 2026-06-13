# THM-473 — The Involution Average: E[det(I+S)] = involutions, E[det(xI+S)] = Hermite

**Status:** PROVED and ADVERSARIALLY VERIFIED (2026-06-11; exhaustive recheck through
n=7 labeled, every identity confirmed) — but READ THE ATTRIBUTION SECTION: parts 1
and 3 are previously published (KMPRS 2025); this file is a corollary-with-
reinterpretation. The new content is part 2 (the involution/SYT evaluation at x=1)
and the Godsil–Gutman genealogy. See 05-knowledge/results/thm473_adversarial_check.out.
**Provenance:** mac-mini-2026-06-10-S2 (renumbered THM-468(avg)->THM-473 per first-come collision; see MSG-099). Companions: THM-468 (setup), THM-472 (ceiling).
Related: HYP-2383/2387, T777, everything-is-the-triangle (Young tableaux mandate).

## Statement

Over the uniform distribution on all 2^C(n,2) labeled tournaments on n vertices
(equivalently: independent fair ±1 coins on the C(n,2) entries of the skew matrix S
above the diagonal):

1. **E[det(xI + S)] = Σ_{k=0}^{⌊n/2⌋} C(n,2k) (2k-1)!! x^(n-2k) = i^(-n) He_n(ix)**,
   where He_n is the probabilist's Hermite polynomial. Equivalently: the expected
   characteristic polynomial of the random tournament skew matrix is the signless
   matching polynomial of K_n.
2. **E[det(I + S)] = I(n)**, the number of involutions on [n] (= the total number of
   standard Young tableaux with n cells): 1, 2, 4, 10, 26, 76, 232, 764, 2620, … (A000085).
3. For even n, **E[det S] = (n-1)!!** = the number of perfect matchings of K_n.

## Proof

det(xI + S) = Σ_K x^(n-|K|) det(S[K]) over principal minors. Odd |K| minors vanish
(skew); even minors are Pf(S[K])². Pf(S[K]) = Σ_μ sgn(μ) Π_{e∈μ} s_e over perfect
matchings μ of K. Hence

  E[Pf(S[K])²] = Σ_{μ,μ'} sgn(μ)sgn(μ') E[Π_{e∈μ} s_e Π_{e∈μ'} s_e].

The entries s_e (e an unordered pair, s_{ji} = -s_{ij}) are independent uniform ±1,
so the expectation of a monomial is 1 if every edge appears an even number of times
(i.e. μ = μ', as a matching uses each edge at most once) and 0 otherwise:
E[Pf²] = #matchings(K) = (|K|-1)!!. Therefore

  E[det(xI+S)] = Σ_{k} C(n,2k) (2k-1)!! x^(n-2k).

He_n(ix) = i^n Σ_k C(n,2k)(2k-1)!! x^(n-2k) gives the Hermite form. Setting x = 1:
Σ_k C(n,2k)(2k-1)!! counts involutions by their 2-cycle support = I(n) (RSK: = total
SYT count). The x-coefficient extraction at |K| = n (even n) gives E[det S] =
(n-1)!!. ∎

## Computational verification (this session)

- E[det(I+S)] = I(n) EXACTLY at n = 3..7 via per-iso-class census with labeled
  weights n!/|Aut| (e.g. n=7: 486539264/2097152 = 232 = I(7)); Monte Carlo at n = 8:
  763.4 vs I(8) = 764 over 200k samples.
- E[det(xI+S)] = matching polynomial of K_n verified exhaustively (all labeled
  tournaments, exact integer char polys) at n = 3..6.
- Script: 04-computation/hadamard_det_census_macmini_s2.py
  Output: 05-knowledge/results/hadamard_det_census_macmini_s2.out

## Notes

- One identity ties the session's three objects together: the same Pfaffian
  expansion gives the FLOOR 2^(n-1) (all matchings forced to ±1, local orders,
  THM-468), the CEILING (n+1)^((n-1)/2) (Pfaffian minors as large as skewness
  allows, DRT/skew-Hadamard, THM-472), and the AVERAGE (matchings counted with
  RSK/Young-tableau weight, this theorem). Floor = circle geometry, ceiling =
  Hadamard geometry, average = tableau combinatorics — all three living on the
  staircase triangle's own objects (matchings, tableaux, involutions).
- The random tournament ensemble has expected characteristic polynomial equal to a
  Hermite polynomial — and this IS in the literature (see Attribution): KMPRS 2025
  for tournaments, Godsil–Gutman 1981 for the symmetric-signing antecedent.
- Engineering hook: E[det] = I(n) gives a cheap sanity invariant for testing
  tournament samplers / enumeration code (any bias shows up against A000085).

## Attribution (adversarial check, 2026-06-11)

Parts 1 and 3 are NOT new. S. Klanderman, M. Montee, A. Piotrowski, A. Rice,
B. Shader, "Determinants of Seidel tournament matrices," Linear Algebra Appl. 707
(2025), 126-151 (arXiv:2406.09697; DOI 10.1016/j.laa.2024.11.011):
- Their Theorem 4.1 = part 3: E(det X) = (n-1)!! over Seidel matrices of
  tournaments of even order n.
- Their Theorem 7.7 = part 1: the expected coefficient of x^(n-2k) in the
  characteristic polynomial is the number of k-edge matchings of K_n; they display
  c(x) = Σ_k C(n,2k)(2k-1)!! x^(n-2k) and themselves remark it is the matching
  polynomial of K_n / probabilist's Hermite up to sign alternation.
- Their Theorem 7.8 extends to any graph G (expected char poly = signless matching
  polynomial of G); they note it also follows from Hou–Lei, EJC 18 (2011) #P156,
  Thm 2.3.

What THM-473 adds: the x = 1 evaluation E[det(I+S)] = Σ_k C(n,2k)(2k-1)!! = I(n)
= A000085 (involutions = total SYT count via RSK) — an immediate corollary not
stated by KMPRS — plus the explicit i^(-n)He_n(ix) closed form and the
staircase/Young-tableau framing.

Antecedent: C. D. Godsil, I. Gutman, "On the matching polynomial of a graph,"
Algebraic Methods in Graph Theory (Szeged 1978), Colloq. Math. Soc. Janos Bolyai 25
(1981), 241-249: the expected characteristic polynomial of a uniformly random
SYMMETRIC ±1 signing of A(G) is the matching polynomial μ(G,x) — the identity
underlying Marcus–Spielman–Srivastava interlacing families. The skew (tournament)
signing yields the signless version. KMPRS do not cite Godsil–Gutman; that
genealogy is this project's observation.

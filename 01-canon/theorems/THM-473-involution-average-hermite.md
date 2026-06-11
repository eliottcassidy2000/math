# THM-473 — The Involution Average: E[det(I+S)] = involutions, E[det(xI+S)] = Hermite

**Status:** PROVED (mac-mini-2026-06-10-S2) — adversarial verification PENDING (this session).
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
  Hermite polynomial — the same answer as the GUE/GOE expected-characteristic-
  polynomial heuristic. The tournament ensemble is "antisymmetric Bernoulli";
  whether this identity is in the random-matrix literature is being checked
  (session literature thread); the proof above is elementary either way.
- Engineering hook: E[det] = I(n) gives a cheap sanity invariant for testing
  tournament samplers / enumeration code (any bias shows up against A000085).

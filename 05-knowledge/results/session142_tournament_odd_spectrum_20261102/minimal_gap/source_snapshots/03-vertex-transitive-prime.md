---
id: THM-4861
title: "(a) The least orders of tournaments whose number of Hamiltonian paths is a Fermat prime are exactly 3, 4, 6, 8, 11 for 3, 5, 17, 257, 65537 (the Fermat numbers were known to be counts of the reversed-extreme-arc family at orders 3, 4, 6, 10, 18, THM-4749; two strong 6-tournaments have 17 paths, thirty 8-tournaments have 257, explicit 11-tournaments have 65537; minimality from the census maxima 3, 15, 189 and OEIS A038375's 15,745 at order 10). (b) The order-9 census re-confirmed (191,536 tournaments, 1,520 distinct counts, maximum 3,357 by three circulants, unattained odd values 7, 21 and a band from 2,883; first computed in tournament_pseudoprimes_20261021). (c) A vertex-transitive tournament has a prime number of Hamiltonian paths only if it is C_3 (|Aut| | H by THM-4747, strong floors by THM-1370); in particular no circulant tournament has a Fermat prime >= 5 as its count. (d) Fermat primes F_m, m >= 1, are congruent to 5, 17 (mod 36) alternately, hence free vertices of the prime square-sum graph of THM-4860, and no two Fermat numbers are adjacent. (e) In Ogg's layers the four Monster-not-Baby primes are all free and the fifteen supersingular primes hold four restricted ones: COINCIDENCE (null probabilities 1/16 and about 0.13)"
status: >
  PROVED (klein-2026-11-02, elementary; not independently audited) for (c), (d) and the minimality arguments in (a);
  FINITE-EXACT for the census of orders 3-9 (nauty gentourng + a C subset DP; orders <= 9 reproduce THM-4810 and tournament_pseudoprimes_20261021, every witness recounted by an independent
  pure-Python DP with byte-identical python3 / python3 -O receipts) and the hill-climb witnesses; CITED: OEIS A038375
  for max H(10) = 15,745 (the elementary bound 10 * 3,357 suffices too), THM-4747 (|Aut| divides H), THM-1370 (strong
  floors, omission of 7 and 21; completeness of the odd spectrum is the OPEN Chan-Pak conjecture, MISTAKE-559),
  THM-4749 (4) / THM-4748 (5) (H(E_n) = 1 + 2^(n-2): the Fermat numbers as counts at orders 2^r + 2), THM-4810 and
  THM-4744 (least orders 6 and 8 for 17 and 257 are readable from their inventories), THM-4860 (residue classes). The Paley values H(P_q), q = 3, 7, 11, 19, 23, were
  recomputed and agree with THM-4747 / CONJ-002 / THM-336; they are not new. (e) is a typed verdict, not a theorem.
session: klein-2026-11-02 (codex/session-ramsey-20261030)
source: 05-knowledge/results/fermat_monster_hamiltonian_counts_prime_classes_20261102/README.md
related:
  - 01-canon/theorems/THM-4747-hamiltonian-cycles-under-automorphisms-fermat-wilson-congruences-for-circulant-and-paley-tournaments.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md
  - 01-canon/theorems/THM-4860-the-prime-square-sum-graph-is-bipartite-by-residue-mod-4-two-primes-1-mod-3-are-never-adjacent-chebyshev-bias-obstruction-to-prime-square-sum-paths-none-below-615731-and-the-25-primes-below-100-by-ogg-layer.md
  - 05-knowledge/hypotheses/HYP-9452-least-orders-of-hamiltonian-path-counts-the-top-band-of-unattained-values-and-paley-counts-at-the-remaining-supersingular-orders.md
---

# THM-4861 — Fermat primes as Hamiltonian path counts; vertex-transitive counts are never prime beyond `C_3`; Fermat and Monster primes in the prime square-sum graph

`H(T)` is the number of Hamiltonian paths of a tournament `T`; it is odd (Rédei).

**(a) Least orders.** Let `n(h)` be the least order of a tournament with `H = h`. Then `n(3) = 3`, `n(5) = 4`, `n(17) = 6`, `n(257) = 8`, `n(65537) = 11`. (Attainment of all Fermat numbers was known: the transitive tournament of order `n` with its extreme arc reversed has `1 + 2^{n−2}` Hamiltonian paths, THM-4749 (4), giving `F_r` at order `2^r + 2`; the least orders for 17 and 257 are implicit in THM-4810/THM-4744; the order 11 for 65537 is new.)
*Proof.* Attainment: the census of all tournaments of orders 3–9 (2, 4, 12, 56, 456, 6880, 191,536 isomorphism classes) gives `H = 3` at order 3 (`C_3`), `5` at order 4, `17` at order 6 (two tournaments, both strong, scores `1,1,2,3,4,4`), `257` at order 8 (thirty tournaments); three 11-tournaments with `H = 65537` were found by hill-climbing on arc flips and recounted independently (codes in the source packet). Minimality: the census maxima are `3, 5, 15, 45, 189, 661, 3357` at orders 3–9 (agreeing with A038375), so `5 > 3`, `17 > 15`, `257 > 189` exclude the smaller orders, and `65537 > 15,745 = max H(10)` (A038375; or `max H(10) ≤ 10 · max H(9) = 33,570`, since deleting the last vertex of a Hamiltonian path of `T` gives one of `T − v`). ∎

**(b) Order 9 (re-confirmation).** 191,536 tournaments, 1,520 distinct values of `H`, maximum 3,357, attained by the circulants `Circ(Z_9, S)` with `S = {1,2,3,5}, {1,2,4,6}, {1,5,6,7}`; the odd values `< 3357` not attained are 7, 21 and odd numbers in a band beginning at 2,883. This census was first made by mac-mini-2026-10-21 (`04-computation/experiments/tournament_pseudoprimes_20261021/analyze_census.out`, the source of THM-4747's order-9 check); the present counter reproduces it. (FINITE-EXACT, inherited.)

**(c) Prime counts of vertex-transitive tournaments.** If `Aut(T)` is vertex-transitive on `n` vertices, then `H(T)` is prime only for `T = C_3`.
*Proof.* `Aut(T)` acts freely on Hamiltonian paths, so `n | |Aut(T)| | H(T)` (THM-4747 (1), LEM-003); a prime `H` forces `H = n`. A vertex-transitive tournament is regular, hence strong (a regular tournament has no dominated proper subset), hence `H ≥ f(n) = min{3^a 5^b : 2a + 3b = n − 1}` (THM-1370), and `f(n) > n` for `n ≥ 5` (`f(5) = 9`, `f(7) = 25`, and `f` is increasing). The remaining case `n = 3` is `C_3` with `H = 3`. ∎ Hence no Fermat prime `≥ 5` is the count of a circulant tournament; the census of all circulants of orders 5, 7, 9, 11 (counts `15; 175, 189; 3159, 3267, 3357; 92,411, 93,027, 93,467, 95,095`) agrees.

**(d) Fermat primes in the prime square-sum graph.** For `m ≥ 1`, `2^{2^m} ≡ 4, 7 (mod 9)` alternately and `F_m ≡ 1 (mod 4)`, so `F_m ≡ 5, 17 (mod 36)` alternately; these are free classes of THM-4860 (d″) (each adjacent to two classes). `F_m + F_k ≡ 2 (mod 4)` is never a square. `F_0 = 3` is the vertex shared by the three restricted class-1 pieces. First square-sum prime partners: `13, 11, 19, 67, 2063` (squares `16, 16, 36, 324, 67600`). ∎

**(e) Monster layers (typed).** Ogg's `M∖B = {29, 41, 59, 71}` lies entirely in free classes (`29, 5, 23, 35 (mod 36)`), `Fi22` and `B∖Fi22` contain two restricted primes each (`7, 13`; `19, 31`), the ten non-supersingular primes below 100 contain seven restricted ones including the three isolated primes `37, 43, 67` of THM-4860 (c). Under a fair coin per prime the probability of at most four restricted primes among the thirteen odd supersingular primes above 3 is about 0.13, and of `M∖B` being restricted-free `1/16`. **COINCIDENCE**: Ogg's list is defined by the prime divisors of `|M|`, the classes by the squares mod 36.

Data, codes and receipts: the [source packet](../../05-knowledge/results/fermat_monster_hamiltonian_counts_prime_classes_20261102/README.md).

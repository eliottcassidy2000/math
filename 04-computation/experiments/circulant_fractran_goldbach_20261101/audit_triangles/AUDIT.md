# Audit: `quadratic_2916_odd_pair_triangles_n4_20261102.md`

**Verdict: SOUND WITH CORRECTIONS.** All numbers and proofs check. The corrections are one INDEX overclaim, missing credits for prior and concurrent work, and minor wording.

Independent audit, 2026-11-02, of commit 917c974679.
- `audit_independent.py` is stdlib only and imports nothing from the audited code. It reads the audited `.out` only to compare the printed triangles.
- 32/32 checks agree (`audit_independent.out`); normal and `-O` runs are identical.
- Both audited scripts replay byte-identically apart from the timing line.
- The `.py` is unchanged since its pre-registration commit 5788ec6eb9.

## Findings

1. **OVERCLAIM.** INDEX.md line 3 says "equal the Goldbach pairs exactly when m-1 or m+1 is prime". Only the "if" direction is proved; the converse is FINITE-EXACT to 20000. The note body is correct. **Fix:** "...when m-1 or m+1 is prime (proved); strictly larger for all 5,820 other even m <= 20000 (FINITE-EXACT; open)".
2. **CREDIT.** §2.3 does not cite `jacobian_eth_prime_wait_20261102/PRIME-MRB.md`, committed 19:39, 52 minutes before pre-registration. It already has:
   - the same p_*(m) = least prime = ±1 (mod m);
   - p_* >= m-1;
   - a census of even m <= 10000;
   - neighbour-prime cases with delta = 1.

   C8 matches its rows. In its notation the proposition is R(m) = G(m) when delta(m) = pi(p_*) - pi(m-2) = 1, which is equivalent to m having a prime neighbour. **Fix:** cite it. The new content is R/G, P < m^2, the 41.8% share and the 16012 record.
3. **CREDIT.** `session146_carrier_primegame_20261103/prime290/README.md` (committed 20:13, before pre-registration) proves the same odd-numerator classification with the same escape inequality. It also treats the jagged triangular odd-pair array. **Fix:** cite it as concurrent. Optionally also cite:
   - `collatz_one_out_edge_systems_20260930.md`, which has the same eight points;
   - Walde-Russo 1994: c = -29/16 is the t = 1 member of their rational 3-cycle family (A11).
4. **MINOR.** §1 and the `.py` docstring date THM-4139 to 2025-08-25; it is 2026-08-25 (git e01cde8b21). Also, "eight quarter-integers plus the exact 3-cycle" reads as 11 points; write "including".
5. **MINOR.** The §2.2 paragraph "What survives" is untyped but GENERIC. Its chain is flip lines then converse, with shapes (2,2,2,2), (2,1,1) and composite (4,2,2) = the f^2 fibres (B6). Equal-shape partition chains are always isomorphic, and this one uses flip lines, not the owner's class nodes.
6. **MINOR.** In §2.1, "units mod 16" is GENERIC (true of any 8 consecutive odd integers). It is also dynamically empty: g(a+16) = g(a)+8 (mod 16) (A8). The §2.2 rebuttal "(Z/16)^× is not GF(2)^3" misses this. B in fact has a canonical cube, b = (a-1)/2 mod 8, in which -x is t -> t+111 (A9, A12). This too is GENERIC.
7. **MINOR.** §2.3, "largest P(m)/m is 17": P(16012) = 272,203 = 17·16012 - 1, ratio 16.99994, the unique maximum. "Having seen exactly the odd primes below m": when P(m) = m-1, it has seen those <= m-3. This is harmless, since the residue m-1 pairs only with 1.
8. **MINOR.** §2.4, "primes modulo increasing even moduli": the inner loop is increasing primes at a fixed m, and m indexes the rows. Rephrase.
9. **MINOR.** "Hourglass atom" occurs nowhere else in the repo. HYP-3053's C(n-2,2) = 1 hourglass cycle is in the full C(n,2)-tile flow, not the 3-tile chart of §2.2. Label the identification a reading.

## Confirmed

1. **The map.**
   - g(a) is odd for odd a. For |a| >= 9, 4(g(a)-|a|) = (|a|-2)^2 - 33 >= 16. With the closure of {±1,...,±7}, this proves bounded iff |a| <= 7 (brute force to 200001).
   - The chain is 8, 4, 3, 3 with fibres (2,2,2,2) = {±x}, then (2,1,1). The depths are as stated.
   - This is THM-4139 §1; the Poonen credit is fair.
2. **The n = 4 metagraph.**
   - The classes have fibres {5,1,1,1}; the merge gives {5,2,1}; there are 4 flip lines.
   - The line multigraph is a star with a loop at the strong class. f on f(B) is a triangle plus a pendant edge. None of the 24 bijections works; the merged level also fails; the composites differ ({5,2,1} vs {4,2,2}).
   - No structural map:
     - Aut(f|B) = Z/2, the swap 3/4 <-> -3/4.
     - The Klein four-groups ⟨flip, transpose⟩ (fixed-point profile (0,0,4)), ⟨-x, Aut(f|B)⟩ (0,2,6) and the units {1,7,9,15} mod 16 (0,0,0) are pairwise non-conjugate.
   - The NUMEROLOGY verdict is fair.
3. **Primes until ±1.**
   - P(m) < m^2; 4179/9999 and 4178/9998.
   - G ⊆ R, and R = G with a prime neighbour. The proof is correct: P(m) >= m-1, and the scanned primes are then the odd primes <= m-3.
   - R ⊋ G for all 5,820 other m.
   - **I cannot prove the converse.** Each case has a witness {p-m, 2m-p} with 2m = p+q, q < m < p and p-m composite (C7). Goldbach for 2m gives such a pair, and it is extra iff p-m is composite. This needs a refined Goldbach statement, so OPEN is right.
4. **{10, 20, 39}.** The formula gives 10, 20, 39 (69 at n = 6). The note's display equals the general formula in the chart_clique_norm README.
5. **The triangles.** Both the odd-pair rows and the PRIMEGAME rows reproduce, and D(n) = 1 iff n is prime. The credits to THM-4856, THM-4853, THM-4858, session65 and THM-4851 are accurate.

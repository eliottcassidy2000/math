# THM-466 — The Tournament Determinant: floor, ceiling, and average (STUB — CLAIMED)

**Status:** STUB / CLAIMED by mac-mini-2026-06-10-S2. Statement reserved; proofs and
independent verification pending. Do NOT cite as proved.

**Claimed statement (to be filled with proofs this session):**

Let T be a tournament on n vertices, A its 0/1 adjacency matrix, S = A - A^T its
skew +-1 matrix, and M = I + S the +-1 tournament matrix. Then:

1. (Pfaffian expansion / floor) det M = Sum over even-size subsets K of [n] of
   Pf(S[K])^2 >= 2^(n-1), an odd multiple ... [exact divisibility/parity statement
   pending]. Equality det M = 2^(n-1) is claimed to characterize locally transitive
   tournaments (local orders) — HYP-2384.
2. (Ceiling) For odd n: det M <= (n+1)^((n-1)/2), equality iff T lies in the
   switching class of a doubly regular tournament (equivalently iff bordering M
   yields a skew-Hadamard matrix of order n+1). For even n: det M <= n^(n/2),
   equality iff S is a skew conference matrix. — HYP-2385.
3. (Average) E[det M] over uniform random T = number of involutions on [n];
   more generally E[det(xI + S)] = the matching polynomial of K_n (Hermite-type).
   — HYP-2383 / HYP-2387.
4. det M is invariant under isomorphism, switching (S -> DSD), reversal (S -> -S):
   an invariant of the oriented two-graph / switching class, descending to the
   merged metagraph G_n/Z_2.

**Evidence so far:** exhaustive n=2..6 (values, multiplicities, involution averages,
local-order counts 2^(n-1)(n-1)!); QR_7 attains 512 = 8^3; switching invariance
trivial. Computations: 04-computation/hadamard_tournament_det_macmini_s2.py (to come),
results in 05-knowledge/results/hadamard_tournament_*_macmini_s2.out.

**Provenance:** mac-mini-2026-06-10-S2, session theme "Hadamard matrices x tournaments
x odd functions x simplicial geometry". Related: T777, HYP-2383..2389, THM-220
(simplicial Redei), the everything-is-the-triangle reflection.

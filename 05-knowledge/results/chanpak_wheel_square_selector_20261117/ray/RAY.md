# Physical 99-slope rails

The first-25-prime exponent wheel of the earlier `recurrence_first25` packet starts from (N_k=99⋅2^k+23). Its least period is 61,209,720; 32,353,278 phases survive, and (k=16) is the first composite survivor at positive heights. The wheel classifies divisibility by these selected primes, not absolute placement. Here every member, **independently of wheel status**, has a literal strong tournament on (k+7) vertices, for every (k≥1). The same deletion tournament recontacts to give the exact companion (N_k-14), and changed sources give (N_k+2) and a 7-multiple in a uniform nearby window at every height.

## One-exit proof

Fix a tournament (R) and a nonempty tagged set (A) in it. Write (J=H(R)). For an ordered two-cover ((X,Y)) by directed paths (possibly empty, all old vertices used exactly once), let (W) count those with (Y) starting in (A); let (U) count those additionally having (X) end in (A). Set (K=J-U).

Create a vertex (q) beating precisely (A) in (R). Add a nonempty tournament (B) such that (R →  B →  q), with the displayed relations uniform. Let (h=H(B)), and let (P_2(B)) count its ordered two-covers with empties allowed. Any maximal run in (B) which does not end the full path must be followed by (q), so there are at most two runs. If all of (B) is immediately before (q), the exterior contributes (W+J): for (X,B,q,Y), either (Y) starts tagged or (Y) is empty. If the sole (B)-run is terminal, in (X,q,Y,B) necessarily (Y) is nonempty, since (q) loses to all of (B). There are (W-U) such exterior choices. For two (B)-runs the order is (X,B_1,q,Y,B_2), and there are (W) exterior choices and (P_2(B)-2h) ordered nonempty covers of (B). Summing gives (W P_2(B)+K h).

All rows below are verified strong, with a nonempty tagged set. Consequently, the assembled graph is strong: (q) reaches one tagged vertex, hence all of (R), then every vertex of (B); all of (B) reaches (q) directly, and all of (R) reaches (q) through (B). For (B=TT_k), (P_2(B)=2^k) by choosing the subset of its unique increasing order assigned to the first path; (h=1). This is the advertised all-height formula.

## Fixtures and exact selection

Rows are zero-based outgoing binary masks; bits of the tag mask specify exactly the vertices beaten by (q). The subset-path DP in `replay.py` independently computes (J,W,U) as incidences of disjoint complementary paths. Each row has (W=99).

| K | rows | tag mask | J | U |
|---:|---|---:|---:|---:|
| 1 | (56,41,43,16,6,24) | 12 | 23 | 22 |
| 7 | (50,8,11,33,14,22) | 6 | 33 | 26 |
| 9 | (42,48,3,6,45,12) | 40 | 37 | 28 |
| 13 | (40,13,25,32,107,6,47) | 40 | 45 | 32 |
| 17 | (114,64,59,51,66,18,44) | 66 | 45 | 28 |
| 19 | (20,21,40,51,36,3) | 33 | 37 | 18 |
| 23 | (42,48,3,6,45,12) | 34 | 37 | 14 |
| 25 | (92,109,24,64,10,93,20) | 24 | 41 | 16 |
| 41 | (118,60,16,53,96,4,46) | 64 | 41 | 0 |

Thus for the two equal-source choices (K=9,23), just redirect (q)'s tag from vertex 3 to vertex 1, keeping vertex 5 tagged. The same exterior has a gap (14H(B)) for every (B). This is **not** equality of its entire continuation profile for every other kind of outside reader.

For a prescribed absolute target (N) in this bank, division by 99 finds (K) as the canonical remainder. Test (K) against the table, then test that the positive quotient is a power of two of height at least one. This recovers the actual source and block, not just a phase. At an arbitrary fixed (k≥1), choose (K=(13,19,17)_{k mod 3}). Exactly that one of these three rails is divisible by 7, since (99≡1 mod 7) and (2^k≡1,2,4 mod 7). It always lies between (N_k-10) and (N_k-4). The quotient is prime at heights 1, 2, 4, 6, 8, 9, 13, 14 (trial division; there is no prime-infinitude claim). The pair with (K=23,25) gives two primes together at k=2,4,6,12,14.

Even this changed-source bank cannot always furnish a survivor of the first-25 wheel. For (k≡0 mod {158340}), all nine rails have an explicit small factor, in the table's order: 5,53,3,7,29,59,61,31,5. Each claim reduces to a finite order dividing 158340, tested in the replay. The wheel still exposes a useful distinction: for the favored (K=23), exact absolute geometry has already been solved for the ray, including rejected heights and composite survivors; there is no need to mistake sieve readiness for a realizer.

The nine rows also provide a complete physical residue selector for arbitrary nonempty B with 7 not dividing h: pick the tag and core from the table so that K ≡ −P₂(B)/h mod 7. Every residue occurs among the K's, and all these core changes use exactly the same slot rule. If 7 divides h, the read is divisible by 7 precisely when 7 also divides P₂(B), independently of the row. This does not force a prime quotient, nor a prescribed absolute quotient.

The simultaneous composite phase generalizes to any finite menu of rails of the form Aᵢ·2^k+Bᵢ which at some common k₀ have odd values greater than one, and each Aᵢ>0. Pick an odd prime divisor pᵢ of each value there and repeat k₀ modulo the least common multiple of the orders of 2 modulo the pᵢ. Every rail keeps its divisor as its value grows. The primes in this general argument need not be among the first 25.

Finally, this is a finite menu of exponentially spaced exact targets, hence (O(log X)) outputs up to (X). It does not imply density one or coverage outside its stated rails.

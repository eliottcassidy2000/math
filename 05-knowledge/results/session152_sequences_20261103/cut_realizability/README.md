# Tournament-specific cut restrictions

Status: exact necessary conditions, a relaxation counterexample, and targeted finite checks. These results do not prove a spectrum induction. They do not claim that the restrictions characterize tournament cut carriers, nor that either H=9 example precludes a different admissible state bank.

Use the convention of `cut_induction`: H is the number of Hamiltonian paths in an n-vertex tournament T, F(S) counts paths after adding x with x→i exactly for i∈S, and K_ij=Q_ij+Q_ji is the nonnegative even two-path-cover coefficient. Every integral carrier c has c_ij+c_ji=K_ij/2 and F(S)=H+2 c(δ⁺S).

## Marked breaks and the one-error statistic

Let B be the number of vertex permutations with exactly one backward adjacent arc. Every ordered cover by two nonempty directed paths determines a permutation with a marked break between the paths, and conversely. A permutation with no backward adjacent arcs contributes n−1 marked breaks; one with exactly one contributes its unique backward break; all other permutations contribute none. Counting by break endpoints gives the exact identity

    Σ_{i<j} K_ij = (n−1)H + B,    Σ_{i≠j} c_ij = ((n−1)H+B)/2.

In particular B≥0 and B≡n−1 (mod 2), using odd H and even K. The carrier total must be at least ceil((n−1)H/2), a condition unavailable from strong support alone. B is a familiar tournament statistic, not claimed new: the repository's `h_drift_black_selfline_scout_codex_S14.out` already counts orders with one backward adjacency as those created by flipping their unique bad arc. The marked-break interpretation here connects it directly to the ear coefficient sum; no exhaustive literature or repository novelty claim is made.

## Every symmetric cut

For a nonempty proper S let C(S) be the sum over old Hamiltonian paths of their numbers of adjacent pairs crossing S, and let B_cross(S) be the number of permutations with exactly one backward adjacent arc, where that unique backward adjacent pair crosses S. Restricting the marked-break bijection to crossing endpoints yields the exact identity

    K(δS) = C(S) + B_cross(S).

Every Hamiltonian path crosses S at least once. Thus K(δS) ≥ H+1, because H is odd but K(δS) is even. More finely, C(S) ≥ H+E_same(S), where E_same counts old Hamiltonian paths whose two endpoints lie on the same side of S. In carrier and response terms,

    c(δ⁺S)+c(δ⁻S) = K(δS)/2 ≥ (H+1)/2,
    F(S)+F(V\S) = 2H+K(δS) ≥ 3H+1.

If T is strong, take a directed Hamiltonian cycle. It crosses S at least twice, and each cycle edge appears in n−1 of its n distinct path rotations. The remaining H−n paths contribute at least one crossing each. Consequently

    C(S) ≥ 2(n−1)+(H−n) = H+n−2,
    K(δS) ≥ 2 ceil((H+n−2)/2).

This supplements the existing one-sided strong-parent growth bound. For example n=5, H=9 forces every symmetric carrier cut to have weight at least 6.

These restrictions are not implied by the previously listed elementary relaxations even together with the global total bound. For n=5 and H=9, take the directed cycle 0→1→2→3→4→0, give 0→1 weight 20 and every other cycle arc weight 2, and add 0→2 of weight 1. The support is strong, every nontrivial outgoing cut is at least 2 (meeting strong-parent growth), and its cut values have gcd one: singleton {2} has cut 2 and {0,1} has cut 3. Its total weight is 29 ≥ 18=(n−1)H/2. Nevertheless the symmetric cut at {3} has weight 4<6, so it cannot carry any strong five-vertex tournament with H=9. This only shows nonrealizability, not failure of a tournament state bank.

## Targeted ternary obstruction

Codes use lexicographic pairs i<j, bit 1 meaning i→j, and S is x's out-neighborhood. Both order-five parents are strong and have H=9. Code 8 has child image {15,17,19,23,25,27,29,33,37,41}; the minimum complementary response sum is 17+23=40; each response 15 pairs with 33, not 15. Target 27 occurs at masks 5 and 11, target 29 at 10 and 21, and 31 never occurs. Code 759 has image {15,17,19,23,29,31,33,37,43}; target 27 never occurs, while 29 occurs at masks 1,5,6,10,21,23,25,26 and 31 at 2,9,22,27. Hence neither parent individually has all 3H, 3H+2, 3H+4 children. H alone cannot determine those choices. This is not an obstruction to a multi-parent bank or to eventual interval coverage.

## Audit and reproduction

`python3 check.py` enumerates all tournaments through n=4, separately enumerates permutations to compute H, B, crossing counts and K, and checks the exact identities and directed-cycle bound whenever a Hamiltonian cycle exists. It independently constructs each ear child and evaluates its Hamiltonian paths by subset DP, checking complementary-response identity. It also verifies the two specified n=5 codes, target masks and images. Output is recorded in `out.txt`. The script checks finite examples; the proofs above are general. Existing session 152 reciprocity and cut-induction notes already supply K's cover meaning and complement formula; this note focuses on the additional marked-break and per-cut restrictions. Search of repository results for the exact phrases 'one backward adjacency' and 'two-path covers' found the cited S14 one-error statistic, but is not a comprehensive novelty search.

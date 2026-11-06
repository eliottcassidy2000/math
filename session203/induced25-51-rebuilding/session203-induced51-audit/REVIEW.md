# Independent audit: induced H25 cannot occur inside H51

6 November 2026. **PASS in the stated scope.** This is an all-order ancestry obstruction and a sharp finite repair, not a missing-count theorem or a spectrum threshold.

## Finding

No finite tournament with exactly 51 directed Hamiltonian paths contains an induced subtournament with exactly 25 directed Hamiltonian paths. No bound on either tournament's order or source family is imposed.

There is a six-vertex H25 tournament which reaches H51 by one old-arc reversal and one birth. Thus the minimum number of old-arc reversals, when arbitrarily many births are allowed, is exactly one for this parent and target.

## All-order reduction

Hamiltonian-path counts multiply over the linearly ordered strong components of a tournament. A nontrivial strong component has count at least 3. Consequently the nontrivial components of an H51 tournament are either one H51 component, or one H3 component and one H17 component; all other components are singletons.

Induced Hamiltonian-path counts cannot increase under deletion: inserting each deleted vertex at the first available compatible position extends each old path injectively. Intersections of an induced subtournament with the ambient strong components remain an ordered join, so their counts multiply, with an empty intersection contributing 1.

In the H3/H17 case, an induced H25 subtournament would have factors a and b with a <= 3, b <= 17 and ab = 25. Since a divides 25, a = 1, which would require b = 25, a contradiction.

In the single-H51-component case, singleton components contribute only 1, so the induced H25 would have to occur inside that strong component. The inherited Busch lower bound H(S) >= 5^((|S|-1)/3) gives |S| <= 8: for |S| >= 9, 5^8 = 390625 > 51^3 = 132651. It remains to check all strong H51 tournaments through order 8.

## Complete finite census

The supplied bank has 7,412 representatives and SHA-256 5d5a5f625b4a85e04c0744aeb847ed36f59456f6a39996479aa82b9eb950afe3. Its completeness is proved by the inherited Session 194 growth checker, not merely by known class totals: starting from the empty tournament, it canonicalizes every one-vertex attachment to each preceding-order representative and verifies exact equality with the next bank. Canonicalization retains full adjacency and exhausts permutations inside fixed outdegree classes. The checker was reread and replayed successfully; its receipt is GROWTH_REPLAY.json. There are 62,423 extension presentations and class counts 1,1,2,4,12,56,456,6880 for orders 1 through 8.

The new independent check.py uses outgoing push transitions for paths ending at a vertex. It validates the bank and computes every induced-subset count. There are exactly 12 strong H51 classes: 8 of order 7 and 4 of order 8. None contains an H25 subset. The full row matrices and induced-count histograms are retained in NORMAL.json.

As an independent count calculation, literal recursive enumeration of every simple directed path gives the same 6,144 subset counts across all 28 H51 representatives, including the 16 weak order-8 presentations. The empty subset has count 1 in both calculations. Positive H25 witnesses exist in the bank (2,5,9 classes at orders 6,7,8), so the test is not vacuous. Normal and optimized Python runs have identical receipts and do not rely on assert statements.

## Sharp physical escape

The parent has outgoing rows (8,13,17,20,35,15), H=25 and two Hamiltonian cycles. Reverse its old arc 0 -> 3, yielding (0,13,17,21,35,15), H=11. Add vertex 6 which beats only old vertex 1. The child rows are (64,13,81,85,99,79,2), H=51 and three Hamiltonian cycles. check_repair.py validates every pair, the exact one-arc edit, preservation of all other old edges, and the birth mask. Literal permutations independently give (H,C)=(25,2),(11,0),(51,3).

The universal obstruction excludes zero old-arc reversals even with arbitrarily many births. The displayed physical construction uses one, so the reversal count is sharp. This does not establish a uniform rebuilding operation for other target counts.

## Reproduction and dependencies

From the shared scratch root:

    g++ -O2 -std=c++17 session194-weak-parent-certificate/verify_completeness.cpp -o /tmp/s203_growth_replay
    /tmp/s203_growth_replay session194-weak-parent-certificate/inherited_tournaments_through8.txt
    python session203-induced51-audit/check.py
    python -O session203-induced51-audit/check.py
    python session203-induced51-audit/check_repair.py

The all-order argument depends on the standard Busch strong-order lower bound, SCC multiplicativity, and induced monotonicity. The bank and growth checker are the already reviewed Session 194 inputs. This audit independently calculates the relevant path counts and replays that completeness proof, rather than claiming a newly implemented isomorphism engine. The author's exploratory final zero counters were unused variables and were not treated as evidence.

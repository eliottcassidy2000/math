# A tempting binary induction loses its root menu

**PROVED scoped identity and explicit failure of closure. FINITE-EXACT census through order 6.**
The fixed-root Z transfer is inherited; the following is a test of a possible sufficient invariant, not a new universal count theorem.

Write a physical rooted state as (D,E,U,d). Suppose E=d and

    U=D+d+c, where c is either -1 or +1.

Then Z sends this state to

    (2D+c, d, 2D+d+2c, d).

The same condition, with the same sign c, survives. Thus a single such root gives the infinite ray

    D_t=2^t(D+c)-c.

This is an immediate specialization of the inherited complete fixed-root Z formula. It supplies one repeated binary branch, not freely chosen binary digits.

A tempting stronger invariant is that one graph has eligible roots of both signs. It would allow a choice between 2D-1 and 2D+1. Even that condition is not preserved by the chosen extension. The four-vertex tournament with outgoing masks

    [8,1,3,6]

has D=5. Root 0 has state (5,1,5,1), of sign -1. Root 3 has state (5,1,7,1), of sign +1. Extending by Z at root 0 gives masks

    [8,17,19,22,1]

and full count 9. Its only eligible sign is -1. Extending instead at root 3 gives

    [24,17,19,6,8]

and full count 11. Its only eligible sign is +1. The complete root profiles in CHECK.json certify both failures, so this is not a failure caused by refusing to reroot the resulting graph.

The independent labelled census through order 6 gives sign-plus counts

    5,9,11,17,19,23,29

and sign-minus counts

    1,3,5,9,13,17,23,25,33.

Only count 5 has a single graph with both signs in this finite universe. The intersections of the two scalar lists are larger: separate graphs with the same count do not supply a common continuation-preserved root menu. No assertion about all graph orders follows from this census.

## Missing sufficient invariant

A valid binary induction would need an explicitly realizable family of full root menus, a finite seed set covering boundary cases, and for each required binary branch a physical extension followed by an allowed rebase that restores one of those menus at the new count. The examples show why a pair of algebraic sign conditions without closure does not suffice. Finding such a closed family remains open.

## Reproduction

Compile probe.cpp with C++17, -O3, -Wall, -Wextra, -Werror and -pedantic. Its complete enumeration is all 2^(n(n-1)/2) labelled tournaments at each order 2 through 6, with exact subset dynamic programming for all root profiles. Compare output to PROBE.txt. verify.py independently counts paths and reserved two-path covers by permutations on every subset for the three explicit graphs; run normally and with -O. This finite hostile is separate from the universal algebraic identity above.

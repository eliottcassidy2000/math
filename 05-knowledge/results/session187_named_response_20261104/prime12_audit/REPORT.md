# Independent audit of the even-order cubic support cutoff

**PASS, including a complete fresh radius-four replay.** This audit does not claim global twelve-vertex exclusion, universal reconstruction, or a full-response collision.

## Proof review and scope

Read the author's proof, C++ search and cubic programs, both Python checkers, and incoming T1–T4 from `two_port_twelve_principles_20261103/PROOF.md`.

For a fixed squarefree support W of (p+q)X, multiplication by U+1 assigns weight equal to the number of odd-total enlargements of W, plus its unexpanded odd indicator, modulo two. With at least two unused names, the enlargement parity is zero. With one unused name and even ambient order, the two contributions cancel; with no unused names, the total support is even and both contributions are zero. Exactly the odd supports of size at most n−3 survive. Disjoint nonempty paths have union size at least two, so the smallest retained size is three. This validates the complete named mixed-cubic cutoff, not just shared-label coefficients.

The coefficient formula includes both choices of the a-start, with actual disjoint supports. When the c-name repeats an a-name, the compatible path is necessarily a singleton, which remains in the DP. The mirrored ACC equality is consistent with the squarefree identity quoted by the author. No ordinary evaluation z=1 is used.

Every converse count difference is even. Integer finite differences therefore may be divided by two before reducing modulo two, as required. At order four, the modular cubic vanishes; the exact integer cubic need not vanish. The integer global-converse polynomial retains all 1885 coordinates, including the 1584 mixed coordinates and zero pure quadratics. Arbitrary ambient module reversal is not treated as cubic.

Inherited T1 applies to the actual port extensions by giving each port the replaced root's exterior incidences. A nonzero intrinsic coefficient modulo four implies some actual intrinsic query differs modulo four, since Boolean finite differences and evaluations are inverse integer transforms. T1 and the odd factor H(Q) therefore rule out this specific P in **every** uniform ambient quotient. This does not rule out every prime twelve-core.

The radius-four search counts source/edit presentations separately and makes no claim about distinct isomorphism classes or disjoint balls. Its endpoint-parity sieve, exact proper-subset primality check, and odd-support path test for the inherited cover-skew obstruction are implemented consistently. Only endpoint-balanced targets undergo the second test, which is sufficient to exclude all candidates passing both necessary tests.

## Independent computation

`check.py` uses literal Hamilton-word coefficient expansion from the sibling collision checker, and independently enumerates all short directed paths and disjoint pairs to implement the support formula. It verifies **5780 mixed coordinates across 81 tournaments**: every labelled order-two and order-four tournament, plus 15 explicitly listed order-six controls. This is distinct from the author's bitset path DP. Normal and optimized audit outputs are byte-identical.

The same checker uses forward Hamilton subset DP, independent of the author's predecessor recurrence, for twenty physical queries on P and its converse. It independently recovers both specified mixed coefficients, −56686 and −20574, and the zero quartic coefficient. The author's Python normal and optimized replays agree byte-for-byte. A freshly compiled full integer C++ calculation reproduces the complete author `cubic.json` byte-for-byte.

An even simpler physical check gives

    F_P(0,4)=263509,
    F_(P^op)(0,4)=206331,
    difference=57178=2 (mod 4).

Thus P was not selected as an example requiring cubic information for its first detection. It is a prime endpoint-balanced full-polynomial control; the all-order support identity is the mathematical increment. The proof does not claim cubic necessity.

The C++ command-line interface originally permitted odd small orders although its shortcut is even-order only. The author corrected it to reject odd orders; the final binary returns the expected error status for order three. The order-twelve computation was unaffected.

## Full radius-four replay

A fresh strict-warning C++ build completed both radius-four balls. All ten cumulative count rows agree exactly with the author receipt after removing elapsed-time fields: 1,537,384 source/edit presentations, 1,407 endpoint-balanced presentations, 1,089 prime endpoint-balanced presentations, and zero combined-sieve survivors. `search.out` retains the raw replay; `search_comparison.json` records the exact comparison. The Python radius-one checker additionally reproduced both sources with exact integer cover convolution rather than the search's odd-path parity shortcut; normal and optimized outputs agree with the author's original receipt.

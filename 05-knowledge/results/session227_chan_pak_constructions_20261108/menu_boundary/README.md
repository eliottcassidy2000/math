# Fixed-menu bounded-intercept extension

**Status: PROVED, independent audit PASS; finite controls PASS.**

This packet extends Session226's bounded-intercept and common-residue-gap
line obstruction from O/Z/J words to arbitrary words in O, Z, and any fixed
finite menu of nonempty rooted grafts. Gadget sizes are unrestricted but
fixed with the menu. It does not claim a new original residue obstruction.

For each gadget, its first and third state rows have a finite, sharp ratio
bound U'≤kappa_B D', computed on the four extreme rays of the inherited
physical cone. Every bounded-intercept history therefore has a core from
a finite list, followed only by O operations. Its closing count is
D+4u·2^k. The explicit core-envelope bound is

    C floor(kappa C)+#{seeds with D≤C}.

Thus bounded intercepts still yield O(log X) distinct counts through X, and
Session226's exact even-g residue-gap argument survives without change to
its essential hypothesis. Full-rank or state-escaping grafts do not by
themselves evade this restricted obstruction.

- `PROOF.md`: all-order theorem, sharp constants, explicit finite thresholds,
  singleton/sink exceptions, and precise limitations.
- `verify.py`: portable standard-library finite checker.
- `CHECK.json`, `CHECK.optimized.json`: matching normal/optimized receipts.
- `SHA256SUMS`: exact-file checksums.

The checker covers all 5,405 rooted labelled gadgets through order five,
obtaining literal permutation/bin statistics and independent subset-DP
states; all 841 rooted source/gadget pairs through order three; all 374,490
words of lengths zero through five over an eight-letter menu from ten
physical seeds; 241,695 exact line-congruence candidates; and explicit
large-tail, stationary-Z, and no-common-residue-gap hostiles. No finite
check is substituted for the arbitrary-order proof.

Reproduce in this directory:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json

Arbitrary nonunit line banks lacking a common residue gap, alternative
closing constructions, and unrestricted scalar coverage remain open.

# Bounded entry grammar

PROVED: in every strong tournament, the sum of root-start counts over any
nonempty m-vertex subset is at least floor(m^2/4)+1. At most4k-1 vertices
can have start count at mostk. The exact SCC deletion decomposition therefore
has a uniformly bounded eligible entry set; its actual entry arity satisfies
|A|^2<4k/q. For start budgets through26, these bounds are103 and10 respectively.

Every same-budget recursive step is a literal J/O/Z chain. Multi-entry and
nontrivial-tail junctions descend to smaller start budgets. This is a complete
physical decomposition, not a finite arithmetic state space or a proof of
F(8,137). The later finite-interface proposal is deliberately not included.

Reproduce author controls:

    g++ -O3 -std=c++17 verify.cpp -o /tmp/bounded-entry-verify
    /tmp/bounded-entry-verify > /tmp/bounded-entry-normal.json
    cmp NORMAL.json /tmp/bounded-entry-normal.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify.cpp -o /tmp/bounded-entry-debug
    /tmp/bounded-entry-debug > /tmp/bounded-entry-debug.json
    cmp DEBUG.json /tmp/bounded-entry-debug.json
    python3 audit/independent.py > /tmp/bounded-entry-independent.json
    cmp audit/INDEPENDENT.json /tmp/bounded-entry-independent.json

The independent permutation script uses Python assertions and is intentionally
run without -O. The author's C++ checks are explicit runtime conditions and
are retained in both builds. audit/REVIEW.md records scope and independent
checks. reviewed-originals/ freezes the exact pre-promotion candidate.

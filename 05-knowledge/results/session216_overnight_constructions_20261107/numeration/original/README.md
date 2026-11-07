# Three-core numeration lane

Read PROOF.md. The all-odd-at-least-65 inverse conjecture remains open.
This packet proves a unique balanced-prefix candidate lemma, implements exact
all-length membership via deterministic alpha stripping, and records a
restricted curvature obstruction and explicit hostile tests.

Run from this directory:

    g++ -O3 -std=c++17 balanced.cpp -o /tmp/session216-balanced
    /tmp/session216-balanced CORE_STATES.txt 0 range 10000
    /tmp/session216-balanced CORE_STATES.txt 1 range 10000
    /tmp/session216-balanced CORE_STATES.txt 2 range 10000
    /tmp/session216-balanced CORE_STATES.txt 1 1000000007 1000000033 1000000000039 1099511627775 1099511627777
    /tmp/session216-balanced CORE_STATES.txt 2 1000000000039 4295971 9948407
    python3 inverse_alpha.py
    python3 verify.py

A balanced result with core=-1 is only a failure of the tested band, including
its explicit b<27 resource ceiling. It is never an all-length absence claim.
Only inverse_alpha.py can return ABSENT_FROM_GRAMMAR after its complete finite
inverse exhaustion. No background run or uncontrolled resource search is
needed for these receipts.

CORE_STATES.txt is copied verbatim from the audited arithmetic three-core
packet. The producer uses the f recurrence; verify.py uses its independent
selected-position pair sum and also performs 320 literal subset-DP graph
counts. LARGE_GRAPHS.json retains all seven displayed large-target graphs
and explicit Hamilton cycles proving their strongness. Their large counts
are certified by the independently evaluated inherited construction formula,
not by infeasible whole-graph subset enumeration.

# Session216 inverse-mask lane

PROVED exact finite inverse-membership theorem; FINITE-EXACT positive
construction certificate. Read PROOF.md for scope, inheritance, and limits.
The executable inverse specializes to two apexes with arc0->1. No eventual
spectrum theorem or enlarged ambient finite range is claimed.

Run from this packet directory with Python3 and a C++17 compiler:

    python3 probe.py
    python3 census.py
    python3 inverse.py
    python3 budget_probe.py
    python3 digits.py
    python3 digits_boundary.py
    python3 variable_digits.py
    python3 rebuild_digits.py
    g++ -O3 -std=c++17 audit.cpp -o /tmp/session216-inverse-audit
    /tmp/session216-inverse-audit < witnesses.txt

probe.py runs341 independent transfer controls. census.py recreates all4998
witnesses on at most12 vertices. audit.cpp independently verifies every saved
graph by subset Hamiltonian-path DP. Timing fields naturally vary.

witnesses.json and witnesses.txt store the same literal mask words: bit i
means apex i beats that chain vertex. Chain arcs point in increasing word
order. A dash in witnesses.txt encodes the empty chain. H(empty tournament)=1.
The positive packet graphs always include both apexes.

Digit probes are bounded hostile tests, not all-order obstructions. The
saturation inverse exhausts reachable states before claiming ABSENT; finite
resource interruption returns INCOMPLETE. PUBLICATION_FILES.txt lists the
release payload, and SHA256SUMS pins it. No dependencies beyond standard
Python and C++ libraries.

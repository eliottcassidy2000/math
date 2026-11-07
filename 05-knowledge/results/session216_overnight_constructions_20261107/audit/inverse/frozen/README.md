# Session216 inverse-mask lane

Read PROOF.md for exact scope, inheritance, proof, and limits.

Reproduction from the containing workspace directory:

    python session216-inverse/probe.py

The original breadth-first exploratory driver in probe.py is memory-intensive
and should not be used to reproduce the certificate. Instead use:

    python session216-inverse/census.py > session216-inverse/census.out
    python session216-inverse/inverse.py > session216-inverse/inverse.out
    python session216-inverse/budget_probe.py > session216-inverse/budget_probe.out
    python session216-inverse/digits.py > session216-inverse/digits.out
    python session216-inverse/digits_boundary.py > session216-inverse/digits_boundary.out
    g++ -O3 session216-inverse/audit.cpp -o /tmp/session216-inverse-audit
    /tmp/session216-inverse-audit < session216-inverse/witnesses.txt

probe.py is imported for its transfer implementation and independent Python
subset counter. witnesses.json and witnesses.txt store identical literal
mask words, with bit i denoting apex i beating that chain vertex. The two
apexes have arc0->1; the chain is increasing. No external dependencies.
Timing fields vary between runs.

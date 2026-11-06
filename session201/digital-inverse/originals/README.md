# Digital gcd inverse

Unreviewed all-order inverse and search-node bound. Universal occupancy is OPEN.

    python verify.py
    python -O verify.py
    python verify.py --target 2134117 --all
    python verify.py --target 100000837

`--all` exhausts the complete finite-template inverse. Without it, a positive
answer stops at its first exact witness. This does not claim optimality.
The code uses exact integer arithmetic, never floating cutoffs.

`PROOF.md` explains the singleton arithmetic boundary and why it is not a
physical singleton-Z operation. The final nontrivial unit indices are the
actual graph certificates.

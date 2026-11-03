# Order-ten finite upper certificate

Read PROOF.md together with the independently certified order-ten lower-bound packet and the previous order-nine lifting proof. The core verifier enumerates all canonical saturated diagrams and checks exact integer annihilation by a 56-dimensional proposed kernel. To replay:

    cat kernel.part00 kernel.part01 > kernel.txt
    sha256sum kernel.txt
    g++ -O3 -std=c++17 -o verify_upper verify_upper.cpp
    ./verify_upper
    python3 verify.py
    python3 -O verify.py

Expected joined SHA-256 is 91c68094d7c5d74467a9de29e949b6bbabe844edd0e2c0193bca9d5b241522fb. Expected output is `kernel candidate rank 56` followed by `checked 32400 zero 17042 result PASS`. Progress is sent to stderr. The saved receipt is verify.log. The kernel matrix has 735 rows and 56 columns, indexed by the canonical representatives from the separately audited order-ten r=4 block.

The included verifier shares a coefficient routine with the discovery producer. The independent `verify.py` builds and runs a separate diagram enumerator, with exact integer checks and hostile controls; `audit.normal.log` and `audit.optimized.log` record identical passes. The preceding lower-bound certificate is essential for exactness.

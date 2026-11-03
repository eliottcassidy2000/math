# Order-nine endpoint-set block certificate

See PROOF.md for the exact theorem, the triangular argument, scope, dependence on the previously certified order-eight lifting, and the higher-order lower bound. The four banks separately certify one diagonal block per endpoint-set size. Replay with a C++17 compiler and Python 3:

    g++ -O3 -std=c++17 -o verify9 verify9.cpp
    for r in 1 2 3 4; do ./verify9 "$r" "block$r.contexts" "block$r.reps"; done
    python3 hostiles.py
    python3 -O hostiles.py

Expected prime-field determinants and rank checks, including deletion and duplication controls, are in verify.log. The separate exact-cover verifier independently reconstructs each signed representative and validates all context endpoints. Each context line is `N k u1 v1 ... uk vk` with old vertices 0,...,8. SHA256SUMS records the files.

`generate_blocks.cpp` is the discovery program; it is not needed for certificate verification. To rerun discovery in a temporary working directory, compile it and run `./generate_blocks r` for each block; output files are written in that directory. Random-shuffle results can vary with the standard library. See the order-seven packet for direct Hamilton coefficient and actual-attachment Walsh controls.

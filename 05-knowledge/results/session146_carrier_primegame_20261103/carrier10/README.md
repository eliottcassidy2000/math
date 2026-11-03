# Order-ten endpoint-set lower-bound certificate

See PROOF.md for the rank argument, dependence on exact order-nine lifting, scope, and all-order lower bound. The r=4 bank is deliberately rectangular: it proves 679 independent rows among 735 columns and does not prove an upper bound on that block or the exact order-ten dimension.

To replay with a C++17 compiler:

    cat block5.contexts.part00 block5.contexts.part01 > block5.contexts
    sha256sum block5.contexts
    g++ -O3 -std=c++17 -o verify10 verify10.cpp
    for r in 1 2 3 4 5; do ./verify10 "$r" "block$r.contexts" "block$r.reps"; done
    python3 hostiles.py
    python3 -O hostiles.py

Expected concatenation SHA-256 is given in BLOCK5_SHA256. Both modular ranks, square-block determinants, and deletion/duplication controls are in verify.log. The verifier independently reconstructs representatives by exact cover and checks all endpoint-saturation conditions. Each context line has `N k u1 v1 ... uk vk`, with old vertices 0,...,9. SHA256SUMS checks the distributed files.

`generate_blocks.cpp` records the discovery program and is not needed for verification. To rerun discovery in a temporary working directory, compile it and run `./generate_blocks r` for each block; output files are written in that directory. Random-shuffle results can vary with the standard library. See the order-seven and order-nine packets for additional direct Hamilton, actual-attachment Walsh, and hostile parser controls.

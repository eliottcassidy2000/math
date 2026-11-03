# Order-seven certificate and higher-order lower bounds

See PROOF.md for the exact theorem, scope, and lifting proof. The two ordered context parts concatenate to the fixed 2668-row certificate. On a system with a C++17 compiler and Python 3, independently replay it with:

    cat contexts7.part00 contexts7.part01 > contexts7.txt
    sha256sum contexts7.txt
    g++ -O3 -std=c++17 -o verify verify.cpp
    ./verify 7 contexts7.txt reps7.txt
    python3 controls.py
    python3 -O controls.py

Expected concatenation SHA-256: ee0c332d4dd580eb3073cf95eba3ca098d005cd1f3f26e14c0aaaa9c6671baef. Expected modular ranks and determinants are 2668/48 modulo 101, and 2668/66 modulo 103; the duplicated-row hostile control has rank 2667. The verifier constructs its own representative census by a different method and compares its representatives to `reps7.txt`; the latter is an additional reproducibility check, not a trusted input to the coefficients.

Optional deletion control: `./verify 7 contexts7.txt reps7.txt --delete-last` should report rank 2667 modulo each prime. The `controls.py` script separately checks direct permutation coefficients, a Walsh extraction from actual fixed attachments, and rejects malformed contexts. `reps6.txt` supports those parser controls. The two control logs document normal and optimized Python runs.

`generate_rank7.cpp` is the deterministic discovery code and is not needed for verification. To rerun the search in a temporary working directory, compile it with `g++ -O3 -std=c++17` and execute the resulting program; it writes `contexts7.generated.txt` and `reps7.txt` in the working directory. Exact pseudorandom shuffle sequences may depend on the standard library implementation, so generated certificate bytes are not the cross-platform correctness criterion.

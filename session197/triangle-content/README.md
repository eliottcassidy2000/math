# Triangle tests recover universal module content

This independently reviewed Session197 packet replaces the general matching-context test bank with repeated applications of one actual triangle wrapper.

- [All-order proof](PROOF.md)
- [Independent audit](audit/REPORT.md)
- [Standalone source verifier](verify.py)
- [Independent physical checker](audit/check.cpp)
- [Independent response-lattice checker](audit/check_algebra.py)

For a nonempty n-vertex input T, let W(T)=C3[T,K1,K1] and d=ceil(n/2). The d counts H(T),H(WT),...,H(W^(d-1)T) generate the integer scalar-response lattice of all whole-module contexts. Their gcd is exactly the universal content G(T). Integer linear combinations may use negative coefficients; this is not a positive target-construction claim.

The possibly weak identity input can be replaced for gcd purposes by a finite bank of exclusively strong W/E outputs, where E(T)=C3[T,Tr2,K1]. This proves that restricting universal whole-module content to strong outputs loses nothing. The strong bank is a gcd replacement, not an identical integral response lattice.

A separate explicit all-order family shows why the more restrictive cyclic content behaves differently under composition. The m-fold ordered join of a ten-vertex witness has contents(3^(2m),3^(2m+1),3^(2m+2)). Its strong triangle wrapper lowers the cyclic content's3-adic valuation by one. An explicit second cyclic stage escapes the original cyclic divisor at every m. At m=1 the actual16-vertex output has H=196004745, congruent to54 modulo81.

This does not establish new prime/7p coverage or infinitude. The weaker fixed-parity content statements and all-order construction scope are distinguished in the proof.

## Reproduce the controls

All three checkers are standalone and have no input-file or external-package dependencies.

    python3 verify.py
    python3 -O verify.py
    g++ -O3 -std=c++17 audit/check.cpp -o independent
    ./independent
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 audit/check.cpp -o independent-debug
    ./independent-debug
    python3 audit/check_algebra.py
    python3 -O audit/check_algebra.py

Source outputs match receipt.json and receipt-optimized.json. Independent C++ outputs match audit/independent-normal.json and independent-debug.json. Independent algebra outputs match audit/algebra-normal.json and algebra-optimized.json.

The physical audit uses literal odd-cycle collections for profiles and arbitrary-walk inclusion-exclusion for actual Hamilton counts, separately checking1099 small profiles,3222 strong outputs and both large graphs. Evaluation-grid algebra checks110 normalized response columns,700 response compilers,700 signed gcd cases and powers1 through24. The all-order proofs are not inferred from these cutoffs.

## Included unchanged mathematical dependencies

The frozen proof retains its session-local source names. Their release aliases are:

- Central-factorial operators and actual D/E wrappers: dependencies/CENTRAL_OPERATORS.md.
- Universal factorial-weighted content: dependencies/UNIVERSAL_MODULE_CONTENT.md.
- Cyclic content and the physical ten-vertex witness: dependencies/CYCLIC_MODULE_CONTENT.md.
- Cyclic Gram kernel: dependencies/CYCLIC_GRAM_FORESTS.md.

All dependencies were independently reviewed before this refinement. This packet adds the one-letter response bank, strong-only refinement and composition counterexample without editing those checkpoints. No literature-priority claim is made. RELEASE_FILES.sha256 binds every other included file.

# A physical rewrite with a shared, arbitrary module

**PROVED and independently audited in the stated scope. No eventual-spectrum threshold is claimed.**

For every tournament module U, the explicit A/B/C rewrites preserve U literally and satisfy

    H(child(U)) = 2 H(parent(U)) + H(U).

B/C replace one core vertex at the same order and change the Hamilton-cycle count from H(U) to 3H(U). A preserves all old vertices, reverses one old arc, and adds one vertex. The empty module is included with H(empty)=1.

For nonempty U, family C also admits a returning count rebase after appending a sink inside the shared module. It preserves the Hamilton-path count only. The audit retains an equal-count/different-cover-profile hostile, so later constructions must not carry the old bin coordinates through that rebase unchanged.

- [Proof and concrete cores](PROOF.md)
- [Independent audit](audit/REVIEW.md)
- [Exact status promotion](PROMOTION.json), with reviewed originals preserved
- Frozen core proof, coefficient receipt, and audit under dependencies/

The module values H(U),B_U(2),B_U(3) remain coupled graph data. These formulas do not make arbitrary correction values or target counts independently available.

## Reproduction

    g++ -O3 -std=c++17 verify.cpp -o /tmp/arbitrary-module-check
    /tmp/arbitrary-module-check > /tmp/arbitrary-module-normal.json
    cmp /tmp/arbitrary-module-normal.json NORMAL.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify.cpp -o /tmp/arbitrary-module-debug
    /tmp/arbitrary-module-debug > /tmp/arbitrary-module-debug.json
    cmp /tmp/arbitrary-module-debug.json DEBUG.json
    python3 audit/check.py > /tmp/arbitrary-module-independent.json
    cmp /tmp/arbitrary-module-independent.json audit/NORMAL.json
    python3 -O audit/check.py > /tmp/arbitrary-module-independent-opt.json
    cmp /tmp/arbitrary-module-independent-opt.json audit/OPTIMIZED.json

The author and independent checks both cover all 76 labelled modules of orders 0 through 4, 228 rewrites, and 684 physical output graphs. The arbitrary-order conclusion comes from the path-partition bijection and exact core coefficients, not these finite controls.

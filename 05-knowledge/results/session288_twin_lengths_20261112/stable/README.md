# Exact fixed-bank free-length classification

See [PROOF.md](PROOF.md). This lane classifies the **union of the original
through-order8 bank's jointly free cylinder families**, not its closure under
adaptive graph operations or newly unlocked stable ports.

Main output:98 irredundant cylinders on88 rigid marked twin quotients, with
exact all-order generating function

    (z³+4z⁵+8z⁶+20z⁷+46z⁸)/(1−z) + (5z⁷+14z⁸)/(1−z)².

Thus the fixed union has19n−49 marked classes at every n>=8. The physical
tail outputs from these79 one-clock and19 two-clock charts number at most
`79L²+19L³` below X, L=floor(log2(X−1)). This applies the inherited fixed-core
sparsity mechanism with the new exact cylinder dimensions, not a new
unrestricted spectrum theorem.

New instruments: anchored support/Stanley–Reisner obstruction complex for
acyclic ports, a canonical marked twin cylinder compiler, and finite
stable-contraction normal forms. Global confluence is conditional.

## Reproduce

From this directory, Python3, NumPy and a C++17 compiler are needed:

    python3 verify.py
    python3 -O verify.py
    cmp replay.json replay.optimized.json
    python3 cones.py
    python3 -O cones.py
    cmp cones.json cones.optimized.json
    python3 normalforms.py > normalforms.json
    python3 explore.py > explore.out
    python3 -O explore.py > explore.optimized.out
    cmp explore.out explore.optimized.out

`verify.py` compiles the local `delta.cpp` into an automatically cleaned
temporary directory. It checks the count kernel against independent permutation
and split enumeration, full symbolic skeleton tensors, literal doubled/mixed
length controls, and canonical marked classes. `explore.py` is only a seeded
12,000-extension hostile search; its m2/m3 screens are not exact stability
certificates. Its largest literal graphs have11 vertices, and every path/cover
accumulator is bounded by2(n+1)n!, well below signed64-bit capacity there.
All normal/optimized checks use explicit exceptions rather than Python asserts.

Input census completeness is inherited from the pinned atlas receipt, whose
SHA256 is in `cones.json` and `normalforms.json`; no claim is made to rerun
that entire original census. The compiler newly verifies all98 retained
all-length defect tensors and every published union coefficient rule.

## New higher-rank births

[CLOCK_BIRTHS.md](CLOCK_BIRTHS.md) gives new strong three-clock and four-clock
cores of orders9 and10, with exact series `z^9/(1-z)^3` and `z^10/(1-z)^4`.
These are disjoint from the old fixed-bank union. The former is the first
three-clock birth using the inherited complete≤8 lower-order bank; no first
four-clock claim is made. Run:

    python3 clocks.py
    python3 -O clocks.py
    cmp clocks.json clocks.optimized.json

These extra literal controls reach order14, still safely below the stated
signed64-bit path/cover bound. The complete320-template immediate next-clock
search is scoped explicitly and does not prove bounded global clock rank.

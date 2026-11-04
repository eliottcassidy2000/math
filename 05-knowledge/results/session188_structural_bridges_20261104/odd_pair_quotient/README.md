# Odd quotient reversal cubes with pair children

**All-order theorem:** for every odd-order quotient tournament Q with at least three labels, replacing each label by a fixed named two-vertex block gives a fully covered laminar family in which quotient converse is visible on actual block-uniform two-port queries. The quotient need not be prime.

The exact mixed-coefficient trace is four modulo eight. Its proof combines the classical odd-cycle expansion with a paired-label cancellation and the incoming singleton-quotient E1 identity. The two-separate-port-cycle contribution has no shared-index cubic and an even total mixed-quadratic coefficient sum. No arbitrary quartic term is assumed absent.

Consequently all 2^(q+1) states obtained by quotient converse and independent pair-interior converses have distinct complete named responses. A fixed coarse battery identifies the quotient bit, then lifted single-module queries identify the child bits. The query bounds in `PROOF.md` are constructive upper bounds and are not claimed optimal.

Even quotient order, larger child blocks and arbitrary legal inversion sequences remain outside the separation theorem. The q=1 boundary is explicitly excluded. The modulo-eight identity has a zero trace at even q and must not be read as a collision claim.

## Reproduction

Run `python verify.py` and `python -O verify.py`; their stdout and generated data must agree. The verifier uses only Python's standard library and an existing C++17 compiler. It compiles temporary local executables, with no installation.

- `trace_mod8.cpp` checks every labelled quotient of orders two through five, 1,098 in total.
- `trace_mod8_samples.cpp` adds eight fixed labelled samples at each of orders six and seven.
- `query_bank.cpp` literally counts actual extended tournaments. `verify.py` uses independent backward recursion for additional physical-query controls.
- `coarse_cases.json` retains all coarse response values and every nonzero coarse contrast coefficient for q=1, q=2, C3, and all three prime-five quotient types. The checker verifies the stronger coefficientwise mixed modulo-eight identity, not just its trace.
- `full_named6.json` retains the entire full named response bank and contrast polynomial for the six-name C3 pair blowup, including every cubic and quartic coefficient. Coarse compression can cancel physical quartics; the full named bank does not discard them.

`receipt.json` binds proof, checker sources and data. `normal.stdout` and `optimized.stdout` are the deterministic execution receipts. Independent review is coordinated separately; its audit packet records the final source pins.

The exploratory code in `../covered_prime/` is not part of this packet and is not needed to reproduce it.

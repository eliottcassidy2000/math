# Cofactor deformations: a second global classification and all two-cycles

The map is F_t(cp)=(c-1)p+2tc-1, where c is the exact power of two in the positive even input and p is odd. For t>1 this is an auxiliary arithmetic deformation, not the aliquot map.

Read in this order:

1. [Finite odd-part reduction](FINITE_ODD_PART_REDUCTION.md): every growing input with oddpart(p+1)>=2t+1 has a smaller second iterate. This reduces the residual work at fixed t to t odd-coefficient lanes with unbounded exponents.
2. [Global t=2 classification](T2_GLOBAL_CLASSIFICATION.md): the next target posed in the reduction note is now closed. Every orbit reaches a fixed point or one of six explicit cycles. Six-step descent holds outside a finite exceptional set, and the example56 shows the bound is sharp.
3. [All two-cycles for every t](TWO_CYCLE_CLASSIFICATION.md): an exact Mersenne-gcd divisibility criterion classifies every two-cycle, with dyadic cores c<=2t-2 and d<=8(t-1)^2. The second bound is attained for infinitely many parameters.

Global periodicity for t>=3 is not proved. The all-t two-cycle theorem does not exclude longer cycles or divergent trajectories. The changing-core mechanism provides a controlled comparison with the t=1 theorem; it is not a proof transfer to Collatz or general aliquot dynamics.

Run `python verify_t2_global.py` with assertions enabled. Its `assert` checks are disabled under `python -O`, so optimized execution is not a verification mode for that script. The all-t checker `verify_two_cycles.py` uses always-active guards and may also be run under `python -O`.

## Independent audit and reproduction

Read the [independent proof review](review/INDEPENDENT_REVIEW.md). The original [t=1 classification](../dyadic_cofactor_global_descent_20261019/README.md) is a separate predecessor. All three results in this package have independent audit PASS.

Use Python 3 with SymPy, from this package directory:

    python3 -B verify_t2_global.py --output replay-t2.json
    python3 -B verify_two_cycles.py
    python3 -B -O verify_two_cycles.py
    python3 -B review/independent_controls.py --output replay-independent.json
    python3 -B -O review/independent_controls.py --output replay-independent-optimized.json

Compare the first generated JSON with `t2-verification.json`, the two-cycle outputs with their normal/optimized receipts, and the independent outputs with the corresponding review receipts. The first command requires assertions enabled and must not be counted under `-O`. The independent verifier uses always-active failures; `--negative-control` must return a nonzero exit status in both modes. Independent receipts include universal symbolic identity, exact-core and strict-negativity certificates, as well as finite regression controls.

`SOURCE_PROVENANCE.json` preserves the original reviewed hashes and documents candidate-to-proved status changes and the package-relative provenance path. Mathematical proof and verifier logic are unchanged. The review manifest hashes the published artifact versions.

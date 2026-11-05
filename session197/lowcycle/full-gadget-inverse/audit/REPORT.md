# Independent review: full gadget-response inverse

Session 197. Result: **PASS**, with no correction required. This review binds the following frozen files in `session197-lowcycle/`:

- `FULL_GADGET_INVERSE.md`: `e1b1e12b59434eb2e3fca80639cb0b9b5ae71f36ab762d612b90a1635e91a261`.
- `verify_full_gadget_inverse.py`: `a2fb1f47f905425e28b76790ea30fd50ba8158a6c5fe04686316ed493a2d2b40`.
- `full_gadget_inverse_certificate.json`: `5be84c45250879ad84d2bf3e032c06909eec79a302c7310501b45b92104d7fae`.

The named physical-gadget and earlier scalar-inverse dependencies are copied and hash-pinned in the review's `full-gadget-inverse/source/` directory. The proof restates all arithmetic identities needed by the new algorithm. The source files were not changed.

## 1. Proof review

The theorem correctly imposes `A>=1, B>=A, 0<=C<=4A, D>=0, 0<=E<=A`. It concerns fixed arithmetic coefficients, not the realizability of arbitrary coefficient vectors or selection of an arbitrary physical source.

The reverse-complement identity uses the minimal width `bit_length(L-1)`. Its singleton case and the prohibition on leading-zero padding are essential and correctly stated. The low/high split and the parity coefficient in the inherited `g` identity have the correct orientation.

The two global cutoffs are valid: `F(L)>=BL`, and `F(L)>=F(2)` for every `L>=2`. For odd `L>=3`, the latter uses `g(L)>=3`, so the negative parity term does not invalidate it. The singleton is tested separately.

For each shell `k>=5`, the selected low width gives `B0>=8H`. On each exact integer high-word interval, the change in the `g` term is at least `-Ax`, the parity term costs at most `Ac`, and the reversal term costs at most `C(H-1)`. Since `u>=c`, the stated strictly positive lower bound follows. This proves binary-search correctness despite the nonmonotonicity of the unsplit functions. The interval endpoints cover every candidate exactly once, including truncated shells and powers of two.

The complexity statement is correctly an arithmetic-operation bound. The geometrically increasing low-word table sizes give `O(sqrt(K) log(2+K))` operations and `O(sqrt(K))` integer words. Table lookup covers every required `g(x+1)` and `g(b+1)`. Once the early rejection permits a nonsingleton search, the coefficient and intermediate bit lengths are bounded by `O(log(2+T))`; the explicit cost caveat for reading potentially larger original coefficients is necessary and present. No polynomial-in-bitlength or global coverage claim is made.

## 2. Source replays

The verifier was copied before running, because it writes its receipt beside itself. Normal Python and Python `-O` both reproduced the frozen receipt byte for byte:

`5be84c45250879ad84d2bf3e032c06909eec79a302c7310501b45b92104d7fae`.

These runs include invalid inputs, singleton/zero-target behavior, both coefficient-cone boundaries, physical gadget controls, complete small inverse scans, and the three larger arithmetic examples.

## 3. Separate checker

`full-gadget-inverse/check.cpp` was written independently and imports no source code. It obtains gadget counts by an assignment enumeration before/after the marked root, checks those values against literal permutations and reserved-path cuts, and evaluates the objective directly. Its inverse chooses the high-word size by the inequality `8H^2<=2^k`, rather than copying the source width formula, and binary-searches the direct evaluator rather than the source's split expression.

The optimized and checked/UBSan builds produced identical deterministic receipts:

- 128 literal physical gadget profiles.
- 32,768 minimal-width identities.
- 8,176 exact split identities.
- 115,200 strict branch comparisons.
- 110,808 complete inverse comparisons with exhaustive candidate scans.
- 162 dyadic-boundary inverse controls.
- All three larger examples, with exactly the stated unique solutions.

Checker SHA256: `e8366ae4d8be60333d209e91dc43fdc35f34fddda723b8bc9ff837cd14358d0b`.

Both independent receipt hashes: `7753e45ec7daa6c6fe748ba0446dbc3b27342bd7f6351b27a417f2bc84edff69`.

Reproduce from the review directory:

    g++ -std=c++17 -O2 full-gadget-inverse/check.cpp -o /tmp/full-gadget-inverse-check
    /tmp/full-gadget-inverse-check
    g++ -std=c++17 -O1 -D_GLIBCXX_ASSERTIONS -fsanitize=undefined full-gadget-inverse/check.cpp -o /tmp/full-gadget-inverse-check-debug
    /tmp/full-gadget-inverse-check-debug

The finite controls support the independently reviewed proof; they are not used to extrapolate the all-order inverse theorem. The later fixed-two-source and two-successive-gadget applications were not silently added to this frozen theorem's scope.

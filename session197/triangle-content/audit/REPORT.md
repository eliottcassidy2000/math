# Independent review: repeated triangle content and cyclic composition

**PASS, no correction.** This review binds the frozen triangle-content proof SHA256242a28612c17288da56ea785341e5189f9b03eea3e70fec16b4110d7fc88fcde, standalone verifier bbd76a9897f30d10bd846e39468ab644098d19d7ed9f36c90241f6bcb357fbcc, and both source receipts5924fa0547fdfec0cc1ec4c6968ffc4b9b1673dfbe6e8f91699071c478379746. The original files were not edited.

## All-order mathematical review

The tridiagonal matrix is used with input basis elements as columns, so its lowering coefficient is beta_k in the asserted orientation. A height-k path contributing to the coefficient of the vacuum must cross every edge i->i-1 at least once. Its weight is therefore divisible by their product. At length k there is only the strictly descending path, and before length k there is none. The exact product is k!gamma_k in both parities, including p_0=1. Bandedness makes every finite power coefficient a finite sum.

Consequently the first d scalar wrapper responses are a lower unit-triangular integer transformation of p_k u_k. Every later response is an integer combination of those coordinates. The universal free-fragment factorial theorem and its Newton transformation also put every whole-module functional in this same integer response lattice. Inverting the unit-triangular matrix proves the stronger response compilation statement, with coefficients allowed to be negative and depending on the input parity/order. This is not a positive physical target construction or a dimension bound independent of input order.

The identity boundary is handled correctly. For the singleton, the first positive triangle responses3 and15 have gcd3, and the initial triangle has universal content3, forcing that same divisor on every subsequent response. Omitting the identity without adding another operation would therefore be incorrect.

The strong-bank argument has a different scope from the integer response-lattice statement. Both wrappers are actual strong outputs for every nonempty input. Their scalar relation3H(WT)-H(ET)=4H(T), together with odd Hamilton counts, recovers H(T) for gcd purposes. It does not generally recover H(T) by an integer linear combination of the two strong responses. The source makes only the correct gcd claim in this section. As an independent hostile, the Tr3 and C3 input rows force recovery coefficients3/4 and-1/4 from D and E, respectively. The graph count and maximum-order bounds in the source include d1 and d2 correctly.

Ordered joins multiply complete labelled-bin polynomials because each bin traverses the ordered blocks in their prescribed order. Therefore ordinary fixed divisors multiply exactly under repeated powers of one polynomial. For the explicit input, the first three evaluations have valuation4 and the fourth/fifth have valuation2 with the stated residues. Integer unit-triangular normalized central evaluation gives the first three weighted coordinates valuation at least4m and the fourth exactly2m. The first factorial with a factor of3 occurs at k3, proving the complete ordered-join contents(3^(2m),3^(2m+1),3^(2m+2)). The prime-support theorem rules out other primes for that family.

For the strong wrapper, evaluations1 and2 have valuation at least4m, evaluation3 has exactly2m+1, and evaluation4 has exactly2m. The lowest terms at4 have coefficient20 times2^m modulo3, so cancellation cannot raise the claimed valuation. These facts prove the three asserted3-adic valuations via the same triangular argument. The source properly avoids excluding other primes in the wrapped contents.

The second cyclic readout4Q(3)+5Q(2)-11Q(1) agrees with the reviewed kernel. Substituting D(P^m) yields precisely12P(4)^m+46P(3)^m-3P(2)^m-21P(1)^m. The first term has strictly smaller3-adic valuation than the others for every m>=1. Thus this actual two-stage cyclic construction escapes the single-stage cyclic divisor at every m. It remains a whole-module context, consistent with universal-content monotonicity. The final iterated-context gcd statements distinguish the identity-including language, the strong W/E bank and the positive one-letter exception correctly.

## Independent computations

Neither independent checker imports the author's code or the path-polynomial helper.

`check.cpp` enumerates literal simple odd directed cycles, anchored at their least vertex, and sums disjoint cycle collections recursively. This independently reconstructs all1099 physical bin profiles through order5 and the complete ten-vertex witness polynomial. A separate inclusion-exclusion calculation over arbitrary directed walks counts every tested actual Hamilton output. It checks3222 strong wrapper graphs, both finite-bank gcds and the D/E identity on the complete small input bank. It also directly counts the12- and16-vertex graphs as164673 and196004745, with the latter54 modulo81. Intermediate walk sums use signed128-bit integers; reported counts are range-checked before conversion. Optimized and assertion-enabled debug outputs are identical.

`check_algebra.py` works by evaluations rather than symbolic polynomial expansion. It computes triangle iterates on an integer grid, verifies110 normalized response columns, and checks700 signed reconstruction/gcd controls. It independently compiles700 arbitrary integral factorial-weighted Newton functionals into the finite response bank by back substitution. It then computes central coordinates of the ordered powers and their wrappers directly from normalized positive evaluations for m1 through24, checking all exact contents, all wrapped valuations and every double-cyclic escape valuation. Normal and Python-optimized outputs agree.

The source verifier was copied into `source/` and replayed normally and under Python optimization. Both receipts exactly reproduce the frozen source receipt. Original source and proof bytes still match their assigned hashes.

## Reproduction and exact artifacts

    g++ -O3 -std=c++17 check.cpp -o /tmp/triangle-content-check
    /tmp/triangle-content-check
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 check.cpp -o /tmp/triangle-content-debug
    /tmp/triangle-content-debug
    python3 check_algebra.py
    python3 -O check_algebra.py

No input files or external packages are required by either independent checker.

- check.cpp:daf75b5577e534b6802bca81a51c2b35106457b8fc7957ca985c766f9b354b6c.
- independent-normal.json and independent-debug.json:e7cef670ada97c1bc22e22b936bb98316a45b195eeeee5598c535c323d2c519f.
- check_algebra.py:75a291568ab64790d58f57eac4f015558eecb8e5efa7d3460eb8303093be6d54.
- algebra-normal.json and algebra-optimized.json:c9301296ddb09d84b78bcaf0f1f96a465bc833980e57b053e66dacfd57b8c4a6.
- source-normal.json and source-optimized.json:5924fa0547fdfec0cc1ec4c6968ffc4b9b1673dfbe6e8f91699071c478379746.

All-order scope follows from the proofs above, not the numerical cutoffs. No correction, unresolved blocker or broader prime-coverage claim is present.

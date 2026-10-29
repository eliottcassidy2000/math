# Independent recurrence and harmonic-transform audit — Session105

## Findings

**PASS, scoped to the proofs and finite controls below.** No arithmetic defect was found in `enumeration/integer_response.py`. The full-factorial and odd-local gauges are compatible with the inherited Session98 projective quotient, including exact negative binary shifts after coefficient multiplication. This is an implementation/arithmetic bridge, not a new Burnside formula or a stricter response quotient.

The actual Session87 comparator is `degree_bucket_dp`, which returns the entire prefix, not its slower direct-partition routine. The stored timing experiment does **not** establish uniform improvement over it: at N=150 and N=200 the new packed variants are slower. The README and proof note state this limitation. Timings are empirical, machine-dependent, and compare different output scopes.

A stale local Session98 path was reported and corrected to `session98_growth`, with the published packet path separately identified. No arithmetic correction was required. The harmonic draft's phrase about “negative coefficients of individual paired terms” was corrected after independent review: whole paired-core Taylor coefficients are nonnegative. For n=m+f the subtracted/original coefficient ratio is f*2^(1-c-f)<=1. Only the second rational summand carries a negative sign; this does not affect the pole theorem.

## Independent arithmetic proof review

For a processed odd-cycle history of mass d=N-r, cycle count t, pair-orbit energy e, and centralizer z:

* `V=(N!/r!)2^(e+rt)/z` is integral because `N!/(r!z)` counts partial permutations on N labels. The displayed transition coefficient counts k disjoint p-cycles among r remaining labels. Direct division of consecutive history weights yields exactly the claimed exponent E. Even when E is negative, every resulting history is integral; histories merged in a projective state share the multiplier, so their sum remains divisible. `exact_shift` explicitly rejects a nonzero remainder, including under `-O`.
* `U=odd(d!)2^(e+rt)/z` is integral for processed lengths at least 3. The odd centralizer divides `odd(d!)`; its power of two is absorbed by `e>=sum v2(m_p!)`. The odd-factorial transfer ratio is the odd part of an integral partial-permutation coefficient. The final unit-tail fold divides only after forming the correct numerator.
* Retained deficits are `D_j=sum(j-gcd(a,j))`, with processed a>j. Hence `0<=D_j<=jt<=d<=N`. The specified field width is sufficient, including powers-of-two boundaries. Ordinary packed addition cannot carry between retained fields. Overflow in discarded higher fields cannot travel downward into retained fields. Masks keep exactly odd coordinates 3<=j<=min(p-2,r').
* All odd lengths are traversed, not only primes. Lengths 9,15,21,25,27 are present in the finite partition controls. The optional three-tail is an exact multiplicity sum, not a deletion of three-cycles. Full centralizer factors are retained.

The checker rederives oracle energies by constructing a representative permutation and counting its orbits on actual unordered pairs. It does not use the target's gcd-energy formula for that oracle. It imports target algorithms only for comparison; the arithmetic-history tests independently reconstruct individual history weights.

## Harmonic interpretation and meromorphic transform

The up/down note's commutator follows by separating deletion of the newly attached vertex from deletion of an old vertex; the latter leaves exactly two choices for the forgotten edge. Normalization gives `[D,E]=I`. Local nilpotence of D gives the finite telescoping right inverse and harmonic projector. Thus `h_n=dim ker D_n=a_n-a_(n-1)`, with `a_-1=0`. The decomposition is a decomposition into signed vector subspaces, not a collection of prime tournament types. On `E^(n-j)H_j`, `ED` acts by n-j, establishing the stated refresh spectrum; the separate up/down audit owns exhaustive matrix validation.

Set `b_n=n!a_n/2^binom(n,2)` and `eta_n=n!h_n/2^binom(n,2)`. Exactly

    eta_n=b_n - n 2^(1-n)b_(n-1),
    H(z)=B(z)-z B(z/2)-(z^2/2)B'(z/2).

The identity is first coefficientwise formal, and analytic on |z|<1. It does not treat the divergent raw ordinary generating function of a_n as analytic. The inherited Session98 meromorphic B extends this expression to the plane.

For each inherited core summand `F_mu(z)=C z^m/(1-qz)^(m+1)`, `q=2^(-(m-c))`, its exact paired harmonic summand is

    C z^m/(1-qz)^(m+1)
      - C(m+1)2^(-m) z^(m+1)/(1-qz/2)^(m+2).

Normal convergence of the inherited core sum on pole-free compacta permits differentiation on slightly larger pole-free neighborhoods. The transformed series consequently converges normally away from its two pole ladders; alternatively its meromorphic definition follows directly from the displayed finite differential transform.

At `rho=4^k`, B alone is singular; at `rho=2*4^k`, only the dilated terms are singular. These sets are disjoint. In each bank the unique maximal core `3^k` supplies the top order. If A_k is the inherited positive leading amplitude at 4^k, the leading amplitude at `2*4^k` is `-(3k+1)(2*4^k)A_k`, nonzero. Orders are respectively `3k+1` and `3k+2`: equivalently every `2^r` is a pole of exact order `floor(3r/2)+1`. The empty core already checks both the simple pole at 1 and double pole at 2. Infinite genuine poles exclude D-finiteness of this normalized analytic germ.

This is a direct consequence of the inherited core expansion and the new harmonic interpretation. It supplies neither an independent prime atlas nor a scalar recurrence that computes h_n without class-count input; no literature-priority claim is justified.

## Positive source-free interpretation — additional audit

**PASS.** A universal source is unique, so it is preserved by every isomorphism. Deleting it and adjoining a new universal source give inverse maps between source-containing n-classes and arbitrary (n-1)-classes. Thus h_n also counts source-free classes; n=0 is handled separately. Iterative deletion uniquely strips an initial transitive sequence and leaves a source-free remainder. This positive dimension interpretation does not place individual source-free basis vectors in ker D: nonzero harmonics in positive degree remain signed vectors.

For an odd permutation with c nontrivial cycles and f fixed vertices, any source must be a fixed vertex. For a specified fixed vertex, its incident edge-orbits are exactly one for each nontrivial cycle and each of the other fixed vertices: c+f-1 independent bits. All must point outward. Source events for distinct fixed vertices are disjoint because their common edge cannot point outward from both endpoints. Consequently the source-free fraction is exactly 1-f*2^(1-c-f), or 1 when f=0. This proves the paired-core nonnegativity with a direct fixed-object meaning. It is not a modular-prime classification.

The independent checker now enumerates all labeled tournament orientation masks for every odd permutation type through n=6, tests actual permutation invariance and the absence of a universal source, and compares the fixed source-free count to the claimed filter. Summing these counts with the centralizer weights gives source-free class counts 1,0,0,1,2,8,44, independently matching the differences a_n-a_(n-1).

## Reproduction and finite scope

Run from this audit directory (the packet's enclosing directory name is irrelevant):

    python3 check_independent.py > normal.json
    python3 -O check_independent.py > optimized.json
    cmp normal.json optimized.json

The checker uses only the standard library and explicit exceptions, never removable assertions. Controls include literal pair-orbit Burnside through n=30, all six target variants through n=30, individual gauge transfers and carry-safe packing through n=22, rejection of deliberately nonexact division, the actual Session87 prefix through 200, both packed gauges at 40/60/100/150/200, and actual Session98 values/state rows at 40/60/100. Independent core coefficient extraction checks the normalized harmonic transform through 30, verifies paired-core coefficient nonnegativity, and checks unique top pole cores through k=10. The initial audit verified baseline byte identity against recovered local sources. The portable checker verifies the bundled copies against their recorded SHA-256 pins, without requiring another checkout or network access. It also hashes each directly audited proof note and the arithmetic implementation. The parent's source provenance pins the source main commit; this audit makes no raw/lazy fetch or publication.

Finite checks corroborate but do not replace the all-order arithmetic, Weyl, and analytic arguments.

## Packaging control

The checker uses its packet root (`Path(__file__).resolve().parents[1]`) and no hardcoded enclosing packet name. Normal and optimized runs passed with identical receipts after this change. A separate copy under a randomly named temporary directory, without the recovered sibling source trees, also passed and produced the identical receipt. Required proof files are hashed unconditionally, so missing packaged dependencies fail rather than being silently skipped.

## Inherited source-free basis and compressed Smith atlas

The recovered 2026-10-19 `PADDING_FROBENIUS_PROOF.md`, §5, already proves the integral source-free induced-count basis, all-order deletion cokernels, and Hasse deletion diagonals. The source-free family, deletion torsion, and finite-field rank classification are therefore **inherited**, not Session105 discoveries. The normalized harmonic generating transform is a consequence combining the interpretation with the later Session98 pole theorem. This audit of the compressed algorithm is conditional on that inherited integral classification; the independent small tournament-matrix checks in the up/down lane separately test its application.

For single deletion, the inherited diagonal has n-r with multiplicity h_r for 0<=r<n. A Smith diagonal is obtained by sorting each prime's valuation list independently and aligning the sorted lists. `compressed_smith.py` performs precisely this operation on multiplicity runs: each valuation run has an inclusive end position; the union of all run endpoints partitions the rank interval. On each half-open interval the test `stop>start` picks the correct valuation for every prime. The interval's length is its multiplicity. Equal neighboring factors are merged. No enormous rank-sized array is allocated by the production algorithm.

Each prime p has at most floor(log_p n) internal valuation breakpoints. Thus at most 1+sum_p floor(log_p n) nonempty intervals occur, independent of the astronomical rank. Empty valuation levels only reduce the bound. Valuations are nondecreasing across intervals, hence consecutive invariant factors divide. Sorting preserves each prime's total valuation, so the cokernel-order prime exponents are exactly sum_r h_r v_p(n-r). The number of unit entries modulo p is sum_(r:p does not divide n-r) h_r, giving the finite-field rank.

The independent stdlib checker verifies every degree 1..200: rank multiplicity, positive interval lengths, divisibility, the independently calculated breakpoint bound, exact cokernel-order prime exponents, and finite-field ranks for every prime through 200. It independently trial-factorizes the small inherited diagonal entries, without using the target's valuation or prime routines. For n<=10 it expands and sorts valuations for direct comparison; for n<=8 it additionally performs elementary gcd/lcm diagonal reductions, an independent integer Smith conversion. Large-degree checks never expand the rank. The `compressed_smith.py` file is included in the receipt hashes. This is a compressed implementation of an inherited classification, not a new Smith theorem.

The stored atlas JSON is also checked, not merely the generating function: normal/optimized atlas bytes agree, its count-prefix digest matches the independently compared Session87 prefix, its factors match recomputation, and its prime-order exponents match the independently factorized inherited diagonal. The atlas JSON is included in the audit hashes.

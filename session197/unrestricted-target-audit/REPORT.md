# Independent review: unrestricted logarithmic cores and numeric-target algorithm

**PASS. No correction or blocker.** The proof establishes an unrestricted finite decision-and-witness algorithm with the stated logarithmic reduced witness bound, exact reduced-order floor and bit complexity. It does not prove the full proposed spectrum or unrestricted binary-input NP membership.

## Frozen source and replays

Reviewed source directory: `session197/unrestricted-target/`.

- `PROOF.md`: 5bea17dabd28f39878de94bd59bff984473bed6e93f2d36beafc4f41929defaa
- `algorithm.py`: 6c5f3191c662c5d41887c30f11901be1e1c9babeee1e57010ed1620825cfe0c9
- `verify.py`: ae7b1033887f1c53b3793eeadff281cec85c6009e8b1b2c966d30a67d4b89219
- `VALIDATION.json`: f8cbc03022162424c8abdd55b09ade79f7034820e3d1e26f16f46ff2dd120343

A copied verifier was run normally and under Python `-O`. Both outputs match the frozen receipt byte-for-byte. Original source hashes remain unchanged.

## Strong-order and reduced-order proof

Busch's exact strong minima and their scope n>=3 are used correctly. In particular order6 has minimum15, not17. The weaker inequality H^3>=5^(n-1) alone is insufficient to establish the proposed component inequality at order4, but the proof explicitly treats that case: a mixed fourth vertex can be inserted in a cyclic triangle, giving four cycle cuts, and odd Hamilton-path parity then gives H>=5. Order3 is equality. The comparison 5^(s-1)>=3^s is valid from order5 onward. Thus 3^s<=H^3 holds for every nontrivial strong component.

The SCC factorization is exact because every cross-component arc points forward in the component order. A Hamilton word cannot return to an earlier SCC, so it is a concatenation of component words. Removing singleton components preserves H. If h>1 the retained graph is nonempty, while h=1 is correctly handled by a separate singleton witness. The implementation's empty-graph count0 and empty component-product1 are not conflated.

Multiplying the component inequality gives 3^N<=h^3, with N<=floor(3 log_3 h). The floor placement is correct, and h=5 genuinely refutes 3 floor(log_3 h). The triangle-product examples show sharpness for the size of a reduced graph, not for the smallest possible alternative realization of that count. The proof explicitly makes this distinction.

The reduced-order minimum formula passes review. The exchanges involving a size5 part strictly improve any partition except the isolated order5 case. Three size4 parts can likewise be replaced by four size3 parts. This leaves zero, one or two size4 parts according to the residue, with the isolated exception g(5)=9. No reduced positive orders1 or2 are asserted. Comparing f(s) with g(s) at3 through8 and then using their respective factor5/factor3 recurrences proves the component replacement argument at every order. The three physical seeds attain the required partition costs.

## Algorithm and factorization edge cases

The capped endpoint recurrence is exact under saturation because every update is a sum of nonnegative values. An overflow value h+1 can never be mistaken for a count at most h. Enumeration of every labelled strong graph through S(h) is complete: every nontrivial component count of an h-witness divides h and is at most h. The displayed exact formula f(n) is strictly increasing, justifying the implementation's first-exceedance loop for S(h).

The divisor search is sound and complete. All factors are greater than one; therefore its product states increase strictly, and repeated factors are allowed. Every partial product in an actual factorization divides h. Conversely every returned factor list is realized by the literal ordered join of its stored strong components. The procedure neither requires prime factorization nor assumes realizability of individual primes.

Special controls cover h=1 with an empty bank, invalid and even targets, an empty ordered join, incomplete factor banks, repeated factors, and products with several nontrivial SCCs. A product of zero factors is never returned as an empty-graph h=1 witness. The code's actual returned h=1 witness is `(0,)`.

## Bit-time and workspace

For b=ceil(log_2(h+1)), both enumeration bounds are O(b). Decoding each labelled graph takes polynomially many bit operations, and its capped subset DP uses O(n^2 2^n) additions of O(b)-bit integers. Consequently total bit time is 2^O(b^2)=h^O(log h).

For the strong-factor implementation, at most h representative counts and at most h divisor states are stored. Each representative has O(b) vertices and O(b^2) adjacency bits. Every factor list has length at most log_3 h, and each entry has O(b) bits. These banks, lists, graph counters and the capped DP all fit within 2^O(b) bits. Even a quadratic scan over the divisor states contributes only 2^O(b) arithmetic work and does not change the stated dominant time bound.

A logarithmic-size graph description does not give a polynomial-in-b Hamilton-count verification algorithm. The proof correctly declines to infer unrestricted NP membership or polynomial binary-input time.

## Independent implementation and results

`check_target.py` regenerates all unlabelled tournaments through order7 by one-vertex extension and exact degree-class canonicalization. The complete counts are1,1,2,4,12,56,456. It counts every graph by literal vertex permutations, independently obtains SCCs using Tarjan's algorithm, and checks the reviewed capped DP against those literal counts.

Results, identical under normal Python and `-O`:

- 532 literal tournament counts;
- 525 nonunit singleton-removal checks;
- 2660 capped-DP comparisons;
- strong minima3,5,9,15,25 at orders3 through7;
- reduced minima3,5,9,9,15 at orders3 through7;
- complete target decisions through23 reproduced by recursive division rather than forward products;
- exact reduced partition minima independently checked through order500 using all published strong component sizes;
- 10012 exact bound checks, including targets next to powers of two through256-bit size.

The independent recursive witnesses sometimes differ from the source's shorter witness, as expected. For example it uses two triangles for9 and a3-by5 component product for15. Literal counts confirm both constructions. The all-order completeness of the decisions through23 follows from f(7)=25, not merely from a finite observation.

Independent checker SHA256: bd51d622bbc598732909c748b9f588c344aa323042233bd1fda05809175fd9e1.

Independent receipt `independent-normal.json` SHA256: 37cf8dd5e59d09bf1c28ee2f0f21155ff54720ae0318315b9b85b17e7124fe0a.

This review verifies the elementary corollary and algorithm. It makes no priority claim and no stronger spectrum claim.

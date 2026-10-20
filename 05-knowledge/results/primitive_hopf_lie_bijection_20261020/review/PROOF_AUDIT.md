# Independent proof audit: primitive Hopf/Lie tournament correspondence

Reviewer: circulant-construction worker. Date: 2026-10-20.

**Verdict: PASS**, with the classical free-Lie/PBW and standard Lyndon leading-word facts imported exactly as declared in the packet. No gap was found in Proposition 4.2 or the other reviewed all-order statements. This is a mathematical audit with bounded independent certificate checks, not formal verification or a new full enumeration.

## Reviewed snapshot

Main proof: `PRIMITIVE_HOPF_LIE_BIJECTION.md`

SHA256: `fd3ed52c18a0bbe924f9a151cc5575f3a7621d28e7655fa4d639f89a9a145e12`.

Supporting producer scripts inspected:

- `primitive_walsh.py`: `99779c7b830c695150dda3c802e92111288cef83126faa0bd61cce0b91b86a69`
- `primitive_enumeration.py`: `17eb60b82cfd0c5ca8277c26b7d2cbc2d166701f38bf064577c50759346d741a`

These hashes identify the reviewed version; later changes need a scoped recheck.

## 1. Fourier normalization and coalgebra types

The tournament coproduct counts actual vertex subsets, with multiplicities. It is coassociative and cocommutative. The ordered-join product is associative, and splitting a join independently in its two factors proves bialgebra compatibility. Connected grading supplies the antipode.

The signed graph quotient is correctly typed as orientation lines. Over Q an odd-sign stabilizer kills its line; a sign-even orbit contributes one line. A chosen graph mask fixes its generator only up to sign. This ambiguity does not affect support matchings or the dyadic determinant claim.

The coefficient 2^(-binom(n,2)) in normalized Fourier is necessary and correct. For each fixed vertex split A,B, crossing orientation bits contribute zero if the graph has a crossing edge. Otherwise their sum is 2^(|A||B|), leaving exactly the two internal normalizations. Thus the result is the component-splitting coproduct with its actual subset multiplicities, without an automorphism factor or a missing binomial coefficient.

The graph coproduct must retain the induced orientation signs rather than simply canonicalize induced graph masks. The statement does so, and `component_coproduct_walsh` deliberately uses the actual induced masks. This is an important correctly handled implementation boundary.

The component coalgebra uses ordinary unscaled monomials: k identical components split with binomial multiplicities. It is not the divided-power normalization. Over characteristic zero its primitive subspace is exactly the connected orientation lines. Consequently the connected-to-primitive restriction of Fourier follows genuinely from coalgebra structure.

## 2. Eulerian idempotent and free-Lie generators

The proof of primitiveness via convolution logarithms is valid. In Hom(H,H⊗H), the maps α(x)=x⊗1 and β(x)=1⊗x commute under convolution because the coproduct is cocommutative, and their product is Δ. Since Δ is multiplicative, it commutes with taking the convolution logarithm. Therefore Δe(x)=e(x)⊗1+1⊗e(x).

For a strong generator s, the k=1 term of the ordered-partition formula is s, and every other term has at least two nonempty ordered-join factors of smaller vertex weights. Replacing every factor by its SCC word preserves that strict increase in SCC-word length. The recursive generator change is therefore invertible and unitriangular on each finite vertex-weight space.

The SCC decomposition identifies H as the free associative algebra on strong types. Replacing them by the primitive generators e(s) identifies its coproduct with the standard primitive-generator tensor-Hopf coproduct. The imported free-Lie/PBW theorem then gives its full primitive subspace. There is no unjustified identification of a raw strong tournament with a primitive element; the triangle subtraction explicitly detects that error.

## 3. Proposition 4.2: why no odd factorial enters the determinant

The determinant-only result is sound. The lattice distinction is important:

1. **Full Fourier over R=Z[1/2].** The labelled Walsh transform is an R-isomorphism, and taking coinvariants preserves an isomorphism. Signed graph orbit coinvariants are free R-modules on the even orientation lines because their only possible stabilizer relation is multiplication by -1, whose 2-torsion vanishes over R. Tournament coinvariants are free on types. Thus the full degree-n Fourier determinant is an R-unit.
2. **Rational triangular changes.** The raw-strong-to-Eulerian generator change has determinant one, even though its entries may have odd denominators. The standard Lyndon bracket has its Lyndon word as leading word with coefficient one. Unique nonincreasing Lyndon factorization makes the ordered PBW monomials unitriangular relative to ordinary words. Their determinant is likewise one, up to harmless basis-order signs. Neither statement asserts that these changes preserve the R-lattice entrywise.
3. **Coradical filtration.** Component count on the graph side and PBW degree on the tournament side are the respective coradical filtrations. Fourier preserves that filtration. On its k-th associated graded piece, apply the k-fold reduced coproduct. For an unscaled monomial it gives the sum of all k! orderings of its primitive factors, including repetition multiplicities. Over Q that symmetrization is injective. Both sides use the same normalization, so the resulting top map is exactly Sym^k of the primitive map. There is no residual k! factor in the diagonal block.
4. **Weight decomposition and induction.** At vertex weight n, the degree-one filtration block is K_n. Every other diagonal block is a direct sum of tensor products of Sym^(r_j)(K_j), with j<n and sum j r_j=n. For a d-dimensional matrix A, det Sym^r(A)=(det A)^binom(d+r−1,r−1); tensor products similarly contribute integral powers of determinants. Therefore the full determinant is det K_n multiplied by integral powers of lower det K_j. Starting at weight one gives det K_n=±2^z for an integer z.

Thus the proof does not confuse a rational PBW basis with an integral primitive lattice. The individual 3-denominators remain fully compatible with a dyadic determinant. For example, a matrix diag(1/3,3) has determinant one but no support matching consisting entirely of 3-adic units; determinant information alone cannot restore an entrywise p-local assertion. The packet correctly declines that inference.

The saved matrices were independently evaluated by recursive Laplace expansion, rather than the producer's Gaussian elimination. Their determinants at weights 1–5 are 1,1,−1/4,1/64,1/2199023255552, agreeing with the stated powers of two. Nonzero support and matching distinctness were also checked.

## 4. The support-certified type bijection

An invertible finite primitive matrix has a nonzero permutation term and therefore a perfect support matching. The ordered feasibility search specifies one matching unambiguously. Composition with the inherited SCC word, nonincreasing Lyndon factorization, and graph component multiset gives mutually inverse type maps. The same declared alphabet order must be retained in both directions.

This is an effective global matching construction, not a direct local atom rule, and not a labelled equivariant bijection. The packet makes these limits explicit. Orientation-line sign changes multiply columns by -1 and do not alter the support condition.

## 5. Converse, edge parity, and signed enumeration

Converse is an algebra anti-automorphism of ordered join. Reversing the block order in the Eulerian formula proves that it commutes with e. On a bracket of letter length r it therefore acts as (-1)^(r−1) times the Lie automorphism that converses the strong letters. On Fourier images, arc reversal multiplies the graph character by (-1)^|E|. This establishes the asserted primitive eigenspace dimensions.

The SCC palindrome formula R(z)=(1+Q(z))/(1−S(z²)) is correct: outer strong letters occur in converse pairs, and the center is absent or self-converse strong. The signed graph component product has factors (1−z^n)^(-a_n)(1+z^n)^(-b_n). Its ghost term at a divisor k is k d_k for an odd repetition and k c_k for an even repetition, exactly as in (5.4).

The length-refined Witt character formula correctly uses trace(theta^j on strong letters), hence Q for odd j and S for even j. The additional factor (-1)^(r−1) belongs to actual converse, not theta itself.

A useful manual weight-six check is: the primitive basis has 35 strong letters, six brackets of the singleton with weight-five strong letters, one length-three bracket of two singletons with the strong weight-four letter, and one length-four bracket of three singletons with the triangle. This gives dimension 43. Their converse traces are respectively 9,−6,+1,−1, totaling 3; hence the parity dimensions are 23 and 20. Ignoring converse-paired strong letters would indeed fail here.

The Burnside restrictions for self-converse tournaments are correct. A vertex cycle of odd length greater than one produces a forbidden odd-length pair orbit without endpoint exchange; an even vertex cycle has an antipodal orbit that forces its length to be 2 modulo 4. Two fixed vertices also fail. Conversely cycles of length 2 modulo 4, with at most one fixed vertex, satisfy every pair-orbit condition. Each pair orbit contributes two choices. The ordinary tournament restriction to odd vertex cycles is also correct.

The independent oriented-graph calculation uses all permutations and the three pair states absent/forward/backward. An inconsistent pair orbit permits only the absent state, whereas a consistent one permits three states. The implementation's 3^good-orbits formula and converse parity modification are therefore correct.

Finally, fixed multisets under an involution include equal multiplicities of each nonfixed component pair. This proves the extra factor at doubled weights and the correction from 727 inverse-Euler atoms to 726 connected self-converse oriented graphs at six vertices. The out-star/in-star pair on three vertices is the correct first witness.

## 6. Bounded independent checks and limits

`check_saved_certificates.py` imports neither producer script. It performs direct Laplace determinant calculations on the saved primitive matrices, verifies their support matchings and odd-denominator hostile, and independently multiplies the unsigned/signed Euler factors through order ten using integer binomial coefficients. All **40 checks pass**; see `certificate_checks.out`.

Reproduction:

    python3 check_saved_certificates.py PATH_TO_SPECIES_PACKET

The full order-sixty enumeration and the graph census were not rerun for this review. Their scripts and saved reports were inspected, and the all-order identities were audited above. No additional hypotheses, hidden factorial normalization, or unsupported integral-lattice assertion was found.

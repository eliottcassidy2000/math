# Independent audit: universal whole-module Hamilton content

**PASS, with no correction.** The gcd over all actual whole-module tournament contexts equals gcd_r (r−1)! r! c_r(T). The shortened physical test bank, central-factorial formula and sharp valuation bound are sound. This is an exact divisor invariant, not an attained-count or prescribed-multiple theorem.

## Frozen source binding

The unchanged source packet is `session197-odd-cycle-path-polynomials/`:

- UNIVERSAL_MODULE_CONTENT.md: `4f86f7ada5eff79e61bb2ef4b8139b76df8efedf6996a71dae327bc66b24d9cb`.
- verify_universal_module_content.py: `c9779056d1a4269d8fa6146a4dc37735c69f6f21b8e379b884794193921e5b78`.
- Helper verify.py: `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`.
- Convention proof PROOF.md: `b88c86ac0afa64d3a1bf399c56fa7f8aad046fb93424c70a9ace8d164708183b`.
- Both supplied receipts and both copied source replays: `775e4a4ae7ce2b2511e42df14b1f0b3f86a91864b9c041902a282e66bbc78e0f`.

The source was copied before execution. No author file was changed.

## All-context factorial and matching-bank proof

A Hamilton word in Q[q:T] splits uniquely into its maximal nonempty T-runs. Their ordering gives r! c_r(T). Every internal exterior fragment starts in O and ends in I and is nonempty. Each contains its own nonempty vertex set, disjoint from all others. Permuting those r−1 literal fragments gives distinct valid exterior words, with the two outer fragments fixed. There is no stabilizer from matching endpoint types or isomorphic paths. Thus (r−1)! divides each actual exterior coefficient for every context, without assuming a formal operator is physical.

The matching bank has exactly 2r−2 exterior vertices. An internal O-to-I fragment uses at least two, so no more than r−1 fragments are possible. At equality each is a single edge, and the only permitted O-to-I edges are the r−1 matching edges. All vertices are used internally, so both outer fragments are empty. There are exactly (r−1)! orders. This gives an integer lower unit-triangular coefficient matrix on w_r=(r−1)!r!c_r. Its integer inverse proves the reverse gcd divisibility using actual contexts. Identity r1 and input order n1 are included.

## Central Newton coordinates and valuation scope

Both parity bases are monic integral and triangular on their parity monomials. Consequently the coordinates u_k are integers, including signed inputs. The normalized q_k and r_k binomial identities give Newton coefficients N_(r,k) with zeros for r<k+1 and unit diagonal at r=k+1. After weighting Newton row r by (r−1)!, the coefficient on k! gamma_k u_k is integral because every nonzero entry has r−1>=k. The first ceil(n/2) rows are therefore again integer unit triangular; all later rows are integer combinations. This proves both central gcd formulas and the shortened physical matching bank, including the even n2 boundary.

For each prime p choose a nonzero gamma_k u_k attaining v_p(g). Since gamma_k is divisible by (k!)² and u_k is integral, v_p(k!)<=floor(v_p(g)/2). Hence v_p(g)<=v_p(G)<=v_p(g)+floor(v_p(g)/2), proving G²|g³. No positivity assumption is used in this algebraic step. It implies equality of prime supports and of prime-square divisibility. The alternative Taylor proof is valid for integer-coefficient parity polynomials, including p2.

The prime/7p scope is correct. If g divides a squarefree target, the valuation bound forces G=g. For target49, g=1 or7 is immediate; when g=49, the identity context gives G|H(T), and H(T)<=49 forces G=49. Thus the refined gcd adds no necessary prime/7p filter once the stated ordinary-content and induced-count bounds already hold. This does not assert realizability of any target.

Composition is also correctly scoped: a previously substituted intact module remains a module after a further whole-module substitution. Every later context is thus one of the original contexts, so G(T) divides G(S). Arbitrary nonuniform births or rebuilding do not inherit this argument.

## Independent numerical checks

The standalone C++ checker imports no author code. It counts ordered path covers by literal permutations and forced cuts, rather than subset path-cover recursion. Exterior coefficients are obtained by literal exterior permutations and every possible set of root gaps. Physical substitution places the input module first, unlike the source layout. A separate walk inclusion-exclusion counter counts Hamilton paths in each literal graph without endpoint subset DP.

Optimized and assertion-enabled unoptimized runs agree on:

- 5,405 labelled rooted contexts through order5;
- 15,925 literal exterior coefficients and their factorial divisibility;
- 855 actual transitive-clone counts by walk inclusion-exclusion;
- 1,099 literal ordered-cover input vectors and their full/shortened bank gcds;
- 147 actual matching-bank substitutions;
- the strict eight-vertex input's full ordered-cover vector;
- the strict input's first four actual matching-bank counts, 459,7155,127467,2621187, whose gcd is27.

The last graphs have orders8,10,12,14. Their ordered-cover vector is exactly (459,5778,31860,99000,187200,216000,141120,40320), giving g9 and G27. Every context contains the input induced, so the gcd27 is not itself attained: the monotone path-insertion injection forces H>=459.

The separate standalone Python checker uses direct signed monomial polynomials, exact binomial inversion for Newton coefficients, and independent central-basis elimination. Normal and optimized runs agree on34 basis columns,629 exact Newton entries and2,200 signed polynomial controls. Every case passes both gcd formulas, the shortened Newton bank, prime support/square tests and G²|g³. The strict example has central coordinates (459,405,48,1), with equality G²=g³=729.

The finite controls corroborate the bijections and algebra; the all-order conclusions come from the preceding arguments, not from a census cutoff.

## Reproduction and artifact hashes

Compile `check.cpp` using C++17 and run it. Run `python3 check_polynomial.py` and its `-O` form. Neither independent checker writes files. The copied source checker needs its adjacent copied `verify.py`; it also writes no files.

- check.cpp: `31c7b772a3d6fa66901dbf4c4d2a0b49a4d44955394ed8cdabe8e39551f029d2`.
- independent-normal.json: `d9a205f0516c6e82610a66c84470dd6105c787b867f31b0d0c14400d662b3c8d`.
- independent-debug.json: `d9a205f0516c6e82610a66c84470dd6105c787b867f31b0d0c14400d662b3c8d`.
- check_polynomial.py: `25a615bf9c7c6dc379e15330e7afcdea5e7af381be5c3af087ccf3560fa3b39b`.
- polynomial-normal.json: `bf4ee40609177aabc7650036c6e3aa42d683a07b168bc97cc23f35cba821867f`.
- polynomial-optimized.json: `bf4ee40609177aabc7650036c6e3aa42d683a07b168bc97cc23f35cba821867f`.

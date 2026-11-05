# Independent audit: response cone, exact arithmetic inverse, and coverage through 10,000,001

**PASS. No correction requested.** The source snapshot was copied before replay and was not altered. The all-order statements are the response cone and exact inversion for a fixed physically reachable prefix. The scalar-coverage statement is finite, through 10,000,001. It is not a cofinite theorem and does not settle the unrestricted {7,21} spectrum.

## Frozen source

Under `session197-lowcycle/`:

- `ARITHMETIC_INVERSE.md`: `4ac1727eb27c85f8a39f5598eb9f0e1b020bc4020cbca0a94e6288cb0676c3b7`.
- `verify_arithmetic_inverse.py`: `a49eca8dea94ce1f177d23ede6abc6fb6dd923a62eb65876268591132ca8a4f0`.
- `arithmetic_coverage.cpp`: `f202e22f53eb11d4f61a7b00b55d0c325325cf0eac5d38e0bd1c0f35d01aa950`.
- `arithmetic-c2-witnesses.json`: `c02dcf59ed38ecda740af73fb2b5755b6ce12c63e2a50a11a4a0e586a380c32c`.
- `arithmetic-hard-witnesses.json`: `04843077aefa963d83b5d400b2189ff47dc596b05d77d3f55a11135e64367707`.
- `arithmetic-inverse-certificate.json`: `d75bc6e2dfdd8e2c18f76cdb29d5e6f535f6e9aa46c1facdeece10a3839dd7c7`.

Copied source normal and Python `-O` runs reproduce the entire frozen certificate byte for byte. Final readback of all six original hashes is unchanged. The certificate also binds the inherited ROOTED/C1 completeness, integer-gadget compiler, and physical two-port compiler/core dependencies. I read the displayed one-seam law and the two-port compiler's physical arc checks; these are inherited results, not newly inferred from the coverage experiment.

## Mathematical review

For an input (h,e,R), the proof correctly uses h>=1, R>=h, and 0<=e<=h. From the gadget recurrences, Q>=g, B>=L, B<=2L-2 for L>=2, and d=2g-Q>=1, with d>=3 for odd L>=3. The latter induction handles even ell and odd ell>=3 separately, with ell=2 as the first odd-output case. Substitution gives both sides of h'<=R'<=2(h'-e'), including L=1 and all endpoint boundaries e=0,h.

The singleton seed (1,1,1) is explicitly excluded from the output cone. Its L=1 output (1,0,2), although not strong, is a nontrivial grammar state and does satisfy the cone. Thus the subsequent 3R/2 growth bound, finite preterminal response bound, and logarithmic number of later steps have no hidden seed exception. The formal inverse hostile (5,1,7) -> (3,2,4) correctly shows that cone arithmetic is not physical inverse closure.

The final-gadget inverse correctly splits L=1 from L=2k+t>=2. Divisibility by gcd(2R,h), residue modulo h/gcd(2R,h), the modulus-one case, parity, and both interval bounds are correct. The final equality test makes this an exact fixed-prefix inverse, rather than a sufficient congruence criterion. The lower interval bound follows from g<=L; the upper follows from g>=1. No prime or odd-h assumption is needed for the arithmetic algorithm itself.

The prefix enumeration is complete within the stated response cap: R'>=h+R is strictly increasing, and L<=B/R is a necessary bound for an output with response at most B. Deduplication is by the full (h,e,R) triple. The finite one-final-gadget image does not certify absence from the unrestricted grammar. The proof explicitly obtains the 29 small C1 exclusions from its separately complete bounded C1 theorem.

## Independent full image regeneration

`independent_coverage.cpp` does not use the source's recursive gadget constants. It reads the literal binary apex pattern and counts possible Hamilton-word left/right base assignments by a closed finite sum. It also generates the entire prefix bank in increasing R buckets rather than the source's growing-list traversal. Its packed image uses 64-bit words and explicit byte serialization.

For completeness of the independent gadget formula: the transitive base vertices on each side of the apex must appear increasingly. An empty right side contributes one. An empty left side contributes the first orientation bit. If the first right index j is below the last left index i, the number of assignments is 2^(i-j-1); if j>i, a nonempty valid split must be adjacent and gives a 0-to-1 bit transition. Summing gives g. The starting two-cover Q and B follow directly from the same pattern.

Both independent optimized and debug full runs reproduce these exact results:

- Response cap 500: 2,717 prefixes and 165,957,887 terminal trials. Packed image SHA256 `bf321764ea259df608e7fafa9d0a9418ad6230925dab3a952221338fd0342f4e`.
- Response cap 1000: 7,497 prefixes and 230,771,132 terminal trials. Packed image SHA256 `39d97d495f7f7a847eec9181642b3e3cbcdcec8b7fd1b1eb8fac6fcf2db4a9cb`.

These match the source hashes for all entries through 10,000,001, not just the number of gaps. The cap-1000 missing odd set is exactly

    7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,
    277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647.

The cap-500 image has precisely four additional gaps: 5,099,741; 5,477,461; 8,231,473; 9,551,177. No claim about larger targets follows from either image.

## Separate arithmetic and physical controls

`check_inverse.py` enumerates all literal left/right base assignments for gadgets L=1,...,2048: 2,796,203 assignments. It checks their g,Q,B,L against both the independent formula and source recurrence. It then checks 599,495 exact modular-inverse instances over 1,495 virtual profiles, including even h, even targets, e=0,h and modulus one; 42,564 output-cone evaluations also pass. Finally it compares every full prefix triple for both caps, rather than just the bank sizes. Normal and Python `-O` outputs agree.

`check_hard_witnesses.py` independently reconstructs all four large graph recipes, verifies every stored profile trace and prefix state, and compares the stored mask/drop word with a fresh physical compilation. Literal replay checks 1,301 unordered arc pairs against the displayed graphs. The inherited, separately audited uncapped 27-coordinate transfer then gives the target H and exactly one Hamilton cycle for each graph:

- H=5,099,741 at order 26, final (h,e,R)=(5,099,741,799,9,281,553).
- H=5,477,461 at order 25, final (h,e,R)=(5,477,461,689,5,524,209).
- H=8,231,473 at order 26, final (h,e,R)=(8,231,473,903,13,847,505).
- H=9,551,177 at order 27, final (h,e,R)=(9,551,177,853,16,249,891).

This is not a claim of an independent order-27 subset-DP count. The graph construction, stored-word arc replay and direct binary constants are independent controls; the final interface evaluation uses the named audited transfer. Normal and Python `-O` outputs agree.

All 26 C2 repair graphs have also passed the separate independent review `session197-hybrid-constructions/audits/C2_REPAIRS_REVIEW.md`, SHA256 `c3d19c035c357d508c2bb9a14188bf8bb69e1f52c552c123d628499bc4a27aa5`. That review reconstructs each graph and checks Hamilton paths/cycles by both subset endpoint DP and a separately written walk inclusion-exclusion counter, plus physical two-port compilation. It binds exactly the recipe file above. The maximum repair order is 13.

Consequently the positive conclusion is supported: every odd N from 3 through 10,000,001 other than 7,21,67 has an actual strong tournament with at most two Hamilton cycles. Absence of those three from the full C<=2 class is a separate inherited complete two-port result, not proved by this finite image. The 29 C1 omissions also use the separately complete bounded C1 result. No extension to all odd N, arbitrary larger prefix caps, or cofinite coverage is made.

## Independent artifact hashes

- `independent_coverage.cpp`: `d6a28eddf404ac76aae69e2a9c11f01109c79189eef0239473ee45c67b7c9cda`.
- `independent-normal.json` and `independent-debug.json`: `ef182db822ff99c3a20cd41ce69da1dc5cbd581c369ea56f0e409fc9dc756c54`.
- `check_inverse.py`: `e4d303863d1d126db4fa51dbde41973096a54d327f796d1cc2d1a41081c938bd`.
- `inverse-normal.json` and `inverse-optimized.json`: `f0d2ffbf97c50393edbbc47fefdc2ecf98288e628d8402a3661d897e4268f508`.
- `check_hard_witnesses.py`: `1e218fbe5c5a99bb383cffc6954d666a8e8fce5e255ab3d365411afd4ff0a478`.
- `hard-normal.json` and `hard-optimized.json`: `57848ff62b9f6b4116ca65f2973e2b78d2280288f5b5382e2323bd4dc9186adb`.

Reproduction: compile the C++ file with C++17 and run with an output-prefix argument and bound 10000001. Run both Python independent checks in this directory, with the copied source tree retained. The source replay also succeeds in that copied tree. This report records no unresolved blocker or mathematical correction.

# Independent review: arbitrary-module unit-offset rigidity

**PASS. No correction or blocker.** The arbitrary-module coefficient identity, the exact B(2)-3H unit characterization, and the resulting pure-unit classification are valid within the specified seven-vertex quotient.

## Frozen source

Reviewed in `session197-odd-cycle-path-polynomials/`:

- `UNIT_OFFSET_RIGIDITY.md`: 3c116086bfd697eaeaf5fcbf46308f1e31cdf653eb2a651053f01bce46ce374b
- `verify_unit_offset_rigidity.py`: 756ea2c09f32d4830a05ebf92bee079812b3b013d9526b173a04ca2f70bbab0f
- imported helper `verify.py`: 10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb
- `unit_offset_normal.stdout`: 566f5064b231081d23adb97d4461ea87b0d1a13957a346a9f679fc80e0852127

The named mathematical dependencies were read: full-profile linearity in `UNIVERSAL_CONTEXTS.md` (8615845e7bf06a6cdb1925fdbd8198a77b1a8864aa9c22704b0459678ff2ed51), and the uniform transitive tensor identity in `session197-hybrid-constructions/positive-clock-family/PROOF.md` (be3c701188a4a354cee62f21b6ee76bf172f67351028d88a872f4c95c9b865c3).

Copied source normal and Python-optimized runs match the frozen receipt byte-for-byte. Original source files remain unchanged.

## Whole-module extension

The required linearity is genuinely on the full bin polynomial, not merely the scalar Hamilton count. Cutting a substituted module into maximal runs gives an unordered path cover. For any fixed number of runs, their possible arrangements with the exterior vertices depend only on that number: every run has the same exterior incidences, and the actual runs are distinct because their vertex sets are disjoint. This supplies a linear functional of the falling-factorial cover coordinates.

For a nonempty transitive module of size n, the bin polynomial is x^n. These monomials, for n>=1, span precisely the polynomials with zero constant term. Every nonempty tournament lies in that space. Thus the transitive identity can be extended successively in each of the five modules, keeping the others fixed. No empty-module specialization is required or inferred.

Evaluation at bases1,2,3 becomes H, B(2), B(3), respectively. In particular a positive-multiplicity term independent of one module's size is its evaluation at1, not the scalar1 on arbitrary inputs. Expanding the inherited formula with this rule gives exactly the stated A,B,C, including all h_U,h_V,h_P,h_R factors. The displayed singleton crossing covers in the inherited proof are valid for the quotient and bound the run degrees uniformly; the extension is not based on an unjustified finite interpolation.

The universal-context operator's stronger central-factorial positivity theorem is not needed here. This audit uses and verifies its whole-module linearity argument only.

## Boundary statistic and unit classification

The author's cut proof is sound: each Hamilton word has n-1 nonempty cuts and two assignments to labelled bins, while a fixed ordered two-path cover has at most two concatenation preimages. Empty-bin contributions are exactly2H. The resulting B(2)>=(n+1)H is valid at n=1 as well.

There is also a direct injective proof. Cut each Hamilton word into d consecutive pieces, allowing empty pieces, and place them in the fixed bin-label order. Concatenating bins recovers both the word and all piece lengths. Consequently, for every integer d>=1,

    B_T(d) >= binom(n+d-1,d-1) H(T).

At d=2 this is the same bound. This independent injection does not replace or invalidate the source's correct two-order argument.

The small-order exception handling is complete. At n>=4 the statistic B(2)-3H is at least2. At n=3 a putative value+1 forces H=1, hence transitivity, but then the statistic is5. A negative value is excluded there. At n=2 the unique tournament has statistic+1, and at n=1 it has statistic-1. The implication H=1 implies transitive follows either from deterministic insertion of triangle Hamilton paths or from the nontrivial-SCC product bound.

Therefore a unit constant C forces h_U=h_V=1 and both boundary factors to be units. U,V must be transitive, while P,R are each K1 or Tr2. Direct substitution gives:

- both K1: C=1 but B=2(b_U+b_V)>0;
- one K1 and one Tr2: C=-1 but B=6b_U or6b_V>0;
- both Tr2: C=1, B=0, A=18*2^(u+v).

Positive exponential amplitude makes the eventual-to-exact argument valid. Thus the stated positive rays72*2^k*3^m+1 are exhaustive in this quotient, and no negative-unit ray or amplitude180 positive-unit ray occurs. No symmetry assumption, transitive-only restriction on the initially arbitrary modules, or claim about other quotients has been inserted into the conclusion.

## Independent checks

`check_modules.cpp` is independent of the author's code. It computes subset Hamilton counts by reverse endpoint recurrence and bin evaluations by explicit subset convolution. It constructs the seven-block graph directly.

Optimized and debug outputs agree on:

- all33,867 labelled tournaments through order6 for the cut inequality and exact +/-1 boundary inventory;
-39 distinct physical (order,H,B(2)) profiles and all2,313,441 four-profile coefficient combinations;
-144 unit-constant cases, of which36 have pure-unit cancellation, exactly as classified;
-363 actual mixed-module graphs, including all32 K1/C3 choices and a broader four-module-type bank with nontransitive strong four-vertex inputs.

An explicit hostile to dropping the absent-variable H factors uses T=P=R=V=K1 and U=C3. Its actual count is219, agreeing with A=48,B=36,C=3. Treating absent-module factors as1 would give201 instead.

`check_cut_fibres.py` separately enumerates every literal ordered two-path cover and every cut incidence on all75 labelled tournaments through order4. It verifies the maximum fibre is2, then directly enumerates labelled bin assignments for d=1,...,5 and checks375 cases of the general injective bound. Normal and Python-optimized outputs agree.

Independent C++ source SHA256: 12acf05e032a2381390f65b001b1d9c71d04d31c581ca63f2a800f4df9f707c2.
Independent C++ receipt SHA256: cc6ec0fbd7876c3f74f4456a6c20d3337d8a0fc13138873b82dba1710e8f69cb.
Literal fibre checker SHA256: 0d7e59c28ce837bd3504e65dbd8ee518ac29f3739787a6580cf13387f138470b.
Literal fibre receipt SHA256: 6bebfe9f46af3f6ebeff37c589b2866cdafa28f34b8eb39555034482de4eee29.

These computations support the proof and explicitly test the difficult module factors. The all-order conclusion rests on the full-profile linear extension and the boundary lemma, not on a finite search.

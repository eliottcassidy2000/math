# Independent review: transitive-triple circulation

Verdict: PASS. No mathematical correction is required. This review binds the final proof SHA256 `df1900581da0ac6755884fd45f64f13d049ab903a0202062e69294f3911807ec`, verifier `4db45f5d111c47fd655852035a9977d8ba949ce368b5aabd3c4e037dd772dc43`, and both producer receipts `318ef31b357ff02a70fa593e459e05e1d9fae391b2aee652e1c85e0cc9daf601` in `../tt3_circulation/`.

## All-order proof review

The squarefree Schur elimination retains physical names and disjoint supports. The effective matrix is the finite Neumann expansion for `B_A + tau J`; local converse transposes it. The complete eight-term disjoint quartic contracts to the four terms in equation (2). The omitted products have at least four A labels and therefore vanish on a three-name module. Independently multiplying this matrix gives all four rows of the displayed coefficient table, including the signs and the `2 tau` terms. Consequently the signed TT3 circulation is exactly the virtual-cycle quartic for arbitrary exterior data, not merely in the finite examples.

The auxiliary `formal.py` verifies this contraction as an exact polynomial identity with thirteen independent commuting exterior symbols and three squarefree A labels. It starts from the complete eight-term quartic, rather than assuming equation (2), and verifies the circulation after full-support extraction. This supplements the direct mathematical review and avoids treating a bounded exterior census as the universal proof.

Grouping X and Y is alias-free: the two distinct marks at each port already saturate the at-most-(2,2) degree bound. The six internal masks and four outer finite-difference corners give exactly the stated 24 actual masks. The strict lower bound follows from the separately proved cyclic ordered-module positivity theorem. The virtual cycle is used only to evaluate the signed coefficient; the actual source remains the named TT3.

For quotient converse, the cited all-cover top-degree identity and exact non-A profile restoration apply. In a fixed ambient, both ports imitate the core's incidences, the augmented core is a common module, and the ordered witness remains a module of the entire exterior after A is deleted. The positive-kernel argument applies with the exterior enlarged; this is not an unjustified generic transfer of a core sign. The proof correctly requires the target distinguished TT3 to retain its physical internal order. Arbitrary original-module families inside that TT3 are not asserted to be covered.

## Independent executable checks

No producer algorithm is imported. `audit.py` reconstructs every retained actual augmented graph and counts Hamilton paths with a separate predecessor-removal dynamic program in `physical.cpp`.

- All 43,434 distinct actual graph counts agree, through order 14.
- All 448 arbitrary-exterior controls agree with a separately enumerated exterior squarefree path matrix and the four exact quartic formulas.
- All 159 positive witness controls and four strong-ambient controls satisfy independently computed positive lower bounds.
- All 841 named polynomial coordinates, including every zero, are recovered by physical finite differences. The contrast has 216 nonzero cubic and 248 nonzero quartic coordinates. All 32 retained full-mask reconstructions agree.
- The independent free-exterior-symbol contraction and circulation both pass.
- Normal and optimized Python modes give byte-identical receipts for both audit programs. Checks use explicit exceptions, not removable assertions.

The source's final wording correction to “bidegree at most (2,2)” is reflected in the pinned proof and source-bound receipts. This review does not infer unrestricted reconstruction, arbitrary prime-quotient closure, or additional family scope from these checks.

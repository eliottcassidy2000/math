# Independent audit: scalar unit-start spectrum through 50,000

Status: PASS, finite exact scope, conditional on the inherited complete physical source grammar.

The Python probe `session205-scalar-unit/probe50.py` was rerun. Its frontier exhausts after order 21, before the hard loop limit. An independently written C++ enumerator uses the nonnegative coordinates (a,f,e,p), their separate transfer formulas, and a `while (!cur.empty())` loop without an order cutoff. It reproduces every frontier cardinality, processes 4,082,826 level-states, and finds the same 29 missing odd values. The largest missing value is 2,647; every odd value from 2,649 through 49,999 is attained. The even upper endpoint 50,000 is not a claimed spectrum value.

## Why the pruning is complete

Use physical responses (h,R,e,d)=(H(U), rooted designated two-path-cover count, root-end Hamilton count, H(U-root)). The rooted edge starts at (1,2,0,1). The audited inherited O/Z/J transfers determine these four coordinates exactly, regardless of omitted higher cover coordinates. Hence merging equal four-coordinate states preserves every scalar continuation. No graph isomorphism or reversal quotient is taken.

On a physical nontrivial unit-start state, R>=d>=e, in fact R>d, so Z and J raise h while O preserves h. A branch with a Z/J child above M cannot later return below M. After at least one O, the response becomes (h,2^k R,0,d); its next Z or J has count h+2^k R. Thus if h+2R>M, all continuations beginning with O either remain at the already recorded h forever, or exceed M at their first non-O letter. Dropping that O child loses no scalar value at or below M. This justifies both pruning tests without a guessed order bound.

All generated states are physical because they start at the literal rooted edge and apply literal operations. The singleton has scalar h=1, already included. It is not used as an input to Z, so the singleton Z/start exception never arises. Completeness for all rooted unit-start sources, rather than only their generated recipes, inherits the independently audited complete grammar in `session203-rooted-cover-release/session203/rooted-cover-closure/dependencies/PHYSICAL_SOURCE_GRAMMAR.md`, Theorem 1. The four-coordinate formulas themselves are also given in its Section 2.

## Independent controls

`independent.cpp` uses the positive-coordinate formulas

- O(a,f,e,p)=(a+e,a+2f+e+2p,0,e+p),
- Z(a,f,e,p)=(2a+f+e+p,f,e+p,0),
- J(a,f,e,p)=(f+p,0,a+e+p,e),

and readout h=a+2e+p. Its output is `INDEPENDENT.out`. All maintained coordinates are nonnegative and bounded by the finite target pruning; unsigned 64-bit arithmetic has ample margin for these M=50,000 states.

`literal.py` separately builds all 364 recipes of length at most five (through order seven), counts Hamilton paths by subset dynamic programming, computes the rooted two-path-cover response by a subset convolution, and verifies the four coordinates and unit-start property. `LITERAL.json` records PASS. Validation uses explicit exceptions, not removable assertions.

Reproduce:

    g++ -O2 -std=c++17 independent.cpp -o independent
    ./independent > INDEPENDENT.out
    python3 literal.py > LITERAL.json

## Required production guard and exact limitation

Replace the exploratory `range(2,24)` with an exhaustion loop, or explicitly reject a nonempty frontier when the limit is reached. The present run actually exhausts at order21, so this issue does not invalidate its stated result. It would be unsafe to reuse the hard cutoff at a larger target without a new completeness argument.

The missing values are:

7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647.

These are absences from the unit-start subclass within the checked finite scalar range, not gaps in the unrestricted tournament spectrum. Exhaustion certifies the finite range across all grammar word lengths. It does not imply any value above50,000 is attained and does not prove eventual coverage. An all-order invariant or constructive induction remains missing.

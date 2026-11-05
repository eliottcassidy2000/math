# An all-order weak-parent obstruction at 41 and 67

Session 194, 5 November 2026. **FINITE-CERTIFICATE RESULT; independently reviewed, PASS.** This packet is separate from the frozen sharp weak-deletion ratio theorem. The unrestricted `{7,21}` Hamiltonian-path spectrum remains OPEN.

## 1. Claim and exact scope

There is no strong tournament T with a vertex v such that S=T-v is weak and either

    H(S)=41, H(T) in {125,127},
    H(S)=67, H(T) in {203,205}.

The old weak parent may have arbitrary order, arbitrary singleton-SCC padding, and arbitrary feedback-vertex number. Only one new vertex is added, and no old arc changes. The conclusion is not asserted for a strong parent, a weak final child, multiple births, or an old-arc reset. It is not an absence theorem for 125,127,203,205: each has a verified strong witness in the bank.

The significance for the prime-atom lane is that 127 is prime and 203=7*29 is a multiplicative atom of the desired spectrum. Neutral padding followed by one strong birth cannot supply these literal ternary children from 41 or67, even after dropping the previous feedback-vertex-two restriction. This does not obstruct other methods of constructing either target.

## 2. All-order reduction to bounded rooted cores

Use the frozen Session 194 exact weak-parent form theorem in `STRONGER_RATIO.md`, section5. For h>=3, a strong child of a weak parent with H(T)=3h+2 or3h+4 has exactly two parent SCCs. Apart from the exceptional C3/H13 and C3+C3/H31 presentations, one SCC is a singleton and the other is a strong core C with H(C)=h. The relevant core-root response must satisfy

    R_a(C)-H(C) in {1,2}, or L_a(C)-H(C) in {1,2},

with the surplus selecting the digit. The values41 and67 are outside both exceptional cases. This theorem excludes all extra singleton padding; the reduction is not a search at a fixed padded order.

For completeness, since41 and67 are prime, SCC multiplicativity by itself already forces each weak parent to have exactly one nontrivial strong component of that count, with only singleton components elsewhere. A nontrivial strong tournament has H>=3 by its directed triangle and OCF, so a component of count1 must be a singleton. The sharper weak-parent theorem reduces the number and placement of those neutral components to one.

The cited Busch minimum B(n), with its standard strong-order domain and monotone three-residue formula, gives

    B(8)=45>41,       B(9)=75>67.

Hence every strong H41 core has order at most7, and every strong H67 core has order at most8. The finite task is now exactly to examine all roots of every strong core of those counts in those order ranges, in both endpoint directions.

## 3. Newly certified completeness, not only inherited deletion closure

The supplied 7,412-row bank is the same frozen bank used in Sessions191--194, with SHA-256

    5d5a5f625b4a85e04c0744aeb847ed36f59456f6a39996479aa82b9eb950afe3.

This packet does not rely on downward deletion closure alone to establish completeness. Its `verify_completeness.cpp` starts from the empty tournament and checks every single-vertex extension of every supplied representative at each preceding order, with every one of the 2^n attachment masks. It proves inductively that the canonical child set equals the next supplied bank at each order1 through8. It also checks that no supplied bank has duplicate isomorphism classes.

Canonicalization retains the complete adjacency matrix and minimizes it over every permutation within each outdegree class, placing those classes in a fixed increasing-score order. Every tournament isomorphism preserves outdegrees, so these permutations cover every possible isomorphism between matrices in that ordering. No score sequence or scalar invariant is substituted for full adjacency. Unsigned64-bit encoding suffices because the largest retained matrix has8*8=64 bits.

Thus equality of the extension-generated set with the supplied set proves both completeness and exact equality to the supplied bank, without invoking an external count of tournament classes. The familiar class counts1,1,2,4,12,56,456,6880 are outputs and consistency checks, not the completeness argument.

## 4. Independent endpoint calculations and hostile controls

`verify_profiles.py` computes H and every rooted R,L,s,e by subset dynamic programming in both endpoint directions. For every H41/H67 core it separately enumerates every vertex permutation. A word with no broken adjacency contributes one to H and to every R_i and L_i. A word with exactly one broken adjacency contributes to R at the vertex immediately after the break and to L at the vertex immediately before the break. More than one break contributes to neither. This is a literal enumeration of the two ordered path bins and is independent of the subset response recurrence.

Every rooted profile agrees. None has R-H or L-H equal to1 or2. The full inventories, including every row matrix and every root response, are retained in `profile_receipt.json`. Together with sections2--3 this gives the claimed all-order exclusion.

Positive controls verify an H5 core with both surplus1 and surplus2, so the code does not exclude all ternary-compatible roots. Literal profiles also verify actual strong witnesses for the two parent counts41,67 and all four target counts125,127,203,205. The theorem is about a particular ancestry operation, not missing target support.

## 5. Reproduction and status boundary

    g++ -O3 -std=c++17 verify_completeness.cpp -o /tmp/s194_completeness
    /tmp/s194_completeness inherited_tournaments_through8.txt > completeness_receipt.json
    python3 verify_profiles.py > profile_receipt.json
    python3 -O verify_profiles.py > profile_receipt.optimized.json
    cmp profile_receipt.json profile_receipt.optimized.json

Run an independent unoptimized C++ build as a compiler-sensitive control and compare its receipt. The all-order consequence also depends on the separately audited weak-parent form theorem and the cited Busch order bound. It must not be promoted unless those dependencies and this certificate pass independent review.

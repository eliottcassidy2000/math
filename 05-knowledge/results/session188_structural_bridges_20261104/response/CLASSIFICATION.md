# Exact order-five classification of the full named quartic tensor

**FINITE-EXACT:** among all 1,024 tournaments on five fixed physical names, the complete quartic tensor has exactly 502 equality classes. Of these, 492 contain exactly a tournament and its whole converse. The other ten classes have four members and are described below. **This is not a full-response collision census.**

For each named three-element set M, orient M as a directed triangle. Of the remaining two names, choose one as a universal source and the other as a universal sink. The four possibilities (two triangle orientations times two exterior orders) have the same named quartic tensor. They form one four-member class. These ten classes, one for each M, are all the exceptions.

Consequently every strong tournament on five fixed names is determined up to whole converse by its quartic tensor alone. There are 544 labeled strong tournaments in the exhaustive census. This statement is bounded to old order five; no extrapolation to order six or arbitrary order is made.

## Exact evidence rather than hash-only comparisons

`classify_five.py` independently regenerates the quartic normal form for all 1,024 tournaments. It groups by literal tuples of **all 100 named bidegree-(2,2) coordinates**, including zeros, not by a hash. It checks that each two-member class consists exactly of complementary edge codes. For every four-member class it checks the unique source, unique sink, directed middle triangle, both exterior orders, and all ten possible named triangle sets. Strongness is tested by reachability from every name.

`classification.json` retains the common named basis and every class's exact vector, codes, and adjacency rows. The underlying `verify.py` census compares every one of these quartics to full Hamilton-word expansion, independently of the cover construction. Its finite controls retain digests for compactness, but the classification itself uses exact vectors and retains them.

## The exceptional quartic collision is not a full-response collision

Take

    T=(10,0,3,6,15),
    U=(6,0,10,3,15).

Both have physical source 4, sink 1, and the directed triangle on names 0,2,3. Only that triangle is reversed. All quartics agree. The scalar Hamilton counts are both 3, and the complete linear parts also agree. Nevertheless the same actual query gives

    F_T(a=1,c=4)=31,
    F_U(a=1,c=4)=29.

The query means a_0=1 and c_2=1, with every other attachment bit zero; x->y remains fixed. Their complete difference has six nonzero quadratic and 36 nonzero cubic coordinates, all retained in the receipt. It has no quartic. Thus high-degree agreement does not imply full-response equality, even for an actual module converse.

This example is compatible with S187 localization: the active quotient SCC is the triangle's singleton contracted root, so an intrinsic separator must survive the ordered source/sink padding. It does not contradict the inherited module-converse obstructions.

## Reproduction

    python classify_five.py > classification.json
    python -O classify_five.py > classification.optimized.json
    cmp classification.json classification.optimized.json

Both runs pass and agree byte for byte. This exact finite classification is an additional response-theoretic result, separate from the all-order quartic normal form and transitive-exterior theorem in `PROOF.md`.

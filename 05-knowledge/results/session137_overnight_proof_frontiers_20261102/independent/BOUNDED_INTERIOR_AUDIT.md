# Independent audit: fixed height and interior normal form

**Verdict: ACCEPT** `constructors/BOUNDED_INTERIOR.md` in the version pinned
by `SOURCE_PINS.json`. This is an all-source classification at fixed (h,L),
not a grammar-restricted census or a global missing-H claim.

The coefficient comparison is correctly normalized: both expressions are for
`C=(F-h)/2`, so `b_uv=R_uv+R_vu`, not twice that quantity. If r singleton SCCs
lie between selected endpoints in distinct SCCs, their even subsets give
exactly `2^(r-1)` even-cardinality directed paths for r>=1. Complementary H
weights are positive. Hence the bound `r<=1+floor(log2 L)` is valid, including
selected intermediate singleton vertices.

The unselected-block map and S0 commute on the entire native four-state
vector. Moving suffix singleton SCCs to the end therefore gives a literal
actual source with the same state; no claim of labelled isomorphism or
compatibility with arbitrary old-vertex re-marking is made. Removing neutral
prefix singletons and the final S0 tail leaves the stated finite base order.

The L=0 boundary is now complete: the base may have g>0, and all unselected
suffix singletons may be removed. A singleton selected SCC is the only
possible retained singleton outside the nonsingleton core, giving order
at most h+1. Empty cut is the separate constant F=h case. This avoids
incorrectly applying the older L=g=0 theorem to L=0,g>0 sources.

For a positive tail k, `q=2^k q0` and g=0, so the defining interior identity
gives `F=2^k q0+h-2L`. At k=0 the original base F0 must be retained separately.
The exact inverse tests this finite bank and a positive integral power of two;
it does not accept arbitrary formal positive states. The O(log X) conclusion
for every finite (h,L) bank follows immediately. Thus bounded height **and**
bounded interior cannot furnish the full spectrum, without asserting that
unbounded resources are sufficient.

At h=L=1, positivity forces exactly two selected vertices in a transitive
parent. Their pair weight is one only when they are adjacent or have one
intervening vertex. The normalized cut words are 11 and 101, yielding native
bases `(1,1,3,1)` and `(1,3,5,1)` and precisely the claimed two affine-dyadic
rays plus their two zero-tail base values.

## Independent replay

`audit_bounded_interior.py` imports no producer module. It builds all
tournament isomorphism classes through order six, takes every labelled cut,
computes endpoint/ordered-cover counts independently, computes SCCs by
reachability, and directly tests the induced base and complete restored
native state. This gives 4,054 marked sources, 3,557 with positive L,
455 nonempty deleted suffixes, 455 nonempty deleted prefixes, and a maximum
deleted suffix of five vertices. Every source passes the appropriate base
bound and exact affine-dyadic state law. The h1,L1 base bank matches exactly.

Normal and optimized outputs are byte-identical. Reproduce:

```
python3 -B audit_bounded_interior.py > /tmp/bounded-interior.json
python3 -B -O audit_bounded_interior.py > /tmp/bounded-interior-O.json
cmp /tmp/bounded-interior.json bounded_interior.normal.json
cmp /tmp/bounded-interior-O.json bounded_interior.normal.json
```

These finite checks support the source-preserving normal-form proof; they do
not replace its universal quantifiers. Source and output hashes are in the
aggregate manifest.

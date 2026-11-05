# Independent H15 classification and cloning audit

2026-11-05. **PASS for the all-order H15 classification, counting formulas,
and complete transitive-cloning responses. No correction was needed.**

Frozen proof SHA256:
`95304b8751cbb6979439afd49b69dd8ca17d702214013bd3db7a0522ab8e23c7`.
Verifier SHA256:
`5c4ba38dbada155f82de452097da3f4b5231bf1b8cb4a517e4cc92dc37633d5a`.
Receipt SHA256:
`164c2e42c0bba417a994ef281830a2409d0022de4fcbe208c4030badde3f3ad5`.

## All-order reasoning

The first-outneighbor insertion is an injection from Hamilton paths of a
nontrivial strong parent to child paths not ending at the added mixed
vertex. A parent path ending at an inneighbor supplies a disjoint extra
path. Strongness supplies that endpoint path; parity then improves strict
growth to an increase of at least two. The proof does not incorrectly
assume that arbitrary deletion preserves a Hamilton path.

The named S195 universal strong-growth theorem was independently read in
the preceding audits. Its mixed-singleton/cyclic-pair alternatives apply
starting from an induced triangle. A cyclic pair after order three has
count at least 29. A pair directly from the triangle gives order five and
count 15, after which strict growth rules out further steps. Without a
pair, the singleton chain reaches the exact order-six minimum 15 before
any order-seven member could reach at least 17. Thus the cutoff at order
six is proved at equality H=15; it does not invoke the strict H<15 module
threshold. Independently, Busch's previously checked bound H^3>=5^(n-1)
also gives the cutoff, since 15^3<5^6.

The complete strong census yields precisely the four stated cores, with
automorphism orders 3,5,1,1. Multiplicativity over the canonical SCC order
then leaves exactly one H15 core, or an ordered H3/H5 pair, together with
singleton padding. These are complete isomorphism invariants. This gives
the all-order unlabelled count (n-3)(n-4), and labelled count
n!*(5n^2-17n-32)/15 for n>=5. The low-order zero/product conventions at
n=5 and n=6 are correct. Automorphisms cannot exchange SCC positions.

The named S196 clone-run and matching-rank results were independently
read in the preceding audits. Every core root has matching rank at most
two, so the literal run coefficients prove the entire functions, rather
than an extrapolation from a few multiplicities. Singleton roots give
the constant function. Roots in the ordered H3/H5 pair give no new
function after multiplication by the untouched component count. Thus
the eight displayed functions are exhaustive for arbitrary, not just
strong, H15 tournaments. At multiplicity two their values are exactly
15,25,27,33,37,41,59. Every nonconstant function is increasing and is at
least 45 from multiplicity three onward. Consequently neither 29 nor 31
occurs at any positive multiplicity, even as a one-off step. No assertion
about other construction operations follows or is made.

## Independent verification

The independent degree-based canonical census from the H67 audit covers
every tournament through order eight, with labelled orbit masses checked
there. It confirms the four strong cores and both counting formulas at
each order. Padding constructions independently reproduce the formulas
through order twelve; their actual Hamilton counts are checked through
order nine. The minimum among all strong order-six classes is exactly 15.

All 2,370 mixed extensions of nontrivial strong representative parents
through order six satisfy strongness, oddness, and the increase by two.
The audit computes clone-run words using a separate memoized word-building
recursion, not the source's permutation-and-gap enumeration, and computes
matching rank by exhaustive edge-subset testing. Every root of all four
cores recovers the stated run support and function. Independently built
clone graphs, with a different vertex layout, give 236 exact controls,
including the two weak product orders and singleton padding. All pass.

Normal and optimized source runs agree byte-for-byte with the frozen
receipt. The independent program is `audit.py`, with results in
`AUDIT.json`; it uses the sibling independent primary and canonical-census
audit helpers. Author source files were not edited. The conclusions
combine proved all-order reductions with complete finite certificates;
no finite-to-infinite extrapolation is used.

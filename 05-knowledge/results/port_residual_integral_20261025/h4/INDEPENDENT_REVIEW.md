# Independent review of the h=4 backward repair

2026-10-25. **REVIEW PASS:** the literal three-black detour in `README.md`
is proved for the stated endpoint-four family and n>=12. No producer
correction was needed. **FINITE-EXACT:** producer replay and a separate
Gaussian-elimination audit described below.

The exact source is a labelled tournament T with directed, closed blue
P=(0,...,n-1), N^-(0)={4,n-3,n-2,n-1}, and 3->1. The affine source
engine imposes precisely these pair-orientation equations. It retains
all unspecified middle arcs and tests every completion, without clone
contraction or an isomorphism filter. Nonempty-source checks are present.

## Literal legality and target

Q1 introduces only0->5 and (n-1)->1. The first uses the exact port row;
the second is the r_P-reflection of (n-2)->0. Q2 introduces2->0 and
1->3, the complements of source arcs omitted by Q1. Q3 restores0->1
and4->5 after two complements. The remaining joins are retained from
the preceding word. All words include every labelled vertex exactly once.

The two edge differences X={01,45,05,1(n-1)} and
Y={12,23,02,13} give E(Q1)=E(P) XOR X,
E(Q2)=E(P) XOR X XOR Y, and E(Q3)=E(P) XOR Y. Three K moves therefore
reach the literal labelled target K_P(T), rather than merely an
isomorphic tournament. This XOR calculation is logically separate from
the legality verification; both have been checked.

## Three anti witnesses audited

For Q1 at T the reversed ordered pair of (0,3) is (5,4).
The former is an outgoing nonport pair and the latter reverses the
original path edge4->5, so their arc values differ.

For Q2 at S=K_Q1(T), the reversed ordered pair of (3,n-1) is (6,0).
The original anti map gives T(3,n-1)=T(0,n-4)=1; 6 is also a nonport,
so T(6,0)=0. Both unordered pairs are omitted by Q1, yielding source
values0 and1 at S. The bound n>=12 guarantees the named nonports and
all witness pairs are distinct from the exceptional port/path pairs.
No unknown middle arc enters this contradiction.

For Q3 at B=K_P(T), the reversed ordered pair of (0,3) is (n-4,n-2).
Both are nonpath pairs of P, so B(0,3)=0 and
B(n-4,n-2)=1-T(n-4,n-2)=1-T(1,3)=1. This uses the specified guard
3->1. Testing at B is legitimate: positional reversal preserves E(Q3),
so complementing all other pairs leaves its anti status unchanged.
Hence the third word is also black at its actual preceding source.

The three comparisons establish the result for every source completion
and every n>=12. In particular, the proof does not infer global symmetry
from a small induced tournament. Existing clones, arbitrary middle arcs,
and literal labels remain in the source and target.

## Hostiles, finite scope, and reproduction

The old0->4 join remains false in this source family. Adding5 to the
incoming-port set makes the new0->5 join false in a nonempty source cell.
Thus the exact port set is a real hypothesis, and the outside-scope
hostile is not discarded as inconsistent. Neither word failure is a
connectivity counterexample.

Producer `check.py` passed normally and with `python -O`: all109 affine
cells n=12..120,327 seeded literal models, and edge/reversal identities
through1000. A separate Gaussian implementation, with no producer import,
certified all completions for n=12..40 and checked witness images and
edge XOR through512. It also independently verified the extra-port
hostile. Its full reproducible code and normal/optimized invocation are
embedded in [the h3 independent review](../h3/INDEPENDENT_REVIEW.md), since
that one audit program checks both families.

This closes the specified h4,3->1 residual and combines with the inherited
h4,1->3 branch. The stronger global black-core connectivity statement
remains outside these proofs.

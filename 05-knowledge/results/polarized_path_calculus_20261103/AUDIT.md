# Proof and computation scope

**PROVED:** the start-budget decomposition and its converse strongness
guards; the exact Fibonacci class generating functions; one-port and
multiport tail formulas; the polarized seam; and sole-new-root birth
operators. The proofs are direct path decompositions. The standard Camion
Hamilton-cycle theorem is used only to ensure that every vertex of a
nontrivial strong component starts a path and that its H is at least three.
The general theory of path covers, the old scalar reflection formula,
four-coordinate reset and golden field convention are inherited, with
links in README. No external priority assertion is made.

**FINITE-EXACT, `explore.py`:**

* 81 rooted input pairs, taken from the nine rooted path charts at orders
  one through three, for the seam at m=1,2,3;
* the same 81 pairs under ordinal sum, with both root choices, at m=1,2,3;
* all nine sole-new-root outdegree-one births;
* 135 one-port funnels: those nine roots times all words of length 0..3,
  at m=1,2,3;
* 176 multiport funnels: all nonempty port sets on the four specified
  literal cores (TT3, C3, and the order-four path charts 1 and 3), times
  every length-two word, at m=1,2;
* all 1,098 fixed-spine tournament charts at orders 3..6, every root,
  filtered to strong sources. The 955 root/chart occurrences with start
  budget one or two satisfy the core-and-tail decomposition. These are
  chart occurrences, not 955 distinct rooted isomorphism classes.

**FINITE-EXACT, `counts.py`:** all fixed-spine charts through order seven
(1+1+2+8+64+1024+32768). Requiring exactly one path starting at vertex zero
leaves exactly one representative of each rooted unique-start class;
strongness is tested independently by reachability. The resulting strong
counts 1,3,8,21,55 and all-source counts 1,1,2,5,13,34,89 agree with the
proved rational generating functions. H histograms are in counts.json.

Independent controls use 123 permutation/weak-composition bin counts
on every rooted fixed-spine chart through order four, at m=1,2,3;
189 direct tail-colouring counts for all words of length 0..5 and m=1,2,3;
and two literal order-eight polynomial constructions at m=1,2,3,4.
Twenty actual Q4 attachments check every one of the five inherited
presentations on four sources, and check q2's F49 transport. Two H31
entry witnesses are expanded by the literal two-birth reset and their
95/97 outputs counted independently from the transfer formula.

The bin DP first counts Hamilton paths of every induced subset, then
convolves disjoint bins. The permutation implementation instead enumerates
each full vertex word and its weak bin-size composition. These are distinct
counting paths. Conditions use explicit exceptions, not removable asserts.

The structural recognizer in `recognize.py` is a third algorithm. It uses
no Hamilton-count DP to recognize the class or construct its unique path.
It independently obtains the class counts from all 33,868 charts through
seven, and agrees with endpoint DP on all 361 rooted chart occurrences
through five. Recognition and canonical labelling have the proved O(n^3)
adjacency-inspection bound; the reference implementation uses recursion.
Rank/unrank round trips cover all 1,596 classes at orders 3..10, including
a reversed-label control for each, and three selected ranks at each of
orders 50 and 100. Large-order checks reconstruct literal tournaments
and use the structural recognizer; they do not enumerate their Hamilton
paths or every class at those orders. The optimal fixed-order bit bound
is proved by the rank bijection, independently of these examples.

Reproduction, from this directory:

    python -B explore.py
    python -B counts.py
    python -B recognize.py

The committed explore.out, counts.json and recognition.json are stdout. Optimized (-O)
runs must agree after newline normalization; verification.json records
that comparison and the source-file hashes. The source pin for the Q4
oracle is golden_inverse_defects_20261103/inverse.py together with its
golden_seam_recurrence_20261103/atom_atlas.json input. No native search or
81/27 computation is silently rerun or claimed here.

## Development corrections retained

An initial interpolation probe assigned E(0)=S(0)=0, confusing the
combinatorial domain m>=1 with polynomial evaluation. The singleton
immediately rejects that assignment; its E=S=1. The corrected compiler
interpolates rooted polynomials on positive m and extrapolates formally.
This is documented in CALCULUS and covered by the independent controls.

An exploratory sentence assumed both H31,F=32/33 examples arose from
small-start parents. Only the first supplied witness has that form; the
second has indegree two and defect two. The final theorem states a
sufficient source construction and retains that distinct entry, rather
than promoting a nonexistent classification of all F-H=2 roots.

The component-sort implementation originally consulted the list being
sorted from inside its sort key. Python temporarily empties that list;
the directed triangle rejected the result. The final code sorts against
an immutable copy. These draft issues were repaired before publication;
no previously published theorem is being retracted.

## Boundaries and open next steps

The polynomials need not classify sources or support arbitrary old-root
switches. Actual source construction and named path data are retained as
sidecars; the 81/27 hostile remains relevant to any wider carrier claim.
No finite-dimensional scalar carrier, arbitrary-complement closure,
complete prime/7p supply, or A000568 extension is proved. The two-start
decomposition is complete, but a closed two-marked recursive compiler
requires propagation of the displayed mixed term and any further joint
data tested by later operations.

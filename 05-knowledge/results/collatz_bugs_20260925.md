# Collatz bugs: admissible motifs, directional distance, and prime overlap

**Status:** PROVED elementary scoped statements below; FINITE-EXACT controls.
Universal connection of bugs is equivalent to the OPEN positive Collatz
conjecture. The square-sum/negative-Collatz and graceful/positive-Collatz
isomorphisms are OPEN proposed connections, with no correspondence established.
No novelty claim. Session: codex-collatz-bugs-20260925.

## Inheritance and concept board

The user's construction takes an even identity I=2^m q, q odd, the spine
I/2 -> I -> 2I -> 4I, and odd branches (I-1)/3 and (4I-1)/3, with
arrows reversed relative to the shortcut Collatz map. It proposes a directional
count of backwards traversals between bugs, and connections to square-sum
Hamiltonian paths, prime encodings, and graceful tree labelings.

Closest proved mechanism: Proposition 11 of
[the inverse-tree note](collatz_procgen_20260924_inverse_tree_mod192.md),
E D^2 = S E for D(x)=2x, E(x)=(2x-1)/3, S(x)=4x+1; its Theorem 1
identifies guarded closure with the known root basin, not all positive integers.
Canonical hostile: the cycle 1<->2 prevents the unreduced graph being a tree;
I=8 fails the integrality guard. Corrected near miss: a residue class decides
only the recorded precision of an inverse branch, not unrestricted ancestry
(same note Proposition 2; 2026-09-24 MISTAKES entry).
Least-used relevant sidecar here: the height of the common ancestor, together
with the choice of clock (individual T edges versus accelerated odd steps).

Anchor: bug distance and global connectivity. Niche: square-sum endpoint
certificates. Wildcard: prime overlaps along the sibling ladder.

| Live object | Invariant / operation | Lost data / cheapest hostile |
|---|---|---|
| Bug motif | legal inverse edges; scale I by 4 | mod-3 guard; I=8 |
| Directed distance | ordered pair of LCA heights | cycle / component; 1,2 |
| Prime encoding | full prime-exponent vector of q | addition carries; 5,21,85 |
| Square-sum graph | existence of a spanning path | endpoint choice; Q23 endpoint22 |
| Graceful tree | vertex interval and distinct edge differences | ancestry metric alone lacks a labeling |

Methods used: META-PATTERNS, "Search the statement before the method" and
"Type every analogy and every implication". No new method card is needed.

## 1. The exact bug and its recursion (PROVED)

Set T(n)=n/2 for even n and T(n)=(3n+1)/2 for odd n. The reversed graph
has x -> 2x always and x -> (2x-1)/3 precisely when x=2 mod3.
Hence the proposed two odd branches both exist **iff I=4 mod6**.
For those identities put a=(I-1)/3. Then

    I/2 -> I -> 2I -> 4I,
    I/2 -> a,
    2I -> 4a+1.

Proof: T(a)=I/2 and T(4a+1)=2I; both leaves are positive odd integers
iff I is even and I=1 mod3. The map I -> 4I preserves this condition.
The two switches are at I/2 and 2I. Thus saying I/2 switches immediately
between a and 4a+1 suppresses the two doubling edges to the second switch.

Writing I=2^m q with q odd imposes q=(-1)^m mod3. For m>=3, I/4 is
another admissible bug. At m=1 or m=2 the upward even recursion reaches q;
q is the first odd leaf of the admissible bug J=3q+1. In the first case
its incoming reversed-graph vertex is T(q)=J/2. In the second case I/4=q;
one has ended the even recursion, not the full graph traversal.

Bugs are overlapping motifs. For example B4 and B16 share vertices 8,16,5
and edges 8->16 and 8->5. The following distance uses identities as marked
vertices in the original graph; it does not assert that motifs partition it.

## 2. A rigorous f and the exact global obligation (PROVED)

Let C={n>=1: T^k(n)=1 for some k}. Removing 1 and 2 from the reversed
graph on C gives a rooted tree at 4. Indeed, any orbit from outside the
root cycle entering it must first pass 4, since the only preimages of 2
are 1 and 4 and the only preimage of 1 is 2. Each nonroot vertex has
one parent, and every vertex has finite height h(n)=min{k:T^k(n)=4}.

For identities X,Y in this tree, define f(X,Y) as the minimum number of
edges traversed against the reversed arrows along a path from X to Y.
Write c=LCA(X,Y). Unique paths in a tree give

    f(X,Y)=h(X)-h(c),   f(Y,X)=h(Y)-h(c).
    d(X,Y)=f(X,Y)+f(Y,X).
    h(X)-h(Y)=f(X,Y)-f(Y,X).

This proves the proposed zero criterion: f(X,Y)=0 iff X is an ancestor
of Y in the inverse tree (equivalently Y's forward T-orbit contains X).
If both values are positive, T^f(X,Y)(X)=T^f(Y,X)(Y)=c is the first
common switch. A genuine branching vertex has c=2 mod3. For X=Y both
values vanish, so the phrase "the nonzero value" requires distinct vertices.

Example: X=10, Y=16 have forward paths 10->5->8 and 16->8.
Thus c=8 and (f(10,16),f(16,10))=(2,1). The meet need not be an
admissible bug identity (8 is not). Compressing to bug identities can lose
the switch and therefore changes this f unless its address is retained.

These numbers count T edges. A distance counting odd jumps requires the
odd-only graph and a redefined unit; counting a change in travel direction
would give at most one on a simple tree path and is a different statistic.

On the unrestricted reversed graph, minimum reversal cost is still defined,
with infinity if there is no finite underlying path. But on its root cycle
f(1,2)=f(2,1)=0 although 1!=2; do not claim tree identities there.

**Equivalence.** Every admissible I has finite f(I,4) iff positive Collatz
holds. The forward implication uses I=3n+1 for each odd n: I is admissible
and T(I)=T(n), so n and I share their eventual behavior. Basin membership
is constant across every underlying graph edge, so a finite connection to
4 implies membership in C. Every even integer halves to an odd one. The
reverse implication follows from finite forward trajectories. It suffices
equivalently to require all admissible bugs lie in one weak component.

This pinpoints the research obligation: prove global connection without
assuming it when defining the domain or asserting a common switch.

## 3. Prime encoding has a concrete local consequence (PROVED)

Every odd q>1 has a unique odd-prime factorization; q=1 has the empty
factorization. Retaining exponents gives a complete encoding of q, but
does not by itself supply a law for ancestry after adding 1 in 3q+1.

There is nevertheless an exact prime statement along a bug ladder:

    a_(j+1)=4a_j+1,
    gcd(a_j,a_(j+1))=1,
    gcd(a_j,a_(j+r)) divides (4^r-1)/3.

Proof: iteration gives a_(j+r)=4^r a_j+(4^r-1)/3; use the Euclidean
algorithm. Thus adjacent odd leaves have disjoint prime supports, while
nonadjacent leaves may share primes: 5,21,85 has gcd(5,85)=5. The
prime-support-only quotient additionally loses multiplicity (3 and9),
which a full prime-exponent vector restores. None of these local coprimality
laws implies eventual arrival at the root.

Also 3a_(j+1)+1=4(3a_j+1), so all ladder members have the same odd
successor under U(a)=oddpart(3a+1). This recovers the existing sibling
ladder theorem; f separates their different locations on the doubling spine.

## 4. The two proposed bridges and their missing predicates

**CITED:** the nontrivial square-sum path sizes are 15,16,17,23 and all
N>=25; Gerbicz records the all-N result in
[OEIS A090461](https://oeis.org/A090461), comment dated 2018-01-21,
checked 2026-09-25. N=1 is a vacuous path under the usual convention, so
"first fourteen impossible" uses a nontrivial-arrangement convention.
For N>=2 the impossible ranges are 2..14,18..22,24.

**PROVED elementary:** negating an integer conjugates positive 3n-1 to
negative 3n+1. Under shortcut 3n-1 the displayed cycles are {1},
{5,7,10}, and {17,25,37,55,82,41,61,91,136,68,34}; they have minima
1,5,17. Exhaustion by these three basins is OPEN, not an eventual
threshold theorem analogous to square-sum Hamiltonicity. The exact
square-sum sequence alone provides no map between the predicates.

For gracefulness, the closest exact bridge is
[THM-4470, pairing ladder](../../01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md),
item 2: the halving edges {i,2i} and odd edges {2i-1,3i-1} each realize
each positive difference i exactly once. This is a genuine shared
difference structure. A graceful labeling of an e-edge tree additionally
requires its e+1 vertex labels to be exactly 0..e and edge differences to
be exactly 1..e. The f pair recovers rooted metric information; it does
not construct those labels. A direct subtree isomorphism also meets a
degree obstruction: the inverse Collatz tree has degree at most3, whereas
the tree K_(1,4) has degree4. A more elaborate representation would have
to specify its map and how it preserves the graceful-labeling predicate.

Connection contracts: bug -> guarded inverse graph preserves legal edges
under identity placement (retain integrality and clock); graph -> f
preserves ancestry and distance on the stated tree (retain component,
root-cycle convention and meet address); ladder -> prime supports
preserves adjacent coprimality (retain exponents to restore labels).
Neither of the proposed problem-to-problem isomorphisms yet has a map.

Board revisit: the meet example shows why bug compression needs switches;
the root-cycle hostile makes global f finiteness the anchor obligation;
the ladder makes the prime wildcard exact but local; the square-sum audit
requires endpoint data, a different missing coordinate from Collatz height;
THM-4470 supplies the graceful mechanism without a convergence implication.

## 5. Reproduction and correction lineage

Run from repository root:

    python3 04-computation/experiments/collatz_bugs_20260925.py

Output: [collatz_bugs_20260925.out](collatz_bugs_20260925.out). Universes:
all even I from2 through5000; 100 odd ladder seeds through199 with7 rungs;
all ordered pairs of admissible identities through298. No pruning filters.
The tree test includes every intermediate trajectory vertex, without a
size cap; a step limit aborts instead of silently dropping a case. Independent
0-1 BFS checks the LCA reversal counts. Controls include valid I=10,
invalid I=8, the root cycle, and nonadjacent shared prime5.

During inheritance, an unrelated proof error was found in the square-sum
note's degree-2 forcing rule: path endpoints may have degree2 in the host
graph. The repaired note and MISTAKES entry retain endpoint information.
This correction does not alter the square-sum existence pattern.

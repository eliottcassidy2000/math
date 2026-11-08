# Exact strongness criterion and a first-birth separation

**PROVED structural criterion; FINITE-EXACT separating example.**
Chan–Pak all-odd completeness remains OPEN. No literature-priority claim is made.

## 1. Inherited operation and exact criterion

Use the rooted graft from [Session225](../../session225_graft_recovery_20261108/PROOF.md).
For nonempty rooted tournaments (A,a),(B,r), retain both induced tournaments,
put B -> a and A-a -> B, and mark r. Write the resulting graph G_B(A,a).
The root r does not affect the following unrooted strongness criterion.

**Theorem.** G_B(A,a) is strongly connected if and only if |A|>=2 and a
reaches every vertex of A. No connectivity hypothesis on B is required.

**Proof.** If A is a singleton, a is a sink in the graft. Otherwise suppose
a reaches all of A. Choose an old nonroot v. The root a reaches v and then
all of B. Every vertex of B reaches a directly. Every old nonroot reaches
all of B directly, then a, then every old vertex. These paths establish all
ordered-pair reachability, independently of the arcs internal to B.

Conversely, a simple directed path starting at a and ending in A cannot
visit B: every return from B to A must enter a, repeating its initial
vertex. Thus any vertex of A reachable from a in the graft was already
reachable in A. Strongness forces all of them to be reachable. This proves
both directions, including arbitrary nonstrong gadgets. ∎

**Corollary.** Start with any rooted tournament of order at least two whose
root reaches all vertices, and apply any finite sequence of nonempty rooted
gadgets. Every nonempty iterate is strong. Its newly marked root reaches
every vertex, so the criterion iterates. The gadgets need not be equal,
strong, or invertible as state operators.

The prior Session225 proof gave a sufficient statement with B strong. The
present criterion removes that unnecessary hypothesis and proves necessity.
It does not say that a strong output determines the source graph or its path
count, or that every strong tournament has such a decomposition.

## 2. A concrete bridge to the inherited H=9 birth collision

Session224 proves that

    A=[14,28,24,16,1], B=[10,28,25,16,1]

are the two nonisomorphic strong classes in the first-birth H=9 fiber.
Rows encode outgoing-neighbor masks on vertices 0,...,4. That classification,
their minimal order, and their self-converse status are inherited, not new
claims here. See [the exact source](../../session224_first_occurrence_20261107/proof/PROOF.md).

Apply the fixed Session225 four-vertex gadget C=[8,1,3,6], rooted at 0,
to each of the five possible source roots. Literal subset dynamic programming
gives the following counts, with roots ordered 0,...,4:

    H(G_C(A,a)) = 423,605,357,185,261,
    H(G_C(B,b)) = 315,541,551,217,207.

The sets are disjoint. Thus this single gadget separates these two inherited
unrooted birth atoms by the resulting Hamilton-path count, regardless of
which source root was chosen. All ten outputs are strong by the theorem.
This is a particular finite separator, not a complete invariant of tournaments,
a minimal separator, or a scalar-spectrum separation from O/Z/J constructions.

For comparison, the source four-coordinate states (D,E,U,d), in root order, are

    A: (9,1,15,1), (9,1,25,5), (9,2,16,5), (9,4,10,5), (9,1,9,1),
    B: (9,1,11,1), (9,1,21,3), (9,1,23,5), (9,3,11,5), (9,3,9,3).

Multiplying each state by the inherited Hamilton-output row (5,-5,27,-22)
gives the same values. The verifier computes the ten graph counts directly;
it does not infer them from those state vectors.

## 3. Verification and scope

Run `python3 verify.py > replay.json` and compare with CHECK.json; repeat
with `python3 -O`. Standard library only. All guards use explicit exceptions.

The exact universe comprises every labelled tournament on orders 1 through 4,
with root label 0, giving 75 sources and 75 gadgets and all 5,625 ordered pairs.
No strongness filter is imposed. Reachability is checked on literal grafts
against the criterion, including singleton sources, sink roots and nonstrong
gadgets. The ten H=9 extensions have nine vertices and are counted by an
independent path subset recurrence. Finite controls support the proof but do
not replace its universal quantifiers.

# Two-sided adaptive peeling is nonconfluent from order six onward

**PROVED:** nonconfluence at every order n≥6 by an explicit critical family.
**FINITE-EXACT:** six is the smallest order; every labelled strong source of
orders 3,4,5 was checked independently. This is a new proposed rule's hostile
probe, not a retraction of the oriented theorem in [ANCESTRY.md](ANCESTRY.md).

The one-sided unique inverse suggests allowing its converse as well. A legal
inverse now deletes either an indegree-one vertex or an outdegree-one vertex,
in each case only when the remaining tournament is strong. Each orientation
has at most one legal choice; their union can have two. That union does not
have a canonical terminal core, even modulo converse.

## Minimal named witness

Use outgoing row-bit notation on names 0,…,5:

    T=(16,9,3,37,14,23).

T is strong. Both named 5 (indegree one) and named 0 (outdegree one) are legal
first deletions.

| Branch | Legal named deletions, with orientations | Terminal names | Terminal rows in that name order |
|---|---|---|---|
| A | 5, indegree one | (0,1,2,3,4) | (16,9,3,5,14) |
| B | 0, outdegree one; 4, indegree one; 5, indegree one | (1,2,3) | (4,1,2) |

Branch A is the irreducible strong five-vertex cyclic composition C3[C3,1,1].
Branch B is C3. Both have no further legal deletion of either orientation.
Their different orders exclude ordinary isomorphism and converse equivalence.
The complete set of terminal ordinary types for this source is precisely these
two types. All vertex names survive the displayed transports; branch B's
second and third numbers are original names, not reindexed positions.

## The critical-pair mechanism

Let x be an indegree-one legal inverse and y an outdegree-one legal inverse
of a strong source of order >3. They are different. Necessarily x→y: otherwise
y's unique outgoing arc is y→x, and deleting x makes y a sink, contradicting
strongness. If U=T−{x,y} is strong, both deletions remain legal and commute.
Therefore failure can only occur when U is not strong.

Write a for x's unique predecessor and b for y's unique successor. They lie
in U. Strongness of U+x forces a into U's terminal SCC; strongness of U+y
forces b into U's initial SCC. This identifies the missing coordinate: the
ordered SCC structure after **both** proposed deletions. One-deletion tests
do not certify the common child.

The smallest hostile has U=b→K with b a singleton and K=C3. One branch
closes this SCC chain into a cyclic composition that blocks both deletions;
the other retains a legal adaptive stack above K.

## An obstruction at every larger order

Let K be any strong tournament of order k≥3, and choose a named a∈K.
Construct T on K plus names x,b,y as follows:

1. Birth x with unique in-neighbor a.
2. Birth b with unique in-neighbor x. Thus b beats every vertex of K.
3. Birth y with unique out-neighbor b. Thus every vertex of K and x beats y.

T is strong. In fact the first two steps are ordinary adaptive births over
strong parents. At the last step y reaches the old source through b and every
old vertex reaches y directly except b, which reaches y through K.

There are two incompatible legal reductions:

* Delete x. The remaining tournament is the cyclic composition
  K→y→b→K, a strong tournament of order k+2.
* Delete y, then b, then x. All three deletions are legal with orientations
  outdegree one, indegree one, indegree one, and recover **exactly** the named K.

The first branch is irreducible under both rules. Every vertex of K has an
internal in-neighbor and an internal out-neighbor, and gains b as an extra
in-neighbor and y as an extra out-neighbor. Thus neither its indegree nor its
outdegree is one. Vertex b has unique in-neighbor y, but deleting b leaves
K→y nonstrong. Vertex y has unique out-neighbor b, but deleting y leaves
b→K nonstrong. No legal choice remains.

The second branch can be continued to some terminal core of order at most k.
It cannot equal the first branch's core of order k+2. Strong tournaments exist
at each k≥3: start with C3 and use an adaptive birth repeatedly. Thus the
construction proves nonconfluence at every source order n=k+3≥6.

For the displayed minimal witness, K has original names (1,2,3), a=3, x=5,
b=4, y=0. The construction is exact, not merely an isomorphic motif analogy.

## The irreducible cores already contain all strong tournaments

The parent lane extracted a further exact consequence of this critical family.
For every strong K of order at least three, define

    W(K)=K→y→b→K.

The proof above shows W(K) is irreducible under both deletion rules. Moreover,
b is its unique indegree-one vertex and y its unique outdegree-one vertex;
all vertices of K have indegree and outdegree at least two. Both ports are
therefore characteristic, and deleting them recovers exactly K. It follows that

    K≅L  iff  W(K)≅W(L),     Aut(W(K))≅Aut(K).

Every automorphism of K extends by fixing b,y; every automorphism of W(K)
restricts to one of K. This proves both directions and the group identification.
Converse exchanges the two ports, and W(K)^op≅W(K^op), so the same injectivity
holds on converse-merged classes. The extension and inverse retain the actual
core vertices, rather than just its Hamilton count.

Thus classification of even the **two-sided** irreducible cores contains the
classification of all strong tournaments via an explicit order-plus-two
embedding. Primitive cores cannot be presumed to form a simpler finite alphabet.
This does not invalidate canonical one-sided ancestry; it locates the full
classification difficulty inside its irreducible objects.

The inherited [Session150 wrapper transfer](../session150_hamiltonian_spectrum_20261103/transfer/TRANSFER.md)
also gives

    H(W(K))=3H(K)+c_2(K),

where c_2 counts **ordered covers by two nonempty paths**. Empty slots are
excluded here. The one-run Hamilton words contribute 3H(K), and words splitting
K into two runs contribute c_2(K). This is not a scalar map of H: that source
already supplies two H15 inputs with c_2=78 and90 and outputs123 and135.
Keeping Aut unchanged does not keep the frame mass H/Aut unchanged. This is
the arithmetic sidecar required when embedding the ordinary core carrier.

## Strongest survivor and repaired state

The oriented unique inverse remains canonical. The two-sided rewrite system
terminates, since each step lowers order, but its local critical pairs need
not join and its terminal objects are not unique. A deterministic tie-breaking
rule would select a normal form only after an explicit gauge is supplied;
it would not make the rewrite system confluent.

An exact two-sided classifier may retain the finite branching reduction DAG
with named deletion ports, or keep the pair of independently computed oriented
ancestries when only converse-merging is needed. The latter pair does not
claim to enumerate all two-sided reduction branches. The next question is
whether critical pairs admit a useful finite set of additional relations or
marked sidecars. The displayed family rules out resolving this simply by
identifying converse terminal cores.

## Reproduction

    python3 two_sided_peeling.py > two_sided_peeling.json

The checker uses the independent literal routines in `ancestry_audit.py`, not
the producer's ordinary-class bank. It examines all 2,24,544 labelled strong
sources of orders 3,4,5, all their legal branches, and full permutation
canonicalization of terminal types. It then checks the complete terminal set
and named branches of the six-vertex hostile. Finally it instantiates the
family at every root of C3, S4, and R5 (12 sources of orders 6,7,8), verifies
the terminal cyclic core has no legal inverse, and replays the other branch
back to exactly K. No extrapolation from these finite controls is used in the
all-order proof.

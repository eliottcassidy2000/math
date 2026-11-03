# Unique adaptive ancestry and the irreducible strong core

**PROVED all-order decomposition and orbit formulas; independent reviews
and small exact controls.** This classifies adaptive ancestry, not the
Hamiltonian-count image or all fixed-order metagraph components. No
literature-priority claim is made.

## Inheritance and definitions

Use [Session150's adaptive birth](../session150_hamiltonian_spectrum_20261103/adaptive/RESULT.md):
from a strong tournament T choose a named root v and add x with v->x and
x->u for every old u other than v. Write this operation B_v(T). It retains
all old arcs and is strong. The inherited exact growth formula proves
H(B_v(T))>H(T). Arbitrary choices of old roots are allowed.

The closest earlier unique-history mechanism is the much narrower
[Session103 five-letter unique-Hamilton-cycle grammar](../session103_words_necklaces_20261029/README.md),
which credits Shelor's classical construction. Its legal root choices
preserve a unique Hamiltonian cycle. No unique-cycle hypothesis is imposed
here. The hostile is a strong tournament with no permissible inverse;
the corrected near miss is conflating unique decoding with surjectivity
onto numerical targets. The least-used coordinate is the identity of the
last inserted vertex as a characteristic vertex of the final graph.

For a strong tournament U of order n>=4, call x a **legal last vertex** if

    indegree_U(x)=1 and U-x is strong.

Its sole in-neighbor v is necessarily the root of its inverse birth.
A strong tournament K of order at least three with no legal last vertex
is called irreducible **for this birth operation**. C3 is an irreducible
seed. This is not indecomposability under homogeneous substitution.

## A1. At most one legal last vertex

**Theorem.** Every strong tournament U of order n>3 has at most one legal
last vertex. Consequently B is injective on rooted ordinary classes:

    [(T,v)] -> [B_v(T)].

**Proof.** Let L be the vertices of indegree one. If |L|>=4, the induced
tournament on L has more than |L| edges, contradicting total indegree at
most |L|. If |L|=3, all three internal indegrees must equal one; L is a
directed triangle and no outside vertex beats any vertex in L. Since
n>3 this contradicts strongness. Thus |L|<=2.

If L={x,y}, suppose x->y. Then x is y's sole in-neighbor. Deleting x makes
y a source in a tournament of at least three vertices, so U-x is not
strong. Only y can possibly be legal. If |L|<=1 uniqueness is immediate.

Every birth from a strong parent has its newborn as a legal last vertex,
so the inverse exists on the image. An isomorphism of two children must
carry the unique last vertex to the unique last vertex and its unique
in-neighbor to the other root. Its restriction is an isomorphism of the
rooted parents. The converse is immediate by extending such an isomorphism
to the newborn. QED.

Existence is not claimed for every strong tournament. At order five, the
regular tournament has no indegree-one vertex. The strong cyclic wrapper
C3[C3,1,1] has one, but deleting it leaves C3->1, which is not strong.
Both have H15. The adaptive language still realizes H15 at order six:
failure to generate a particular class is not a missing numerical value.

## A2. Canonical core and chronological root word

Repeatedly remove the unique legal last vertex while one exists. The
process terminates, is invariant under isomorphism, and leaves a unique
induced irreducible core K. Record each deleted vertex and its unique
in-neighbor in the remaining graph. Reading this list backwards reconstructs
U exactly from K by adaptive births.

For a fixed named core K of size q, a root word of length m has at position
i, starting at zero, one of the q original core vertices or one of the i
previously born vertices. Its alphabet size is q+i. Birth indices are
chronological, not interchangeable anonymous vertices.

**Theorem.** Ordinary descendants of an irreducible core K correspond
exactly to Aut(K)-orbits of these root words. Aut(K) acts on the original
core letters and fixes every chronological birth letter. Moreover

    Aut(U) = {g in Aut(K): g fixes every core vertex used as a root}.

**Proof.** Unique peeling makes every newborn characteristic, successively
in each remaining induced graph. Thus every final automorphism fixes all
newborns and restricts to an automorphism of K. It must fix the root of
each birth; a root that was itself born is already fixed, while an original
root imposes the displayed pointwise condition. Conversely any such core
automorphism extends by fixing all newborns and preserves every birth
incidence. The same argument for an isomorphism between two words gives
exactly the stated orbit equivalence. QED.

This yields a terminating canonical description of any strong tournament:
peel, canonically identify the terminal core, and minimize the recorded
root word over the core's automorphisms. The difficult classification of
irreducible cores is retained explicitly. It is not erased by naming them.

## A3. Exact descendant counts and the triangle-seed family

Let G=Aut(K). For g in G write f(g) for its number of fixed core vertices.
Burnside's lemma and A2 give the number d_K(m) of ordinary descendants with
m births:

    d_K(m) = (1/|G|) sum_(g in G) product_(i=0)^(m-1) (f(g)+i).   (A3)

The empty product is one. A fixed word chooses a fixed core vertex or an
already born vertex independently at each position, proving the product.
Equivalently, as a formal exponential generating function,

    sum_(m>=0) d_K(m) z^m/m! = (1/|G|) sum_(g in G) (1-z)^(-f(g)).

For K=C3, nonidentity automorphisms fix no vertex. After the first birth,
its selected core root kills the remaining cyclic symmetry. Therefore
every triangle-seed descendant of order n>=4 is rigid and

    A_n=(n-1)!/6, n>=4;       A_3=1.                       (A4)

The ordinary class counts at orders4,5,6,7 are1,4,20,120. The unrestricted
root histories from one fixed labelled C3 have exactly three presentations
per ordinary class: the three initial rotations. Fixing the first root to
vertex0 removes precisely that duplication. Later old-root choices may
not be merged by their H values or by equal rooted profiles.

Two useful nonfree controls are the regular five-core, with Aut=C5 and
no nonidentity fixed vertices, giving d(1)=1,d(2)=6; and C3[C3,1,1], whose
C3 action fixes two core vertices, giving d(1)=3,d(2)=14. The latter
demonstrates why merely dividing every word count by |G| is wrong.

## A4. Labelled core counts

Let s_n count all strong tournaments on a specified n-element set and
c_n count the irreducible ones. For n>=4, A1 gives the exact disjoint
decomposition into irreducible graphs and graphs with a legal last vertex:

    c_n=s_n-n(n-1)s_(n-1).                                (A5)

Indeed choose the newborn in n ways, its root in n-1 ways, and the strong
parent on the remaining labels in s_(n-1) ways. No two choices give the
same labelled child, by unique inversion. At n=3, c_3=s_3=2. Iterating yields

    s_n=n!(n-1)! sum_(q=3)^n c_q/[q!(q-1)!].              (A6)

The small controls s_3=2,s_4=24,s_5=544 give c_4=0,c_5=64. The two order-five
cores contribute24 and40 labelled copies, from automorphism orders5 and3.
These formulas count classes/objects, never Hamiltonian paths. Their
factorials give no prime realization theorem.

## A5. A retained chronology can be a simple tree

For the triangle-seed family choose initial labels0,1,2 with
0->1->2->0. At every later chronological vertex i, retain its parent
r_i<i. The directed edges from an older vertex to a newer vertex are
exactly the tree edges

    0->1, 1->2, and r_i->i for i>=3;

every other pair points from the newer vertex to the older one. Thus a
recursive rooted tree with the specified initial chain is a faithful
labelled construction representation. It retains actual root identities,
birth times and all arcs, rather than just the shape of an unlabelled tree.
Different linear extensions or forgotten root names cannot be identified
without checking the resulting tournament. A4 specifies exactly the
remaining threefold initial gauge.

For an original triangle vertex v, or any newborn once created, its
indegree at time t is

    1 + number of later births whose chosen root was not v.       (A7)

The selected root keeps its indegree; every other old vertex gains one.
This supplies a cheap structural measure of the cost of returning to an
older root. The rooted defect in [ADAPTIVE.md](ADAPTIVE.md) is zero exactly
at indegree-one roots and is bounded below by twice the number of pairs
of in-neighbors. New-root resets erase that defect; general old-root
switching reintroduces it.

## What this changes in the proof attempt

Unique ancestry supplies a proof of injectivity and a way to avoid repeated
ordinary classes during forward generation. It supplies an inverse from
an **actual graph** to its history. It does not supply an inverse from a
target integer p to a graph. Those are distinct obligations.

The construction's H values strictly increase, so exact target-bounded
search terminates. Its failure at one target would refute coverage by this
particular seed language; success at every tested target is not an
induction. A prime/7p proof must choose an actual parent and root with the
required next count and retain the compatible root-response data for the
next stage. The sharp collision in ADAPTIVE.md prevents replacing that
data by the four numbers of just one root.

The all-strong core decomposition and the native move compiler in
[METAGRAPH.md](METAGRAPH.md) address different classification problems:
one describes birth ancestry; the other transports orientations using
native moves and an auxiliary vertex. Neither implies fixed-order native
connectivity, nor that all allowed H values occur in one ancestry tree.

# Arbitrary old-port births: exact SCC reduction and bounded all-order failures

**PROVED reduction; FINITE-EXACT complete source classification for H<125.**
The unrestricted odd spectrum is OPEN. These failures concern one declared
constructor, not the existence of tournaments with the displayed counts.

## Actual operation, not a scalar update

For a tournament T with a selected vertex x, let B_x(T) adjoin one vertex z,
with z->x and u->z for every old u!=x. Old labels and arcs are unchanged.
Let F=H(T) and define the integer port charge

    I_x(T)=H(B_x(T))-2H(T).

A port supplies a binary digit precisely when I_x=-1 or+1. After each birth
we allow **any** old or new vertex to be chosen as the next port. This is
strictly larger than Session142's keep-old/mark-fresh S,R alphabet. That
packet's all-depth no-go theorem does not cover this larger language.

Use the inherited marked-parent coordinates for A=T-x and S=N+(x):
h=H(A), q counts ordered covers of A whose second path is nonempty and
starts in S, and g counts A paths ending in S. Then

    I_x=q+g-F,       q+h>=F+g.                            (1)

The first equation is the exact inherited singleton-birth formula. The
second is its nonnegative interior guard. All states below come from actual
adjacency, not from solving(1) numerically.

## Theorem: padding never supplies a new binary port

Suppose F>=3. If I_x is+1 or-1, then T has the form

    T=P -> C,

where P is transitive, C is a nontrivial strong tournament, and x lies in C.
Moreover I_x(T)=I_x(C). Conversely every such prefix preserves C's charges.
Births at its eligible vertices preserve P as a neutral prefix, so deleting
P preserves the entire subsequent binary history with arbitrary re-marking.

**Proof.** If a nonempty SCC suffix R follows the component containing x,
write U for all components up to that component. Treat U as the old marked
source and append R, every vertex of which is beaten by x. Write h0,F0,q0,g0
for U's source coordinates, H_R=H(R), and E_R for its ordered two-path
cover count allowing empty paths. The inherited block action gives

    F=H_R F0,
    q=E_R q0+(E_R-H_R)h0,    g=H_R h0.

Since E_R>=2H_R and q0+h0>=F0, we obtain

    I_x=E_R(q0+h0)-H_R F0 >= H_R F0=F>=3.

Thus x must lie in the last component. If it is a singleton sink, I_x=-F,
also impossible. Otherwise write T=P->C with x in the last strong component.
In every cover counted by q, the second path starts in C and cannot return
to P; all of P is in the first path. Endpoint counts factor likewise. Hence
(F,q,g,h) is H(P) times the corresponding state in C, and

    I_x(T)=H(P) I_x(C).

The values+1,-1 force H(P)=1. A tournament has one Hamiltonian path iff it
is transitive: a directed triangle gives at least three by induced path
monotonicity; a transitive tournament has its unique increasing path.
Finally a birth in C leaves every vertex of P before both C and z. All P
ports remain ineligible by the suffix argument. This proves the whole-word
statement, not just a one-step count equality. QED.

This reduction preserves the actual strong source and the selected ports.
Prefix deletion subtracts its fixed number of vertices from the clock; it
does not identify different strong witnesses or change future birth costs.

## Complete finite test with an all-order consequence

The inherited Busch floor has f(10)=125 and increases with strong order.
Therefore **every** binary port on a source with 3<=H<125 is represented
by a strong tournament of order at most9. There is no hidden padding bound.

All classes at orders3..9 were read from the pinned nauty banks: counts
2,4,12,56,456,6880,191536. Exactly687 of their strong members have H<125.
Their complete labelled port vectors are in ports.json. This yields all
missing signs in that range, recorded in BINARY_SUMMARY.json.

Examples of global constructor obstructions:

- Every H13 source lacks a+1 port:27 is not obtained by a binary B birth
  from any H13 source, at any order.
- Every H19 source lacks a-1 port: the same holds for37 from19.
- Every H31 source lacks **both** signs: neither61 nor63 can be obtained
  by this operation from any H31 source, even after arbitrary re-marking
  or neutral SCC padding.

All27,37,61,63 are attained tournament counts in the inherited census. These
are constructor failures, not additional forbidden Hamiltonian-path values.
The H3/+ and H11/- failures recover the genuine forbidden values7 and21.

## A stronger finite hostile for recursive two-sided sources

Call T two-sided when it has both signs. Among all classes through9, the
numbers of two-sided classes at orders4..9 are

    1,2,3,8,16,30.

For **none** of them can one choose a minus birth and a plus birth such
that both children are again two-sided. This checks all ports, not one
representative of each sign. Source labels, new vertices, and full child
adjacency remain in the test. It is a complete finite obstruction, not an
all-order theorem about two-sided sources.

There are actual positive controls: the H5 diamond has a minus birth to a
two-sided H9 source; its plus birth gives H11, which cannot have a minus
port because21 is forbidden. Higher two-sided counts37,57,65,... occur,
but the recorded branch test still fails. Thus the larger language really
adds choices while the most direct closed-family proposal does not survive.

The next repair must permit more general ears, richer multi-vertex attachments,
or count-preserving changes to the strong source—not merely a new name
for the same port. A set of witnesses per count is more flexible than one
two-sided witness, but the all-order H13/H19/H31 failures show that even
this set-valued repair is insufficient for the present birth alone.

Dependencies: Session142 `constructions/BINARY_NO_GO.md` and
`minimal_gap/PROOF.md`; Session134 `constructor/PROOF.md`; the current
`spectrum_descent_jacobian_20261102/DESCENT.md`; source hashes are retained.

The later TWO_BIRTH_THEOREM.md strengthens the word boundary analytically:
arbitrarily many B births cannot repair any missing binary transition.
INDUCED_OBSTRUCTIONS.md then closes selected branches even for arbitrary
old-edge-preserving extensions. Its47->73->117 witness also proves that
count-level arrows cannot be composed without a compatible middle source.

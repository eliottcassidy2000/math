# Endpoint-four repairs through movable cuts and ordered edge syndromes

21 October 2026. **PROVED structural reductions, with a FINITE-EXACT
24-system boundary certificate.** Independent audit is required before
canon promotion. No priority claim, full endpoint-four theorem, or global
black-core connectivity claim is made.

## Result and remaining boundary

Let P=(a,b,x,...,y,c,z) be a closed blue Hamiltonian word of a tournament
T on n vertices. A word is blue when its actual positional reversal is
an anti-automorphism, and black otherwise. Put K_Q(T)=T^op xor E(Q),
where E(Q) is the set of unordered consecutive pairs. Assume

    n>=15,    d_T(a)=n-5,    d_T(z)=4.

Use P positions 0,...,n-1 as labels, so a=0, b=1, x=2, y=n-3,
c=n-2, z=n-1. There is a three-element middle incoming set I with

    N^-_T(a)={z} union I,
    N^+_T(z)={a} union J,    J={n-1-i : i in I}.

The path condition excludes b from I. Overlap between I and J is allowed.
The middle is open when b->c, and closed when c->b.

**Theorem.** There is an exact three-black K detour from T to the literal
labelled target K_P(T) in every one of the following cases:

1. The middle is open, with no further restriction on I.
2. The middle is closed and c is not in I, with no further restriction.
3. The middle is closed, c is in I and y is not in I.
4. The middle is closed, I={y,c,h}, x->c, and h is neither x nor n-4.

Thus the only part of the n>=15 endpoint-four layer left open by this
packet is the closed middle with I={y,c,h}, where either c->x, or x->c
and h is x or n-4. A failed template is not asserted to be a bridge.
Orders below 15 are not classified here. The companion order-ten packet
closes the endpoint-three n=10 boundary; the inherited n>=11 theorem
and that finite certificate concern a different endpoint layer.

The result is unfiltered: there is no strong-target, mixedness, module,
or isomorphism-class condition. Every proof retains the physical ports,
the actual word at each intermediate tournament, and its path-edge XOR.

## Inheritance and live concepts

The closest proved mechanisms are endpoint_three_repairs_20261020h:
root.md R1--R4, theory.md sections 2--9, lift.md, and the open end-pair
ledger in probe_new_ledgers.md. The current endpoint-three theorem is
its README, not an exploratory sidecar. In particular, its degree-free
movable-cut corollary already handles case 3 above; that is inherited,
not a new theorem of this packet.

The canonical hostile is the order-six endpoint-three source whose target
is C3->C3 and has no black words. The corrected near miss is lifting a
singleton's positional reversal through a larger block without recomputing
it. The least-used sidecar here is the affine parity syndrome of two
specified anti-reversals, projected onto the actual endpoint row.

The concept board is: ordered marked graphs; physical exceptional ports;
movable cuts; exact path-edge XOR; permutation-orbit capacity; affine
endpoint-row projection. Anchor: the next endpoint layer. Niche: turn
three exceptional ports into a small parity-budget contradiction.
Wildcard: whether the marking-aware trijection supplies a genuine legal
operation compiler rather than an unlabelled connectivity analogy.

Sources are pinned in SOURCE-MANIFEST.json to repository commit
c4a7875652a5ea516e93a31c3c6410d2204773b2. No external literature claim is
used. The metadata/correction, quotient, and actual-consequence checks of
the repository protocol were applied before building these arguments.

## 1. The ordered-graph trijection as a legal operation compiler

In P positions, extend the backward-bit array x_uv to all unordered pairs
u<v by setting x_(i,i+1)=0. The free bit x_uv for v>=u+2 is the graph edge
{u,v-1} on ordered {0,...,n-2}. This is the inherited exact trijection:
marked tournaments, ordered graphs, and ordered tournaments of one lower
order. P-blue is reflection invariance of this graph. K_P complements
its entire free host. The order and marked word are indispensable.

For a candidate word q and a source reached after t K moves, let
D be the XOR of the preceding t path-edge sets. For i<j, the backward
bit in q coordinates is exactly

    x_{min(q_i,q_j),max(q_i,q_j)}
      xor (t mod 2) xor [ {q_i,q_j} in D ] xor [q_i>q_j].       (1)

The same formula at j=i+1 must vanish for q to be directed. At every
nonadjacent pair it is the new graph coordinate. Equality of coordinates
at (i,j) and (n-1-j,n-1-i) is exactly the actual blue predicate for q.
All assertions here follow pair by pair, not from quotient connectivity.

Consequently, a proposed ledger supplies affine equations for directedness
and for the failure of each blackness claim. The endpoint condition
|I|=3 is a cardinality condition on one marked row, not merely a score
multiset. Equation (1) retains the bit needed to distinguish a chart
change from a tournament operation.

For any three directed words, the equality

    E(Q1) xor E(Q2) xor E(Q3)=E(P)                              (2)

implies that three K moves end exactly at K_P(T): three full complements
leave one, and the retained path pairs add by XOR. This is a labelled
identity. It does not identify, quotient, or lift any abstract graph path.
Blue/black status is the same at the two endpoints of a word's own K
move, because both full converse and E(Q) respect its reversal.

## 2. A four-port obstruction for the relocated-b final word

Let

    F=(a,x,3,4,...,y,b,c,z).

Suppose it is directed at B=K_P(T). P is still blue at B. If F were
blue, g=r_F r_P would be an automorphism of B fixing a and z. Its
middle cycles have lengths n-2 for odd n, and two cycles of length
(n-2)/2 for even n. This follows from the explicit action

    g(1)=2, g(i)=i+2 for 2<=i<=n-5,
    g(n-4)=1, g(n-3)=n-2, g(n-2)=3,
    g(a)=a, g(z)=z.

The middle outneighbours of a in B are exactly {b} union I, of size
four. Since g fixes a, this set would be a union of its middle orbits.
At odd n>=7, the single orbit has odd size n-2>=5. At even n>=8,
write m=(n-2)/2>=3. If m is even, g cannot be a tournament automorphism:
an even-order permutation has a nontrivial involution power, and an
involution reverses the arc between a transposed pair. If m is odd,
the available nonzero union sizes are m and 2m, neither four. Hence
F is black for every n>=7 in this endpoint-four layer whenever directed.

The bound and proof are specific to the four-element retained port set;
they do not erase the marked endpoint or assert a generic score criterion.

## 3. A new movable cut for a second-position port

Assume the middle is closed, x->c, x in I, and y not in I. Since y not
in I, x->z. Since x is in I, z->y. Choose k to be the first middle
outneighbour of z after x. Then k>=3, k<=y, and its predecessor beats
z. Put L=(x,...,k-1), R=(k,...,y), both nonempty. Use

    Q1=(a,b,L,z,R,c),
    Q2=(a,L,z,R,b,c),
    Q3=(a,L,R,b,c,z).                                         (3)

Q1 follows the path except for predecessor->z and z->k. At the first
target, a->x is the complement of x->a, y->b is the complement of
b->y (the latter equals x->c by P-blue), and b->c is the complement
of c->b. These are exactly the new joins of Q2. The omitted pairs
predecessor,k and c,z are complemented twice and restored for Q3.
All remaining joins are retained. The added pairs predecessor,z; z,k;
a,x; y,b; and b,c cancel in pairs, leaving (2).

Q1 reversal sends c to a and b to y. It would equate true c->b with
false y->a. Q2 reversal sends a to c and b to x. It would equate the
retained true a->b with the complemented false x->c. Both words are
black. Q3 is the word F of section 2, which is black. This proves a
multi-port movable-cut theorem rather than a fixed-position guard.

The directedness and first two blackness proofs do not use |I|=3.
For other endpoint layers, only the final four-port argument must be
replaced. This pinpoints the remaining degree-dependent obligation.

### When both x and y are exceptional

Now let I={x,y,h}, still with c not in I and x->c. If h!=3, use (3)
with k=y. Its predecessor beats z because 3 is not in I. Q1 reversal
would equate a->z (false) and x->c (true), so it is black. Q2 and Q3
have the same proofs as before.

If h=3, choose k=n-4. The two relevant ports of z are n-4 and y;
its predecessor n-5 beats z at n>=9. The Q1 permutation moves the
last vertex of P into the start of its last four positions. If both
P and Q1 were blue, their reflection product would consist of two
disjoint four-cycles on the first and last four positions, with all
remaining vertices fixed. This has even order and is impossible for
a tournament automorphism. Q2 and Q3 are again black as before.

This completes every closed-middle case with c not in I, x in I,
and x->c, without dropping the consecutive-port boundary.

## 4. Repairing an open adjacent end-run

Suppose the middle is open and I={y,c,h}. Let t=y unless h=y-1, in
which case t=y-1. Put L=(x,...,t-1), R=(t,...,y). Thus L is nonempty,
a beats its last vertex, and the first vertex of R beats a. Use

    Q1=(a,b,c,z,L,R),
    Q2=(b,c,z,L,a,R),
    Q3=(b,L,a,R,c,z).                                         (4)

The new joins of Q2 are supplied by complementing the two stated
incidences with a. The original pairs b,x and y,c are complemented
twice and restored for Q3. The pairs b,c; z,x; last(L),a; and a,first(R)
cancel in pairs; all other surviving pairs are exactly E(P). This
proves directedness and the exact target.

The first word is the inherited penultimate block word. Its reflection
product with P has one n-cycle when n is odd; two n/2-cycles when
n=2 mod 4; and four n/4-cycles when n=0 mod 4. In the odd case, blue
would force regularity and n=9. For n=2 mod 4, the a-orbit excludes
z and y and contains at most the two other incoming ports; regularity
forces n<=10. For n=0 mod 4, the a-orbit is
{0,2,6,10,...,n-6}, excluding z,y,c, so it has at most one incoming
port and regularity forces n<=12. Thus Q1 is black at n>=15.

Q2 retains a->b. If |R|=1, its reversal would equate this with y->c,
which was complemented by Q1 and is false. If |R|=2, its reversal
instead equates it with y->z. Here h=y-1, hence x is not in I; y->z
was true originally and was complemented by Q1, so is false. Q2 is
black in both cases. Q3 reversal would equate b->c, retained twice
and true, with x->z, retained twice and false. It is black as well.

The moving insertion point is essential precisely at the new port
h=y-1. Keeping the old endpoint-three insertion would make its new
join point in the wrong direction.

## 5. Two endpoint-row budget lemmas

There are two inherited fixed ledgers whose second words need stronger
blackness proofs when a third exceptional port is allowed. Write
W=(x,...,y), D=(3,...,y).

Open penultimate branch, y in I and c not in I:

    Q1=(a,b,c,z,x,D),
    Q2=(b,c,a,z,x,D),
    Q3=(b,x,D,c,a,z).                                         (5)

Closed branch, c,x not in I and x->c:

    Q1=(c,b,z,a,W),
    Q2=(z,c,b,a,W),
    Q3=(a,z,b,W,c).                                           (6)

Their directedness and XOR proofs are unchanged by additional ports:
(5) needs only b->c, z->x, a->c; (6) needs c->b, b->z, z->a,
a->x. The reflected incidences supply precisely these joins.

Put X_v=T(a,v), a forward arc bit. At S=K_Q1(T), suppose Q2 is blue.
Combining the P and Q2 anti equations gives the following constraints
for every n>=27:

- For (5): X_4 xor X_(n-6)=1, X_8 xor X_(n-10)=1,
  X_12 xor X_(n-14)=1.
- For (6): X_5 xor X_(n-3)=1, X_9 xor X_(n-7)=1,
  X_13 xor X_(n-11)=1.

These are not extrapolations from a search. In the ordered-pair bit
coordinates, use the following chains, alternating P-blue and
Q2-blue equations; the total right-hand side of each chain is one.
The first and last pairs are the two endpoint-row variables.

For (5):

    (0,4), (n-5,n-1), (0,n-6)
    (0,8), (n-9,n-1), (4,n-6), (5,n-5), (0,n-10)
    (0,12), (n-13,n-1), (8,n-6), (5,n-9),
      (4,n-10), (9,n-5), (0,n-14).

For (6):

    (0,5), (n-6,n-1), (0,n-3)
    (0,9), (n-10,n-1), (5,n-3), (2,n-6), (0,n-7)
    (0,13), (n-14,n-1), (9,n-3), (2,n-10),
      (5,n-7), (6,n-6), (0,n-11).

Equation (1) checks every link: the first Q2 link in each chain has
right-hand side one and all other links zero. All displayed pairs are
properly ordered and distinct where required at n>=27. In each branch
the three endpoint pairs are disjoint, so each consumes at least one
incoming port.

For (5), y is an additional incoming port outside the six endpoints.
For (6), vertex 4 is an additional forced incoming port: Q2-blue maps
the true S(y,b) to S(n-5,z). But N^-_S(z)={b} union J, and n-5!=b;
hence n-5 belongs to J, or 4 belongs to I. Again it lies outside the
six endpoints. Either branch would therefore require at least four
middle incoming ports, contradicting |I|=3. Q2 is black.

### The exact finite boundary 15 through 26

Only 24 affine systems are needed, one for each n and each ledger.
Use one bit per unordered pair, all P-blue equations, the P path,
z->a, the stated branch incidences, and the proposed Q2-blue equations
at S. Do not fix I and do not enumerate tournaments.

Each equation is x_u xor x_v=s, with a fixed zero node for constants.
In a consistent component, each variable is a component-root bit plus
a known parity. Project only the row T(a,v), v=2,...,n-2. A component
fixed by the zero node contributes exactly its forced zeros to |I|.
A free component with p projected variables of parity zero and q of
parity one contributes at least min(p,q). Sum across components. This
is the exact minimum |I| in that affine family, since root choices
are independent. A contradictory component has no tournament at all.

The checker implements this both with parity union-find and independently
with breadth-first propagation. Every one of the 24 systems is either
inconsistent or has minimum |I|>=4. In fact the consistent systems have
minimum 6,8, or 10, as recorded verbatim in receipt-normal.json.
This is a complete FINITE-EXACT boundary proof, not an all-n inference.
The three explicit chains supply the unbounded tail.

## 6. Completing the stated partition

Here are the remaining first/last-word arguments, so no check is replaced
by a path-edge identity alone.

Open middle with y not in I: use inherited root.md R2 at k=4. Its
cyclic first word is black for n>4k-2=14; the second has endpoint
score sum n-2; the last has a nontrivial involutory reflection product.
This gives all three black words and the literal target.

For (5), Q1 has the orbit structure used in section 4. The odd case
forces n=9. At n=2 mod 4 the a-orbit excludes z,y and has at most
two incoming ports, forcing n<=10. At n=0 mod 4 it has size n/4 and
at most two incoming ports, so n<=20. The only surviving n>=15
would be n=20, since n=16 gives an even orbit. But Q1-blue also
requires n-5 in I: its reversal compares a->(n-5) with the false
path-opposite c->y. At n=20 this extra port lies outside the a-orbit,
leaving at most one incoming port there, a contradiction. Q1 is black.
The last word of (5), tested at B, would equate true a->b with false
z->x. Thus it is black. Section 5 handles Q2.

Closed middle with c not in I and c->x: use

    Q1=(c,W,b,z,a), Q2=(z,c,W,b,a), Q3=(a,z,b,W,c).              (7)

If Q1 were blue, its endpoint equality forces d_T(c)=4. The mandatory
outneighbours b,x,z then force x,3 into I under Q1 reversal. Its
reflection product with P has one n-cycle at odd n, forcing n=9,
or two n/2-cycles at even n. The a-orbit is {0,1,2,4,...,n-4};
only x and possibly the third port are incoming there, so its size is
at most five and n<=10. Hence Q1 is black for n>=15. At its target,
z and a have degree n-5 and 3, summing to n-2, so Q2 is black.
For Q3 use inherited theory.md section 2 at k=4: the common final word
is black whenever n>4k-2. Its directedness and XOR are those of (7).

For (6), the first-word reflection product has even order at every
even n. At odd n it is an n-cycle, which would force n=9. Thus Q1
is black at n>=15. The common final word has the same bound as (7).
Section 5 supplies Q2. Section 3 supplies all remaining c-not-in-I
cases with x->c.

Closed middle with c in I, y not in I: use the inherited degree-free
movable-cut theorem in endpoint_three theory.md section 9. Choose
the first outgoing port of z after b, not a fixed chosen port. Its
three words are (a,L,z,R,c), (L,z,R,a,c), (L,R,a,c,z), with the cut
before that port. The arc witnesses and the n-2 score sum do not depend
on the number of extra ports. All hypotheses hold in the present case.

Finally, if I={y,c,h}, x->c, h not in {x,n-4}, the inherited two-block
repair works with the third port retained. With D=(3,...,y), its K words
are

    (a,x,D,c,z,b),
    reverse(a,b,x,c,reverse(D),z),
    (b,z,D,c,x,a).                                           (8)

The guards are exactly a->x and 3->z, in addition to x->c. The
cyclic first word is black above 14, the middle word has endpoint
sum n-2 (equivalently use its R-form), and the final anti test compares
unchanged z->a with the reversed path arc b->x at R_P(T). The two-block
XOR proof gives the literal endpoint. This is a guarded extension,
not permission to ignore h=x or h=n-4.

## 7. Verification and hostile boundaries

Run from any directory:

    python3 verify.py > receipt-normal.json
    python3 -O verify.py > receipt-optimized.json
    cmp receipt-normal.json receipt-optimized.json

The checker imports no inherited producer. It constructs complete
reflected tournaments directly as sets of directed arcs, checks words
in their actual sources, computes every K operation, and compares the
literal final arc set. The independent graph-coordinate checks retain
the marked path and validate reflection/complement rather than graph
isomorphism. The finite boundary is checked by two different component
algorithms. Explicit port chains and permutation formulas are controlled
through n=200; their written proofs, not this finite range, give the
all-order statements.

Small exhaustive inputs are every marked endpoint-four source at n7..9;
only the stated new-template subfamilies with n>=8 are replayed there.
Other small templates are explicitly counted as below threshold, not
claimed. Larger inputs are deterministic samples at every n15..60,
including adjacent-port, low-position-port, and tail-run hostiles. Every
claimed larger template must be directed and black at all three steps.
The residual closed end-pair boundary is counted as unclaimed rather
than silently filtered out of the theorem.

Two explicit negative controls are included: adding h=y-1 to the old
open end-pair ledger makes its last(L)->a join illegal, and adding a
third consecutive port to I={2,3} makes the old fixed cut point at an
actual outgoing port of z. The repaired movable cuts are replayed on
these same labelled sources; their complete rows are in the receipt.
These are failures of fixed templates, not connectivity counterexamples.

The source/correction frontier was read before derivation; the repository
and remote were not modified by this research lane. SOURCE-MANIFEST.json
pins the inherited mathematical sources. MANIFEST.sha256 pins the
publication payload. The ordered-graph map proves a legal-state compiler
and finite affine certificates; no quotient metagraph path is assumed
to lift, and no global black-component connectivity conclusion follows.

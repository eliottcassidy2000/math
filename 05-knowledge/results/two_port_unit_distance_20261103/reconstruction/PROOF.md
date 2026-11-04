# Deletion contraction and a half-order reconstruction bound

**Status: PROVED** by the exact identities below, subject to the packet's
independent audit. Universal reconstruction from TWO ports is **OPEN** here;
no exact two-port collision is asserted. No historical priority or efficient
evaluation of Hamiltonian counts is claimed.

## 1. Recovered target and retained information

Let T be a tournament on a fixed named vertex set V, |V|=n. Attach m ordered
base vertices z_1->...->z_m, with every core/base edge chosen independently.
The complete m-port response R_m(T) is the resulting Hamiltonian vertex-word
count for each such named attachment. Equality always aligns every old label,
every port, and every attachment bit. It does not mean equality after separately
forgetting the labels of different observations.

Session27 already establishes reconstruction with n-1 ports and an exact
two-port tomography of H, endpoint counts, and the ordered two-path-cover
matrix C. The latter matrix alone is not the entire two-port response:
ordered three-path covers can also contribute. Sources:
[faithfulness](../../session27_repunits_lattices_20261022/session27-response/faithfulness-proof.md),
[modular refinement](../../session27_repunits_lattices_20261022/session27-response/MODULAR_FAITHFULNESS.md),
and [response proof sections 7-8](../../session27_repunits_lattices_20261022/session27-response/proof.md).

The immediately inherited [positive-circulation theorem C4](../../proof_principles_20261103/continuation/PROOF.md)
says that equality of the complete ONE-port response forces the changed arcs,
oriented as in the first tournament, to carry a positive balanced circulation.
Its proof uses Session176's positive cut weights. The known opposite-C3
collision, the 81/87 strong ambient witness, and the regular-seven 403/407
hostile are retained; none supplies a two-port collision.

The new exact map below extracts a full one-port response of **every named
vertex deletion** from the full two-port response. Iterating the same map
yields a new all-order reconstruction bound and stringent support conditions
on any two-port counterexample. These deductions preserve all endpoint data;
they do not replace it by degree scores or a formal signed flow.

## 2. Two varied edges delete one old vertex

Fix v in V, and two adjacent ports x,y with x->y. Give x and y identical
relations to every old vertex other than v. They also have identical relations
to all other ports because they are adjacent in the transitive base order.
Choose signs

    alpha=+1 iff x->v,  beta=+1 iff v->y,

and independently vary these two edges through their four assignments.
Let T_(alpha,beta) denote the resulting tournament. Contract the formal block
x,v,y to a single port z, using x/y's common relations to the remaining
vertices. The contracted tournament U has core T-v and m-1 ordered ports.

**Theorem R1 (exact deletion contraction).** For every such attachment,

    sum_(alpha,beta in {-1,+1}) alpha*beta*H(T_(alpha,beta))
      = 2 H(U).                                          (1)

In particular R_m(T) determines R_(m-1)(T-v) for every named v whenever m>=2.

**Proof.** Expand each H into Hamiltonian words. A word omitting one of
the two varied undirected edges cancels on summation over that unused sign.
Every surviving word contains both edges xv and vy. The old vertex v
therefore has its two path incidences occupied, and the consecutive block is
either x,v,y or y,v,x. The first order requires alpha=beta=+1; the second
requires alpha=beta=-1. In either case the character alpha*beta equals +1.
Contracting the block gives a Hamiltonian word of U: its first/last port
vertex has exactly the same exterior adjacency as z. Conversely each
Hamiltonian word of U expands uniquely in each of the two block orders,
contributing twice. Old edges incident to v, and the edge xy, are unused in
every surviving word. There are no other words. QED.

The factor two is exact over integers; no division in a residue ring is
being assumed. Equation (1) uses four actual tournament queries, not a
claimed positively weighted realization of their signed sum.

### Exact depth, padding and repeated deletion

Knowing exactly m ports also supplies every smaller depth: prepend a chain
of universal-source ports, whose forced initial Hamiltonian segment can
be deleted without changing H. Thus if R_m(T)=R_m(T'), repeated application
of (1) proves

    R_(m-|S|)(T-S)=R_(m-|S|)(T'-S)
    for every S subset V with |S|<=m-1.                (2)

The required ordered ports can always be chosen adjacent. After each
contraction, the remaining ports retain their transitive order. Equation
(2) is a statement about full response functions, not merely the unextended
Hamiltonian counts of the deletion cards.

## 3. Counterexamples require many changed neighbors

For T,T' on V let D be the directed graph consisting of exactly the
unordered pairs on which they differ, oriented as in T. It is a partial
orientation of a simple graph: there are no loops or opposite arcs.

**Theorem R2.** If R_m(T)=R_m(T') and T!=T', then every vertex incident to
D has at least m outgoing arcs and at least m incoming arcs in D.
Every nontrivial weak component of D therefore has at least 2m+1 vertices.

**Proof.** Pad down to one port. The inherited positive-circulation theorem
makes every incident vertex of D have both an incoming and an outgoing arc.
Suppose a vertex u has d<m outgoing neighbors in D. Delete their entire
set S. Equation (2), then padding, gives equal one-port responses of T-S
and T'-S. In their difference graph D-S, u has no outgoing arcs, but retains
all its incoming arcs, since no vertex can be both an incoming and outgoing
neighbor of u. The one-port positive-circulation theorem now requires a
sum of positive incoming weights to equal zero, a contradiction. Reverse
all directions for the incoming bound. The m incoming and m outgoing
neighbors are disjoint, yielding 2m+1 vertices in u's component. QED.

A stronger version retains the whole circulation obstruction: after deleting
any at most m-1 old vertices, every remaining changed arc must lie on a
directed cycle in D-S. It is not enough that D initially be strongly
connected or have balanced ordinary degree; the actual positive weights are
different on different deletion cards.

### Odd-order converse cannot survive ordered-cover tomography

**Lemma R3a.** At every odd order n>=3, the ordered two-path-cover matrix
C of a tournament cannot be symmetric. Consequently T and T^op have
different two-port responses at every such order.

**Proof.** Every Hamilton path has n-1 arcs, and a cover counted by C_uv
concatenates to a Hamilton path precisely when u->v is an actual arc.
Thus sum_(u->v) C_uv=(n-1)H. If C were symmetric, the total number q_2
of ordered nonempty two-path covers would be

    q_2=sum_(u!=v) C_uv=2(n-1)H.

Put H(empty)=1. Summing over every subset S gives

    sum_(S subset V) H(T[S])H(T[V\\S])=2H+q_2=2nH.       (3b)

Each nonempty tournament has odd H by the inherited tournament parity
theorem, and H(empty)=1 is odd too. Complementation groups the left side
into 2^(n-1) equal-product pairs, each contributing 2 modulo4. For n>=2
that is an even number of pairs, so the left side is divisible by4.
For odd n the right side is 2 modulo4, a contradiction. Finally converse
transposes C, and Session27's actual two-port tomography recovers C.
QED. This parity argument was independently supplied by the packet's
rank-audit lane; it is not inferred from finite reconstruction data.

**Corollary R3 (half-order reconstruction).** For n>=4,

    R_floor(n/2)(T)=R_floor(n/2)(T') iff T=T'.          (3)

For even n=2m, R2 would require a difference component with at least
2m+1>n vertices. For odd n=2m+1 with m>=2, R2 forces every changed vertex
to be joined by changed arcs to all other vertices; hence D is the entire
tournament, and T'=T^op. Since m>=2, R3a rules this out. At n=3 two ports
suffice by R2; at n=2 one suffices. Orders zero and one have a unique
labelled tournament and require no queries. This improves the recovered
n-1-port upper bound for n>=4; it does not claim the optimal depth.

### A fixed query battery and an exact decoder

For n>=4 set m=floor(n/2). Choose every S subset V with |S|<=m-1. On T-S,
recover the one-port values at masks of weight zero, one and two. The
one-port response is quadratic, so these values determine it completely.
For each desired value, apply (1) backwards, splitting a port around each
deleted vertex in a predetermined order. Each split uses four signed
queries and adds one port. The resulting 4^|S| queries use |S|+1<=m
ports; pad to exactly m with universal sources if desired.

Add Session27's fixed two-port tomography battery recovering C, of size
1+2n+3binom(n,2), padded to m ports. (It may repeat queries from the deletion
battery; the bound need not remove duplicates.) Thus a sufficient fixed
battery has at most

    sum_(s=0)^(m-1) binom(n,s) 4^s
        [1+(n-s)+binom(n-s,2)]
      + 1+2n+3binom(n,2)                               (3a)

queries. This overcounts repeated contexts. All masks, split orders, and
padding conventions can be chosen from the names and n alone, independently
of every unknown core arc. The signed sums divide by 2^s exactly because
they arise from actual tournament observations. Equality of the recovered
one-port deletion tables is enough for R2's argument; the ordered-cover
tomography eliminates the remaining odd-order converse case. Hence this
battery separates every pair of labelled cores at the stated m. For even
n, the added tomography may be omitted because R2 alone suffices.

An exact decoder is therefore available: enumerate all 2^binom(n,2)
candidate cores and retain those matching the recovered tables. There is
exactly one. This is a constructive information theorem, not an efficient
enumeration algorithm; evaluating H, large integers, and the exponential
candidate search all retain their actual costs.

For the requested TWO-port problem, R2 says every changed vertex needs at
least two changed inneighbors and two changed outneighbors, and every
changed arc must remain on a directed cycle after deleting any one vertex.
In particular an isolated triangle reversal is impossible at every ambient
order. R3 proves two-port reconstruction through order five without an
exhaustive census. Combining with Session27's finite complete two-port check
through order six, a global counterexample must have n>=7. The all-order
bound (3) does not settle whether two ports already suffice for every n.

**Current finite extension:** the parent [complete sidecar census](../census/RESULT.md)
now proves labelled two-port faithfulness through order9, including every
sidecar automorphism. The current counterexample lower order is therefore10;
the order-six comparison above identifies the inherited starting point.

### The graph obstruction is not sufficient

Let C5 be the cyclic tournament with outgoing masks (6,12,24,17,3), and
compare its converse. Their difference graph has two incoming and two
outgoing neighbors at every vertex; every one-vertex deletion is strong.
They have identical complete one-port responses. Yet attaching x with
outmask 1 and then y with outmask 2 (with x->y) gives **89 versus 91**.
The finite control checks this directly. Thus satisfying R2's graph
conditions does not certify an actual two-port collision.

## 4. Exact finite compression of the full TWO-port response

Write the two attachment masks as a,b on V, with x->y fixed. Every
Hamiltonian word uses at most two edges incident to x and at most two
incident to y. Its indicator, and hence the response polynomial F_T(a,b),
has Boolean multilinear degree at most two in each mask separately:

    F_T(a,b)=sum_(|A|<=2, |B|<=2) c_(A,B) a_A b_B.     (4)

Here a_A=product_(u in A) a_u and likewise for b. The formula includes
the endpoint, adjacent-port, and separated-port words, including every
orphan after deletion. The degree bound is not a truncation assumption.

Let L_n={A subset V: |A|<=2}. All coefficients are determined by

    c_(A,B)=sum_(I subset A, J subset B)
      (-1)^(|A|-|I|+|B|-|J|) F_T(1_I,1_J).            (5)

Therefore the entire full two-port response is exactly determined by
(1+n+binom(n,2))^2 actual queries. Equality of these entries is equivalent
to equality on all 4^n mask pairs. This is an exact oracle-query upper
bound, not a polynomial-time algorithm for computing the queried H values.

The old-old pair weights after adding the first probe contain degree-two
functions of a, hence genuine three-path-cover information. Restricting
to the smaller tomography of H,S,E,C is not licensed by (4). Conversely,
the weights involving the new probe yield only its already known ordered
two-path-cover current:

    2w_(x,u)(b)=B_u+E_u+sum_v C_(v,u)
                +sum_v b_v (C_(u,v)-C_(v,u)),          (6)

with b_v=1 iff x->v. To prove (6), in a cover ending at x the first path
is either the singleton x or an old path ending at v followed by x;
in a cover starting at x the dual cases apply. Adding both ordered cover
counts gives (6). It explains why simply rereading these new-old weights
cannot recover the missing orientation by itself.

## 5. Scope and next decisive obligation

Source: full named attachment response. Exact map: a two-edge character
contracts an old vertex between adjacent exterior ports. Preserved
predicate: every Hamiltonian count of the resulting smaller response.
Lost information: one unit of available probe depth. Restoring the depth
requires an additional port; it cannot be assumed free.

R1-R3 are proved without a new exhaustive census. The remaining universal
two-port obligation is precise: rule out or exhibit a pair with the same
complete bidegree-(2,2) response and the robust positive-circulation
constraints on every single deletion. Neither ordinary degree balance,
the first derivative matrix C alone, nor the n-1-port theorem closes it.

Reproduction and the exact finite universes are recorded in `check.py` and
its normal/optimized receipts. Those controls corroborate the identities;
they do not replace the all-order proofs or assert global two-port injectivity.

Run `python3 check.py` and `python3 -O check.py`; outputs are byte-identical.
The standard-library checker covers all 75 labelled cores of orders 1-4:
16,932 full two-port entries and their bidegree reconstructions, 75 full
ordered-cover tomographies, 2,153 deletion contractions, and 4,306 new-old
weight identities. It adds all 1,604 adjacent-pair contractions for three
ports on all cores of orders 1-3, 120 repeated two-deletion contractions
using 1,920 actual three-port queries, and the parity identities for every
labelled core of orders 2-5. Independent literal-permutation checks cover
548 two-port and 272 three-port counts. The C5 hostile includes actual
one-deletion strong-connectivity checks. Checks use explicit exceptions,
remain active under `-O`, and bind the executable source hash before and
after computation. No randomized universe or floating-point step is used.

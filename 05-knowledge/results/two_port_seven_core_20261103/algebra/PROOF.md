# Full-response content and ambient converse descent

**PROVED** all-order statements B1-B5 below. **FINITE-EXACT** controls are
separate. Universal two-port reconstruction remains **OPEN**; no exact
collision or literature-priority claim is made here.

The [independent B1-B5 review](../recovery/ALGEBRA_REVIEW.md) accepts these
proofs and replays the complete checker; its separate receipt pins that
review's executable dependencies. The later complete modular sieve and
its [independent review](../residual/INDEPENDENT_REVIEW.md) exclude the
combined B3-B4 residual for every even core of size at most10. Together
with the all-order odd-module theorem, reversing one nontrivial module
of size at most11 is therefore always visible in every ambient quotient.
Any remaining module-converse candidate has even, nonselfconverse core
size at least12. This does not assert arbitrary order-eleven reconstruction.

## Inheritance, observer, and live concepts

The closest mechanisms are [A1-A4, odd-module current and full three-path
transfer](../../two_port_beyond_nine_20261103/algebra/PROOF.md), and
[P1-P4, squarefree parity and even self-converse separation](../../two_port_beyond_nine_20261103/audit/PARITY.md).
The hostile is the eight-vertex nonselfconverse tournament with equal start
and end parity. The corrected near miss is that full response functions
must have relatively prime integer values. The least-used sidecars are
pointwise automorphism stabilizers and fixed-slot endpoint path covers.

The concept board is full attachment polynomial, stabilizer arithmetic,
uniform module runs, endpoint parity, ordered cover skew, and physical
multi-in-port defect. The anchor is universal two-port reconstruction; the
niche is arithmetic content of the entire response; the wildcard is the
exact zero-ready boundary of the recent multi-in-port grammar.

All tournaments have named vertices, with the actual arc u->v as the
intrinsic binary relation, no ties, and no quotient of names. For m ports,
fix their transitive orientation x_1->...->x_m. Let F_T^(m)(a) be the number
of Hamilton paths after attaching them, where a_(i,v)=1 means v->x_i.
All old-port bits vary independently. Set H(empty)=1. Two ports also allow
the opposite exterior arc by exchanging the two port names.

The complete two-port response is retained, not replaced by its H,S,E,C
sidecar. In the squarefree old-vertex ring, with G=(I-ZA)^(-1)Z, e=1^t G,
s=G1, a_v=[v->x], c_v=[v->y], b=1-a and d=1-c, A4 gives

    F_T(a,c)=[V]{(1+e a)(1+b^t G c)(1+d^t s)
                  +(1+e c)(d^t G a)(1+b^t s)}.       (1)

In particular D_(i,j,k,l)=[V]e_i G_jk s_l retains the three nonempty paths
with end(P)=i, start(Q)=j, end(Q)=k, start(R)=l. The case j=k is allowed;
all supports must be disjoint. None of our arguments discards this tensor
or asserts that its individual entries are all separately observable.

## B1. A coefficient sees only finitely many symmetry-breaking labels

**Theorem B1.** The unique Boolean multilinear polynomial F_T^(m) has
total degree at most 2m. Let G_0 be any subgroup of Aut(T). A coefficient
whose old-port variables mention the set S of old labels is divisible by
the order of the pointwise stabilizer (G_0)_(S).

**Proof.** Expand by literal Hamilton words on all old and new names. Each
word contributes either zero or a product of factors a_(i,v) and
1-a_(i,v), one for each old-port edge it uses. Each new vertex has at most
two neighbors in a word, so the degree is at most 2m. A word never uses
the same old-port edge twice, so these terms are already multilinear.

The stabilizer fixes all ports and every label in S, preserves the old
arcs, and permutes the words. The coefficient of the selected monomial
in each word's contribution is invariant under this action. This includes
zero contributions, because the stabilizer permutes only variables outside
the selected set. Its action on words is free: a word contains every old
name, so an automorphism fixing that word is the identity. Summing over
orbits proves divisibility. QED.

Consequently the gcd of the stabilizer orders over all |S|<=2m divides
every actual response value. This is an exact statement about the complete
response polynomial, including its highest-degree path-cover information.
It loses integer scale if the observer subsequently normalizes the table;
it does not itself identify any two normalized response functions.

## B2. Exact content of arbitrarily many triangle components

Let T_s be the ordered join of s directed triangles, with all arcs from
earlier to later triangles. Let cont_m(T_s) be the positive gcd of all
values of F_(T_s)^(m), including every old-port mask.

**Theorem B2.** For every s,m>=0,

    cont_m(T_s) = 3^max(s-2m,0).                     (2)

**Lower bound.** The independent rotations of the triangles give the
automorphism group (C_3)^s. A monomial mentions at most 2m distinct old
labels, so at least s-2m triangles avoid all its labels. Their independent
rotations fix those labels. B1 makes every coefficient, and hence every
query value, divisible by the right side of (2).

**Sharpness.** Taking every port to be a universal source gives the value
H(T_s)=3^s. Thus the content has no prime factor other than three. Put
t=min(s,2m). Partition the first t triangles into consecutive blocks of
two, with one singleton block if necessary. Assign one port to each
block. Freeze its attachment to earlier triangles as old->port, and to
later triangles as port->old; ports are in block order. Unused ports are
universal sinks. The full extension is now an ordered join of the local
blocks with their respective ports and of the untouched triangles.

For one C3 block, setting one named attachment bit from zero to one
changes H from3 to5, giving coefficient2. For two consecutive C3 blocks,
select one bit in each triangle. The four local values, at bits
(0,0),(1,0),(0,1),(1,1), are respectively9,15,45,41, so the mixed
coefficient is -10. These small constants can be obtained by inserting
one port into the nine Hamilton words and counting the remaining words;
the standalone check also counts all7! permutations independently.

Factorization under ordered join makes the coefficient on the t chosen
bits equal to

    3^(s-t) (-10)^floor(t/2) 2^(t mod2).             (3)

A Boolean coefficient is an integer alternating sum of actual query
values. Therefore the content divides (3) as well as3^s. Their gcd is
exactly the right side of (2), proving sharpness. The case t=0 includes
m=0 and s=0, with the empty product equal to one. QED.

The same content holds if the observer also varies every port-port arc:
B1's lower bound applies to each fixed exterior orientation, while the
transitive exterior orientation already supplies the sharpness queries.

For two ports this yields arbitrarily large common factors: T_s has
content3^(s-4) for s>=4. Thus the proposed universal primitive-table
shortcut is **REFUTED**, already by the fifteen-vertex T_5, whose content
is3. This is not a response collision: no different named tournament with
a proportional table is supplied. It records the arithmetic normalization
that a proposed projective-faithfulness proof must retain.

The sharpness coefficient uses up to four distinct old labels at m=2.
It belongs to the full three-path quartic layer in (1). For T_4, selected
labels0,3 at x and6,9 at y give quartic coefficient100 while H=81, so the
complete table has content one. Appending untouched triangle components
multiplies this local certificate by three for each appended triangle.

## B3. Endpoint parity of a module descends through every ambient quotient

Let M be a nonempty module substituted at z in an arbitrary quotient Q;
write T=Q[M] and T'=Q[M^op], retaining all labels and cross arcs.
Let S_u,E_u count Hamilton paths starting/ending at u. Let C_uv count
ordered covers by two nonempty paths, first ending at u, second starting
at v, including absent-arc gradient directions.

**Theorem B3.** If T and T' have equal H,S,E,C, then

    S_u(M) = E_u(M) mod2 for every u.                (4)

For k>=1 let B_k(u),E_k(u) count ordered nonempty k-path covers of M,
with the first fixed slot starting/ending at u. For k>=3 their difference
is even: the other k-1 slots can be freely permuted, so each count is
divisible by (k-1)!. For k=2 the same parity follows from P2. Indeed in
the squarefree ring over F2, let P=1^t G1=U-1, U=product_v(1+z_v),
e=Us and U^2=1. Then

    e_u P = (Us_u)(U-1) = s_u+e_u = s_u P.          (5)

Taking the full old support gives E_2(u)=B_2(u) modulo2. Thus
E_k(u)-B_k(u) is even for all k>=2.

Decompose an ambient Hamilton word into its nonempty M-runs and exterior
vertices. There are nonnegative integer coefficients d_k, independent of
u and of M's internal arcs, with

    E_u(T)=sum_k d_k E_k(u),     d_1=E_z(Q).         (6)

The distinguished M-run is the final run of the word; its endpoint is u.
All other M-runs occupy ordered slots. Reversing M replaces E_k by B_k
with the same coefficients, since fixed slots may be permuted.

For each exterior vertex w there is also an expansion

    C_uw(T)=sum_k d_(k,w) E_k(u),
    d_(1,w)=C_zw(Q).                                (7)

Here differentiate the potential cross edge u->w before evaluating the
Hamilton polynomial. Its removed factor marks the M-run immediately
before w; its endpoint is u. The remaining exterior skeleton imposes all
other cross arcs. The edge z->w need not exist for this derivative. At
one M-run this is exactly the quotient derivative C_zw(Q). These are
integer skeleton counts, with no falling-factorial division.

At least one of d_1 and the d_(1,w) for actual z->w is odd, because

    E_z(Q)+sum_(z->w) C_zw(Q)=H(Q),                 (8)

and H(Q) is odd. Every quotient Hamilton word either ends at z or uses
one actual outgoing edge at z, proving (8). If the ambient sidecars agree,
subtract the corresponding equation (6) or (7) for M and M^op, and reduce
modulo2. Terms k>=2 vanish by (5); an odd leading coefficient forces (4).
QED.

In particular, **reversing a positive even self-converse module is visible
in H,S,E,C in every ambient quotient**, by inherited P4. This strengthens
the whole-core even self-converse obstruction: the ambient tournament
itself need not be self-converse. The conclusion also covers every even
nonselfconverse M with unequal endpoint parity.

## B4. Ordered-cover skew modulo four also descends

**Theorem B4.** Under the same ambient sidecar equality,

    C_uv(M)=C_vu(M) mod4 for all u,v.                (9)

**Proof.** Retain the directed-arc variables X before evaluation. Write
q_k(M) for the ordered k-path-cover polynomial. Its unordered partition
expansion is

    q_k = k! sum_(partitions Pi into k blocks)
                      product_(B in Pi) H(M[B]).   (10)

For a directed derivative with endpoints u,v, only the block containing
both contributes. The difference of the uv and vu derivatives of that
block's H is its ordered-cover skew, which is even by P3, including every
induced subtournament. Consequently

    partial_uv q_k - partial_vu q_k is divisible by 2 k!

at every tournament adjacency point. For k>=2 it is divisible by4.
The uniform-run identity from A1 is H(T)=sum_k c_k q_k(M), with
c_1=H(Q) odd. Therefore

    C_uv(T)-C_vu(T)
       = H(Q)(C_uv(M)-C_vu(M)) mod4.                (11)

Reversing M transposes its internal gradient entries in T, as in A2.
Equality of the ambient matrices makes the left side zero. Multiplication
by the odd H(Q) is invertible modulo4, proving (9). QED.

Together B3 and B4 require both endpoint parity balance and ordered-cover
symmetry modulo4 for a hidden module converse. These are necessary, not
sufficient, conditions. Odd modules of size>=3 already violate (9) by A1.
Even self-converse modules violate (4). The previous even nonselfconverse
parity hostile has rows

    (4,17,42,131,45,75,31,119).

It satisfies (4) but C_02=162 and C_20=12, whose difference150 is2 modulo4.
Its converse reverses the difference. Hence neither of the two known
order-eight endpoint-parity hostile classes can hide its converse in
any ambient quotient. We do not claim no larger even core satisfies both
congruences; this is an explicit remaining algebraic obstruction.

## B5. The multi-in-port zero-ready boundary is exactly one physical port

In Session181's [multi-in-port grammar](../../session181_multiport_defect_20261104/multiport/PROOF.md),
take nonempty prefix A, suffix B and nonempty named W subset V(B). Its
rooted extension satisfies

    E_r(M,2)-H(M)=H(A) J_W,
    J_W=sum_(u,v in W) C_uv(B).

The inherited positive cut weights are w_uv=(C_uv+C_vu)/2, positive
integers for every distinct pair. Thus

    J_W=2 sum_({u,v} subset W) w_uv >= |W|(|W|-1).   (12)

**Theorem B5.** The root is zero-ready if and only if |W|=1. With at
least two physical in-ports its defect is at least
H(A)|W|(|W|-1)>0. Equality holds exactly when every pair weight inside W
is one. This applies at every suffix and prefix size; no strongness or
finite-menu restriction is required. It closes this specific proposed
zero-ready escape, not the general return grammar or two-port problem.

## Reproduction and scope

Run `python check.py` and `python -O check.py` in this directory. The
standalone standard-library script pins its own source at start and end,
uses explicit exceptions, retains actual response values, and emits
deterministic JSON. The normal and optimized receipts are frozen beside
the script. They test ambient endpoint and cross-channel expansions,
internal skew transfer, full-response coefficients and content, and the
physical-port defect. The finite universe and hostile controls are stated
in the receipts; none of these samples replaces the all-order proofs.

The receipts agree byte-for-byte. They include 889 ambient module pairs,
9,564 exact fixed-slot endpoint/cross-channel identities, and14,612 internal
skew entries, with contexts permitting up to four runs. Every labelled
module of orders1-4 is tested in every labelled quotient of orders1-3;
the receipts name the additional even hostiles and deeper contexts.
Eleven complete response tables contain4,750 queries. Their full Boolean
coefficients are checked against the actual pointwise stabilizer orders,
not merely the weaker uniform lower bound. The local constants also have
literal permutation counts. Actual queries certify quartic coefficients
100 at old order12 and300 at old order15. Finally126 physical multi-in-port
extensions test both inherited identities and the singleton defect boundary.

The source/target maps now separate three issues: full response to integer
content loses scale but has exact stabilizer control; an ambient converse
collision descends to two intrinsic module congruences; a set of physical
in-ports produces a strictly positive pairwise defect once its size is two.
No universal injectivity assertion, complete order-eleven census, unit
distance realization, or six-colour obstruction follows from this packet.

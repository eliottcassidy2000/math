# Odd modules cannot hide from two ports; the complete three-path transfer

**Status: PROVED** exact algebra and an all-order family obstruction;
**FINITE-EXACT** controls are stated separately. Universal two-port
reconstruction remains **OPEN**. No novelty priority is claimed.

## 1. Inheritance and the precise advance

The earlier [two-port packet](../../two_port_unit_distance_20261103/README.md)
retains H, both endpoint vectors, the ordered two-path-cover matrix C,
every one-port deletion response, and the complete bidegree-(2,2) response.
It proves reconstruction with floor(n/2) ports and exact two-port faithfulness
through order nine. Neither its positive-circulation constraints nor its
ordered-cover tomography alone is asserted globally sufficient.

The closest new mechanism is Session180's [squarefree named endpoint transfer](../../session180_170807_returns_20261104/ports/PROOF.md).
We use its full matrix G, not separate endpoint margins, to derive the
complete two-port formula below. The hostile is a pair of opposite balanced
modules: [the prior contextual theorem C1-C3](../../proof_principles_20261103/continuation/PROOF.md)
makes every one-probe response identical in every uniform ambient quotient.
The least-used sidecar is the path-cover run count k in an internal gradient.

The new theorem closes that entire proposed hiding mechanism at the next
probe: **reversing any odd module of size at least three changes ordered C,
in every ambient tournament**. It also separates every selection of reversals
of pairwise disjoint odd modules, at every order. This goes beyond the older
whole-core odd-converse lemma and the single-triangle gradient theorem.
It does not classify all one-probe collisions or all possible disagreements.

## 2. Differentiate the full path-cover polynomial

For a tournament M on r named vertices, q_k(M) counts ordered k-tuples of
nonempty directed paths partitioning its vertices. Set q_0=0 for r>0 and
q_k=0 for k>r. Before evaluation, regard q_k as a polynomial in independent
directed-arc variables X_uv for u!=v. It is homogeneous of degree r-k.
Permuting its k ordered slots acts freely and preserves their product weight;
therefore every coefficient of q_k is divisible by k!, and in particular
q_k(M) is divisible by k!.

Two exact derivative identities are

    sum_(u->v in M) partial_(uv) q_k = (r-k) q_k,       (1)
    sum_(u!=v) partial_(uv) q_k = k q_(k+1).            (2)

Equation (1) is homogeneous Euler at a 0/1 tournament adjacency point.
For (2), differentiate an edge in a path of an ordered k-cover and cut
that edge. This leaves an unordered (k+1)-cover, an ordered choice of two
of its paths to concatenate, and an ordering of the resulting k paths.
Each unordered (k+1)-cover has k(k+1) possible ordered joins and k!
orders after joining. Its contribution is therefore k(k+1)k!; dividing
by (k+1)! gives exactly k times the ordered-cover count. No actual edge
is required for the join, since that directed-arc variable was removed by
the derivative. This is why absent arcs remain in (2).

Let Q be any tournament with a distinguished vertex z, and replace z by
the module M, inheriting all uniform cross edges. There are nonnegative
integer exterior skeleton counts c_k, independent of M's internal arcs,
such that the polynomial identity is

    H(Q[M])=sum_(k>=1) c_k q_k(M).                     (3)

This is the recovered uniform-run decomposition, now retained before
evaluation so it can be differentiated. Its first coefficient is

    c_1=H(Q),                                        (4)

because one M-run contracts to the vertex z in an ordinary Hamilton word
of Q. In particular c_1 is odd by tournament Hamilton-path parity.

## 3. A mod-four internal current in every odd module

For T=Q[M], let C_uv(T) count ordered spanning two-path covers of T whose
first path ends at u and second starts at v. This is the directed gradient
partial_(uv)H(T), including derivatives in absent arc directions. Define
the internal current using M's actual arcs:

    J_M(T)=sum_(u->v inside M) [C_uv(T)-C_vu(T)].       (5)

**Theorem A1.** If r=|M| is odd and r>=3, then, for every Q,

    J_M(Q[M]) = 2 mod 4.                              (6)

**Proof.** Apply (1) and (2) to (3). Each internal unordered pair has one
actual orientation, so

    J_M(T)=sum_k c_k [2(r-k)q_k(M)-kq_(k+1)(M)].       (7)

For k>=2, the bracket is divisible by4: 2q_k is divisible by4, and
k(k+1)! is divisible by4. For k=1, put H(empty)=1 and observe

    P_M(2):=sum_(S subset V(M)) H(M[S])H(M[V(M)\\S])
           = 2H(M)+q_2(M).

Complementation pairs the 2^r subsets into 2^(r-1) pairs, each contributing
2 modulo4, since every induced Hamilton count is odd (including the empty
one). For r>=2 the number of pairs is even, so P_M(2)=0 modulo4. Thus the
k=1 bracket is

    2(r-1)H(M)-q_2(M)=2rH(M)-P_M(2)=2 modulo4

when r is odd. Equation (4) makes its multiplier c_1 odd. All other terms
vanish modulo4, proving (6). QED.

The nonzero residue is independent of the size, strong connectivity,
orientation, or symmetries of the exterior quotient. M need not be regular,
vertex-transitive, or endpoint-balanced. The module condition is essential
to the fixed coefficients in (3); an arbitrary odd subset is not covered.
For a concrete hostile, the transitive five-vertex tournament with outgoing
masks (0,1,3,7,15) has internal current -12 on the nonmodule subset
{0,1,4}, hence zero rather than two modulo4. Singletons also have zero
internal current, explaining the r>=3 boundary.

## 4. Reversing odd modules is always visible at two ports

**Theorem A2.** Let M be an odd module of size at least three in T. Reverse
exactly its internal arcs to obtain T'. Then C(T')!=C(T), and consequently
their complete two-port responses differ.

**Proof.** Reversing every path of a cover proves q_k(X^t)=q_k(X) as a
polynomial identity. The same exterior coefficients occur in (3), so for
u,v in M,

    C_uv(T')=C_vu(T).                                 (8)

If the two full matrices agreed, the internal matrix would be symmetric,
contradicting (6). Session27's actual two-port tomography recovers ordered
C, so at least one actual query in that battery distinguishes T,T'. QED.

There is an integer contrast, not just a nonzero abstract coefficient.
Using the actual arc directions of the first M as fixed reference, (8)
and (6) give

    sum_(u->v in M) [(C_uv-C_vu)(T)-(C_uv-C_vu)(T')]
       = 2J_M(T) = 4 mod 8.                           (9)

The ordered current itself is an integer combination of actual zero-,
one-, and two-port query values. In Session27's convention, a mask bit
means old vertex->port, and

    D_uv=F(e_u,e_v)-F(e_v,e_u),
    P_u=F(e_u,0)-F(e_u),
    C_uv-C_vu=P_u-P_v-D_uv.                            (10)

Thus the contrast (9) is detectable even when each queried H is supplied
only modulo8. No division by two in a residue ring is used. If M is
regular, summing the potential terms P_u-P_v over its arcs gives zero,
so only the singleton-pair two-port contrasts are needed.

**Corollary A3 (many modules).** Suppose M_1,...,M_s are pairwise disjoint
odd modules of size at least three. The complete two-port response separates
all 2^s labelled tournaments obtained by independently reversing them.

To justify the independence, collapse all s modules to their quotient.
The joint skeleton expansion of H is a sum of products
c_(k_1,...,k_s) product_i q_(k_i)(M_i). Reversing any other module leaves
each of its q-polynomials unchanged, so it leaves every internal gradient
entry of M_i unchanged. Reversing M_i transposes those entries. Therefore
two distinct reversal selections are distinguished by (9) on any module
where they differ. This is a labelled statement, not a claim that these
2^s tournaments are pairwise nonisomorphic.

Every nontrivial module in the earlier endpoint-balanced hiding theorem
has odd size: its k=r-1 endpoint balance forces ordinary regular degree.
Hence A3 applies to that entire exponential family of one-probe collisions.
One probe loses those reversal bits; two probes recover them at every
ambient order. A counterexample to universal two-port reconstruction must
use a mechanism outside this specified family.

### An eleven-vertex control beyond the earlier finite frontier

Inflate one vertex of a directed triangle by the cyclic tournament on nine
vertices, and compare its internal converse. Outgoing masks are

    T =(542,572,632,752,992,961,899,775,527,1024,511),
    T'=(992,961,899,775,527,542,572,632,752,1024,511).

Both have H=52875 and identical complete one-port responses by the inherited
balanced-module theorem. Add x with x->0 and all other old vertices->x;
then add y with y->1 and every other vertex, including x, pointing to y.
The first counts are both133903; the two-port counts are **414055 and417121**.
These exact values are finite controls. The all-order conclusion is A2,
not extrapolation from this witness or an order-eleven exhaustive census.

## 5. Full squarefree two-port transfer, retaining three-path endpoints

Work in R_V=Z[z_v:v in V]/(z_v^2:v in V). The inherited matrix

    G=(I-ZA)^(-1)Z

enumerates every simple directed path with its full vertex support. Let
e=1^t G (column endpoint sums) and s=G1 (row endpoint sums). For vectors
a,c with a_i=[i->x], c_i=[i->y], put b=1-a, d=1-c; fix x->y. Coefficient
extraction [V] always means the full squarefree monomial of the OLD vertices.

**Theorem A4 (complete two-port formula).**

    F_T(a,c)=[V]{(1+e a)(1+b^t G c)(1+d^t s)
                +(1+e c)(d^t G a)(1+b^t s)}.           (11)

**Proof.** Every Hamilton word has x before y or y before x. In the first
case the old vertices form three successive possibly empty runs: before x,
between x and y, after y. The first and last run give 1+e a and 1+d^t s.
The middle run is either empty (the actual arc x->y) or a path counted by
b^t G c. In the second port order, the middle run must be nonempty, since
y->x is absent; the three factors are 1+e c, d^t G a, 1+b^t s. Multiplication
in R_V kills all overlapping supports and retains every valid disjoint
assignment. Taking [V] forces coverage. The construction is bijective,
including single-vertex middle runs and all paths that become orphans when
either port is deleted. QED.

For an explicit joint-endpoint expansion define

    M_jk=[V]G_jk,
    U_ijk=[V]e_i G_jk,
    W_jkl=[V]G_jk s_l,
    D_ijkl=[V]e_i G_jk s_l.                            (12)

Here D counts three ordered nonempty paths with end(P)=i, start(Q)=j,
end(Q)=k, start(R)=l. The middle singleton j=k is allowed. The outer paths
and the middle path must be disjoint, as enforced by the ring. U and W
are joint two-path endpoint tensors, not products of endpoint margins.
For nonempty V, the adjacent x,y contribution is

    sum_i E_i a_i + sum_l S_l d_l + sum_(i,l) C_il a_i d_l.

The nonadjacent x-before-y contribution is

    sum_(j,k) M_jk b_j c_k
      +sum_(i,j,k) U_ijk a_i b_j c_k
      +sum_(j,k,l) W_jkl b_j c_k d_l
      +sum_(i,j,k,l) D_ijkl a_i b_j c_k d_l.

The y-before-x contribution is the same displayed expression after the
simultaneous replacements (a,b,c,d)->(c,d,a,b). For empty V the adjacent
contribution is one. This is the complete response, with no lower-order
proxy substituted for its last tensor. Equality of full responses means
equality of the resulting polynomial; we do not claim it individually
recovers every entry of D.

## 6. A concrete part of the three-path tensor visible in the quartic layer

Let K_(A,B) denote the coefficient of a_A c_B when |A|=|B|=2 in the unique
Boolean multilinear response polynomial. Its entire bidegree-(2,2) part is

    [V]{(e a)(a^t G c)(c^t s)
         +(e c)(c^t G a)(a^t s)}.                     (13)

This follows by expanding b=1-a and d=1-c in (11). Repeated a- or c-labels
within these products already vanish by repeated vertex support, so no
higher-degree term can fold down into this layer. Formula (13) shows
K_(A,B)=K_(B,A)>=0 and, more precisely,

    K_({u,v},{u,v})=0,
    K_({u,v},{v,w})=C_uw(T-v)+C_wu(T-v)>0             (14)

for three distinct labels u,v,w. For the second identity, the only
surviving middle path is the singleton v, because it carries both a_v
and c_v. The two port orders leave exactly the two ordered covers of
T-v with endpoint pairs (u,w) and (w,u). Positivity is the inherited
positive cut-weight theorem. The first identity also follows directly:
the four selected edges would force a cycle x-u-y-v-x in a Hamilton path.

For disjoint A,B, (13) retains genuinely joint three-path information.
When n=4, every such coefficient is exactly4: each of the four cross pairs
between A and B contributes once in its actual direction, across the two
port orders. These identities give exact checks on any full-response
implementation and show explicitly which three-path terms a pair-marginal
compression would lose. They do not establish injectivity of K alone.

## 7. What remains

The source/target map for A1-A3 is an odd uniform module in an arbitrary
realized tournament, mapped to its actual ordered-cover current. The
preserved discriminator is a nonzero residue2 modulo4; the observer reads
it through genuine two-port queries modulo8. Labels, module incidence and
path-cover multiplicity are all retained. The destroyed information in
the one-probe observer is exactly recovered for this reversal family.

A4 instead retains the complete old support and all joint endpoints in a
squarefree algebra. It is an exact sufficient evaluation state, not a
small autonomous compression or a proof that its partial observations
recover every arc. Remaining candidates must survive the full formula
(11), not merely match H,S,E,C. Universal two-port reconstruction beyond
the finite census and outside A3 is still OPEN.

## Reproduction and exact scope

Run `python3 check.py` and `python3 -O check.py` from this directory. The
normal and optimized JSON receipts agree byte-for-byte. The standalone
standard-library engine hashes its source before and after execution;
every check uses an explicit exception and remains active under `-O`.

The module controls include 4,280 ambient pairs and 97,144 internal
transpose entries: every labelled three-vertex module in every quotient
of orders1-4; every labelled five-vertex module in every quotient of
orders1-2; eight specified five-vertex modules in all quotients of orders3-4;
and 32 deeper alternating contexts exposing higher run coefficients.
All 1,032 labelled module cores at orders3 and5 have their profile and
derivative normalization independently counted by permutations/cuts and
subset dynamic programming. Twenty-four quotient/module families check
128 independent-reversal variants, including modules of different sizes.

For every labelled core of orders0-4, the engine checks all16,933 full
two-port entries against the squarefree formula; 549 are additionally
counted by literal Hamiltonian permutations. It checks all17,065 joint
three-path tensor entries in those cores, including1,200 nonzero entries
with a singleton middle path, and2,378 quartic coefficients. The eleven-
vertex witness uses exact subset dynamic programming; it is not a full
order-eleven enumeration. The odd nonmodule hostile is retained explicitly.

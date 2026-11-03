# Run-content filters and actual local transport

**PROVED elementary consequences of pinned inputs; FINITE-EXACT base
atlases; same-agent audit using independent counting algorithms.** Named
statements below are local to this file, not reserved global theorem IDs.
The full Hamiltonian-path spectrum conjecture remains **OPEN**.

## T1. Arithmetic content of a homogeneous module

Write C_j(B) for ordered nonempty spanning j-path covers, C_1=H, and

    d_k(B)=gcd(C_1,...,C_k),   d(B)=d_|B|(B).

For a marked quotient (Q,v), the inherited maximal-run bijection gives

    H(Q[v<-B]) = sum_j b_j(Q,v) C_j(B),   b_1=H(Q),          (1)

where b_j counts quotient words using v exactly j times and every other
vertex once. Successive letters must be valid quotient arcs. Thus d(B)
divides every homogeneous contextual count. If the context allows at most
k visits, the stronger divisor d_k applies. This is a necessary test;
arbitrary coefficient vectors b_j are not asserted realizable.

Let E_B(r) count ordered r-slot spanning path covers, allowing empty slots.
Then

    E_B(X)=sum_j C_j(B) binom(X,j),
    d(B)=gcd(E_B(1),...,E_B(|B|))=gcd_(r in Z) E_B(r).       (2)

The finite evaluation matrix is triangular integral with diagonal1; its
inverse is the finite-difference matrix. Generalized binomial coefficients
are integral on all integers, proving the last equality. Applying THM-002
to every colour class of an r-colouring gives the stronger monomial form

    E_B(X)=sum_D 2^|D| X^(|B|-sum_(cycle in D)(|cycle|-1)),  (3)

over collections D of disjoint odd directed cycles. Consequently E_B is
monic in Z[X], has zero constant term for nonempty B, and has parity |B|.
For a prime ell the following are equivalent:

    ell | d(B);
    E_B(a)=0 for every a in F_ell;
    X^ell-X divides E_B(X) over F_ell.                      (4)

The last equivalence uses the distinct roots of X^ell-X. This tests prime
divisibility only; it does not recover prime-power precision. Also
|Aut(B)| divides every C_j, because an automorphism fixing an ordered tuple
of directed paths fixes every vertex. Equality with d is not implied.

### Prime-count lemma and the two lane filters

If H(B)=p is a prime greater than3, then d(B)=1. Indeed SCC
multiplicativity leaves one nontrivial strong component S with H(S)=p;
all other SCCs are singletons. Write |S|=s. Busch's cited exact strong
minimum gives H(S)>s for s>=4; S cannot be a triangle. Hence s<p.
No directed cycle crosses an SCC, so (3) gives E_B=X^a E_S. If p divided
d(B), E_S would vanish on every nonzero element of F_p, and also at0.
This contradicts its monic degree s<p. Since d(B) divides H(B)=p,
the conclusion follows.

Now let B be a nonempty proper homogeneous module of a strong T. Its
quotient Q is strong, has at least three vertices, and H(Q)>=3. Equation(1)
therefore gives H(B)<=H(T)/3.

* If H(T)=p is prime, d(B)|p and d(B)<=H(B)<p, so d(B)=1.
* If H(T)=7p with prime p>3, then d(B)|7p and d(B)<7p. For p!=7 the
  possibilities are1,7,p. If d(B)=p, the positive odd multiple H(B) of p
  is at most7p/3, forcing H(B)=p. The prime-count lemma contradicts d=p.
  For p=7 the divisors below49 are already1 and7. Thus d(B) is1 or7.

These filters apply to all proper homogeneous modules, with no assumption
about the automorphism group or number of module vertices. They do not
assert sufficiency or restrict nonhomogeneous cuts.

### Sharp run-depth example

The five-vertex B13 in transport.py has

    C=(13,78,198,240,120),  (d_1,d_2,d_3,d_4,d_5)=(13,13,1,1,1).

Any strong context with at most two visits remains divisible by13 and has
count at least39, so cannot be prime. Three visits require at least four
exterior vertices: between two visits to a homogeneous module there must
be at least two exterior vertices, since tournaments have no digons.
Three visits are genuinely enough to permit a prime: the two-arm family
below with B13 and m=3 has H=9497. Its m=2 predecessor has H=2261,
which is composite. The vanishing divisor is not a primality guarantee.

There is a sharp additional **FINITE-EXACT boundary**. Exhaust all1024
labelled five-vertex quotients, with root0. The544 strong quotients give
exactly the following counts after substituting B13:

    351,429,507,611,637,663,897,1293,1345,1449,1501,1605,2261.

All are composite. Independent direct nine-vertex subset counts and
quotient-word contraction agree on all544 cases. Together with the
two-visit obstruction for at most three exterior vertices, this proves
that **five is the minimum number of exterior vertices** in a strong
prime-count tournament containing B13 as a proper homogeneous module.
The H9497 witness supplies equality. Arbitrary nonhomogeneous attachments
and nonstrong singleton padding are outside this statement.

## T2. The native double-triangle compiler

For a Hamilton path P in the actual labelled tournament T, tau_P reverses
every edge outside P. The re-marking used here moves one selected vertex
within P and is legal only if the resulting word is a Hamilton path of the
current tournament. The inherited cancellation law is

    tau_Q tau_P(T)=T xor (E(P) symmetric_difference E(Q)).  (5)

Choose interior indices i,j with |i-j|>=2. Move P_i to position j, obtaining
Q. The three commands C, R_(i,j), C are legal exactly when Q is a Hamilton
path of tau_P(T). If legal, they reverse exactly six edges, forming two
directed triangles with the moved vertex as their only common vertex.
Every named out-degree, and both path endpoints, are preserved.

Proof for i<j: put v=P_i, u=P_(i-1), w=P_(i+1), a=P_j, b=P_(j+1).
The removed path edges are uv,vw,ab and the inserted edges uw,av,vb.
Legality in tau_P says the latter three point oppositely in T. Hence
u->v->w->u and v->a->b->v are directed cycles of T. Equation(5) reverses
both cycles. The i>j case is the inverse operation. Each vertex loses and
gains one outgoing edge per incident cycle.

Conversely, on a five-vertex union of two such cycles, a legal frame always
exists. Write them u->v->w->u and v->a->b->v. If w->a, choose
P=(u,v,w,a,b); otherwise use P=(b,v,a,w,u). Relocate position1 to3.
In larger tournaments existence of an extending Hamilton frame is an
additional guard, not a consequence of the local bowtie.

If these five or seven named vertices instead form a homogeneous module B,
any selected Hamilton path of the quotient can be expanded by a certified
module path P. The relocated module path Q has the same endpoints. Thus
the same three commands are legal on the full source and change only the
six internal edges. **All exterior and cross relations survive literally.**
The initial selected frame is prescribed. No algorithm or bound for
reaching it from an arbitrary marking is included in this theorem.

This is a source-preserving realization of the previously numerical module
exchange, not a free reset to another isomorphism class. Its inverse is the
reverse three-command word, with the retained physical path Q.

## T3. Uniform strong constructors and their transport law

Put T_m^(b)(B)=C3[B,TT_m,TT_b], with cyclic arrows B->TT_m->TT_b->B,
m>=1 and b=1 or2. These are strong, even if B itself is not strong.
The diagonal law of THM-3121 and the transitive cover counts

    C_1(TT_m)=1, C_2=2^m-2, C_3=3^m-3*2^m+3

give, with C_j=C_j(B),

    H(T_m^(1))=(2^m+1)C_1+(2^m-1)C_2;                    (6)
    H(T_m^(2))=(3*2^m-1)C_1+(2*3^m+2^m-5)C_2
                         +(2*3^m-4*2^m+2)C_3.           (7)

There are no terms beyond C_3 because the last block has only two vertices.
For (7), enumerate balanced positive run contents with the last coordinate
1 or2 in the C3 kernel. This yields the displayed three coefficients.

The mixed five-vertex modules have C=(h,78,198,240,120) for h=11,13,15.
The regular seven-vertex modules have

    C=(h,1638,6678,15120,20160,15120,5040), h=171,175,189.

Their higher-run agreement is inherited, and independently replayed here.
The new compiler implements the following certified base exchanges:

    B11 -> B13, B11 -> B15, P7 -> R7, R7 -> U7, U7 -> P7,

up to the explicit labelled target rows in verification.json. No unspecified
isomorphism is used when expanding a certificate into a larger source.
For any homogeneous context the change is

    H(new)-H(old)=(h_new-h_old) H(Q).                     (8)

For T_m^(2), H(Q)=3*2^m-1. Use global path P followed by both transitive
blocks to obtain the literal three-command certificate at every m.

In particular, the prime and7p candidate formulas are

    B11: 552*3^m-681*2^m-5;
    B13: 552*3^m-675*2^m-7;
    R7:  7*(2376*3^m-3507*2^m+713).                       (9)

All are actual tournament counts for every positive m. Primality of their
values or quotients is a separate predicate. Exact examples include
B11,m=10 giving prime31897699 and R7,m=8 giving7*14691857=102842999;
both order17 counts were also verified by full subset DP. No infinitude
of prime outputs, coverage of all primes, or interval extension follows.

The one-arm R7 family is7*(259*2^m-209). The two-arm scalar sequence
satisfies a recurrence with polynomial (X-1)(X-2)(X-3), by (7). These
scalar clocks suffice for these supplied sources, not arbitrary attachments.

## T4. Arc-transitive one-vertex attachments and a sharp Paley escape

Let B be an arc-transitive tournament of order n>=3. Add x which beats a
chosen subset S of size s and loses to its complement. Then

    H(B+x_S)=H(B)+[C_2(B)/(n(n-1))] s(n-s).               (10)

Proof. Arc transitivity implies vertex transitivity. Paths with x at an
endpoint contribute H(B): uniform start and end counts give
s H/n+(n-s)H/n. For internal x, let M_uv count ordered two-run covers
whose first run ends at u and whose second starts at v. Arc transitivity
makes M_uv+M_vu constant on unordered vertex pairs, namely
2 C_2/[n(n-1)]. Vertex transitivity makes the row and column sums of M
equal. Thus M-M^t has zero divergence, so its sum over any directed cut
is zero. The oriented cut from the complement of S to S is half the
symmetric cut sum, proving (10). No claim of uniform individual M_uv is
needed, and the coefficient is permitted to be rational in this derivation.

For Paley7, C_2=1638 and the coefficient is39. Its128 one-vertex cuts have
counts189,423,579,657 with multiplicities2,14,42,70. None is prime.
There is a nine-vertex prime witness: append a vertex beating ports{0,1},
then a vertex with cut mask251 in the resulting eight-vertex source.
Its H is1259 by two independent exact algorithms. Consequently exactly
two added vertices are necessary and sufficient for a prime extension
that keeps every original Paley7 edge. This is a minimum number of added
vertices, not the minimum prime value across all two-vertex extensions.

Alternatively, T2 first rewrites Paley7 to a labelled R7. Attach one vertex
beating exactly port0. Uniform endpoint counts in vertex-transitive R7
give H=175+1638/7=409. The verifier checks the literal source, both selected
paths, the changed six edges and the new cut. The edge-rewriting route and
the edge-preserving route are different operation languages.

## Remaining obligation

THM-4094 reduces completeness to all odd primes p!=7 and7p with p!=3.
The constructions here supply an all-size legal exchange mechanism and
arithmetically filtered source families. They do not choose a realizable
quotient for every prescribed p. The missing step is a uniform realization
theorem for the joint coefficients b_j and the compatible selected frame.
An integer solution to (1), a Boolean inverse for source laws, or an
infinite sequence of candidate counts is not that theorem.

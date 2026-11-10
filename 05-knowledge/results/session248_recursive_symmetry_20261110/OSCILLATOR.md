# The tournament growth oscillator and its primitive invariant modes

**Status: PROVED below from the session's exact birth/death identity and
signed-graph support decomposition.** This is an explicit representation
of the oscillator algebra on tournament-type functions. No claim of
priority for abstract creation/annihilation algebra or graded-graph theory.
No claim that tournament automorphism groups are oscillator groups.

## 1. The operators have a concrete combinatorial meaning

Let H_n be the finite-dimensional real Hilbert space of functions on
n-vertex tournament isomorphism classes, with inner product weighted by

    pi_n(T)=n!/(2^binom(n,2) |Aut(T)|).

Write t_n=dim H_n and t_0=1. The birth lane defines two stochastic kernels:

- B_n(T,U): append one vertex with uniformly chosen incidence;
- D_(n+1)(U,T): delete a uniformly chosen vertex.

As operators on FUNCTIONS their directions reverse the state arrows:

    (B_n f)(T)=sum_U B_n(T,U)f(U),       B_n:H_(n+1)->H_n,
    (D_(n+1)g)(U)=sum_T D_(n+1)(U,T)g(T), D_(n+1):H_n->H_(n+1).

Detailed balance makes these adjoints. The exact combinatorial identity is

    (n+1)B_n D_(n+1) - n D_n B_(n-1) = I_(H_n).       (G1)

Indeed the deletion after a birth removes the newborn with probability
1/(n+1); otherwise it is fair replacement of a uniformly chosen old
vertex. This proof retains multiplicity and the actual incidence law.

On the algebraic direct sum H=direct_sum_(n>=0) H_n, define

    a f_n = sqrt(n) B_(n-1) f_n,                     n>=1,
    a f_0 = 0,
    a* f_n = sqrt(n+1) D_(n+1) f_n.                 (G2)

The star denotes the adjoint pairing on this finite-support domain. These
are generally unbounded operators on the Hilbert-space completion; no
closure or exponentiation statement is needed here.

**Theorem G1 (PROVED).** On the finite-support domain,

    [a,a*]=a a* - a* a=I.                             (G3)

Proof: on H_n the two products are exactly the two terms of (G1).

Thus combinatorial growth gives a literal creation/annihilation relation.
The square-root factors are forced by orbit-volume detailed balance;
uniformly weighting unlabelled classes would not give this adjoint model.

## 2. Every newly born invariant supplies a vacuum

For v>=0 put

    V_v=ker(a:H_v->H_(v-1)),

with V_0=H_0. Call V_v the primitive invariant space at size v. The bridge
lane's rank theorem says that B_(v-1) has rank t_(v-1), so

    dim V_v=c_v=t_v-t_(v-1),   v>=1;   c_0=1.         (G4)

Equivalently, V_v is spanned by the inherited signed-even graph Fourier
modes with support on all v vertices. The graph types are invariant
coordinate labels, not a bijection preserving tournament automorphisms.

Choose an orthonormal basis {h_(v,j)} of V_v. For r>=0 define

    h_(v,j;r)=(a*)^r h_(v,j)/sqrt(r!).                (G5)

**Theorem G2 (PROVED).** These vectors form an orthonormal basis of H,
with h_(v,j;r) lying at tournament size n=v+r, and

    a h_(v,j;r)=sqrt(r) h_(v,j;r-1),
    a* h_(v,j;r)=sqrt(r+1) h_(v,j;r+1).              (G6)

Consequently H is an orthogonal sum of c_v oscillator ladders starting at
each size v. There is one constant ladder at v=0, no new ladder at v=1 or
2, one at v=3, two at v=4, and eight at v=5.

Proof: the commutation relation gives
`a(a*)^r=(a*)^r a+r(a*)^(r-1)`. Since a annihilates h_(v,j), induction proves
(G6) and unit norms. For vectors at the same total grade but born at
v<w, move all raising operators from the w-born vector across the inner
product. The resulting inner product contains
`< (a*)^(w-v)h_(v,j), h_(w,l)>`, which is zero because a annihilates the
second vector. Equal-v vectors retain their initial orthogonality. At
size n there are sum_(v<=n)c_v=t_n such vectors, so they are complete.

The “birth of an invariant” is therefore not just a difference of counts:
it produces an exact primitive vector and an entire normalized growth
ladder. Appending isolated vertices to its Fourier partner is the
coordinate-side realization of this ladder; averaging tournament deletions
is the physical operator-side realization.

## 3. Two notions of size, and the exact memory generator

Let L act on H_n by multiplication by n, and let N=a* a. Then on the
vector h_(v,j;r),

    L=v+r,     N=r,     S:=L-N=v.                    (G7)

S measures primitive SUPPORT size, while N counts subsequent growth steps.
They are not the same number operator. Moreover S commutes with both a
and a*, so the birth-support label is preserved along the whole ladder.

At a fixed size n, the delete-then-birth kernel is K_(n,n-1)=D_n B_(n-1).
Therefore

    n(K_(n,n-1)-I)=N-L=-S.                            (G8)

If each of the n vertices is refreshed at rate1, the Poissonized evolution
is exactly exp(-tS): every support-v invariant decays at rate v, independent
of the ambient size once this clock is used. The discrete one-refresh
factor is 1-v/n. The first nonconstant vacuum has v=3, explaining the
exact type-level gap3/n. Labelled edge functions already have support2,
explaining why their gap is only2/n.

This is the clean recursive mechanism: growth adds descendants to each
primitive invariant without changing its support label; refresh removes
information according to that label; new size-n vacua are exactly the
modes destroyed by one discrete refresh at size n.

The separate edge-noise operator has eigenvalue determined by the number
of Fourier edges e. It commutes with S but is not determined by S. The
repository's old edge-degree calculus and the new support-size calculus
therefore supply two independent exact coordinates.

## 4. The lowest nonconstant ladder is triangle memory

Let c_3(T) count cyclic triangles, and put

    Z_n(T)=c_3(T)-binom(n,3)/4.

For any retained k-subset, each triangle not entirely retained has at
most one fixed arc; averaging its remaining fair arcs gives cycle
probability1/4. Hence

    K_(n,k) Z_n = [binom(k,3)/binom(n,3)] Z_n.         (G9)

This identifies the unique support3 ladder explicitly. Under fair birth,

    E[c_3(U)|T]=c_3(T)+binom(n,2)/4,
    B_n Z_(n+1)=Z_n.                                 (G10)

Under uniform deletion,

    D_(n+1)Z_n=[(n-2)/(n+1)] Z_(n+1).                (G11)

These are closed FIRST-MOMENT identities, not a Markov quotient on the
scalar triangle count. Different tournaments with the same c_3 can have
different next-step distributions and different automorphism groups.

A normalized version is

    h_n=4Z_n/sqrt(3 binom(n,3)), n>=3.

Its norm is1: cyclic-triangle indicators have variance3/16 and distinct
triangles are pairwise uncorrelated, because disjoint edge sets are
independent and, when they share one edge, each conditional cycle
probability remains1/4. The resulting raising law is

    a* h_n=sqrt(n-2)h_(n+1),

exactly the vacuum-v=3 oscillator ladder. Thus the most familiar
intransitivity statistic is the lowest primitive mode of the full
recursive invariant space, not merely a convenient empirical coordinate.

## 5. Connection contract and boundaries

Source: actual random birth/deletion on tournament classes.
Target: oscillator ladders in the invariant Hilbert space.
Map: (G2), using reciprocal-automorphism weights.
Preserved: exact transition averages, adjointness, grading, multiplicities,
and the support spectrum.
Lost if only ladder eigenvalues are kept: values of each invariant on
individual tournament classes, actual automorphism groups, and specific
birth incidences.
Restoration sidecar: the signed-graph basis values and pointed birth data.

The finite 3-adic Heisenberg group inside the module lane and the real
oscillator relation here are distinct constructions on distinct spaces.
No group identification follows from the common word “Heisenberg.” The
former acts on actual tournament vertices; the latter acts linearly on
observables across all tournament sizes. Both now have explicit maps, so
their similarities and differences can be studied without conflating them.

No extra historical audit or independent census is used. The proof rests
on the new session identities in [birth/PROOF.md](birth/PROOF.md) and
[bridges/PROOF.md](bridges/PROOF.md), whose small probes concern those new
kernels. The oscillator normalization and ladder decomposition are
algebraic consequences proved here.

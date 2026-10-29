# Three runs: an index-54 guard and a native complement reconstruction

2026-10-29, 09:00 UTC continuation. **PROVED integral identities /
FINITE-EXACT actual tournament results.** Same-agent alternate-algorithm
audit. The named-port weight language is specified below; no claim of a
smallest carrier for every unnamed re-marking language is made.

## T1. The ordered three-run observer

For a tournament T on n named vertices, index rows and columns by ordered
triples of distinct vertices. Let M_T[s,t] count spanning covers by three
ordered nonempty directed paths, the ith having start s_i and end t_i.
Singletons are allowed and different paths are vertex-disjoint. This is
the k=3 instance of [TWO_RUN R1--R2](TWO_RUN.md), with all joint incidence
retained in the count. It remains a marginal of Session69's full
complemented-relation response.

The closest proved mechanism is the two-run parity fiber product.
The canonical hostile is the81/27 pair, whose two-run certificate is
retained unchanged. The near miss is treating a rational symmetry split
as an integral product. The less-used sidecar is the actual integer
source lattice, now for the nonabelian run group S3.

The live concepts are run permutations / integer realization / precision /
chain relations / actual complement / named versus unnamed ports.

## T2. Exact integral Fourier coordinates for S3

Order the six permutations by their image tuples:

    012,021,102,120,201,210.

For integer coefficients x_g, retain the trivial coordinate H=sum x_g,
the sign coordinate S=sum sign(g)x_g, and the four entries of
sum x_g rho(g). Here rho is the standard integral representation on
the basis e0-e2,e1-e2 of the sum-zero lattice. In row-major order for
that two-by-two matrix, this map has matrix

```text
F = [1  1  1  1  1  1
     1 -1 -1  1  1 -1
     1  1  0 -1  0 -1
     0  0  1 -1  1 -1
     0 -1  1  1 -1  0
     1 -1  0  0 -1  1].
```

The representation multiplication is checked directly for all36 products.
Thus F is an algebra homomorphism from Z[S3] into
Z direct-sum Z direct-sum Mat2(Z), and an isomorphism over Q.
Its exact Smith factors are

\[
\boxed{1,1,1,3,3,6;\qquad |\det F|=54.} \tag{T1}
\]

For a split candidate y=(H,S,a,b,c,d), the necessary and sufficient
integer-source guard is F^(-1)y integral. It is equivalent to the four
independent congruences

\[
\begin{aligned}
H-S&=0\pmod2,\\
H+S+2a+2d&=0\pmod3,\\
H-S+2a-2b-2d&=0\pmod3,\\
H-S+2b+2c&=0\pmod3.
\end{aligned} \tag{T2}
\]

Proof: six times the inverse is the integer matrix stored in the receipt.
All its rows reduce to H+S modulo2. Its first three rows modulo3 are
independent and span the remaining three. Hence (T2) is equivalent to
every inverse numerator being divisible by6. There is one binary and
three ternary constraints; omitting them invents source coefficients.
The exact quotient lattice is (Z/3)^2 direct-sum Z/6.

At modulus m, the image count is consequently

\[
\boxed{|\operatorname{im}(F\bmod m)|=
{m^6\over\gcd(m,3)^2\gcd(m,6)}.} \tag{T3}
\]

For powers2^k and3^k with k>=1 this is 2^(6k-1) and3^(6k-3).
These count images of the coefficient map. They are not the number of
all solutions to a differently reduced lattice-membership problem.
The checker additionally exhausts all source vectors modulo2,3,4,6;
at modulus6 it also checks every target against the full guard (T2).

## T3. Apply the guard to actual three-run response matrices

Every unordered triple supports a free S3 orbit of ordered triples.
The joint response commutes with simultaneous run permutation, so each
six-by-six orbit block is an integral group-convolution operator.
Applying F separately to every orbit gives an integer ambient split
operator N=F_big M F_big^(-1), but the legal lattice is

    F_big Z^(n(n-1)(n-2)),

with index54^binom(n,3). The inverse guard must be checked on each orbit.
The ambient integer sector coordinates generally contain extra states.
The intertwiner F_big M=N F_big is checked exactly, block by block.

Actual complete ordered-cover counts and rational ranks are:

| Tournament | Matrix order | Three-run covers | Rational rank |
|---|---:|---:|---:|
|C3|6|6|6|
|T5, steps1,2|60|210|56|
|Paley7|210|6678|210|

For Paley7, the full native and ambient-split primary exponent
multiplicities are:

| Prime | Native counts at exponents0,1,... | Ambient split counts |
|---:|---|---|
|2|165,31,2,3,9|160,35,6,3,6|
|3|201,2,6,0,1|197,9,3,1|

Their total exponent sums agree, respectively80 and18, as required by
the equal nonzero determinants under rational conjugation. Yet the native
images have32 times as many elements modulo2 and81 times as many modulo3.
Thus both predicted gluing primes are visible in an actual strong atom.

For any displayed exponent distribution b_e, the full precision recurrence
is log_p |im(M mod p^k)|=sum_e b_e max(k-e,0). In particular the native
Paley7 formulas are

    f2(k)=165k+31(k-1)_++2(k-2)_++3(k-3)_++9(k-4)_+,
    f3(k)=201k+2(k-1)_++6(k-2)_++(k-4)_+.

Both eventually gain210 prime digits per precision step. A rational
representation decomposition alone misses their finite initial layers.

## T4. Four apparently free directions are tetrahedron boundaries

The T5 three-run matrix has a genuine four-dimensional rational kernel.
Its entire kernel consists of run-symmetric coefficient laws generated by
oriented tetrahedron boundaries on the five named vertices:

\[
\partial[abcd]=[bcd]-[acd]+[abd]-[abc]. \tag{T4}
\]

Each triangle coefficient is copied to all six orders of its vertices.
The five tetrahedron boundaries span a four-dimensional space, and the
checker verifies M_T5 times every lifted boundary is zero. Conversely a
nonzero56-minor modulo5 proves rank at least56, so these four independent
relations exhaust the kernel. This is a literal small chain-complex
connection, not an analogy inferred from matching dimensions.

Global converse transposes the joint response matrix. All four directions
remain in the kernel after this operation, so observing both T5 and its
global converse still has rank56. This does not justify dropping them
before actual path complement.

## T5. One actual complement exposes every hidden direction

The declared weight state is v in Q^60, with coordinates indexed by
named ordered endpoint triples. Its observation is M_T v. A legal named
Hamilton-path complement changes T to U=T^op xor E(P), retaining all
vertex names and the same coefficient vector v. This is the actual
literal relation update, not global converse. We study this specific
linear weighted observer; fresh unnamed path multisets are a different
language.

**Finite-exact theorem.** For every one of the15 Hamilton paths P of T5,
the target U is strong, and

\[
\boxed{\operatorname{rank}_{\mathbb Q}
\begin{pmatrix}M_{T5}\\M_U\end{pmatrix}=60.} \tag{T5}
\]

Each rank is certified already modulo5. Thus no nonzero linear direction
of this port-weight space may be discarded while preserving the original
response and even one prescribed legal complement response. Depth zero
has four hidden directions; depth one has none. This establishes a sharp
one-step reconstruction result for this source and this query language.

An explicit equal-mass nonnegative witness makes the loss concrete.
Let v+ put weight1 on all six orders of each of the triples124 and014;
let v- do the same for024 and012. Both have total weight12. Their
difference is the lifted boundary of0124, hence all original observations
agree. For P=(0,1,2,3,4), the target outneighbor masks are

    U=[26,20,9,18,4].

At start triple(0,1,2), the before counts are1 and1, while the after
counts are2 and0. This is legal complement observing a correlation that
the original joint-response quotient erased. The witness is symmetric
under permutation of the three runs, so merely forgetting run order does
not repair it.

Vertex names are a separate distinction. If weight laws are required to
be invariant under the full order-five vertex automorphism group of T5,
there are12 ordered-port orbits and the original observer has rank12 on
that invariant space. None of the four hidden directions lies there.
This prevents extending the named counterexample to fully vertex-unnamed
states without specifying how a chosen complement path supplies a phase.

## Inheritance, controls and reproduction

The source-preserving coefficient arithmetic complements the incoming
[C3 coefficient guard](../overnight_research_20261029/CYCLIC_COEFFICIENT_GLUE.md).
S3 adds the sign character and a matrix-valued standard representation;
its index54 is not obtained by ignoring the interaction between those
sectors. The two-run parity guard is recovered by restricting attention
to exchanging two run positions, but it is not the full three-run guard.

The81/27 two-run separator is replayed in its original scope. This packet
does not compute its1320-by1320 three-run matrix. The six-vertex revisit
hostile also remains: a joint count matrix is not the full complemented
relation with all port marks. The sharp reconstruction above concerns
port weights with T and its literal complement relations known.

The checker independently enumerates every cover on C3,T5,Paley7 by
permutations with cuts and by three disjoint subset-path convolutions.
It verifies every vertex transporter and all run permutations. Modular
unit elimination computes every nonzero2- and3-primary layer; good-prime
ranks and explicit null vectors establish the rational ranks. T5 receives
an additional full integer Smith computation. The15 legal complement
targets are enumerated completely and checked for strong connectivity.

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/three_run.py --write
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/three_run.py
```

Next questions: which native response kernels survive *all* legal
complements; whether an analogous small chain kernel occurs in other
atoms; and how the fully unnamed operation averages or retains the phase
of the chosen path. The number56 here is a matrix rank; no identification
with the56 unlabelled six-vertex tournaments has been established.

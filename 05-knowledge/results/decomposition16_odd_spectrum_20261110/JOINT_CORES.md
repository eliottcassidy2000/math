# Primitive named deletion cubes have finite strong kernels

2026-11-10. **PROVED scoped reduction; FINITE-EXACT small controls.**
This extends the one-mark factorization in [induction/PROOF.md](induction/PROOF.md).
It retains all simultaneous deletions of named roots. These counts are
not the full joint attachment response, and no reconstruction claim is made.

## J1. Factorization and an exact size bound

Let T be a tournament with r>=1 distinct named vertices x_1,...,x_r.
For S contained in [r], put h_S=H(T-{x_i:i in S}), using H(empty)=1.
Write N=h_empty and d=gcd_S h_S. A strongly connected component (SCC)
is *marked* when it contains at least one of these named vertices.

The product of H(C) over all unmarked SCCs divides d. Indeed, deletion
of any subset of named vertices leaves every unmarked SCC intact, with
all its exterior arcs still uniformly forward. Ordered-join
multiplicativity factors this same product out of every h_S.

Consequently, if d=1, every unmarked SCC is a singleton: a tournament
with H=1 is transitive, and its only strongly connected instances are
singletons. Remove these unmarked singleton components. The remaining
ordered list C_1,...,C_k consists of k<=r marked SCCs, including any
marked singletons. It retains the whole deletion cube exactly, and

    |C_1|+...+|C_k| <= k + floor(log_5(N^3))
                         <= r + floor(log_5(N^3)).                 (J1)

To prove the bound, Busch's minimum theorem gives
H(C)>=5^((|C|-1)/3) for each nontrivial strong C; for a singleton the
same inequality is 1>=1. Multiply over the k components and use
N=product H(C_i). This is an effective finite bound for a fixed target
cube, not a bound uniform in N.

The external input is [Busch, Theorem 1 and Corollary 1 (2006)](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/),
checked at the primary PDF on 2026-11-10. SCC multiplicativity and the
single-count rational series were already used in the repository; see
[the inherited fixed-count spectrum grammar](../spectrum_basins_20261020e/spectrum.md).

## J2. Exact rational class series with names retained

Fix a positive odd vector h=(h_S) with gcd 1. Isomorphisms must fix
each of the r distinguished names individually. For each k, let P_(h,k)(z)
count the finite ordered marked-core lists of J1 by their total order.
The list must carry every named root exactly once and have deletion
cube h. Then the ordinary generating series for all marked tournament
classes with cube h is

    G_h(z) = sum_(k=1)^r P_(h,k)(z)/(1-z)^(k+1).                  (J2)

Proof: a unique transitive gap of arbitrary length lies before, between,
and after the k marked SCCs. There are k+1 gaps. The ordered SCC
decomposition and the retained names recover all those lengths, so this
construction neither omits classes nor double-counts them. The numerator
is a finite polynomial by J1. Therefore the number at order n is eventually
a polynomial in n of degree at most r. This counts marked classes; it
is not the number of labelled matrices or Hamiltonian paths.

The primitive hypothesis is essential for this particular denominator.
For example an unmarked C3 multiplied into every coordinate has gcd
divisible by3. More general fixed-count grammars are still finite, but
must retain those additional nontrivial unmarked components.

There is a sharp degree hostile: for the constant cube h_S=3 with r
named roots, every graph is TT_a -> C3 -> TT_b, and all roots must lie
outside the triangle. Its exact series is

    G_h(z)=(r+1)! z^(r+3)/(1-z)^(r+2).

Indeed the r named singleton components and one triangle have (r+1)!
ordered skeletons, each with r+2 independent unmarked singleton gaps.
The eventual degree is r+1, exceeding the primitive bound r.

## J3. A two-root test for compulsory fusion

For two named roots x,y, write the cube in order
(N,h_x,h_y,h_xy), and set

    Delta = N*h_xy - h_x*h_y.

If x and y belong to different SCCs, Delta=0. To see this, factor the
contributions from their two components and everything outside them:
the four entries are (Rac,Rbc,Rad,Rbd). Their 2x2 determinant vanishes.
This necessary condition does not require gcd 1.

If also gcd(N,h_x,h_y,h_xy)=1 and Delta!=0, the primitive kernel is
a **single** marked strong component. Thus J2 reduces to

    G_h(z)=P_h(z)/(1-z)^2,

and the kernel order bound sharpens to 1+floor(log_5(N^3)), regardless
of having two marked roots. This is a concrete joint coordinate that
prevents separate-component explanations of a proposed attachment.

The converse is false: a zero determinant does not guarantee separate
components. The checker records a strong example with cube (15,3,5,1),
with out-row masks (16,9,3,5,14) and named roots0,1.
Nor does a primitive cube identify the graph or determine its next
attachment. The lost coordinates are the actual named edges, endpoint
incidences, and mixed path-cover responses.

## Verification and interpretation

`joint_check.py` enumerates every labelled tournament of orders2..5
with roots0,1, counts paths by literal permutations, and computes SCCs
by separate directed searches. It checks factor divisibility, J1, J3,
the marked-class rational series through order5, sixteen larger padding
controls, and the strong zero-determinant hostile. Both normal and
optimized runs agree.

The [independent joint audit](induction/JOINT_AUDIT.md) additionally
checks three named roots, all 1,096 labelled tournaments of orders3..5,
and all 109 resulting primitive marked-cube series through order5,
using a different path-counting algorithm.

Connection contract: the map removes only unmarked singleton SCCs;
it preserves every h_S and all root names; it discards gap lengths,
restored by k+1 nonnegative integers. It does not discard vertices
inside a marked strong component. For target existence this gives a
finite kernel problem; for recursive construction the actual marked
core and its attachment sidecars still have to be produced.

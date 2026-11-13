# The integral mutation invariant and the triangle's factor three

**PROVED, scoped elementary matrix statements. FINITE-EXACT:** the
displayed tournaments and bounded controls. The congruence mechanism is
standard integer linear algebra applied to the cited Fomin–Zelevinsky
mutation formula; no literature-priority claim is made. All-length
Hamilton-count invariance under ambient mutation and Chan–Pak remain OPEN.

## M1. Every ordinary mutation is an integral congruence

Use the convention in the inherited
[short-return theorem](../../triples_ports_sieve_20261112/ports/PROOF.md):
B is skew-symmetric, with b_ij=1 when i beats j in a tournament.
The ambient states allow arbitrary integer entries. The mutation formula
is the one from [Fomin–Zelevinsky, Cluster algebras I](https://arxiv.org/abs/math/0104151).

For a pivot k define E=I except

    E_kk=-1,      E_ik=max(b_ik,0) for i != k.

Then det E=-1, E squared is I, and direct multiplication gives

    mu_k(B)=E B E^t.                                    (1)

For i,j distinct from k, the extra term is
[b_ik]_+ b_kj + b_ik[-b_kj]_+, equal to
(|b_ik|b_kj+b_ik|b_kj|)/2. Incident entries reverse sign.
This proves (1) for every integer skew matrix, including all intermediate
non-tournaments. It follows that **every** finite ordinary mutation word,
of any length, preserves the integer Smith factors, the determinant,
and ranks over every field. For every modulus q the kernel modules are
isomorphic, via E^(-t); no primeness assumption is needed.

These are necessary conditions for an ordinary mutation connection,
not sufficient conditions and not formulas for H. They also do not apply
to graph birth/deletion, arbitrary arc reversal, or a sign projection of
an ambient quiver. For even matrix order, the Pfaffian changes sign at
each mutation, since Pf(E B E^t)=det(E)Pf(B).

## M2. The three-cycle has an exact integral factor-three interface

Let T=C3 -> R, with triangle a -> b -> c -> a, and let B_R be the
skew matrix of the arbitrary residual tournament R. In the lattice with
original basis e_a,e_b,e_c,e_r, use the new basis

    f1=e_a, f2=e_b, z=e_a+e_b+e_c,
    f_r=e_r+e_a-e_b.

The basis matrix P is integral with determinant1. Pairing these vectors
under B_T gives

    P^t B_T P = J2 direct-sum D3(R),
    J2=[[0,1],[-1,0]],
    D3(R)=[[0, 3*1^t],[-3*1, B_R]].                    (2)

Indeed f1,f2 pair to1, z pairs to neither, and f_r was chosen to
cancel both cross pairings with f1,f2. The residual mutual pairings
stay B_R. Each of the three triangle vertices pairs to an exterior
vertex with value1, so z pairs to it with value3. This is an integer
basis identity, valid also for empty R, not a statement inferred
from a small-order spectrum.

Two useful consequences are:

    nullity_F3 B_T = 1 + nullity_F3 B_R;
    det B_T = 9 det B_(1 -> R).                         (3)

For the first, D3 modulo3 is an isolated zero coordinate plus B_R.
For the second, D3 is obtained from the skew matrix of a new source
followed by R by multiplying its first row and column by3.
The determinant identity holds at every order, but is only a
nonzero obstruction when |R| is odd. For even |R| both determinants
vanish by skew symmetry in odd order.

Separately, ordinal concatenation of Hamilton paths gives
H(T)=3H(R). Equations(2)-(3) explain another appearance of3 in the
same physical construction, but these two invariants must not be
identified: one is an alternating-lattice interface, the other a
path concatenation count.

## M3. Equal H and divisibility by three still do not supply a reset

The following two actual tournaments have six vertices and H=9.
Rows below are outgoing-neighbor bitmasks in the stated labels:

    U = (52,33,34,39,46,0),
    V = (58,60,57,16,32,8).

U is the five-vertex strong tournament with edge mask10 in lexicographic
pair order, followed by a global sink. V=C3 -> C3. Their
skew determinants and Pfaffians are

    det B_U=25,  |Pf B_U|=5;
    det B_V=81,  |Pf B_V|=9.                            (4)

Therefore **no ordinary mutation word of any length connects U to V**,
even if arbitrary weighted quivers occur between them. Moreover U cannot
reach any six-vertex tournament with an ordinal C3 factor.
Determinant alone proves this: all six-vertex tournaments with an ordinal
C3 factor have determinant9 or81. If the factor is extreme, (3) applies
and the three-vertex residual is transitive or cyclic. If it is interior,
the remaining three vertices split into nonempty transitive blocks before
and after it. Moving the preceding singleton sources to the back makes
the triangle extreme, preserving the determinant by(1).
Neither possible determinant is25. Thus no H-invariance assumption
is needed to exclude every such six-vertex endpoint.

This proves a realization obstruction stronger than a short-word bound
for this source. It does not prohibit changing graph order by an
additional legal operation. The three-divisibility of H_U is therefore
not a certificate that triangle-transport reset is accessible.

The checker independently counts both H values by subset dynamic
programming and all720 permutations, and computes determinants by
integer elimination and squared Pfaffians. Their modulo3 ranks are6
and4, respectively, making the obstruction visible already over F3.

## M4. A finite residue projection needs the sign information

Although each fixed mutation is integrally invertible, mutation is not
an autonomous operation on entry residues modulo an odd q. Its matrix E
depends on the signs of the integer entries.

For q=3 take two skew matrices with b_02=b_12=1 and respectively
b_01=1 and b_01=-2. The matrices agree modulo3. At pivot1 their new
b_02 entries are respectively2 and1. Thus they no longer agree modulo3.
Both are legitimate ambient integer quivers; no claim that they lie
in one supplied tournament's mutation class is required for this
failure of a global residue-only update rule.

This resembles the Collatz table's need for a branch and a division
receipt, at a precisely specified level: a finite observed residue
does not determine the integer operation. It is not a conjugacy
between the two dynamics. Retaining the literal mutation word and
integer intermediate matrices resolves the ambiguity. A proposed
finite controller must prove an adequate sign/size guard.

## M5. An actual prime17 return is separated from a ready17 chart

The [two-root birth analysis](../returns/PROOF.md) supplies a more direct
prime-lane test. Start from H5 at rows(10,12,1,4), with root1 of readiness
defect1. Two successive sole-in-neighbour births give

    A=(10,28,1,4,45,15),         H(A)=17, det B_A=49.

Its only indegree-one roots are4 and5, with respective defects0 and3.
The separate six-vertex two-root chart

    B=(58,60,1,52,36,4),         H(B)=17, det B_B=1

has indegree-one roots0 and1 with defects0 and1. Thus B is ready to give
3*17+2=53 at its second root, while the two indegree-one roots of A
do not have the required positive defect1 or2.

No ordinary FZ mutation word, of any length, takes A to B. The distinction
is visible in rank modulo7 as well as the integral determinant. This is
an obstruction to this specified same-count source repair, not to every
possible ready H17 presentation or to a repair that changes graph order.
The two-root classification in the return note identifies B's small
chart; it does not classify all rooted H17 tournaments.

For these checks readiness means F-H, where F counts ordered spanning
two-path covers with the first path ending at the named root and the
second possibly empty. The endpoint condition is included in the
certificate. Counting roots alone would miss the difference.

## Connection ledger and reproduction

| Source -> target | Exact map and preserved predicate | Lost data or extra obligation |
|---|---|---|
| Ordinary integer quiver -> integral lattice | B -> EBE^t; Smith factors and residue kernels | Adjacency signs and H are not recovered from the lattice |
| Extreme triangle -> residual interface | Unimodular basis(2); complete alternating pairing | Basis vectors are not three newly realized graph vertices |
| Equal-H sources -> proposed neutral reset | U,V are actual H9 graphs; A,B actual H17 graphs | Equal scalar count loses determinant25/81 or49/1 |
| Integer quiver -> residues | Entry reduction | Sign and lift required for mutation; q=3 hostile |
| Triangle -> scalar count | Ordinal path concatenation | H/3 must itself be realized; no division operation is constructed |

Run `python3 -B check.py` here, and repeat with `-O`. The explicit finite
universe is all residual tournaments through order4 for(2)-(3), all
tournament pivots through order4 followed by two specified additional
pivots for(1), and the four literal six-vertex source controls. The all-order
claims follow from the displayed identities. No long mutation-ball
search or full spectrum census is used.

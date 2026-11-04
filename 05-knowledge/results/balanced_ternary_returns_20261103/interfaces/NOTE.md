# Old ports stay nonready under every induced extension

**PROVED scoped obstruction and fresh-root construction; FINITE-EXACT
controls.** Let C_m, m>=2, be the actual balanced output of the preceding
session, with H(C_m)=15*2^m-1. No vertex already in C_m can become one-ready
or two-ready in any tournament containing C_m as an induced named source.
This allows arbitrarily many interior insertions and arbitrary arcs to
the new vertices. It excludes a particular return interface, not the
existence of a target Hamilton count by another operation.

A fresh ready root can be built in one native birth. Its count is
26*3^m+23*2^m-27, already larger than both requested ternary targets.
Thus readiness can be rebuilt, but this particular repair fails the
prescribed target budget.

## Inheritance and exact interface

For a rooted tournament (T,r), use

    h=H(T), s=start_r(T), q=H(T-r), F=E_r(2).

F counts ordered spanning path covers (P,Q), P nonempty ending at r,
Q possibly empty. A root is j-ready when F-h=j. The inherited native
two-birth reset NR gives h+2F, so one-ready and two-ready inputs produce
3h+2 and 3h+4. A balanced root has F=h and supplies only 3h.

The nearest prior theorem is [READY section 4](../../two_entry_funnels_20261103/READY.md#4-re-marking-cannot-make-the-standard-reset-reusable):
the output of a reset of a strong source of order at least four has no
one-ready or two-ready root. That theorem concerns the output itself.
The present result proves persistence under every induced extension for
the particular C_m family. Its mechanism combines the inherited rooted
defect with a fixed-start insertion lemma. No numerical coverage result
or new obstruction to an H value is claimed.

The relevant inherited formulas are in
[ADAPTIVE sections 2 and 3](../../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md):

    delta_r=F-h-q+s=2 sum_{{a,b} subset N^-(r)} b_ab(T-r),
    b_ab(U)>=H(U-{a,b}),    delta_r>=d^-(r)(d^-(r)-1).

Also q>=s. If r has unique in-neighbor a, then delta_r=0 and

    F_r(T)-H(T)=start_a(T-r).                          (I0)

This identity holds for weak as well as strong tournaments: a path
starting at r is exactly a path of T-r whose first vertex is not a.

## I1 A conditional endpoint monotonicity lemma

Let K be an induced subtournament of U containing the named vertex a.
Then

    start_a(U)=0  or  start_a(U)>=start_a(K).          (I1)

There is no strongness hypothesis. The zero alternative is essential:
adjoining a universal source kills the old starting count.

Proof. Suppose start_a(U)>0. A Hamilton path starting at a shows that a
reaches every vertex of U. Start with any Hamilton path P of K beginning
at a. If P is not spanning U, there is an arc from its vertex set to a
vertex z outside it; otherwise a could not reach the outside. Select the
least named such z and insert it immediately after the last vertex of P
that beats it. The next vertex, if present, loses to z, so this is again
a directed path, with the same first vertex. Iterate until spanning.
Deleting the added vertices recovers the original path of K. Thus this
deterministic procedure is injective, proving I1.

This does not say that every chosen one-vertex insertion preserves a
positive starting count. The reachability condition chooses an insertion
order in the full target U. It specializes the fixed-start insertion
argument already used in ADAPTIVE section 3 to an injection of whole
families of paths.

## I2 A finite named certificate for permanent nonreadiness

Call r an old vertex of a fixed tournament T. Each of the following is a
sufficient certificate that r is never one-ready or two-ready in any
induced extension U of T:

1. r has at least three in-neighbors in T.
2. r has exactly two old in-neighbors a,b, and T-{r,a,b} contains a
   directed triangle.
3. r has a unique old in-neighbor a; T-{r,a} contains a directed triangle;
   and there is an induced K contained in T-r, with a in K and
   start_a(K)>=3.

Proof. In the first case the indegree and defect bounds persist, so
F_r(U)-H(U)>=6. In the second case either the indegree rises to at least
three, or the same in-neighbors remain. In the latter case the retained
triangle and pair-deletion bound give delta_r(U)>=2*3=6.

For the third case, if r still has only a as in-neighbor, I0 and I1 give
gap zero or at least three. If it has exactly one additional in-neighbor
b, that b is new. Therefore U-{r,a,b} retains the certified old triangle,
and the gap is at least six. Larger indegree is covered by the first case.
These alternatives are exhaustive. No strongness assumption on U is
used, and the statement holds for all extension orders.

This is a sufficient obstruction certificate, not a classification of all
permanently nonready vertices. The sidecar is literal: the in-port name,
the three triangle names, and the marked subgraph K with its starting
count. A bare scalar H or primitive pair does not specify those objects.

## I3 Every old port of the balanced family has that certificate

Recall the prior construction. Start with B_m given by the two-apex word

    2 0^(m-1) 32,

whose physical feedback pair is (c_m,x). Its other names, in word order,
are y,c_1,...,c_(m-1),birth_a,birth_b. Add p beaten only by x, then q
beaten only by p. This is C_m; the names p and q are reset_root and
reset_child in the receipts. All names and old arcs are retained.

The prior proof gives H(B_m)=3*2^m+5 and H(C_m)=15*2^m-1.
Deleting x from B_m leaves

    TT_m -> C3,

with the initial chain y,c_1,...,c_(m-1) and final triangle
c_m -> birth_a -> birth_b -> c_m. Since x beats y, every Hamilton path
of B_m-x can be prefixed by x, and all of those paths start at y. Hence

    start_x(B_m)=H(B_m-x)=3.                         (I2)

Every vertex of B_m other than x receives both new in-neighbors p,q.
The old source B_m is strong. Vertex x already has m>=2 in-neighbors
and receives the additional in-neighbor q. Thus all vertices of B_m have
indegree at least three in C_m, satisfying certificate 1.

For p, the unique in-port is x. Take K=B_m. Equation I2 supplies three
starts, and the displayed triangle avoids p and x. Certificate 3 applies.

For q, the unique in-port is p. Take K=B_m+p, the first native child.
The inherited newborn formula gives

    start_p(K)=H(B_m)-start_x(B_m)=3*2^m+2>=14.

Deleting q,p leaves B_m and its displayed triangle. Again certificate 3
applies. Therefore **every old vertex of C_m remains nonready under every
induced extension**. In particular arbitrary root switching among those
names, even after arbitrary interior births, cannot make the standard
ready reset reusable at an old port.

The condition “induced extension” is essential. Edge changes or deletion
of an old vertex can destroy the certificate. The theorem does not
forbid ready roots among the newly added names, nor other ways to hit a
ternary target without the ready reset.

## I4 A fresh ready root exists but its count grows on the wrong scale

The old vertex c_m has exactly one Hamilton path starting there in C_m.
Its endpoint-cover count is

    F_(c_m)(C_m)=26*3^m+8*2^m-25.                    (I3)

Here is an all-size finite-skeleton proof, rather than a recurrence fit.
Let k=m-1>=1 and compress the old transitive module
U={c_1,...,c_(m-1)}. Its exterior has fixed size seven. Its only incoming
exterior vertices are y,p,q. In a two-bin cover there can be at most
five nonempty U-runs: at most two begin bins, and at most three enter
after one of y,p,q. After compressing each run to a placeholder, a fixed
exterior skeleton with r runs expands in exactly r! S(k,r) ways. Every
assigned run is an increasing subsequence of the old transitive chain.
The coefficients are nonnegative integers independent of k.

The exact values at k=1,...,5 determine these five coefficients by a
triangular integer system; the certificate gives

    (d_1,d_2,d_3,d_4,d_5)=(241,250,78,0,0).

Consequently F=241+250*(2^k-2)+78*(3^k-3*2^k+3), which is I3. For a
Hamilton path starting at exterior c_m, no U-run starts the path, so at
most three runs occur. The exact first three start counts are all one;
the same triangular argument gives coefficients (1,0,0), proving the
unique-start assertion for all k>=1.

Adjoin z beaten only by c_m. The native newborn formulas give

    H(D_m)=H(C_m)+F_(c_m)(C_m)-1
          =26*3^m+23*2^m-27,
    F_z(D_m)-H(D_m)=start_(c_m)(C_m)=1.              (I4)

Thus D_m has a genuine fresh one-ready root with all old names intact.
But

    H(D_m)-(3H(C_m)+4)=26*3^m-22*2^m-28>0

for every m>=2. At m=2 the counts are 59 -> 299, while the desired
branches are 179 and 181. The first member of this fresh-root family
already overshoots. No minimum over all possible
fresh-root constructions is asserted.

The new 3^m term has a combinatorial meaning: a readiness repair observes
three possible clone runs, while the old H family used only the dyadic
Hamilton count. This identifies the hidden growth cost behind this route.

## Interfaces and remaining bridge

The seventeen-query two-apex suffix observer is complete for its stated
word operation. The general seventeen-coordinate two-port birth observer
is a different object, retaining two named active ports under its birth
operations. Neither assertion allows an unrecorded change of old root
or attachment cut. The present proof uses the full named graph and a
small permanent obstruction certificate, not an identification of those
two observers.

The source map C_m -> D_m preserves the full induced old tournament and
the old port names; its additional sidecar is the new root z, its sole
in-port c_m, and the actual attachment cut. It preserves availability of
a one-ready reset but destroys the desired ternary count budget. I2--I3
show why searching only old-root choices cannot repair that loss. The
remaining possibilities include changing the source relation while
transporting the ports, or a different target-budgeted fresh-root
construction. No arbitrary-target selector is supplied here.

## Reproduction and scope

Run `python3 -B check.py` and `python3 -B -O check.py` in this directory.
Both outputs must equal `check.json`. The checker includes all labelled
tournaments through order five and every induced endpoint subset for I1;
all one-birth cuts of C_2 and C_3 for the old-port obstruction; exact
five-run skeleton evaluations; literal fixed-start controls through the
ten-vertex source; and a literal eight-vertex two-bin count. Full named
source and fresh-child rows are serialized. None of these finite checks
substitutes for the quantified proofs I1--I3.

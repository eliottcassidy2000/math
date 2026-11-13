# A named tournament lift of the powers-of-two Collatz orbit

**PROVED scoped graph construction. FINITE-EXACT implementation controls.**
This is an actual lift of one elementary infinite orbit family. It does
not encode arbitrary Collatz inputs or prove the Chan–Pak conjecture.

## Recover the almost-transitive family

Let G_j, j>=1, have the ordered names

    a,b,c,r1,...,rj.

Its arcs are a->b->c->a, both a and b dominate every ri, the ri form
a transitive chain, and every ri beats c. This is the first low-ready
family in [PROOF, R3](PROOF.md), with every residual vertex in I.
In the order a,b,r1,...,rj,c, every edge points forward except c->a.
It is exactly the inherited almost-transitive tournament R_(j+3), not a
new family. **THM-1370**,
[h-spectrum-omits-7-21-all-n](../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md),
records its count:

    H(G_j)=1+2^(j+1).

At the named b, the complete scalar response needed for two births is

    H=1+2^(j+1), F=2+2^(j+1), s=2^j, q=1+2^j.

Thus b is1-ready. The formulas are also the j-specialization of the
independent root classification in PROOF.

## The exact deletion staircase

Delete b. In the remaining graph, read the names in the order

    (a,r1,c,r2,...,rj).

Taking the old r1 as the new b gives exactly G_(j-1), with G_0=C3.
All arcs agree directly with the displayed construction. Consequently

    H(G_j-b)=1+2^j=(H(G_j)+1)/2.                     (L1)

The old a and c survive every such deletion. The active root b is
explicitly replaced by the next residual name; it is not silently
identified with an old vertex that was deleted.

For j>=2, M=H(G_j) is1 modulo8. Under M=4x+1 it corresponds to
x=2^(j-1), so L1 is exactly the even branch of the plus moving-modulus
Collatz map. At j=1, M=5 and x=1; the odd guard applies instead.
Although deletion L1 still exists there as a graph operation, the
guarded arithmetic lift does not use it.

## The exact cycle at5

G_1 has rows(10,12,1,4), in names(a,b,c,r1). Add x beaten only by b,
then y beaten only by x. These actual births give the states

| State | Rows | H |
|---|---|---:|
| G_1 | (10,12,1,4) | 5 |
| L_9 | (10,28,1,4,13) | 9 |
| L_17 | (10,28,1,4,45,15) | 17 |

The two births together implement the odd arithmetic step5->17=3*5+2.
Delete y, recovering exactly L_9; delete x, recovering exactly G_1.
Their counts give the next two arithmetic steps17->9->5. Thus the
four physical operations are two births followed by their reverse
deletions, while the three arithmetic edges are

    5 -> 17 -> 9 -> 5.

Every one of the four G_1 vertex names and arcs survives the entire
cycle. It can be repeated indefinitely with fresh temporary names.
The first arithmetic edge is a two-operation word; there is no claim
that each arithmetic edge is a single physical mutation.

Combining L1 with this cycle gives a constructive lift for every
starting label1+2^k with k>=2, equivalently every positive power-of-two
Collatz input. It uses an unbounded length parameter for the staircase
and three fixed named charts for the cycle.

## History is essential even in this small lift

The staircase's canonical G_3 has H17, but is not the transient L_17.
In the same six-vertex convention,

    G_3  = (58,60,1,52,36,4),
    L_17 = (10,28,1,4,45,15).

Their skew-adjacency determinants are1 and49, respectively. Likewise
the canonical H9 chart G_2 is1-ready at b, whereas L_9 has a2-ready
newborn. A scalar count does not determine the next port or deletion.
The saved receipt retains the graph, its chart and its named operation.

The determinant contrast also prevents an ordinary FZ mutation word
from identifying these particular H17 charts, since such mutation is
an integral congruence with determinant-one squared change of basis.
The all-length determinant argument is in the
[matrix proof, M5](../matrix/PROOF.md); the two integer determinants are checked here as
an independent control. No conclusion about every possible H17 ready
source is asserted.

These constructions are births and deletions. They are not the
four-step ordinary FZ return from the user's premise. They prove
exact guarded count equations for the specified family, not arbitrary
integer realization or convergence. The larger missing induction is
still an exact-target source switch that works beyond these special
addresses and retains the required joint ports.

## Reproduction

The shared [check.py](check.py) verifies all labelled deletion identities
for j=2,...,8, their exact counts, the three cycle charts, the literal
birth/deletion identities and the two determinants. The all-j graph
identity above supplies the unbounded quantifier. Both normal and -O
receipts retain the literal adjacency rows and agree exactly.

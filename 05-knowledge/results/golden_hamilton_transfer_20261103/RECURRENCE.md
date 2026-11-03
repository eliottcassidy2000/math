# Exact single-bridge composition and the golden matrices

**PROVED all orders, under the operations defined here.**

## G1. The retained source and operation

A marked tournament is a Boolean orientation of every vertex pair. Its
chosen source path is0->1->...->n-1. Actual path complement K retains these
arcs and reverses all other arcs. This differs from global converse.

The golden source is an ordered strong-component word in singleton s,
C3 c, and Q4. With nonadjacent pairs in lexicographic order and bit1
meaning i->j for i<j, the five Q4 charts are q0,q1,q2,q4,q5. Their
costs are1,3,4,4,4,4,4. This seven-symbol language and its literal K image
are inherited from [COMPLEMENTS C2-C3](../fibonacci_parity_precision_20261103/COMPLEMENTS.md).

We process K's output in reversed vertex order. The displayed word is
therefore in the reverse of source component order; within each block
its chart also undergoes the indicated reversal of physical labels.
This convention is checked against literal source charts in verify.py.

Between consecutive processed blocks all cross arcs go from old to new,
except one: new first vertex b -> old last vertex a. Write this join
T star B. The port of the result is the last vertex of B. Its associativity
means that nested expressions specify the same actual arcs. There are no
ties and no relaxed weights. The full word preserves the actual labelled
orientation; the numerical observer below deliberately forgets it. Keep
the vertex clock separately as n'=n+cost(letter). The minimum below is
for future H observations, not for simultaneously outputting the word or n.

## G2. A general formula, not restricted to golden blocks

For any nonempty tournament T with a marked vertex a, let

- H(T) count directed Hamilton paths;
- E(T,a) count such paths ending at a;
- C(T,a) count Hamilton paths starting at a, plus spanning covers by two
  nonempty directed paths, exactly one of which starts at a.

There is no factor2 in the last term: the marked component fixes which
path is which. Equivalently, with H(empty)=1,

    C(T,a) = sum_(a in S subseteq V(T)) start(T[S],a) H(T[V(T)\S]).

Define D(B,b) by reversing every arc in the definition of C; its marked
path ends at b. Also let S(B,b) count Hamilton paths starting at b. Then

    H(T star B) = H(T)H(B) - E(T,a)S(B,b) + C(T,a)D(B,b).       (G1)

**Proof.** A path not using b->a has a T-run followed by a B-run. The
missing a->b edge removes exactly E(T,a)S(B,b) concatenations. A path
using b->a has unique runs T1 B1 T2 B2, with T1,B2 allowed empty,
B1 ending at b, and T2 starting at a. The T-runs form one of the covers
counted by C and the B-runs one counted by D. Both other cross-edges,
when present, exist automatically: a and b already belong to T2 and B1.
This is a bijection, proving the positive product in (G1).

This theorem permits arbitrary internal tournaments B. It does not permit
an arbitrary cut at each old vertex; its single exceptional cross-edge is
essential.

## G3. The coarsest observer and its linear dimension

Set v(T,a)=(A,B,C)=(H-E,E,C). The empty source has formal unit (1,0,0).
For a singleton attachment, (G1) and direct path-cover decomposition give

    v(T star s)=(C,A,A+B+C),
    H(T)=A+B,  H(T star s)=A+C,
    H(T star s star s)=H(T)+2C.                              (G2)

For completeness, paths ending at the new singleton are exactly the old
Hamilton paths not ending at a, so B'=A and (G1) gives A'=C. For C', a
marked path beginning at the singleton either stays isolated, contributing
H(T), or starts singleton->a. The latter paths and two-path covers are
exactly those counted by C(T,a). Hence C'=H+C, proving (G2) directly.

Consequently the three future-H responses determine the coordinates:

    C=(H_ss-H)/2,  A=H_s-C,  B=H-A.                         (G3)

Conversely (G1) shows that every future single-bridge block response
depends only on these coordinates. Associativity allows any future word
to be treated as one block. Hence equality of v is exactly equality of
all future H responses: it is the coarsest deterministic observer for
this operation language. It is a countably infinite quotient, not a
finite three-state automaton. The legal word must be retained separately
if actual orientation, phase provenance, or another output is required.

For a fully explicit block compiler, let B_j=B star s^j, j=0,1,2,
and keep b as its original first vertex. Put

    g_j=(H(B_j), H(B_j)-S(B_j,b), D(B_j,b)).

Then H(T star B_j)=v(T,a) dot g_j. The columns of its transfer matrix are

    C-column=(g_2-g_0)/2,
    A-column=g_1-C-column,
    B-column=g_0-A-column.                                 (G4)

This proves exact rational linear closure for arbitrary supplied blocks;
compiling a large arbitrary block still requires its path counts. There
is no general polynomial-time counting claim. For the golden and strong5
alphabets below all entries are nonnegative integers.

Minimum rational linear dimension is3: the singleton-only response
Hankel minor is

    [[1,1,1],[1,1,3],[1,3,5]],  determinant -4.

Any d-dimensional linear representation factors this matrix through a
d-dimensional space, so d>=3. The matrices below attain this bound.
The same rank occurs after restricting to nonempty source prefixes.

## G4. The uniform golden recurrence

Starting at v0=(1,0,0), multiply row vectors on the right by M_l and read
H=A+B. All five marked Q4 charts remain separate letters:

```text
M_s  = [[0,1,1], [0,0,1], [1,0,1]]
M_c  = [[1,0,4], [1,0,3], [3,1,6]]
M_q0 = [[1,0,8], [1,0,7], [7,1,18]]
M_q1 = [[3,0,6], [3,0,5], [9,3,14]]
M_q2 = [[4,1,7], [3,1,6], [6,1,14]]
M_q4 = [[3,0,12],[3,0,9], [5,1,14]]
M_q5 = [[4,1,9], [4,0,7], [7,2,11]]
```

These are derived in two ways: (G4) with literal paths and covers in small
blocks, and the directed-path-forest transfer in transfer.py. The latter
retains endpoints and component pairing, rejects cycles, and rejects two
heads which can never receive a future incoming edge. In reversed order
all backward edges have length at most3. At a block boundary only the
most recent vertex can still receive a future edge. This leaves seven
forest states: empty; singleton; one path with port internal, ending, or
starting; two paths with port starting a nonsingleton or isolated path.
Anonymous older tails remain eligible for outgoing edges; multiplicity
and endpoint pairing are retained, so this is not a cycle-only invariant.

The exact projection of these seven states to (A,B,C) is

    R=((1,0,0),(0,1,1),(1,0,0),(0,1,0),(1,0,1),(0,0,1),(0,0,1)).

For every letter, the finite transfer matrices satisfy M7_l R=R M_l.
The seven-state matrices and the checked equality are in transfer.json.
Together with (G1)-(G4), this proves the recurrence for every word length;
the check through order10 is corroboration rather than the induction.

Small contexts distinguish every Q4 phase:

| Chart | H(q) | H(sq) | H(qs) | H(sqs) |
|---|---:|---:|---:|---:|
|q0|1|9|9|33|
|q1|3|15|9|31|
|q2|5|11|11|29|
|q4|3|9|15|31|
|q5|5|13|13|29|

The singleton-only sequence1,1,1,3,5,9,17,31,57,... is the inherited
tribonacci complement family; see [THM-4731, tribonacci partner](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md)
and [future-boundary section5](../tournament_future_boundary_20261020/README.md).
It is not claimed as a new recurrence.

## G5. An aggregate rational series

Let W_n=sum H(K(source)) over all marked golden source words of vertex
cost n. Different source presentations contribute separately, even when
their outputs are isomorphic. With f=(1,1,0)^t,

    sum W_n x^n = (1,0,0) [I-x M_s-x^3 M_c-x^4 sum_q M_q]^-1 f
                = N(x)/D(x),
    N=1-x^2-6x^3-76x^4-72x^5+x^6+10x^7-7x^8,
    D=1-x-x^2-8x^3-95x^4-87x^5-64x^6-157x^7
        -629x^8-202x^9-27x^10-227x^11-585x^12.

The polynomials are coprime. Thus for n>=12,

    W_n=W_(n-1)+W_(n-2)+8W_(n-3)+95W_(n-4)+87W_(n-5)
       +64W_(n-6)+157W_(n-7)+629W_(n-8)+202W_(n-9)
       +27W_(n-10)+227W_(n-11)+585W_(n-12).

Initial W_0,...,W_11 are
1,1,1,4,32,154,465,1573,7501,31452,112644,411846.
This is a different observable from both the inherited marked-source
count 1/(1-x-x^3-5x^4) and A000568's ordinary class count.

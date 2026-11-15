# Clean guards, selective opening, and the ternary iteration tax

A follow-up to the
[ternary emptying packet](../chanpak_ternary_emptying_20261115/README.md).
It separates a symmetric guard construction from a two-contact
deformation and computes exactly what repeated symmetric guards do
modulo 9. No broader spectrum coverage claim follows.

## The exact guard operator

Let P_B(t) count ordered path covers of B by t slots, allowing empty
slots, as an integer polynomial via the
[odd-cycle collection expansion](https://arxiv.org/html/2307.05569).
For two new vertices with b→a and a→B→b, call the result G(B).
A cycle touching either guard must touch both, and opens to a nonempty
odd-order path in B. If no such cycle is chosen, the guards are
singletons. If such a cycle is chosen, it contributes weight 2t.
Selecting odd order is half the difference between the ordinary
and alternating optional-path convolution. In the squarefree path
ring, the ordinary convolution multiplies by P, and the alternating
one multiplies by P^(-1), by
[rooted slot reciprocity](../chanpak_ternary_emptying_20261115/ROOTED_SLOT_RECIPROCITY.md).
Thus exactly
```
P_G(t) = t² P_B(t) + t [P_B(t+1) − P_B(t−1)].
```
For nonempty B this gives H(G(B))=H(B)+P_B(2).
If B is odd-order, the odd-slot powers t^c contribute factors
1+2^c, divisible by 3 with valuation 1+v_3(c); that is a
termwise statement, not a valuation of the whole sum.

## Repetition does not keep multiplying the ternary valuation

Suppose B has odd order. Write f=P_B, h=f(1), d=f'(1), c=f'(0),
all integers. Then for every k>=3 the exact mod-9 law is
```
H(G^k(B))/3 ≡ (-1)^(k−3) (d+c) mod 3.
```
In particular, from a singleton (d=c=1) every positive
guard iterate has valuation exactly one. Counts start
1,3,15,123,1479; continued cyclic wrapping is not ordinal ×3.

Proof: f is odd and integral, so f(2)≡−h+3d and f(3)≡3c mod 9,
and f(4)≡h mod 3. The first transformed values at indices 1,2,3,
divided by 3 and reduced mod 3, are
`(d, h+d+2c, 2h)`. The operator sends every odd polynomial to
values divisible by 3 at all integer points: just read residues
0,1,2. Its value at index 3 is 3(f(4)−f(2)) mod 9 and therefore
becomes zero modulo 9 after the second iterate. At 1 and 2
the normalized updates are (A,B)→(A+B,A+B+2C), C the normalized
value at 3. Thus the third iterate has equal normalized values
d+c at 1 and 2. From then on both double modulo 3 at each step,
proving the law. The first two singleton iterates (3,15)
have valuation one as well.

## Why clean guards alone cannot reach 63

H(G(B))>=1+2^|B| by the empty collection in the two-slot expansion.
For |B|=1,2,3,4,5 the exact formulas are respectively
`3, 5, 9+6c3, 17+10c3, 33+18c3+6c5`, where ci is the number
of directed cycles of length i. At order five c3>=2 exceeds
63; if c3<=1 then c5=0. Indeed a 5-cycle makes B strong,
and each vertex of a strong tournament belongs to a triangle
(otherwise its predecessors dominate all its successors,
preventing a route back to it). This yields exactly
`3,5,9,15,17,27,33,37,51` up to 63. The replay verifies
all labelled hosts of these orders. Closing this grammar
under ordinal products still excludes 7,21,63, since any
nontrivial factorization there would require 7 or 21.
It does admit 5,15,45.

## A selective opening creates 61 and 63

Keep a a clean guard but allow b to beat a marked subset S of
a transitive B, in addition to a. Put I=B−S. Every cycle uses
b. Cycles through both guards open to odd-order B-paths ending
in I; b-only cycles open to even-order B-paths from S to I.
With weighted path counts C and J, respectively, the ordinary
count is 1+2C+2J. At general t the weights of these terms are
2t and 2t² times the same path sums weighted by P of the
unused B vertices (the second still contains singleton a).

At middle order six, opening positions 2 and 3 (indexed from
zero in descending order) gives C=26,J=5 and H=63, the known
witness. Opening only position 3 gives H=61. This is a different
channel from a matched end-only repair; it moves by one odd
half-step. Minimality of two contacts holds in this transitive
middle subclass: with zero or one contact the last remaining
inlet has index at least n−2, giving at least 1+2^(n−1)
for n>=7 from cycles through both guards. At n<=6 the formula
for one opened position i is
`1+2^n−2 f(i)+2 sum_(j>i) f(j−i−1)`, where
f(0)=1 and f(r)=2^(r−1) for r>0. This is at most 47 for
n<=5. At n=6 it gives 95,79,69,61,51,33 for the six positions.
The only two-contact subset giving 63 at that size is {2,3},
as the simple subset replay checks.

The operator and the shifted witness provide a concrete test
for a constructive strategy: exact symmetry of the full slot
polynomial is not enough, and a selective contact changes
two different opened-path sums. Those sums, not just a scalar
valuation, are the data an absolute inverse needs.

Run `replay_guards.py` for literal count and iteration controls.

# The half-index boundary and a transitive-spine no-go for J

In this note h(i)=2i+1 is the affine half-index chart, not the
notation for a tournament's count. The two holes are h(3)=7
and h(10)=21. The requested gluing h(i)+h(j)−1=h(i+j)
would add half-indices. The allowed prime 11 is h(5); thus
the proposed J for primes at least 11 literally fails at
(11,11), as 5+5=10 lands at the forbidden count. The
natural repaired target excludes that pair. Also h(2)=5,
so the roles of 2,5,10 include three different positions
in this chart, not a count-scaling identity. Ordinal
multiplication comes from the different operation
i⋆j=2ij+i+j, since h(i⋆j)=h(i)h(j).

## An all-order obstruction when the two inputs share a transitive spine

Take a transitive spine W with two apices A,B. Let p and q
be the counts on W+A and W+B; assume both exceed 1. Then the
whole count equals p+q−1 **only if p=q=3**. Indeed the odd-cycle
packs contained in just one input already contribute exactly
p+q−1, as W is transitive. A mixed odd cycle through both apices
or two disjoint triangles, one per apex, would give positive
additional weight.

Assume A→B by naming the apices. Regard the letters 0,1,2,3
as 2[A wins]+[B wins]. A letter 1 gives the mixed triangle
A→B→site→A, so it cannot occur. Let r be the last letter 0,
which exists because A has a directed triangle. No letter
3 can occur at an index i≤r−2: the arc A→B, that letter, any
one intervening site, and r give a directed mixed 5-cycle.
Let s be the first letter 3 (B has a triangle). If s>r,
the two apices' triangles occupy disjoint spine intervals.
Thus necessarily s=r−1.

The pair (s,r) now gives a directed triangle through either
apex. If p>3, A has a triangle other than (s,r) (there
would otherwise be only that adjacent spine pair and no
longer odd cycle). Its early out-neighbor can be chosen at
some i<s; its letter is 2. Then A→i→B→s→r→A is mixed.
If q>3, B has another triangle. If its start is after r
it is disjoint from A's triangle (s,r). Otherwise its
start is s and its end y is after r, with letter 2.
Then A→y→B→s→r→A is mixed. Both give contradictions.
The tape 30 itself with A→B attains p=q=3 and H=5.

Thus enlarging the overlap is not enough if the shared
reservoir remains transitive. The actual [cyclic overlap
examples](PRIME_GLUE.md) evade this test. Separately,
if p=11 and q>243, keeping a q-count input induced in the
output would require a positive increase ≥2⌈log_3 q⌉>10
by the [induced-gap bound](../triangular_backbone_chanpak_20261116/CORE.md).
So a general J must sometimes change the larger input,
not merely cover it with more overlapping vertices.

## Arithmetic preparation is a separate issue

If exact repeatable children 3H+2d for d=0,1,2 on a ready
class were available, the canonical inverse for an odd H
has its unique d determined modulo 3. The initial roots
needed for the admissible arithmetic tree are
1,23,25,63,65,67: the last five have forbidden parents
7 or 21, and 1 is the base. These are scalar counts, not
yet a proof that they have the same reusable ready state.
For a centered alphabet H↦3H−2,3H,3H+2 the corresponding
roots are 1,19,23,61,63,65. Choosing a shortest arithmetic
descent only pays when the ready-state inverse exists.

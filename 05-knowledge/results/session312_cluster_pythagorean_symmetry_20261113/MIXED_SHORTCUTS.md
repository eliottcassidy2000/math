# The mixed Pythagorean trees have a shortcut at every vertex

Use the named branches and Euclid-parameter convention of [the arithmetic note](ARITHMETIC.md), with words read chronologically. The Berggren and Price trees are individually unique-parent enumerations. Their overlap as graphs is tiny, yet their union has a short relation valid on the whole chamber.

## Exactly two common edges

Set $x=m/n\in(0,1)$. For each pair of child maps, cross-multiply the two output slopes. The resulting numerators, with rows $L,M,R$ and columns $A,B,C$, are:
 
| | $A$ | $B$ | $C$ |
|---|---|---|---|
| $L$ | $(x-1)(2x-1)$ | $-x(x-3)$ | $x(x-1)$ |
| $M$ | $-2x^2-3x+1$ | $x(x+1)$ | $-x(x+3)$ |
| $R$ | $-x(3x+1)$ | $(x+1)(2x-1)$ | $-2x^2-x-1$ |

Only $L=A$ and $R=B$ have a rational solution in the open chamber, both at $x=1/2$; the positive root of the middle-row first entry has discriminant $17$. Primitivity forces the source to be $(1,2)$. Thus the only shared parent edges are $(3,4,5)\to(5,12,13)$ and $(3,4,5)\to(15,8,17)$.

## A mixed relation, and exponentially many native paths

For every legal $(m,n)$, 
$$
L(A(m,n))=L(2m,n+m)=(n+m,2n)=C(m,n).
$$
So the mixed chronological word $AL$ and the single Price step $C$ always have exactly the same endpoint, with no exceptional seeds or rescaling. They are genuinely different paths: every native edge strictly increases the second coordinate, and the intermediate $A$-child is strictly between the source and the common endpoint in that coordinate. In the undirected union these edges form a triangle at every source. At the root there are extra label coincidences: its $C$-child $(3,4)$ is also reached by two Berggren $L$ steps, and its $A$-child is its Berggren $L$-child. This gives the small loop on triples $(3,4,5),(5,12,13),(7,24,25)$; the triangle on $(3,4,5),(5,12,13),(21,20,29)$ is a separate root coincidence, $BA=M$.

More generally, from any seed to $C^k(m,n)=(m+(2^k-1)n,2^k n)$, choose $C$ or $AL$ independently at each of the $k$ stages. These give $2^k$ distinct paths even after forgetting branch labels: the optional intermediate vertices occur in disjoint open intervals of the second coordinate. Among these paths exactly $\binom{k}{j}$ have length $k+j$. This is a count of the constructed paths, not an assertion that no others exist.

The identity is a useful contrast with a [frieze ear](../frieze_collatz_reset2_20261007.md): an ear can preserve a transfer matrix without preserving the arithmetic guard, whereas $AL=C$ is both an operator identity and a legal rewrite on every primitive positive seed. It also preserves the dyadic memory cost, since both routes spend exactly one Price step. What it does not preserve is path length or the native tree's parent. This provides a concrete target for a Chan–Pak repair: prove the short and long routes are both legal, agree on the physical endpoint data, and keep the resource used by the next step. The identity does not itself yield a tournament repair.

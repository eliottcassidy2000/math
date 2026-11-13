# The determinant as a symmetry clock

This extends [Session 311](../session311_recursive_symmetry_20261113/WREATH_ATLAS.md). Each leaf color is represented by a different rigid tournament of the same order and full path-slot polynomial. The nonconstant cyclic-join construction of that session gives arbitrarily large rigid equal-polynomial palettes. At fixed depth and palette all the examples below have the same polynomial; they are strong at positive depth.

## Two Pythagorean trees with different memory

Write Euclid parameters as a column $(m,n)$ with $0<m<n$, coprime and of opposite parity. The primitive triple is $(n^2-m^2,2mn,n^2+m^2)$. In the convention of [THM-3357](../../../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md), the Berggren children are
$$
(n,2n-m),\quad(n,2n+m),\quad(m,2m+n).
$$
The children in [Price's second tree, §5](https://arxiv.org/pdf/0809.4324) in his first-column order are
$$
(2m,n+m),\quad(n-m,2n),\quad(n+m,2n).
$$
The respective determinants are $(1,-1,1)$ and $(2,-2,2)$. All preserve coprimality and opposite parity: for Price no odd prime divides the matrix determinant, and the output has an odd coordinate.

Two primitive mod-two columns agree projectively modulo $2^k$ precisely if their determinant is zero modulo $2^k$: the same coordinate is odd and may be normalized to one. There are $2^k$ classes, half in each parity chart. For any shared address, each Price step adds exactly one to the 2-adic valuation of the determinant between two seeds; Berggren adds zero. Thus after $q$ Price letters, projective agreement at precision $k$ is equivalent to seed agreement at precision $\max(k-q,0)$.

The inverse branches illustrate different admissibility rules. Berggren separates the three output slopes $n/m$ into $(1,2),(2,3),(3,\infty)$. For Price, even $m$ identifies the first branch; even $n$ requires the extra test $m<n/2$ versus $m>n/2$. The inverses are respectively $(m/2,n-m/2),(n/2-m,n/2),(m-n/2,n/2)$. Equality at the boundary gives the root by primitivity. Parity or real slope alone is insufficient for Price, but together they give an exact decoder.

## The mixed-schedule theorem

At each level choose one family globally. Color the depth-$d$ leaves injectively by their projective residues modulo $2^k$, $k\ge1$, then use the directed triangle at every join. Count Price levels backward from the leaves. If the $k$-th one is at chronological level $t$ (root step numbered one), the group is exactly $\Gamma_{t-1}$, the full cyclic-rotation group of a height-$(t-1)$ ternary tree. If there is no such level the group is trivial. Pure Berggren is always rigid; pure Price gives $\Gamma_{\max(d-k,0)}$.

Proof: two same-schedule subtrees are isomorphic precisely when their seeds agree at precision $\max(k-q,0)$, with $q$ the number of Price levels below their root. Induct from the leaves. While the child precision is positive, parity of the first coordinate is visible. Berggren gives child parity word $(1-x,1-x,x)$ for parent parity $x$; different parent bits give different weights, and equal bits force identity alignment of the cyclic triple. Price always gives $(0,1,1)$ and again forces identity alignment. Once aligned, equality of the residue of any child is equivalent to parent agreement at the same precision for Berggren and one fewer bit for Price, by the determinant law. At zero child precision all three children are isomorphic. Each layer with positive child precision contributes no rotation; those with zero precision form the full wreath tower above identical rigid blocks. This proves the exact group, not merely its cardinality.

At height three and precision one, move a single Price step from top to middle to bottom: the group orders are $1,3,3^4$, with the same full path-slot polynomial. The more general marking-budget theorem applies to the compressed identical rigid blocks. The mixed schedules need not produce every intermediate group order; Session 311 obtains those by other leaf arrangements.

## The raw spinor geometry is not the address geometry

[THM-3357](../../../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md) gives Berggren sibling minors $b,c,a$ for the parent triple $(a,b,c)$. Modulo three, $c$ is a unit and exactly one of $a,b$ vanishes. Price's three minors are $2b-a,-(n-m)^2,-2b$. All are units exactly when $m+n=0\bmod3$; otherwise exactly one vanishes. At a parent satisfying that congruence the three Price children have respectively $n'=0,\ m'=n',\ m'=0\bmod3$ and all fail it. Common unimodular changes of coordinates and unit column rescalings cannot change those valuations. Thus neither raw spinor construction yields the uniform ternary address metric at every level, even though Price does so at its first split.

The lesson for Chan–Pak is a candidate design principle, not coverage: use a joint arithmetic and port decoder to decide legal inverse branches, and prove that its resets and decreases actually act on the objects. The construction here changes symmetry inside a fixed count fiber and does not realize new targets.

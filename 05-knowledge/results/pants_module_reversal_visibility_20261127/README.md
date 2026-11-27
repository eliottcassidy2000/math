# Which module reversals are visible to switched path counts?

2026-11-27. The [same-square inverse](../pants_three_path_square_inverse_20261127/README.md)
left independent internal reversals in a bank of regular tournament
modules. Here some reversals are invisible even to the isomorphism
type at each cut; others are already exposed by one local six-path
coefficient.

## 1. Two colors exhaust the small carousels

A module \(M\) has uniform arrows to each outside vertex. In the
carousel of odd size \(p=2m+1\), label vertices modulo \(p\) and set
\(i\to j\) when \(1\le j-i\le m\) cyclically. The reflection
\(x\mapsto t-x\) gives an isomorphism to the reversed carousel.
For \(p=3,5\), every cut has a minority color class of size at most
two. For a singleton \(\{a\}\) take \(t=2a\); for a pair \(\{a,b\}\)
take \(t=a+b\); for an empty class take any reflection. The reflection
preserves the cut signs, fixes the exterior, and respects all
cross arrows. Thus reversing this module leaves the switched
tournament isomorphic to the original switched tournament at each
cut. The isomorphism may depend on the cut. All path and cycle
counts, including Hamiltonian counts, agree on every labelled cut.

Taking \(k\ge3\) triangle modules gives \(2^k\) distinct labelled
orientations with identical complete responses, and no two are
equal modulo labelled switching and global converse. Preserving all
cross arrows forces the switch signs constant globally. If a converse
were also taken, the switch signs on any three modules would need
three pairwise products equal to \(-1\), which is impossible.
The unswitched tournaments in this example are ordinarily isomorphic.

## 2. The quartet calculation for any regular module

Let \(C\) be the signed adjacency of a regular module. Reversing
it preserves \(G=A^2\): the internal matrix has zero row and
column sums and annihilates the uniform cross matrix on either
side. Thus [all responses through five vertices](../pants_three_path_square_inverse_20261127/README.md)
are unchanged. For internal vertices \(a,b,c\) and an exterior
vertex \(d\), let \(\tau=A_{ad}=A_{bd}=A_{cd}\) and write \(s=\pm1\)
for the internal sign. Define
\[
 H_{abc}=\sum_x C_{ax}C_{bx}C_{cx},\qquad
 J_{abc}=(C^2)_{ab}(C_{ca}+C_{cb})
 +(C^2)_{ac}(C_{ba}+C_{bc})
 +(C^2)_{bc}(C_{ab}+C_{ac}).
\]
The part depending on \(s\) of the coefficient of the cut character
\(e_a e_b e_c e_d\) in the raw six-path count is
\[
 \frac{\tau s}{4}(3H_{abc}-J_{abc}).
\]

To see this, expand the five arc factors in a six-vertex path.
Reversing an odd-length selected run cancels its contribution.
A term with four distinct cut endpoints must consist of two
disjoint two-edge runs. For the pairing \(ab|cd\) of the quartet
\(S\), write \(U_{ab}=\sum_{x\in S-\{a,b\}} A_{ax}A_{xb}\)
and \(Q_S=\sum_x\prod_{i\in S}A_{ix}\). Summing over their distinct
centers outside \(S\) gives
\[
 (G_{ab}-U_{ab})(G_{cd}-U_{cd})-Q_S.
\]
The coefficient is one fourth the sum over the three pairings
(two run orders, four endpoint orientations, and denominator 32).
For an explicit cancellation of the products \(U_{ab}U_{cd}\),
write the upper triangle on \(a,b,c,d\) as \(x,y,z,u,v,w\)
in row order. Their three products are \(2uz+2yv\),
\(-2uz+2xw\), and \(-2yv-2xw\). They sum to zero.
The sign-dependent part of \(Q_S\) is \(-\tau s H_{abc}\).
That of the linear \(U\) corrections is \(-\tau s J_{abc}\);
the external constant in \(G_{ab}=(C^2)_{ab}-|V-M|\)
cancels when summed against the six internal terms of \(C\).
This gives the formula. A fourth flip difference uses the
16 subsets of \(S\) and reads 16 times that coefficient.

## 3. Every larger carousel and every larger doubly regular module

For the carousel with \(p\ge7\), take \(a,b,c=0,1,3\).
For \(1\le j-i\le m\) cyclically, \((C^2)_{ij}=-p+4(j-i)\),
by the row disagreements off \(i,j\). Therefore \(J_{013}=8\).
The products in \(H_{013}\), in the successive ranges
\(2;\ 4,\ldots,m;\ m+1;\ m+2,m+3;\ m+4,\ldots,2m\),
have signs \(-,+,-,+,-\), giving \(H_{013}=0\).
So for any exterior vertex, flipping the module changes this
coefficient by exactly four in absolute value.

**Doubly regular modules admit an all-order detection theorem.**
For a doubly regular tournament of order \(p\),
\(C^2=J_p-pI\), with \(J_p\) the all-ones matrix. Then
\(J_{abc}=0\) since \((C^2)_{ab}=1\) off the diagonal and the
six internal entries cancel. The exact energy identity is
\[
 \sum_{a<b<c} H_{abc}^2
   =\frac{p(p-1)(p-3)(p+1)}{6}.
\]
Indeed, first sum \(H_{abc}^2\) over all ordered triples,
including repeats, and exchange the sums. It equals
\(\sum_{x,y}((C^\top C)_{xy})^3=p(p-1)^3-p(p-1)\).
The fully repeated terms vanish by regularity.
The triples with exactly two indices equal contribute
\(3p(p-1)\), since \(H_{aab}=-C_{ba}\) for \(a\ne b\).
Divide the remaining sum by six. Consequently, for every
doubly regular \(p>3\), reversal is detected in every host
with nonempty exterior. The sum of squared coefficient
changes over all internal triples for a fixed exterior
contact is \(3p(p-1)(p-3)(p+1)/8\).

There is also a triple-type selector. Write \(t=1\) for
a transitive triple and \(0\) for a cyclic one. For columns
outside the triple, expand the sign product
\((1-2b_a)(1-2b_b)(1-2b_c)\) modulo eight, where \(b_i=1\)
denotes a negative entry. The sum of single negative
indicators is \(3(p-1)/2-3\). The sum of pair indicators
is \(3(p-3)/4-t\): each pair has \((p-3)/4\) common
predecessors, and the source of a transitive triple accounts
for the one removed contribution. Thus
\(H_{abc}\equiv p-3-4t\pmod8\). Every cyclic triple detects
reversal for \(p\equiv7\bmod8\), and every transitive triple
does so for \(p\equiv3\bmod8\).

For **arbitrary** regular modules, the replay finds a detector
for every order-seven labelled matrix, and for every class at
orders 9, 11 and 13 in [McKay's complete regular-tournament lists](https://users.cecs.anu.edu.au/~bdm/data/digraphs.html).
The minimum numbers of detecting triples there are 26, 80 and
130, respectively. Each minimizer is unique up to isomorphism
and has no nontrivial proper module, checked by all subsets of
the representative. These are finite checks of the same local
formula, not a claim for all regular modules of larger order.

## 4. Exact profile count in the block bank

Use the [same-square ordered bank](../pants_three_path_square_inverse_20261127/README.md)
of at least three odd carousel blocks, with internal signs \(s_i\)
and a common cross sign \(s_0\). For each block of size at
least seven, the quartet changes only with \(s_i s_0\),
independently of the other internal signs. In the quartet formula
$G$ is fixed and $U$ uses only local edges. In $Q_S$ the sum
of terms with a center in the exterior contact's own block is zero
by that block's zero row sum; other exterior terms use only cross
arrows. Thus no other internal sign contributes. These relative signs
are all recovered from the six-vertex response. Conversely,
if they agree, a global converse first aligns the cross signs,
then only small carousel reversals are needed. Section 1 shows
these preserve every path and cycle response. Thus for \(b\)
large blocks there are exactly \(2^b\) complete profiles in the
bank, and no path or cycle count further separates them.

In particular the 16 same-square orientations of the
\((3,5,7)\) bank fall into exactly two profiles of eight.
For the quartet consisting of the outside vertex 0 and
carousel labels \(0,1,3\) in the seven-block, the two
coefficients are \(-15.5,-19.5\). The replay counts actual
six-paths on the 16 quartet flips for each orientation.

## Replay

Compile with `g++ -O3 -std=c++20 replay.cpp -o replay`
(add `-mpopcnt` when available). Run `./replay` for the
order-15 bank, all labelled regular order-seven tournaments,
and Paley controls of the doubly regular identity. Feed
McKay's upper-triangle tables as follows:

    ./replay 9 < rt9.txt
    ./replay 11 < rt11.txt
    gzip -dc rt13.txt.gz | ./replay 13

The tables and exact outputs are linked above and recorded
in [replay.out](replay.out).

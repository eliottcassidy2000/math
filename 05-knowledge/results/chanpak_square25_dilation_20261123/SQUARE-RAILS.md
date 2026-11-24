# Two positions below each dyadic odd square

Read a displayed odd label \(x\) as the requested height \(H=x-2\). The hole 7 appears on square 9; the hole 21 appears on label 23, two below square 25. Thus the two positions to track around a new odd square \(r^2\) are **ordinary heights** \(r^2-2\) and \(r^2-4\). The square itself is a third useful control.

**PROVED.** For \(r=2^m+1\), \(m\ge1\), all three heights \(r^2-4,r^2-2,r^2\) are attained, with exactly the original exceptions \((r,H)=(3,7),(5,21)\). The \(r^2-2\) host is strong at every \(m\ge2\) and has \(2m+2\) vertices. The \(r^2\) hosts are also strong. For the leftmost height the general construction permits reducible tournaments.

| Root | \(r^2-4\) | \(r^2-2\) | \(r^2\) |
|---:|---:|---:|---:|
| 3 | 5 | absent | 9 |
| 5 | absent | 23 | 25 |
| 9 | 77 | 79 | 81 |
| 17 | 285 | 287 | 289 |
| 33 | 1085 | 1087 | 1089 |

## A six-vertex norm formula

Take the strong tournament with outgoing masks \((2,16,11,35,13,23)\). Replace vertices 0 and 5 by arbitrary nonempty tournaments \(A,B\), with cross arcs fixed by the original tournament. Write \(h_A=H(A)\) and \(q_A=Q_2(A)\) for its ordinary count and its number of ordered two-path covers with empty slots allowed, and use the corresponding notation for \(B\). Then
\[
                    H=4q_Aq_B+q_Ah_B+3h_Aq_B-h_Ah_B.
\]
Indeed the only outneighbor of port 0 is the ordinary vertex 1 and the only inneighbor of port 5 is the ordinary vertex 3. Each port therefore occurs at most twice in a compressed Hamilton word, even after arbitrary inflation. For one or two visits to each port, the numbers of compressed words in the four cases \((1,1),(2,1),(1,2),(2,2)\) are \(23,9,11,4\). These enumerate just the four ordinary vertices in their feasible placements around the two specified symbols. For one visit to \(A\) the weight is \(h_A\); for two it is \(q_A-2h_A\), since both parts must be nonempty. Expanding gives the formula. Strength follows by lifting the outer directed paths and a directed cycle through each port. [Skeleton and literal replay](square_rails.py)

For identical fillers the formula is the norm
\[
 H=(2q_A+h_A)^2-2h_A^2
    =N_{\mathbb Q(\sqrt2)/\mathbb Q}(2q_A+h_A+h_A\sqrt2).
\]
This is also the field governing the [six-block shifted-square inverse](BANK.md). The product grammar constrains prime factors of the norm. The present construction realizes it directly by path placements. Exchanging two *different* fillers changes the height by twice the determinant of their \((h,q)\) vectors up to sign. Every nonempty filler has even \(q\), by swapping the path slots, so this exchange changes the height by a multiple of four.

For transitive chains of lengths \(a,b\ge1\), put \(X=2^a,Y=2^b\). The respective slot pairs are \((1,X),(1,Y)\), giving
\[
 H=4XY+X+3Y-1,\qquad 4H+7=(4X+3)(4Y+1).
\]
For a prescribed integer \(H\), its factor pairs give **every** pair of chain lengths in this family: test whether \((u-3)/4\) and \((v-1)/4\) are powers of two at least two for \(uv=4H+7\). In particular at \(a=b=m-1\) the factors are the consecutive odd integers \(2r+1,2r-1\), giving \(H=r^2-2\). The formal substitution \(a=b=0\) would give 7. Actually deleting both marked vertices leaves the four-vertex tournament of height 5: nonempty-run factors cease to count placements at width zero. This explains the boundary of the formula. Global absence of 7 and 21 follows from [THM-1370, spectrum omissions](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).

The [parallel unequal-chain construction](../chanpak_h23_golay_chains_20261124/UNEQUAL-CHAINS.md) obtains the same formula with the two ports exchanged, proves injectivity and gives a three-test integer inverse. Its [Golay contact theorem](../chanpak_h23_golay_chains_20261124/SQUARE-FIBRE.md) determines which of 4096 choices of 24 marked contacts retain the square at a fixed width.

There is a rigidity range in the inverse. If \(c=\min(a,b)\ge |a-b|+6\) then \(H+2\) is a square **only if \(a=b\)**. To see this, put \(K=2^{|a-b|}\). A square root \(s\) satisfies
\[
 s^2=4\cdot 2^{2c} K+2^c L+1,\qquad L=K+3\text{ or }3K+1 .
\]
Modulo \(2^c\), write \(s=2^{c-1}u+\epsilon\), \(\epsilon=\pm1\); this follows by factoring \((s-1)(s+1)\), one factor having valuation exactly one. Now \(2^{c-2}(u^2-16K)=L-\epsilon u\). Since \(2^c\ge64K\) and \(L\le4K\), the square-root equation gives \(0<u\le4\sqrt K+2\). Thus \(|L-\epsilon u|\le10K<2^{c-2}\). Both sides vanish, so \(L=4\sqrt K\). The possibilities from the two linear forms in \(K\) are \(\sqrt K=1,3,1/3\); only 1 is a square root of a nonnegative integral power of two here.

## The neighboring positions and their physical boundaries

For \(r^2\), take one apex over an ordered spine with contacts \(O I^m O^m I\), where \(O\) denotes an outgoing apex contact. The weighted inversion formula gives \(1+2+(2^m-1)+2^{2m}+(2^m-1)=r^2\): there are two adjacent \(OI\) pairs and the three indicated ordinary weighted sums. The first and last contacts close a spanning cycle. [One-apex word formula](../chanpak_heegner_absolute_targets_20261123/MERSENNE-RAIL.md)

For \(r^2-4\), the factorization is \((2^m-1)(2^m+3)\). At \(m\ge11\) the offset \(-1\) row for odd part one in [the first-25 family](SQUARE25.md) gives the first factor, and the \((G,\beta)=(32,3)\) row of the [signed dyadic profiles](../chanpak_absolute_liouville_leech_20261122/TRIANGULAR-OFFSETS.md) gives the second. Joining the two tournaments multiplies their counts. At \(4\le m\le10\) both factors are at most 1027, and the [small-height selector, THM-5271](../../../01-canon/theorems/THM-5271-the-count-is-a-potential-on-the-flip-metagraph-the-seam-ladder-implies-chan-pak-the-microcosm-3-1-6-1-21-the-square-sum-cycle-at-32-and-the-order-10-near-join-reach.md) covers all of them. At \(m=3\) the factor 7 fails, but the strong six-bank seed of height \(77=7\cdot11\) supplies the target intrinsically. At \(m=1\) the height is 5, and at \(m=2\) it is the genuine hole 21.

The two attainable heights \(r^2-2,r^2\) admit a single physical flip in the displayed lower host for \(r=5,9\). At \(r=17\) this particular strong host of height 287 has **no** incident edge of gain or loss two, as direct enumeration of all its 45 one-arc reversals shows. This prevents an induction on the marked arc of this host. Other representatives give minimum-order physical flips at all eight dyadic roots from 5 through 513, with a [unique reversed two-scale inclusion](PHYSICAL.md) at the first roots; the [minimum-order inverse](BANK.md) and [three macroscopic prime windows](BRACKETS.md) give actual single-flip witnesses.

# Matching a tournament of height h to a fork of height 3h

For each row and every integer \(s\ge1\), a tournament of height \(h_s\) and a tournament of height \(3h_s\) can be placed on the same vertex set with the stated number of differing arcs. The latter tournament has two marked arcs. Reversing one gives \(3h_s-2\) and reversing the other gives \(3h_s+2\).

| Source height \(h_s\) | Common order | Differing arcs in the displayed alignment | First \((h,3h-2,3h,3h+2)\) |
|---|---:|---:|---|
| \(42\cdot2^s-13\) | \(s+10\) | 7 | 71, 211, 213, 215 |
| \(280\cdot2^s+7\) | \(s+10\) | 5 | 567, 1699, 1701, 1703 |
| \(1019\cdot2^s+219\) | \(s+11\) | 8 | 2257, 6769, 6771, 6773 |

**PROVED with explicit hosts.** The second row supplies infinitely many multiples of the forbidden value 7 with all three children. The second and third rows admit paths of the indicated lengths from \(h_s\) to \(3h_s\) whose intermediate heights stay in that closed interval. A different choice of parent gives a single reversal for the first row even at order \(s+9\), as part of the longer family starting at 29. [Direct fork](DIRECT-FORK.md)

## How to match two infinite constructions

Start with a tournament \(R\), an apex and its mask of outneighbors in \(R\). Append a transitive right flank with binary address \(y\), in the orientation specified by the [star construction](../chanpak_gap_star_ports_20261122/replay.py). The inherited height formula with empty left flank is
\[
       K+cy+e[y\text{ odd}]+h f(y).
\]
Here \(h=H(R)\), \(K\) is the Hamilton count with apex and no flanks, and \(c\) sums over nonempty vertex subsets the paths starting at an apex outneighbor times the height of the complement. The signed version \(e\) counts the same terms positively at odd sizes and negatively at even sizes. The empty complement has height one. The binary term is specified by \(f(0)=0,\ f(2v+1)=f(v),\ f(2v)=f(v)+v+[v\text{ odd}]\). [Tagged path formula](../chanpak_decimal_kernel_20261122/fermat_receiver.py)

For odd \(w>0\) and \(s\ge1\), two classes of right flanks make this expression affine in \(2^s\). Iterating the recurrence on the trailing digits gives the following coefficients.

| Address | Coefficient of \(2^s\) | Constant term |
|---|---:|---:|
| \(2^s w\) | \(w(c+h)\) | \(K-h(w-1-f(w))\) |
| \(2^sw-1\) | \(wc\) | \(K+e-c+hf(w-1)\) |

Keeping a distinct arc for each sign, ask which edits preserve the coefficient and change the constant by exactly two. To match a third tournament to a resulting middle chart \(A2^s+B\), its constant must be \(B/3\) and its slope multiplied by three must have the same odd part as \(A\). This permits an integer shift in the padding length.

**FINITE-EXACT.** For underlying orders 1 through 8, all ordinary tournament representatives, all apex masks, both address types with odd \(1\le w\le31\), and all single edits within the core or between apex and core, the search gives **62** distinct positive-slope middle charts with both signs. Exactly **8** of those charts match a parent under a shift in the run length within the same address and core-size bank (the parent's address type and word can differ). Ordinary charts with the same two coefficients count once. None of these finite totals imposes a bound on other constructions.

For the multiple-of-seven row the source has outgoing adjacency rows \((8,5,9,18,7)\), numbered from zero, with apex mask 11 and address \(5\cdot2^s\). The center has rows \((16,9,67,21,38,15,59)\), apex mask 68 and address \(2^{s+1}-1\). Its lower and upper arcs have endpoint pairs \((2,4)\) and \((0,2)\), respectively. Substitution in the table gives the three counts \(840\cdot2^s+19,\ 840\cdot2^s+21,\ 840\cdot2^s+23\). The [literal vertex construction](physical_navigation.py) specifies all three rows of the theorem.

## A finite interpolation proves whole edit paths

The first \(s\) right-flank vertices form an identically oriented transitive block with uniform external incidences. Every outside vertex either beats the block or loses to it. Adding ordinal sources or sinks preserves the Hamilton count and makes the sizes of these two parts match between the two hosts. Relabel each part independently and order the growing block identically. The only arcs to edit are now in a fixed finite exterior. The alignment search in the replay minimizes the number of differences for this particular pair of parts. The resulting 7, 5 and 8 are upper bounds for distance in the full isomorphism-class metagraph.

Remove the growing block from a Hamilton path in any intermediate tournament. The exterior permutation and the positions of the occupied slots recover its pattern. For \(k\) slots, filling them takes \(k!\,{s\brace k}\) choices, the ordered partitions of the block into \(k\) nonempty pieces. The transitive order fixes every piece. Each slot except a possible final one needs a different exterior successor beaten by the block. There are at most three such vertices in these alignments, so \(k\le4\). Hence every intermediate height lies in the span of \(1^s,2^s,3^s,4^s\) for positive \(s\). Four values from direct integer subset counts determine it by the Vandermonde determinant. This argument makes the finite interpolation conclusive for every run length.

For the second family write \(u=2^s\). One five-edit route has heights
\[
 280u+7,\quad 780u+45,\quad 468u+45,\quad
 312u+67,\quad 840u+19,\quad 840u+21 .
\]
Each belongs to \([280u+7,840u+21]\) for \(u\ge2\). The penultimate host is itself the lower child. In a fixed alignment these steps reverse \((5,6),(8,9),(5,7),(6,9),(3,5)\), in order. The replay computes the analogous eight-edit path.

## What survives replacing a transitive block by any tournament?

**PROVED.** Fix any nonempty block and its uniform incidence to a finite exterior. For a tournament \(M\) in the block write \(P_k(M)\) for its number of ordered covers by \(k\) nonempty directed paths. With \(b_k\) the number of legal exterior permutations with \(k\) specified slots, the total is \(\sum_{k\ge1} b_k P_k(M)\). Deleting the occupied paths and reinserting them are inverse. The coefficients depend only on the exterior, and \(P_1(M)=H(M)\). [Uniform slot interface](../chanpak_continuation_observability_20261123/INTERFACE.md) · [Two-slot strong templates](../chanpak_collatz_edge_buckets_20261123/STRONG-RAILS.md)

Testing a linear identity between exteriors on transitive blocks of sizes \(1,\ldots,K\) proves it on *every* tournament block if none of the exteriors can occupy more than \(K\) slots. Indeed the matrix \(P_k(T_s)=k!{s\brace k}\), \(1\le s,k\le K\), is triangular with diagonal \(k!\). A constant discrepancy of two on all transitive inputs becomes precisely \(2H(M)\) on arbitrary inputs. In particular, for the [direct fork](DIRECT-FORK.md) with source \(21\cdot2^s-13\), writing \(D=H(M)\) and \(Q_2=P_2(M)+2D\) gives
\[
 h=21Q_2-13D,\qquad h\longrightarrow 3h\longrightarrow 3h\pm2D.
\]
Here \(Q_2(T_s)=2^s\). The same three exterior arcs work for all \(M\), including tournaments with internal cycles. With any nonempty \(M\) in the family, \(D=1\) holds exactly for a transitive tournament. This identifies the role of the all-transitive input in the unit correction.

Run from this directory.

~~~sh
c++ -std=c++17 -O3 -Wall -Wextra -Werror physical_forks.cpp -o /tmp/cp-ternary-forks
/tmp/cp-ternary-forks ../../../session212/cycle-levels 8 31
python3 -B physical_navigation.py
python3 -B fork_routes.py /tmp/cp-ternary-forks
~~~

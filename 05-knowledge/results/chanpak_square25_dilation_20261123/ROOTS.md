# Every residue at every order can coexist with no shifted square at any order

**PROVED.** In the [fixed family of twelve marked eleven-vertex tournaments](EXTENSIONS.md), let the marked chain have any positive width. Every odd residue modulo 24 occurs at every order at least eleven. Across **all** of these widths and all twelve tournaments, there is no path count of the form \(s^2-2\). This contrasts with the [six ordered-join generators](BANK.md), which give precisely the one root 113, and the [two independently expandable vertices](SQUARE-RAILS.md), which give infinitely many such roots. All assertions here concern the specified families.

## Twelve short square certificates

Index the marked tournaments by their width-one residues \(1,3,\ldots,23\). From their counts \((c_1,c_2,c_3)\) of one, two and three separated visits, put \(c_3=0\) when absent. At width \(w\ge1\) the shifted height is
\[
 H_w+2=2+c_1-2c_2+3c_3+(c_2-3c_3)2^w+c_3 3^w.
\]
Here are all possible residues in the indicated ranges. Every entry in the last column is a nonsquare modulo its listed modulus.

| Width-one residue | Widths | Modulus | Possible residues of \(H_w+2\) |
|---|---|---:|---|
| 1, 5, 9, 13, 17, 21 | all | 4 | 3 |
| 3 | all | 15 | 7, 11, 13, 14 |
| 7, 15 | at least 2 | 8 | 5 |
| 11 | all | 19 | 18 |
| 19 | all | 9 | 3 |
| 23 | all | 15 | 5, 7, 13 |

The two exceptions to the width ranges are \(7783+2\) and \(9327+2\), strictly between \(88^2,89^2\) and \(96^2,97^2\), respectively. For a reproducible finite proof of each infinite row, iterate the joint state \((2^w,3^w)\bmod m\) from its first listed width until it repeats; the form above depends only on this state. [Exact obstruction and inverse replay](roots.py)

There is also a sharp **single modulus** for this entire family. The least positive integer \(m\ge2\) for which every one of the twelve infinite sequences avoids squares modulo \(m\) is
\[
             m=1368=8\cdot9\cdot19 .
\]
For each \(2\le m\le1368\), iterate the same joint power state from width one. At each state test the twelve linear forms against the set of square residues. Every smaller modulus accepts a square at some width, while 1368 rejects all twelve at every state through the first repetition. This minimum is a finite calculation with a proof of exhaustiveness over widths. [Modulus replay](roots.py)

The width-one residue 11 gives the especially simple family
\[
 H_w=1026\cdot2^w-193.
\]
At every \(w\equiv3\pmod{10}\), its shifted label is 1 modulo 8, 1 modulo 3, and 9 modulo 11. Thus it is a square modulo *every* power of 2, 3 and 11 simultaneously, by unit square lifting and CRT. But \(19\mid1026\) gives \(H_w+2\equiv-1\pmod{19}\), a nonsquare. The first case is the actual height 8015, with shifted label 8017. Compare the [first false local square at sum 97](../chanpak_tiling_viewpoint_20261123/SQUARE-WINDOW-CRT.md). Here the obstruction persists on an entire sequence of strongly connected, rigid tournaments inside the sharp residue selector. The missing prime 19 also occupies one of the five prime positions in the interval \((9,25)\); the specific reason it detects this family is its divisibility of the slope 1026.

## Both translated hole positions and the square itself

**PROVED.** Across all twelve templates and all positive widths, there are no heights \(s^2-4\) either. There is exactly one square height, \(5625=75^2\), at width one on the template of residue 9. In particular the twelve-template rule meets the three square positions \(s^2-4,s^2-2,s^2\) precisely once. The two-port construction meets all three at every dyadic root \(s=2^m+1\) from \(m=3\) onward.

The following moduli certify the two further assertions. In the indicated order of the templates, iterate \((2^w,3^w)\) from width five until repetition, and check that \(H_w+d\) never meets a square residue modulo the corresponding entry. Direct evaluation at widths one through four supplies just the single collision above.

| Shift \(d\) | Moduli for width-one residues \(1,3,5,\ldots,23\), in order |
|---|---|
| 0 | 65, 4, 3, 4, 8, 3, 7, 4, 3, 4, 315, 4 |
| 4 | 8, 4, 8, 4, 1353, 4, 255, 4, 41, 3, 8, 4 |

## Exact inverse and its effective stopping point

Use the [marked-reader interpolation bound](../chanpak_moonshine_prime_blocks_20261118/EXACT_READER_INTERPOLATION.md). For any marked strongly connected tournament on at least three vertices, take a directed cycle through the mark. Remove its other vertices, choose any Hamilton path of the remaining tournament, and insert the cycle at the mark. The resulting full word has two separated visits, so \(c_2\ge1\). The number of ordered surjections onto \(k\) slots is nondecreasing with the width by appending the last symbol, and for two slots is \(2^w-2\). The full height therefore satisfies \(H_{w+1}-H_w\ge2^w\). The exact inverse on a fixed template has at most one solution: for target \(N\), binary-search the positive widths through \(\lfloor\log_2 N\rfloor+2\). Thus any finite collection of such marked templates gives only \(O(\log X)\) distinct integer heights up to \(X\).

For \(E=2^{60}(2^{61}-1)\) the closest height of these twelve infinite families is at width 111 on the width-one residue 11 family. Its signed distance is
\[
 H_{111}-E=2^{112}+2^{60}-193=5{,}192{,}296{,}858{,}534{,}828{,}781{,}452{,}000{,}936{,}066{,}879 .
\]
The search compares the neighboring widths bracketing \(E\) on each of the twelve monotone sequences. Neither perfect-number neighbor belongs to this width selector. Strong tournaments elsewhere attain many shifted squares and the general perfect-number crossing remains open.

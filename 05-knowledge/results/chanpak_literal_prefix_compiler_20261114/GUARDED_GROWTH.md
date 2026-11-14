# A growing bank of literal prefixes that the target selects back

The square-root inverse can be made to select many *different*
actual combinatorial sources, without searching for them by
order. This is a positive subfamily, inside the terrace and
independent of the arbitrary-core bridge.

Use the \(q=1,2\) guarded source language:
\[
\begin{array}{c|c|l}
q&(a',b')&\text{required}\\\hline
1&(a+2b,\ 2^{r+1}(a+b)) & a+2b\not\equiv0\pmod3\\
2&(3(a+b),\ 2^r(3a+4b)) & b\not\equiv0\pmod3.
\end{array}
\]
The seed is \((1,2^s)\) for \(s\ge1\).
The first coordinate is odd, the second even, and their mod-3
pair is nonzero. The last \(q\) is recovered by divisibility
of \(a'\) by 3; the last width follows from \(v_2(b')\) with
offset 1 for \(q=1\) and offset 0 for \(q=2\).
The two displayed linear forms then invert exactly.
Strict descent in the first coordinate terminates at the
single positive seed or rejects. This guarded inversion is the
same valuation mechanism as in the earlier
[terrace pressure lower bound](https://github.com/eliottcassidy2000/math/blob/edf094e6bc983f5e8c8599aecdd2e0c00783773c/05-knowledge/results/chanpak_frame_decoder_pressure_20261114/PRESSURE.md).

What is new is compatible anchoring at the **maximal**
square-root inverse depth. Take any guarded word and choose
a legal last \(q\) whose new first coordinate \(a\le A_k/2\),
where \(A_k=(47,56)M^k(1,2)^t\).
Its unscaled second coordinate is at most \(4a/3\).
Choose the width of that last seam to give \(3A_k\le b<6A_k\).
Then \(p=(47,56)M^k(a,b)^t\) has \(A_k^2\le p<A_{k+1}^2\),
and \(0<a<2v_k\), exactly as in the bridge construction.
So the inverse on **the integer \(p\) itself** recovers
the chosen \(k\), the least line point \((a,b)\), and the
literal guarded source. This also proves that the coding is
injective across the selected words and depths.

To see that the supply really grows with depth, start at \((1,2)\)
and use successive pairs of \(q=1\) steps, the first with width
1 or 2 and the second with width 2. The mod-3 slope \(a/b=2\)
goes to infinity and back to 2, so the guard always holds.
One pair grows \(a+b\) by at most \(100\). Optionally first
prepend the two legal steps \((q,r)=(2,1),(1,1)\), which
also return the slope to 2 and grow the count by at most 66.
For \(k\ge2\) one can use \(\lfloor(k-2)/3\rfloor\) binary
choices in this way: \(100^{1/3}<5\), \(A_k\ge159\cdot5^k\),
and the extra first-coordinate factor in the last step is
at most 2. Thus the required \(a\le A_k/2\) holds.
The optional step changes the corridor types as well as widths.
This gives exponentially many distinct actual inner words with
the exact target-selected maximal depth, not merely many congruences.

For example the replay reaches prime cofactors 1,081,163
and 1,313,652,139 with respective real seed coordinates
\((197,3232)\), \((53,126976)\) at depths 1 and 3.
The latter has a genuine \(q=2\) inner seam.
The interval-packing replay counts the final ordinary
tournaments, and no claim of infinitely many primes in these
subfamilies is made.

Run `python physical_prefix.py`.

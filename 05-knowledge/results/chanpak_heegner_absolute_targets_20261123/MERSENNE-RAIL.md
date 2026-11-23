# Strong tournament rails with Mersenne offsets

Take a transitive spine, orienting every arc forwards, and one apex.
Write O if the apex beats a spine vertex and I if that vertex beats
the apex. Every directed cycle uses the apex, so there is at most
one cycle in an odd-cycle packing. Fix the first and last spine
vertices of an odd cycle at positions \(i<j\). Its contacts must be
O and I, and its internal spine vertices can be any even subset
of the \(j-i-1\) intervening vertices. Thus the odd-cycle formula gives
\[
 H(w)=1+\sum_{i<j:\ w_i=O,w_j=I}2^{j-i-1}
         +\#\{i:w_iw_{i+1}=OI\}.
\]
The final term handles the empty intervening set at adjacent
positions. Any word beginning O and ending I gives a strong
tournament via the full directed cycle.

With \(a,b,c\ge1\) and \((x,y,z)=(2^a,2^b,2^c)\), take the word
\(OI^a OI^b OI^c\). There are three adjacent OI pairs, and the
geometric sums from the three groups of I's give
\[
 H=4+(x-1)+(2x+1)(y-1)+(4xy+2y+1)(z-1)
   =4xyz+(2y+1)(z-x)+1-y.
\]
In particular the symmetric word \(c=a\) gives the exact rail
\[
 H=2^{2a+b+2}-(2^b-1).
\]
Its inverse within this rail is exact: \(b=v_2(H-1)\) and
\((H-1)/2^b+1=4^{a+1}\).

At \(b=3\), every odd exponent \(q=2a+5\ge7\) has the literal
strong realization \(H=2^q-7\), on \(q+2\) vertices. By
[Nagell's theorem](https://archive.ymsc.tsinghua.edu.cn/pacm_download/116/6967-11512_2007_Article_BF02592006.pdf),
the only square counts on this whole rail are \(11^2=121\) at
\(q=7\) and \(181^2=32761\) at \(q=15\).
The smaller Nagell cases have separate physical words:
\(q=4\) has \(OIII\) with height 9, \(q=5\) has \(OIIOOI\) with
height 25, and \(q=3\) gives height 1, necessarily transitive.

For completeness, among positive squarefree Heegner numbers
\(1,2,3,7,11,19,43,67,163\), only 7 can occur in
\(X^2+d=2^q\) with \(q\ge3\). An odd \(d\) forces odd \(X\),
hence \(d\equiv7\pmod8\); the sole even case \(d=2\) is impossible
modulo 8. This is a dyadic norm statement, not an identification
of each Heegner number with a tournament count. See
[Stark's class-number-one determination](https://doi.org/10.1307/mmj/1028999653).
The other Mersenne offsets in the Heegner list are 1 and 3,
which give no squares on the displayed rail.

More generally, for every even \(b=2e\), the symmetric word has
square height exactly at \(a=e-2\ge1\), with root \(2^{b-1}-1\).
Indeed put \(Y=2^{a+e+1}\). If \(a>e-2\), the offset \(2^b-1\)
is too small to reach the next square below \(Y^2\). If
\(a<e-2\), then \(X^2\equiv1\pmod{2^b}\) implies
\(X\equiv\pm1\pmod{2^{b-1}}\), since consecutive even integers
\(X-1,X+1\) have smaller 2-adic valuation exactly one. Thus
\(X>1\) forces \(X\ge2^{b-1}-1>Y\), a contradiction.
At equality, \(Y=2^{b-1}\) and the offset is exactly \(2Y-1\).

Replay: `python3 heegner-rail-check.py` uses literal adjacency
and subset DP for the two square-height strong gadgets, as well
as a separate weighted-contact evaluator.

## One signed shell realizes all exponents at least ten

The one-apex rail is not needed for even exponents. Take seven
ordinary vertices with outgoing masks \((50,84,25,19,32,78,29)\)
and mark the subset with mask 35. Call this core \(A\), the marked
set \(I\), and add a singleton \(b\) and a nonempty filler \(X\).
The cross arcs are \(b\to X\to A\) and \(I\to b\to A\setminus I\).

Write \(D,Q\) for the filler's one- and two-slot counts, with
empty slots allowed in \(Q\). Let \(G\) count ordered splits of
the core into two directed paths \((R,S)\), with \(R\) nonempty
and ending in \(I\), and \(S\) arbitrary. Let \(L\) count such
splits with \(S\) nonempty and starting in \(I\).
A path with two filler runs must have the pattern
\(X_1,R,b,X_2,S\): every bridge must use the unique vertex \(b\)
that beats the filler. Thus the coefficient of \(Q-2D\) is
\(G\). For one filler run, placing \(b\) before it contributes
\(H(A)+G\), according as the path starts with \(b\) or has
a core prefix ending in \(I\). Placing \(b\) after the filler
contributes \(G-L\): the path begins in the filler and reaches
\(b\) from \(R\), then continues to \(S\) only when its start
is outside \(I\) or it is empty. The resulting height is
\[
 GQ+(H(A)-L)D .
\]
The literal subset-DP readings are \((H(A),G,L)=(121,512,128)\).
With transitive filler size \(s\ge1\), \(D=1,Q=2^s\), so this
one literal strong shell gives \(H=2^{s+9}-7\) on \(s+8\)
vertices. Strength follows from the strong core (verified
in the replay), \(I\ne\varnothing\), and the cycle of parts
through \(b,X,A\). The sign is an actual count, not a formal
subtraction of tournaments. At exponent 15 this gives a second
realization of \(181^2\), now on 14 vertices.

For the two small even exponents omitted by the shell, exponent
6 has a strong tournament with outgoing masks
\((102,4,8,19,39,78,30)\) and height 57; exponent 8 has the
strong apex word \(OI^4OI^3\) with height 249. Combined with
the earlier small words and odd rail, this gives a **strong
realization of \(2^q-7\) for every \(q\ge4\)**. The height
at \(q=3\) is 1, attainable only by transitive tournaments.

Replay the signed shell with `python3 heegner-shell-check.py`;
it uses only the standard library and checks both the coefficient
derivation's inputs and literal counts for six filler lengths.

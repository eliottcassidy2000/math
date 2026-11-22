# One ordinary tournament count can encode an entire slot profile

**PROVED.** For every order \(n\ge1\) there is a common tournament construction \(C_n(T)\), using the input \(T\) once, such that
\[
 H(C_n(T))=H(C_n(T'))\quad\Longleftrightarrow\quad P_T=P_{T'}
 \qquad (|T|=|T'|=n).
\]
The output is an ordinary Hamiltonian path count. A marked arc in the input stays one physical marked arc in the two outputs. The construction deliberately separates counts widely. The same separation gives an obstruction when a **unit** gap is prescribed.

## A recoverable transitive return

For \(m\ge1\) let \(R_m(S)=C_3[1,\mathrm{TT}_m,S]\), replacing three blocks of a cyclic triangle by a singleton, a transitive order of size \(m\), and \(S\). Put \(M=2^m-1\). The earlier changing-reader operator specializes at arity one to
\[
 H(R_m(S))=2H(S)+M\bigl(Q_2(S)-H(S)\bigr).
\]
In fact \(P_{\mathrm{TT}_m}(j)=j^m\), since each assignment of vertices to slots supplies one increasing path per slot. Substitution in the cited cyclic transfer at \(j=1\) proves the identity. At arity \(j\) the operator reads only \(P_S(j-1),P_S(j),P_S(j+1)\), with upward coefficient \(j((j+1)^m-j^m)>0\). [Changing-reader transfer and factorial divisibility](../chanpak_changing_sources_readers_20261114/FACTORIAL_LOCK.md)

For an input of order \(N\), choose any \(m\) with \(M>N!\). This choice depends only on order. Because \(M\) is odd, reducing the output modulo \(M\) recovers \(H(S)\) uniquely from \(2H(S)\), and the displayed equality then recovers \(Q_2(S)\). This radix can be sharpened when an upper bound on attainable heights is known.

Write \(r=\lfloor(n-1)/2\rfloor\). For \(r\le1\) the initial ordinary count already determines the profile by its parity, known leading term and zero constant. For \(r>1\), make \(r-1\) successive returns with widths fixed by the preceding order. Reading the final height backwards gives all the preceding ordinary heights. The height after \(i\) returns reads initial arities through \(i+1\). Its coefficient on \(Q_{i+1}(T)\) is a product of \(i\) strictly positive upward coefficients. Thus the heights determine \(Q_1(T),\ldots,Q_r(T)\) by triangular inversion. After removing the common leading term, \(P_T\) has \(r\) possible powers, with the parity of \(n\), starting at 1 or 2. Divide by the lowest power and interpolate in \(j^2\) at the \(r\) distinct squares. This proves injectivity on slot profiles. Equal profiles propagate through every operator and give equal heights. The earlier [finite list of ordinary readings](../chanpak_arc_fiber_grammar_20261122/SQUARE-EARS.md) used one recording at each depth. The increasing moduli here retain every preceding recording inside the final integer.

## Why a wide return misses the perfect-number unit

**PROVED.** Suppose reversing an arc **inside** a common \(n\)-vertex block \(S\) changes the count on \(R_m(S)\) from \(N\) to \(N+2\). If \(\delta\) denotes the signed change in \(H(S)\), reduction modulo \(M\) gives \(\delta\equiv1\pmod M\). Each ordinary count is odd, so \(\delta\) is even and \(|\delta|\ge M-1\). Consequently the larger of the two inner counts is at least \(M\). Cutting each full path at its \(n+1\) positions injects into the ordered two-slot covers, so \(Q_2(S)\ge (n+1)H(S)\). Apply this to the orientation with larger inner count to obtain the target-only bounds
\[
                         M\le n!,\qquad M(2+nM)\le N+2 .
\]
In particular a radix chosen as in the preceding theorem excludes a unit change inside the common block, at any target.

At \(N=E_{61}-1\), **all orders** require \(m\le58\). For \(m\ge59\) the second inequality gives \(n\le8\), while \(8!<2^{59}-1\). At **total order 123**, \(m=58\) gives \(n=64\), contrary to \(64(2^{58}-1)^2>E_{61}+1\). Thus \(m\le57\) there. For general perfect-shaped centers and \(p\to\infty\) the bound is
\[
 m\le p-\tfrac12\log_2 p+\tfrac12\log_2\log_2 p+O(1).
\]
It suffices to consider \(m\ge p/2\). The factorial bound gives \(n\gg p/\log p\), because otherwise \(\log_2 n!\le n\log_2 n<m-1\). Insert this bound and \(M\ge2^{m-1}\) into the second inequality. These exclusions apply to the placement of the marked arc inside the growing module of this precise cyclic return. The sharper congruences modulo 12 and 24 for the width-one return remain available. [Square-return residue theorem](../chanpak_arc_fiber_grammar_20261122/SQUARE-EARS.md)

## Exact small-order sizes

The strong width bound is sufficient for every input. Exhaustive small catalogs give a smaller cost.

| Initial order | Distinct slot polynomials | Distinct ordinary counts | Distinct counts after \(R_1\) | After \(R_2\) |
|---:|---:|---:|---:|---:|
| 5 | 8 | 7 | 8 | 8 |
| 6 | 24 | 19 | 24 | 24 |
| 7 | 156 | 77 | 139 | 156 |
| 8 | 1,728 | 320 | 1,442 | 1,726 |

At order eight there are 1,727 distinct pairs \((H,Q_2)\). The single width-three return separates all 1,727. The sole pair of distinct polynomials with equal first two arities is the archive’s \((345,5640)\) example. [First three-slot separation](../chanpak_crossing_scales_20261115/SLOT_SPECTRAL.md) Every width leaves exactly the same two profiles inseparable **whenever it separates all attainable pairs**; the cyclic formula depends only on those two readings. A further arity can tell these two profiles apart. Run `python3 -B physical.py` to repeat the table through order seven or `python3 -B physical.py --order-eight` for all 6,880 order-eight isomorphism classes.


## The smallest physical flip that stays hidden through one full return

**PROVED using finite-exact input.** Ten is the first order at which a physical flip can have the same ordinary count before and after one common width-one return, but differing counts at a later return. The distinction already appears after two returns, attaining the even-order divisibility bound 1,440.

Use the nine-vertex rows \((62,508,472,496,480,452,129,257,65)\), numbered from zero. Add vertex 9 beating the vertices in mask 17 or 81. These differ at bit 6, so the two strong hosts differ only at the arc between 6 and 9. Their first four slot readings are

| Mask | \(H\) | \(Q_2\) | \(Q_3\) | \(Q_4\) |
|---:|---:|---:|---:|---:|
| 17 | 2,637 | 75,216 | 928,341 | 7,121,472 |
| 81 | 2,637 | 75,216 | 927,621 | 7,115,712 |

For the first-minus-second order of subtraction, parity and the common degree-10 leading coefficient give the exact difference \(2q^2(q^2-1)(q^2-4)\), determined by the four readings. The first return gives 77,853 in both orientations. The second gives 2,230,125 and 2,228,685. Direct path counting on the two 10-, 12- and 14-vertex hosts checks these claims. At the analytic reading of this packet the unshifted difference is \(2/7-2+8/3=20/21\). [Earlier 32-vertex delayed example](../chanpak_arc_fiber_grammar_20261122/HIDDEN-FLIP.md) · [Physical replay](physical.py)

For minimality, root a marked flip at one of its endpoints. Deleting that endpoint gives an ordinary tournament on at most eight vertices. Take one representative of every such class, every contact mask, and every zero-to-one contact edit. This covers all physical flips up to simultaneous isomorphism, with repetitions allowed. To compute all four readings, let \(w_{ab}(j)\) sum \(2j P_{T-V(P)}(j)\) over even-order directed paths \(P:a	o b\) in the parent. An odd packing uses the new vertex either as a free point of color chosen in \(j\) ways or in one odd cycle whose removal leaves precisely one such path. Thus each new slot reading on mask \(I\) is exactly
\[
       jP_T(j)+\sum_{a\in I,\ b
otin I}w_{ab}(j).
\]
All 7,261,277 parent-mask-contact choices through order nine (including order two) have zero occurrences of equal \(H,Q_2\) but unequal \(Q_3\) or \(Q_4\). Through order nine four arities determine the full polynomial by parity. This rules out delayed separation at *any* later depth below order ten. Reproduce the finite lower bound from the repository root:

~~~sh
c++ -std=c++17 -O3 05-knowledge/results/chanpak_decimal_kernel_20261122/hidden_min.cpp -o /tmp/chanpak-slot-min
for n in 1 2 3 4 5 6 7 8
do
    /tmp/chanpak-slot-min "$n" < "session212/cycle-levels/classes$n.rows" || break
done
~~~

# Variable tails give an exact two-parameter +2 selector

This pass turns the old guard run into a free second parameter. It gives three rather than four physical readouts when a modulus is specified, and a complete absolute membership test for the resulting two-dimensional family. It does not yet give all odd absolute inputs.

## 1. An arbitrary transitive tail keeps the literal +2 arc

Under a single common apex use ordinal blocks \(X\to F\to V\). Here \(X\) is a transitive chain with binary contact word \(w\) ending in 0, \(F\) is one of the two six-site middle forks, and \(V\) is a transitive tail with contact word \(t\) starting in 1. Put
\[
 B=\operatorname{bin}(w),\qquad j=\operatorname{bin}(\overline{\operatorname{rev}(t)}).
\]
The endpoint constraints mean \(h_X=e_X=1\) and \(h_V=s_V=1\), in the [ordinal transfer notation](ENDPOINT-AND-SIGNED-DECK.md). The two-lane coordinates are \(R_X=B+1\) and \(J_V=j\): choose the last vertex of a nonempty left lane (incontact) and then freely allocate the earlier vertices in that lane; the reversed formula applies to \(R_X\). Thus the general ordinal formula gives
\[
 \begin{split}
 C_{X F}&=C_F+29(h(w)-1)+B J_F,\qquad R_{XF}=320B+R_F,\\
 H&=(J_F+320j)B+29(h(w)-1)+(R_F-29)j+29(h(t)-1)+C_F.
 \end{split}
\]
Between the two endpoints, only \(C_F\) in this final formula changes, from 111 to 113. This is an actual reversal in \(F\), not two independently constructed arbitrary counts. Every nonnegative even \(B\) and \(j\) is permitted: take the binary word or its reverse complement with the appropriate first/last contact, using the one-bit exceptional words for zero.

Let \(B=2v,\ j=2u\). Reflection preserves the one-apex height and \(h(2v)=2v+1-2q(v)\). Thus with the [two six-site signatures](GEWIRTZ-AND-DUAL-GUARDS.md) the lower endpoints are
\[
 F_A(v,u)=1280vu+382v+316u-58(q(v)+q(u))+111,\qquad
 F_B(v,u)=F_A(u,v).
\]
In particular transposing the two scalar flanks exchanges the two guard **readings**. This is not an assertion that arbitrary unmarked tournaments in the two chart names are distinct or identical. The short-guard charts are only the slice \(u=1\). The replay counts both full orientations of both modules with every tail of length at most five starting in 1, with \(w=1010\), for **62 pairs**. For tail 100 (\(u=3\)) the actual counts are \(22053\to22055\) and \(21921\to21923\).

## 2. Exact inverse for an arbitrary prescribed absolute count in this enlarged family

Given desired lower endpoint \(H\), set \(s=\min(v,u)\). The product term and \(q(z)\le z/2\) give \(s\le\sqrt{(H-111)/1280}\). Enumerate \(s\) and either ordering; for the other coordinate \(y\ge s\) the equation is now
\[
 H=D_s+L_s y-58q(y),
 \qquad L_s=1280s+316\text{ or }1280s+382 .
\]
Thus only the integer interval
\[
 \max\left(s,\left\lceil\frac{H-D_s}{L_s}\right\rceil\right)
 \le y\le\left\lfloor\frac{H-D_s}{L_s-29}\right\rfloor
\]
need be considered (reject negative numerator). For fixed \(s\) it meets at most two positive binary widths, and has size \(O(1+H/(s+1)^2)\).

Within an \(\ell\)-bit width choose a suffix of \(b\) bits, taking the largest \(b\ge0\) such that \(2^{\ell-1}\ge (L_s+29)2^{2b}\), or \(b=0\) if none. In a fixed high-prefix bucket \(y=P2^b+w\), use the [reverse-bit formula](TWO-COORDINATE-RIGIDITY.md#3-a-sublinear-exact-integer-inverse-for-that-short-guard-family). If \(b>0\), the height strictly decreases as the reversal of \(w\) increases: an increase by one in reversal contributes at least \(58\lfloor P/2\rfloor\) downward, while the entire possible adverse change from \(L_s w-58q(w)\) is at most \((L_s+29)(2^b-1)\). The extra parity term is nondecreasing before the minus sign. If \(b=0\), simply test each singleton bucket.

For large \(H\), \(D_s=O(s+1)\) and \(y\) is of order \(H/(s+1)\) wherever candidates exist, so \(2^b\) is of order \(\max(1,\sqrt H/(s+1))\) up to absolute constants. Consequently there are \(O(1+\sqrt H/(s+1))\) relevant high prefixes at fixed \(s\). The total is \(O(\sqrt H\log H)\); run binary search in each, and test exact equality and bounds on any result. With elementary bit scans this is \(O(\sqrt H\,\log^3 H)\) bit-sized arithmetic/scan steps (extra factors may be used for integer multiplication). The routine returns **all** exact coordinate pairs for this chart. To hit the upper endpoint of a physical pair, run it on \(H-2\). To restrict prefixes or widths, apply those tests to the returned coordinates. This is broader than a single-output short-guard inverse, but nonempty output is not guaranteed.

The replay includes 900 small exact inverse fixtures. It also tested \(E_p\pm1\), \(E_p=2^{p-1}(2^p-1)\), for \(p=7,13,17\). None occurs as either endpoint in this entire variable-tail family. This is only a finite chart exclusion, not a statement about general tournaments or about \(p=61\).

## 3. Three physical channels when the modulus is specified

Keep chart B at \(j=2\), chart A at \(j=2+2a\), and keep the [one-apex flank channel \(W\)](FANO-AND-FOUR-SEAM-DECODER.md). The normalized readings after removing each known constant obey
\[
 Y=798v-29q,\qquad X=(831+640a)v-29q,\qquad
 W=4v^2+\lfloor9v/2\rfloor-2q.
\]
For a prescribed modulus \(M\), strip powers of 2 and 5 from it to form \(m\). Choose \(a\equiv-32\,640^{-1}\pmod m\) (use zero if \(m=1\)). Then \(33+640a\) is a unit modulo \(2M\): it equals 1 modulo \(m\) and remains a unit at 2 and 5. From the three physical heights modulo \(4M\), normalize modulo \(2M\), recover \(v\) by dividing \(X-Y\) by this unit, then compute \(q\bmod M\) as
\[
 (798v-Y)-14(4v^2+\lfloor9v/2\rfloor-W).
\]
All states are reachable after any prefix by the preceding routing theorem. The extra tail has \(O(\log M)\) sites using the least nonnegative \(a\). The earlier four-channel minimality remains correct for **that fixed four-channel list**. Here the selected chart depends on the requested modulus.

Run `node two-tail-inverse.js` for the exact orientations, inverse fixtures and perfect-number chart exclusions. An optional comma-separated argument replaces the exponent list.

## 4. A complementary Moore sign obstruction

In the full Gewirtz approach let \(C=D-I\) be the [mod-3 square-zero sign matrix](GEWIRTZ-AND-DUAL-GUARDS.md) on the 57 branches. It has diagonal \(-1\) and off-diagonal entries \(\pm1\). No sign pattern invariant under a 57-cycle, **even up to independent type changes at the branches**, can be valid. If a cycle shifts \(C_{ij}\) with factors \(s_i s_j\), set \(w_{i+1}=\epsilon s_i w_i\) with \(\epsilon=\prod_i s_i\); the recursion closes because 57 is odd. Conjugating by the diagonal signs \(w_i\) makes \(C\) circulant. A sign change can complement the swap count for a branch pair, but preserves its forbidden extremes 1 and 55.

Over \(\mathbf F_3[x]/(x^{57}-1)\), write the circulant as \(f(x)\). The square-zero condition and \(x^{57}-1=(x^{19}-1)^3\) force \((x^{19}-1)^2=x^{38}+x^{19}+1\) to divide \(f\), as \(x^{19}-1\) is square-free. Therefore the coefficient rows repeat at distance 19. In any equal pair of rows, the direct route swaps type because of the diagonal \(-1\), but all other two-hop products are \(+1\). This gives exactly one type swap, impossible. Unlike a restriction on automorphisms of a full Moore graph, this closes a sign-template ansatz *before* filling in any group elements. It does not exclude noncyclic sign patterns.

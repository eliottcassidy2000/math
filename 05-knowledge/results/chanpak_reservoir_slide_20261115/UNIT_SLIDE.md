# A single, strong unit-slide pump

The [two-core carry pump](../chanpak_window_descent_20261115/GROWING_WINDOWS.md)
needed a pair of independently sliding core counts 13 and 29 and a bounded
Bézout sum of two relatively large steps. One marked **eight-vertex** core
now gives a literal unit half-step on its own. For every \(k\ge0\) it gives
\(k+1\) consecutive odd outputs of the same full marked profile and at the
same order \(11+12k\). All the output tournaments are strong.

## Fixed core, words, and matching response

Use outgoing rows
\[
  C=(4,1,130,7,143,159,63,75),\qquad
  \text{mark mask }35=\{0,1,5\}.
\]
A bit \(2^v\) is vertex \(v\). The ordinary spine count is \(h=45\), and
the two-slot count is \(Q=1440\). Marked start and end unsigned weights
are \(s=517,r=818\). The corresponding signed marked-start and
unmarked-end weights are \(d=31,f=25\), and the local apex count is 265.
Here a vertex is marked when the apex beats it, and the unmarked end
weight equals \(Q-h-r=577\).

The two 12-bit transitive words are
\[
 z=000100011001,\qquad w=010111000111.
\]
They have values \(281,1479\) and reversals \(2440,3642\).
The [two-ended slide identities](../chanpak_window_descent_20261115/MATCHED_SLIDE.md)
are literal:
\[
 4095\cdot517+45\cdot281=1440\cdot1479,\qquad
 1440\cdot2440=4095\cdot818+45\cdot3642.
\]
Since \(45\mid Q\) and \(\gcd(s,45)=1\), every same-core, same-mask slide
length \(d\) satisfies \(2^d\equiv1\bmod45\). The order of 2 modulo 45 is
12 (its orders modulo 9 and modulo 5 are 6 and 4). Thus the present tile
is as short as possible **for this return type**, not for other cores or
a return that changes the mask.

Give each base a left guard marked 1 and a right guard marked 0:
\[
  B(k,j)=1\;w^j\;C\;z^{k-j}\;0,\qquad 0\le j\le k.
\]
The guards are exactly the last mark of \(w\) and the first mark of \(z\).
The ordinal slide proof therefore matches both unsigned response
coordinates and both signed coordinates for each \(k\). Moving once
drops the full apex count by **two**, since the exact one-tile counts of
\(1Cz0\) and \(1wC0\) are 11,273,425 and 11,273,423, respectively.
One can also check this switch without a transfer formula by running
the literal 23-vertex subset DP in `replay_unit_slide_literal.cpp`.

The apex reaches the leftmost vertex and all spine vertices; all spine
vertices reach it through the last vertex. This proves strongness. To
move to any other order \(n\ge11\), use \(k=\lfloor(n-11)/12\rfloor\)
and append \(n-11-12k\) common unmarked singletons after the original last
guard, keeping a last vertex that returns to the apex. The common
extension shifts all counts by the same amount. This yields a consecutive
interval of \(\lfloor(n-11)/12\rfloor+1\) odds at **every** such order.

## Exact integer profile and prescribed local output

Let \(b=4096\), \(D_i=(b^i-1)/(b-1)\), and put
\[
 S_i = 1957 b^i+12645D_i,\qquad
 T_k = 5760 b^k-1726-2880\cdot2440D_k.
\]
The full common profile \((h,Q,S,d,T,f)\) is
\[
  (45,\ 5760b^k,\ 2S_k,\ 0,\ T_k,\ 0).
\]
The \((k,j)\) apex count is precisely
\[
 898+S_k+1655\sum_{i=0}^{k-1}S_i+5175k-2j.
\]
Thus a prescribed local output is decoded by one division by two
and the guard test \(0\le j\le k\).

For the formulas, `1 C` has \((N,S,d)=(867,1957,31)\).
One copy of \(z\) has ordinary count one, two-slot count 4096, apex
count 115, marked start weight 281, and unmarked end weight 1655
with signed unmarked end 1. The singleton mark equations give these
values directly. Appending \(z\) adds \(45\cdot114+1655S+d\) to \(N\),
updates \(S\mapsto4096S+12645\), and resets \(d\) to 45. The guard 0
adds \(S+d\), doubles \(S\), and resets \(d\) to zero. The accumulated
initial and following signed terms are \(31+45k\). The marked end
weight is \(2\cdot818+45+2880\cdot2440D_k\) by the reversed-word
formula, giving \(T_k\). In particular \(T_k\equiv29\bmod45\) at
every height.

## Plugging a changing left reservoir into this repair

Write \(B_+\) for the displayed count at \(j=0\). A marked left
reservoir with ordinary count \(g\), apex count \(n\), and total
marked start weight \(s_L\) combines with the single pump to give
exactly
\[
   45n+T_k s_L+g(B_+-45)-2g j,\qquad 0\le j\le k.
\]
The terminal signed-end channel of the pump is zero, so **the left
reservoir need not match, or control, its signed start channel**.
For a prescribed odd target \(N\) the literal repair succeeds iff
\[
\begin{split}
&45n+T_k s_L\equiv N\pmod g,\\
&B_+-2k\ \le\ 45+\frac{N-45n-T_k s_L}{g}\ \le\ B_+.
\end{split}
\]
The quotient automatically is odd: \(45n\) and \(g,N\) are odd and
\(T_k\) is even. Thus the ordinary count of the left reservoir scales
the repair by \(2gk\) but imposes exactly one modular admission
condition. This works for any actual left \(n,s_L\), with no unknown
additive error.

There is a simple literal way to program that residue, distinct from
placing the absolute window. Let the left reservoir end with any
core \(C_L\), of counts \(g,Q_L\), **all unmarked**. Precede \(C_L\)
with a nonempty transitive marked word whose binary value is \(x\),
with leading zeroes permitted. Its start weight is \(x\). The
combined left weight is \(s_L=Q_L x\), and the local left count is
\(n\equiv Q_L x\bmod g\): its other terms in the count transfer are
multiples of \(g\). Hence
\[
 45n+T_k s_L\equiv(45+T_k)Q_L x\pmod g.
\]
Every residue is programmable if \(\gcd(g,(45+T_k)Q_L)=1\), for
example \(g=5,Q_L=32\) at every \(k\); or more generally any prime
ordinary reservoir count that does not divide \((45+T_k)Q_L\).
For a strong glued tournament one would also require the first
vertex to be marked, choosing a longer representative if needed.
Neither the congruence construction nor leading-bit choices control
the required absolute interval. That remains a separate placement
inequality.

## A mixed four-tile return and its limitation

The complete compatible width-12 tile set for this *specific* marked
core with guard bits \((1,0)\), from the two exact readout equations,
is
\[
\begin{array}{c|c|r}
z&w&\Delta\\\hline
000100011001&010111000111&2\\
001100011001&010111010111&-48\\
010100011001&010111100111&34\\
011100011001&010111110111&-16
\end{array}
\]
A word of tile choices \(a_1,\ldots,a_k\) permits a mixed return
\[
  1\;w_{a_1}\cdots w_{a_j}\;C\;z_{a_{j+1}}\cdots z_{a_k}\;0.
\]
Within each such word the entire response remains fixed as \(j\)
varies, and each move drops the apex count by the corresponding
signed \(\Delta\). A consecutive run of \(\ell\) unit tiles thus
gives a certified repair radius of \(\ell\) half-steps in a carrier
whose other tiles can change its profile and baseline. This is an
exact nonuniform-family operation, not a new binary branching
assertion.

In fact this entire four-tile carrier by itself remains thin.
The \(j=0\) base has start count before its last guard at least
\(1957\cdot4096^k\); that guard adds at least this to its apex count.
Sliding can lower it by at most \(48k\). There are at most
\(4^k(k+1)\) variants. Hence through \(X\) these carry at most
\(O(X^{1/6}\log X)\) distinct counts. The placement work must either
use a genuinely larger left reservoir supply or a richer developing
core, not just all four tapes.

## A two-sided hyperbolic address equation

Keep one pump height and write \(Q,S,T\) for its common unsigned
response. The signed channels at both ends vanish. Prepend an
arbitrary transitive tape with mark-word integer \(x\); append an
arbitrary transitive tape of reverse unmark address \(y\). The
[exact binary surcharge](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md)
\(p(t)\) is given by \(p(0)=0\),
\(p(2u+1)=p(u)\), and \(p(2u)=p(u)+u+(u\bmod2)\).
The left tape's local apex count is \(1+p(x)\): take its converse,
which reverses and complements the tape, and apply the same decoder.
The right tape's local apex count is \(1+p(y)\). Thus the two ordinal
transfers give the literal count
\[
 B_+-2j+Qxy+Tx+Sy+45\bigl(p(x)+p(y)\bigr).
\]
Equivalently a physical target placement within this family is a factor
pair satisfying
\[
 Q\bigl(N-B_++2j-45[p(x)+p(y)]\bigr)+ST=(Qx+S)(Qy+T),
 \quad x,y\ge0,\quad 0\le j\le k.
\]
The correction is attached to the actual factors, so factoring the
uncorrected target or matching only a residue does not suffice.
Both factors can also be divided by two, giving nontrivial unit
residues modulo \(Q/2\). Indeed \(v_2(S)=v_2(T)=1\),
\(\gcd(Q,45)=45\), \(S\bmod45=3914\cdot4096^k\bmod45\) and
\(T\bmod45=29\); hence \(\gcd(Q,S)=\gcd(Q,T)=2\).
The left side divided by four must therefore be composite (both
normalized factors exceed one), although the correction is not
fixed independently of \(x,y,j\).

For a prescribed target, choose \(x\), put \(A=Qx+S\) and
\(Y=N-B_+-Tx-45p(x)\), and test only integers
\[
 \max(0,\lceil Y/(A+45)\rceil)
 \ \le y\le\lfloor(Y+2k)/A\rfloor .
\]
This follows from \(0\le p(y)\le y\); the final exact test is
that \(Ay+45p(y)-Y\) is even and between 0 and \(2k\).
When searching both directions, at least one of \(x,y\) is at most
\(\lfloor\sqrt{(N-(B_+-2k))/Q}\rfloor\), because all the additional
terms are nonnegative. One may therefore enumerate only this short
address on either side and solve for the other. This is an exact
placement algorithm in a fixed carrier, not a general existence result.

Even allowing every **one-sided** transitive suffix on the four-tile
bank does not fix its capacity. For each word of \(k\) tiles the common
\(S\) is at least \(3914\cdot4096^k\). An appended address \(y\) adds
at least \(Sy\), so the upper density through \(X\) of *all* its
one-sided suffixes, summing over all tile words and all slide
positions, is at most
\[
 \frac{1}{3914}\sum_{k\ge0}(k+1)\left(\frac4{4096}\right)^k
   =\frac{1}{3914(1-1/1024)^2}.
\]
The additive \(O(X^{1/6}\log X)\) term from suffix addresses \(y=0\)
and rounding follows from the same height bound as above.
This does not apply to the two-sided hyperbolic equation.

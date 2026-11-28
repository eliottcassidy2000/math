# A literal one-arc +2 pair in every odd square shell

2026-11-28. **Theorem.** For every odd \(h\ge5\) there are two strong tournaments, on the same labelled vertices, differing by one arc and with counts \(L,L+2\), such that
\[
h^2<L<L+2,\qquad L+3<(h+1)^2.
\]
For \(h\ge2^{21}\) there is a uniform algebraic decoder on a fixed core. The remaining 1,048,574 cases are covered by a finite replay. This upgrades the [quadratic Mills tree](../chanpak_branching_mills_shells_20261128): at every node its two children can be an actual one-arc pair. The [one-apex obstruction](../chanpak_square_one_apex_flip_density_20261128) remains intact; the repaired chart has a nontransitive deletion core.

The same fixed chart also proves a smaller absolute target window: **for every integer \(N\ge10^{10}\) there is a literal +2 pair in \((N,N+512N^{2/5})\)**. Neither statement is exact prescribed-integer coverage or a prime claim.

## 1. Reuse the core, not the other Farey blocks

Take the six ordinary vertices \(D\) from the [exponent-15 Farey witness](../chanpak_decimal_kernel_20261122/E15-FAREY.md), with outgoing masks \((16,5,9,19,38,15)\) and apex tags 35 (bit \(i\) means the apex beats vertex \(i\)). Reverse the pair \((1,4)\) for the lower host. The exact profile \((h_D,q_D,c_D,t_D,e_D,r_D,k_D)\) is
\[
(29,320,162,129,18,13,111)\quad\text{and}\quad
(29,320,162,129,16,11,113).
\]
Here \(h_D\) and \(k_D\) count full paths without and with the apex. Sum \(H(S)H(D-S)\) over subsets \(S\) for \(q_D\) (empty slots allowed); for \(c_D,t_D\), impose respectively that the \(S\)-path starts tagged or ends untagged. The full-path versions of ending tagged and starting untagged are \(e_D,r_D\). The replay recomputes every entry by subset DP, not by fitting samples of the polynomial.

Use a left ordered flank with tags given by the usual binary word of \(x\) and a right ordered flank with tags given by the reversed complement of the binary word of \(y\). Address zero means empty. All ordinary cross-arcs go left to \(D\) to right. For nonnegative even \(x,y\) the exact lower height is
\[
L(x,y)=111+129x+162y+320xy+29(f(x)+f(y)),
\]
where \(f(0)=0\) and \(f(2n+d)=f(n)+(1-d)(n+[n\text{ odd}])\). The upper is \(L+2\).

For clarity the needed [ordinal transfer](../chanpak_decimal_kernel_20261122/LOCAL-WORDS.md) can be seen directly. For ordinary blocks \(X\to Y\) with a shared apex,
\[
k'=h_Y k_X+h_X(k_Y-h_Y)+c_X t_Y+e_X r_Y,\qquad
c'=q_Yc_X+h_Xc_Y.
\]
For the first identity separate full paths with no \(Y\) before the apex; then paths with some \(Y\) before but no \(X\) after the apex; then the rest. The first case gives \(h_Yk_X-(h_X-e_X)r_Y\). The second gives \(h_X(k_Y-h_Y+r_Y)\). The last independently picks a post-apex \(X\)-slot starting tagged and a pre-apex \(Y\)-slot ending untagged. For \(c'\), the selected slot starts in \(X\) or uses only \(Y\). The left flank has \(c_X=x,\ e_X=0,\ k_X=1+f(x)\), and the right flank has \(t_Y=y,\ r_Y=0,\ k_Y=1+f(y)\). The readings \(c_X,t_Y\) count arbitrary subsets after the first tag or before the last nontag, respectively, giving positional binary weights. Apply the identities first on the left and then on the right; empty flanks can just be omitted. The two varying readings \(e_D,r_D\) cannot enter on this even chart. Both versions of \(D\) are strong and have apex contacts in both directions. A nonempty left flank is reached from the apex at its first site and exits through \(D\); a right flank is reached through \(D\) and its last site beats the apex. This makes both full hosts strong.

At the formerly blocked parent \(503\), take \(x=20,\ y=38\). The two heights are \(253{,}265,253{,}267\) on 18 vertices, checked by full graph subset DP.

## 2. The balanced decoder

For a power of two \(Y\ge2\) set \(x=Ya, y=Yb\). Telescoping the bit recurrence gives \(f(Yz)=(Y-1)z+f(z)+[z\text{ odd}]\); also \(0\le f(z)+[z\text{ odd}]\le z\), by induction. Write
\[
A=320Y^2,\quad B=158Y-29,\quad C=191Y-29,\quad
p=C/A,\quad q=B/A,\quad E=111-BC/A.
\]
Then
\[
L= A(a+p)(b+q)+E+\epsilon,\qquad 0\le \epsilon\le29(a+b).
\]
For a real center \(T\), let \(r=\sqrt{(T-E)/A}\), \(m=\lceil 2r-p-q\rceil\), and \(s=m+p+q\). On the line \(a+b=m\), the smaller real root producing exact quadratic part \(T\) is
\[
a_0=s/2-p-\sqrt{s^2/4-r^2}.
\]
Round \(a_0\) to the nearest integer \(a\), then set \(b=m-a\). Since \(2r\le s<2r+1\), the radical is at most \(\sqrt{r+1/4}\). Changing \(a_0\) by at most \(1/2\) changes the quadratic part by at most
\[
A\big(\sqrt{r+1/4}+1/4\big).
\]
Both coordinates after rounding are positive whenever \(r\ge4\), using \(p,q<0.6\), \(s/2\ge r\) and the radical bound. No floating-point oracle is required. When \(T=t/d\), put \(M=A(t-111d)+BCd,\ S=Am+B+C\); compute \(m=\lceil(\sqrt{4M/d}-B-C)/A\rceil\) and \(a=\lfloor(S-2C+A-\sqrt{S^2-4M/d})/(2A)\rfloor\). In both formulas a square root under the ceiling, or subtracted under the floor, can first be rounded up to an integer with no change to the rounded outcome. The replay implements these integer square-root comparisons.

**All large square shells.** Set \(Y=4\) and \(T=h^2+3h/5\). For \(h\ge2^{21}\) one has \(r<(h+1)/71\). Thus the rounding error is \(<0.45h\) and \(29m<0.82h\). For example, the rounding coefficient bound follows by plugging \(2^{21}\) into \(5120(\sqrt{(h+1)/71+1/4}+1/4)/h\), a decreasing function for positive \(h\). It follows that
\[
h^2+0.15h<L<h^2+1.87h.
\]
This leaves room for \(L+3\) strictly below \((h+1)^2\). The constructed hosts have \(O(\log h)\) vertices.

**Sub-square-root target window.** For integer \(N\ge10^{10}\), choose \(Y=2^{\lfloor(\log_2 N)/10\rfloor}\) and \(T=N+256\lceil N^{2/5}\rceil<2N\), implemented using bit length and an integer fifth-root comparison. Then \(N^{1/10}/2<Y\le N^{1/10}\), \(0<E\le111\), and
\[
4<r<0.08\sqrt N/Y .
\]
The rounding error is at most \(91Y^{3/2}N^{1/4}+240Y^2\le93.4N^{2/5}\). Here \(\sqrt{r+1/4}\le\sqrt r+1/2\), and \(240 N^{1/5}\le2.4N^{2/5}\). The nonlinear correction is \(<29(2r+1)\le9.28N^{2/5}+29\). Thus \(|L-T|<105N^{2/5}\), which proves \(N<L<L+2<N+512N^{2/5}\). The balance comes from making the two errors \(Y^{3/2}N^{1/4}\) and \(\sqrt N/Y\) have the same power of \(N\); it retains the literal arc, unlike the original [free reflected seam](../chanpak_wheel_square_selector_20261117/FREE_SEAM.md) that exchanges addresses.

## 3. Finite coverage and the square tree

The deterministic replay processes every odd \(5\le h<2^{21}\). It accepts a lower count only on the exact integral test \(h^2<L,\ L+3<(h+1)^2\). The counts by first witness are:

| First witness | Parents |
|---|---:|
| One-apex four-parity catalogue through \(p=741461\) | 1,026,067 |
| Balanced core decoder with two positive flanks | 22,484 |
| Core chart on the \(x=0\) axis | 21 |
| Explicit small witnesses | 2 |

The finite balance uses \(Y=4,2,8,16\), the six centers \(h^2+ch\) with \(c=0.6,0.3,0.9,1.2,1.5,0.1\), the nominal sum \(m\) or its next two integers, and the two nearest roots with offsets \(0,-1,1,-2,2\). Approximate roots only suggest finite candidates; the acceptance inequalities and counts are evaluated exactly with safe-sized integers. On the axis, increasing even \(y\) up to the linear bound suffices. For \(h=5\) use order-6 masks 80 and 336, yielding \(27,29\). For \(h=25\) use order-8 masks 98750095 and 98749967, yielding \(657,659\). These are lexicographic forward masks and the replay checks both full literal graphs, strong connectivity and the differing arc.

Starting at \(h=5\), choose either successor indefinitely. The closed intervals \([h_j^{1/2^j},(h_j+1)^{1/2^j}]\) strictly nest and shrink, with disjoint options at each branch. Their unique common point lies in the interior of every interval, so the floors are genuine. The one-arc statement is between siblings at a given node, not between a parent tournament and a child.

Run node replay.js. The large decoder and window checks include targets with 100 digits, in addition to the complete small sweep and literal graph checks.

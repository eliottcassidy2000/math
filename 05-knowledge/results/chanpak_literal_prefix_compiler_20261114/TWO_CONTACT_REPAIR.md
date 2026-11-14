# Two changing contacts give actual bivariate unit repairs and new prefix sources

Extend the [arbitrary-root bridge](ARBITRARY_ROOTED_BRIDGE.md)
by allowing its first new terminal \(w\) to beat a specified
subset \(U\) of the old vertices, not just one marked vertex.
Keep \(r\ge1\) forward transitive followers. The second, gearing
bridge has exactly one contact (at \(w\)) and \(s\ge1\) followers.
All widths are positive and all outputs are literal tournaments.

## 1. Exact endpoint correction and inversion

Write \(h=H(T)\), \(R=\sum_{u\in U}R_+(T,u)\), and let \(E\)
sum, over old even-size subsets \(S\), the number of Hamiltonian
paths of \(T[S]\) whose two endpoints are in \(U\), weighted by
\(H(T-S)\). Let \(x=2^r,\ y=2^s,\ S_0=h+xR\).
The first bridge satisfies
\[
                 H'=S_0-2E,\qquad R_+'(w)=S_0.
\]
Indeed the odd-cycle expansion through \(w\) uses a path starting
in \(U\), followed by an ascending chosen subset of followers.
Half of the follower choices have the right cycle parity.
The only invalid counted closure has **no** followers and the
old path ends in \(U\), where the arc points from \(w\) instead
of back to \(w\). Such a cycle would be odd precisely at an
even old path size. Removing those terms is exactly \(2E\).
The \(w\)-start subset-path sum has no closing arc and therefore
no subtraction. This proves that after the second bridge its
ready deletion profile and ordinary count are
\[
 (a,b)=(S_0-2E,yS_0),\qquad H=(1+y)S_0-2E.
\]
For strong \(T\) and nonempty \(U\), the first bridge is strong
(the new terminal reaches \(T\) and every follower, and they can
reach it). Then the second bridge preserves strongness.

For a fixed actual pair \((T,U)\), a prescribed tuple \((a,b)\)
is realized exactly when
\[
       (a+2E-h)/R \quad\text{and}\quad b/(a+2E)
\]
are integer powers of two at least 2. A prescribed ordinary
total \(N\) is tested by factoring
\[
                         N+2E=(1+y)(h+xR)
\]
and checking the two positive power-of-two conditions. Thus both
widths vary independently and the test is on the absolute target.

## 2. A true adjacent prime and 7-times-prime rail

Use the same strong five-vertex \(T\) with outgoing masks
\((24,17,3,6,12)\). The two contact sets
\(U_4=\{0,2\}\) and \(U_5=\{0,1\}\) have respective profiles
\((h,R,E)=(15,42,4)\) and \((15,42,5)\).
They change only two arcs and give **exact** ready difference
\((2,0)\), hence ordinary totals differing by 2 for every pair
of positive widths. This is a two-dimensional physical unit
repair with a nonunit source count, unlike a free width on
a single transitive stem.

At \(x=16,y=2\) the states are \((679,1374),(677,1374)\).
The literal 12-vertex totals are **2053 prime and
2051 = 7 times 293**. The \(E=5\) rail is divisible by 7
exactly when \(y\equiv2\pmod7\), i.e. \(s\equiv1\pmod3\),
while the adjacent \(E=4\) rail is never divisible by 7.
Both are units modulo 3 at all widths. The congruences
follow from \(15+42x\equiv1\bmod7\) and \(0\bmod3\).
They do not show that the cofactors on the 7 side or
counts on the adjacent side are prime infinitely often.
The replay also checks a second core giving the smaller
adjacent outputs 677 prime and 679 = 7 times 97.

## 3. Bypass the terrace boundary on a real prescribed target line

Set \(y=2\). The resulting physical state has exact
near-boundary defect
\[
                    d=(b-2a)/2=2E.
\]
For either contact pair above, \(d\) is fixed even
(\(8\) or \(10\)) and \(a\) grows arbitrarily. In the prior
[physical terrace boundary](PHYSICAL_BOUNDARY.md),
\(a>6d>0\) forces a single final \(q=1,r=1\) from a physical
prestate whose first coordinate equals \(d\) and is odd.
Our states therefore **cannot** be terrace prefixes for any
sufficiently large width, even though they are ready for exactly
the same future transfer. This is a different escape mechanism
from the dyadic-slope gap example in the original bridge.

It also works at the target-selected maximal depth. With the
previous \((u_k,v_k)\) and \(A_k=u_k+2v_k\), take the least \(r\)
giving \(a=15+42\cdot2^r-2E\ge A_k\). For these cores
\(A_k\le a<2A_k\) and
\[
    p=u_ka+v_k b=A_k a+4v_k E,\qquad
    A_k^2\le p<3A_k^2<A_{k+1}^2.
\]
The first inequality uses \(a\ge A_k\); the upper bound uses
\(v_k\le A_k/2\) and \(A_k>2E\). The inverse therefore assigns
the correct \(k\), but necessarily has to move past the least
progression point because \(a\ge A_k>2v_k\).

For example, at \(E=5,r=5,s=1,k=1\), the prime cofactor is
\(p=1,229,903\). The least relaxed point is \((77,3802)\).
The new physical 13-vertex prefix is \((1349,2718)\), with
even defect 10 and exact literal count 4067. Apply the one
outer \(M\) and the length-7 terminal to reach
\(H=8,609,321\) on 26 vertices. The prime cofactor and prefix
are checked in the literal replay. We make no claim that all
other points on this example's target line are impossible,
or that every line gets a positive hit.

An ordinary \(\pm2\) difference at this ready source does
**not** remain a unit difference after a shared terminal-seven
continuation: its cofactor difference becomes \(2u_k\).
Keep that distinction when using this as a possible
prime/cofactor coupling.

Run two_contact.py for the direct subset DP, the exact inverse,
the new target line, and every two-contact choice on labelled
tournaments through four vertices. The original two-contact
identities and the adjacent prime pair were found in background
research; the even-defect escape and maximal-depth construction
are consequences developed here.

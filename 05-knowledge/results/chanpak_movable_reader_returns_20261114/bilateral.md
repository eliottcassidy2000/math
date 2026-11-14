# A repeatable exact return by moving four contacts across a core

A carry in the reader can cancel a change of central tournament, keeping
the same absolute Hamiltonian count while changing its rooted response.
Stronger, a four-bit exchange between the wings can be repeated arbitrarily
many times while keeping the full rooted profile. Both mechanisms have
exact target inverses in their chart.

## 1. The two-sided interface

The binary-spine one-apex function \(g\) has \(g(0)=1\) and
\[
g(2v+1)=g(v),\qquad g(2v)=g(v)+v+(v\bmod2).
\]
It counts paths for the binary contact word. Leading 0s and trailing 1s
do not change the count. Inductively \(1\le g(x)\le x+1\): if an even
step has odd predecessor \(v\), use \(g(v)\le(v+1)/2\).

Take a nonempty central tournament \(B\) and apex contact mask \(U\subseteq B\).
Put transitive left and right wings around \(B\), orient all inter-block
arcs forward, and add an apex \(v\). A contact bit is 1 iff the apex
beats that vertex. Use the binary bits of \(x\) on the left, \(U\) on \(B\),
and the reverse complement of the bits of \(y\) on the right. An empty
wing is zero. Leading left 0s and trailing right 1s were discarded as
neutral sources and sinks. If \(x,y>0\), the first left and last right
vertices give \(v\to\mathrm{first}\to\mathrm{last}\to v\); every vertex
reaches and is reached from the apex, so the graph is strong.

Write \(h=H(B)\) and \(f=F_2(B)\), where \(F_2\) counts ordered
two-path covers (allowing empty paths). Let \(P\) sum the rooted
positive responses \(R_+(B,u)\) at \(u\in U\), and \(Q\) sum the same
responses in the converse at \(u\notin U\). A rooted response sums
the number of paths on \(W\) starting at the root, times \(H(B-W)\),
over \(W\) containing the root. Let \(S\) count full \(B\)-paths
ending in \(U\), and \(T\) count those starting outside \(U\).
These are the negatives of the corresponding signed responses by the
[signed endpoint identity](https://github.com/eliottcassidy2000/math/blob/2d70eb9cd31437bbfa8cf418c6285f79fc8afaaa/05-knowledge/results/chanpak_11_13_26_ready_barrier_20261114/READY_CLASSIFICATION.md).
Finally \(K=H(B+v)\) with no wings. Put \(p(b)=b\bmod2\).
The full count is
\[
 N=K+h(g(x)+g(y)-2)+Qx+Py+fxy+T p(x)+S p(y). \tag{1}
\]
Use the expansion over packings of disjoint directed odd cycles.
An inter-block cycle must pass through the apex, so there is at most
one. Apex cycles staying within one old block account for \(K\) and
the two wing increments weighted by \(h\). A crossing cycle consists
of paths in the visited old blocks in block order. The selected old
vertices must number an even amount; the parity factor one half
cancels the cycle weight two. The unsigned choice for a left start
is \(x\), and for a right end \(y\), since the extremal selected vertex
fixes one positional bit and later or earlier vertices in that wing
are free. Thus the unsigned adjacent terms are \(xQ\) and \(Py\);
the left-to-right term is \(xfy\), permitting every two-path split
in the center. The signed split of a nonempty interior old block
vanishes (the two colors have weights 1 and -1), so only adjacent
wing/center pairs contribute a signed term. Within a transitive wing,
any optional vertex pairs off signed choices; only the adjacent
boundary singleton may survive, with sign \(-1\) and weight \(p(x)\)
or \(p(y)\). Its product with the signed central endpoint gives
\(T p(x)\) or \(S p(y)\). Packings omitting an apex cycle give the
base \(h\), included once in \(K\).

Put \(m=\operatorname{bitlen}(y)\), and let \(\operatorname{rev}_m\)
reverse exactly \(m\) bits (\(\operatorname{rev}_0(0)=0\)).
The apex-root observables are
\[
 R_+=(fx+P+h)2^m-h\operatorname{rev}_m(y),
 \qquad R_-=-h\quad(y>0). \tag{2}
\]
In \(R_+\) the apex singleton contributes \(h\); paths entering the
left or center first contribute \(fx\,2^m\) or \(P\,2^m\).
Paths starting in the right wing contribute \(h\) times its contact
word, \(2^m-1-\operatorname{rev}_m(y)\). For the signed endpoint,
if \(y>0\) every Hamilton path in the old part has its last right
vertex terminal, and that vertex points to the apex. Exactly \(h\)
full paths therefore end at the apex.

## 2. Matched physical centers and neutral carries

Central lists are outgoing vertex bit masks; \(U\) is also a mask.
Each pair shares all six coordinates \((h,f,P,Q,S,T)\).

| center masks | \(U\) | \(K\) | common profile |
|---|---:|---:|---|
| 18,20,9,3,12 | 3 | 41 | 15,120,42,63,6,9 |
| 18,20,9,3,12 | 5 | 43 | same |
| 56,57,51,4,40,8 | 17 | 65 | 17,224,72,135,6,11 |
| 54,32,26,51,2,20 | 10 | 69 | same |
| 100,13,104,1,47,74,26 | 29 | 373 | 95,1276,694,445,62,39 |
| 100,13,104,1,47,74,26 | 57 | 375 | same |

The first pair gives a persistent \(+2\) switch for arbitrary choices
of both wings, with \(R_\pm\) unchanged when \(y>0\).
For a neutral return, compare \(y=2q\) with \(y+1=2q+1\).
By (1), this reader change alone shifts \(N\) by \(P+fx+S-h(q+p(q))\).
To cancel a center increment \(\Delta\), solve
\(fx\equiv-(P+S+\Delta)\pmod{2h}\), then take the even quotient
\((fx+P+S+\Delta)/h\) as \(q\), or one less. Examples are
\[
\begin{array}{c|c|c}
 K\to K+\Delta &x&y\to y+1\\ \hline
65\to69 &1+17r&36+448r\to37+448r\\
373\to375&7+95r&204+2552r\to205+2552r
\end{array}\qquad(r\ge0).
\]
The first exchanges the actual central tournament. The right change
flips exactly the apex arc at its first right vertex. Its width \(m\)
is unchanged, so (2) shows \(R_+\) changes by \(-h2^{m-1}\).
At \(r=0\), the first row gives \(N=11411\), with \(R_+\) changing
from \(19879\) to \(19335\); the second gives \(N=1982811\) and
\(R_+\) changing from \(2483731\) to \(2471571\). Both rails are
unbounded, but they do not cover every target.

## 3. An arbitrarily repeatable full-profile return

Keep either \(h=15\) center. For \(X,Y>0\) compare wing choices
\[
 (x,y)=(16X+7,Y)\qquad\hbox{and}\qquad (x,y)=(X,16Y+8).
\]
The first ends the left wing with binary 0111; the second adds 1000
to the reverse-complement code on the right. The responses \(R_+\)
are equal identically: reversing the appended right nibble adds
\(2^{\operatorname{bitlen}(Y)}\) to \(\operatorname{rev}(Y)\), and
\(120\cdot7=15(42+15)-15\) cancels the offset in (2). Also
\(R_-=-15\) on both sides. The binary recurrence gives
\[
 g(16X+7)=g(X)+X+p(X),\qquad
 g(16Y+8)=g(Y)+8+14Y.
\]
Substitute into (1). The first absolute ordinary count minus the
second is exactly \(6(p(X)+p(Y)-1)\); all \(XY,X,Y\) terms cancel.
Thus opposite parity of \(X,Y\) suffices for a full-profile return.

For every \(d\ge0\), odd \(X>0\) and positive even \(Y\), define
\[
 x_i=16^{d-i}X+7(16^{d-i}-1)/15,\qquad
 y_i=16^iY+8(16^i-1)/15,\qquad 0\le i\le d.
\]
Each step has an odd left prefix and even right word, so repeatedly
changes both wings while preserving \((H,R_+,R_-)\) exactly.
The \(d+1\) rooted outputs are pairwise nonisomorphic: deleting the
apex leaves a unique nonsingleton strongly connected component
(the regular 5-vertex core), whose position in the ordered component
list moves four places at each step. In fact they are pairwise
nonisomorphic even up to converse. In the forward direction only
four transferred vertices and five core vertices change their degree
scores; the apex beats three of the transferred vertices before and
after. Subtract \(m=\operatorname{bitlen}(Y)\) from these nine scores.
The local multisets before and after a step are
\[
 (9,7,6,5,2,2,3,3,3),\qquad (6,6,7,7,7,3,2,1,1).
\]
The means agree and the squared sums are 226 and 234, so every step
raises the entire squared outdegree sum by 8. Converse preserves
the squared sum since it replaces each outdegree \(a\) by \(n-1-a\)
and the total outdegree is \(n(n-1)/2\).
At \(d=1,X=1,Y=2,K=41\) the two 13-vertex graphs have profile
\((7163,11253,-15)\).
The \(K=43\) center gives 7165 with the same other channels.
At \(X=1,Y=20\), the four-way construction yields \(57839\)
and \(57841=7\cdot8263\), with 57839 and 8263 prime.

An arbitrarily long bank of equivalent sources preserves subsequent
one-contact bridges because the full rooted profile, not just \(H\),
is fixed. The exchange itself cannot change the target: the two free
wing words instead choose the baseline.

## 4. Absolute target inverses

**Prescribed pair \((N,R_+)\).** For a center \(P>0\), enumerate
\(0\le m\le\operatorname{bitlen}(\lfloor N/P\rfloor)\) by (1).
Since \(h\) is odd, (2) uniquely determines
\(\operatorname{rev}_m(y)\equiv -R_+/h\pmod{2^m}\).
Reverse those \(m\) bits, require \(\operatorname{bitlen}(y)=m\),
solve the exact division (2) for a nonnegative integer \(x\), then
test (1). This costs \(O(\log^2 N)\) arithmetic operations for
Euclidean inversion and word replay across the widths (for
physically possible \(R_+=O(N^2)\) with fixed center). A subsequent
positive bridge target \((a,b)\) is tested at \((a,b/2^r)\) for each
\(1\le r\le v_2(b)\). This is a complete physical-source test in the chart.

**Count \(N\) alone.** For \(P,Q>0\), try the smaller wing
\(u\le\sqrt{N/f}\), including 0. Try both orientations (swapping
\(P,Q\) and \(S,T\)). The larger coordinate \(z\ge u\) must solve
\[
 A z+h g(z)+b p(z)=N-B+h,\quad A=P+fu,\quad b=S,\quad
 B=K+h(g(u)-1)+Qu+T p(u).
\]
We have \(z\le Y=\lfloor(N-B)/A\rfloor\). Large carries prevent
direct bisection on arbitrarily unbalanced wings. Instead peel
\(r\) low bits of the larger variable, with \(D=2^r\) the least
power satisfying \(AD^2>2hY\) and \(AD>2h\). Try each suffix \(j<D\)
that can occur. The substitutions \(z=2w\), \(z=2w+1\) turn
\((A,b)\) into \((2A+h,h)\), \((2A,0)\) respectively, subtracting
\(A+b\) from the target in the odd case. After peeling, the linear
coefficient is at least \(AD\) and the prefix is at most \(Y/D\).
The remaining down-step of the nonlinear terms is bounded by
\(hY/D+h<AD\) using \(1\le g(w)\le w+1\) and \(0\le b\le h\);
bisection is now exact. At fixed center \(D=O(1+\sqrt N/(u+1))\).
There are \(O(\sqrt N\log N)\) suffix trials, and
\(O(\sqrt N\log^3 N)\) arithmetic cost with direct word replay.
To request a depth-\(d\) return, filter the decoded witnesses for
an even right word and the repeated suffix 0111 on the left, leaving
an odd positive prefix. These inverses certify chart failures as
well as returning physical witnesses. A bank guaranteed to pass
for all admissible absolute targets is still missing.

Run python replay_bilateral.py. It verifies central profiles,
small physical wings by Hamiltonian subset DP, bounded inverse
tables, carry rails, and the repeated four-bit identities. The
odd-cycle expansion, signed endpoint identity, one-apex recurrence
and empty-center half-inverse are inherited. The new ingredients
are the matched intransitive centers, compensating reader carry,
four-bit exchange, and exact absolute inverses with no balance
restriction.

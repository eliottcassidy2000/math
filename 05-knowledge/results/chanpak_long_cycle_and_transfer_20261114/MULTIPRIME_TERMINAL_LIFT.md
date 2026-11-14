# Four-prime local freedom in a two-clock physical tournament grammar

**Proved.** For every odd integer \(n\) and arbitrary nonnegative
\(A,B,C,D\), a **strong tournament** in a grammar with eight rooted base
cores, at most nine initial followers and two variable terminal clocks
has Hamiltonian path count \(H\) divisible by seven and
\[
                       H/7\equiv n\pmod{2^A3^B5^C7^D}.
\]
This includes cases where \(n\) is divisible by 3, 5, or 7. It is
not an absolute realization theorem: the grammar supplies at most
\(O((\log X)^2)\) distinct counts up to \(X\).
Separately, the absolute target 180131, for which **every** length-seven
terminal on a ready source fails, has a literal length-three terminal
on a different ready source.

## 1. The literal construction

Use the inherited arbitrary-root positive bridge and Tribonacci
kernel. A width-\(r\) bridge adds \(r\) transitive followers and
then a vertex \(v\) beating only the old root; all other arcs point
forward between old and new vertices. Width \(t\ge1\) turns a core of
count \(a\) and positive split profile \(P\) into a ready pair
\((a,b)=(a,2^tP)\). Then attach \(q\) width-zero bridges, each
at the current root, and a terminal width-\(r\) bridge. Put
\[
 C_q=\begin{pmatrix}e_q&f_q\\f_q&g_q\end{pmatrix},\qquad
 (e,f,g)_{0,1,2}=(1,1,1),(1,2,2),(3,3,4),
\]
with each sequence satisfying the Tribonacci recurrence. The count is
\[
 H=aU(q,r)+bV(q,r),\qquad F(q,r)=(U,V)^t=C_q(1,2^r)^t.
\]
The cores are strong or singleton. Each positive bridge produces a
strong tournament, since the old root reaches every follower directly,
the last vertex returns to the old root, and every follower beats
the last vertex. Each width-zero bridge also produces a strong
tournament when applied to a strong source of size at least two.

## 2. Two-dimensional local clocks at 7 and 5

For a Tribonacci state column the step matrix is
\(M=((0,1,0),(0,0,1),(1,1,1))\). Suppose \(M^L=I\bmod p\)
and \(2^R=1\bmod p\) for an odd prime \(p\).
On a base clock \(q=q_0\bmod L,r=r_0\bmod R\), the two columns
\[
 J_q=\frac{F(q_0+L,r_0)-F(q_0,r_0)}p,\qquad
 J_r=\frac{F(q_0,r_0+R)-F(q_0,r_0)}p
\quad\bmod p
\]
give the digit-lifting map. Indeed at stage \(k\ge1\), changing
clocks by \((iLp^{k-1},jRp^{k-1})\) changes \(F\) by
\(p^k(iJ_q+jJ_r)\) modulo \(p^{k+1}\). This is the binomial
expansion of \((I+pE)^{p^{k-1}}\equiv I+p^kE\bmod p^{k+1}\)
and its scalar counterpart, with their mixed term zero. The
columns depend only on the base clock. If the Jacobian is invertible,
every target vector with the appropriate base residue lifts uniquely
digit by digit.

At 7, \(L=48,R=3\). The full list of base clocks with
\(F\equiv(0,0)\) is
\[
            (7,2),(12,0),(23,2),(28,0),(39,2),(44,0).
\]
Their Jacobian determinants are \(5,4,3,1,6,2\bmod7\),
respectively. Thus on each clock the normalized coefficients
\((U/7,V/7)\) range independently over all residues at every
precision (equivalently, each defines a bijection on two
7-adic parameters). Use only \(q_0=12,28,44,r_0=0\) for the
simultaneous result. Our choices of \(a\) are 7-adic units, so
we can prescribe \(V/7=0\) and \(U/7=n/a\bmod7^D\).

At 5, \(L=31,R=4\). All **25** coefficient vectors occur at
nonsingular base clocks with \(r_0\in\{1,3\}\).
The replay gives a representative for each vector, so
\(H\equiv7n\bmod5^C\) can always be arranged with \(r\) odd
provided \(a,b\) are not both divisible by 5.

The initial three state vectors have determinant 1, so verifying the
three coordinate recurrences at offsets \(L,L+1,L+2\) proves the full
matrix identities \(M^L\equiv I\bmod p\). The finite mod-5 and
mod-7 clock tests need only 62 and 144 base clocks.

## 3. Compatibility with 2 and 3

For the eight source choices use the short literal roots and profiles
in replay_terminals.py. Their \(a\)-values in residue order modulo 16
are \(1,19,5,23,25,11,13,31\). Choose \(a\equiv7n\bmod16\);
take \(b=2^tP\) with \(4\le t\le9\). Both \(a,P\) are 3-adic
units, and \(a,b\) are not both divisible by 5.

**At 2.** Take \(r\ge\max(5,A)\), making the width term vanish
at the required precision. For \(q_0=12,28,44\),
\(e_{q_0}\equiv1\bmod16\) and
\[
                    (e,f,g)_{q_0+48}-(e,f,g)_{q_0}
                           \equiv(16,8,24)\pmod{32}.
\]
These differences remain the same across each base clock, since
\(M^{48}-I\) is divisible by 8 and its square vanishes modulo 32.
More generally
\[
 M^{48\cdot2^j}\equiv I+2^j(M^{48}-I)\pmod{2^{j+5}}.
\]
The nonlinear binomial terms vanish: the \(k\)-th term for \(k\ge2\)
has valuation at least \(j-v_2(k)+3k\ge j+5\).
Thus changing \(q\) by \(48\cdot2^j\) adds exactly
\(2^{j+4}\) to \(ae_q+bf_q\bmod2^{j+5}\), since \(16\mid b\)
and \(a\) is odd. Every residue congruent to \(a\bmod16\)
therefore lifts on this clock, by fixing \(q\) modulo
\(16\cdot2^{\max(A-4,0)}\) after its residue modulo 3 is fixed.

**At 3.** \(M^{13}\equiv I\bmod3\) and \(2^2\equiv1\bmod3\).
Start with \(r=3\bmod6\). If \(7n\) is a 3-adic unit, take \(t=4\):
the finite test finds \(q'_0\bmod13\) giving the correct scalar
modulo 3 and a unit \(q\)-derivative. Lift \(q\) alone, with \(r\)
fixed to 3 modulo the required scalar period.
If \(7n\equiv0\bmod3\), choose \(t\in\{4,\ldots,9\}\) and
\(q'_0\bmod13\) to get the correct scalar modulo 9 with a
unit \(r\)-derivative; the tests certify this for 0, 3, and 6
for every source. Lift \(r\) alone from modulus 9 onward,
preserving \(r\equiv3\bmod6\), and keep \(q=q'_0\) modulo the
required Tribonacci period. A unit \(q\)-derivative alone
cannot handle the zero class.

Choose \(q_0\in\{12,28,44\}\) with the \(q\bmod3\) selected by
the 3-lane. The 7-lane period for \(q\) shares only the factor 3
with the 3-lane period, and their residues agree. The 5-lane
period for \(q\) is coprime to both; the 2-lane shares only
the 2-part of the base 7-clock and uses that same clock.
For \(r\), the 3-lane shares only 3 with the 7-lane and only
2 with the 5-lane. The residues agree since \(3\mid r\) and
\(r\) is odd. Apply CRT, then add a common multiple to \(r\)
for the width lower bound. This is a physical construction
at every finite precision, not just a formal limit.

## 4. Exact escape from the length-seven obstruction

Use the strong core with outgoing masks
\((6,28,104,113,165,195,147,15)\) and root 5. It has
\((a,P)=(563,1066)\). Take \(t=3\), so \(b=8528\),
and terminal clocks \(q=3,r=1\). The coefficients are
\((17,20)\), giving \(17\cdot563+20\cdot8528=180131\).
The replay counts paths independently on all 17 literal vertices.
So changing the core and terminal family can escape the exact
length-seven obstruction without departing from ready geometry.

Positive Tribonacci coefficients grow at least exponentially with
\(q\), and the width term grows at least as \(2^r\). Every \(H\le X\)
in the eight-source grammar forces \(q,r=O(\log X)\), leaving
only \(O((\log X)^2)\) possible outputs. Hence even the complete
four-prime local property does not provide the absolute target.
A growing source bank, deeper branching, or a target-decreasing
physical recursion remains necessary in this approach.

Run python replay_terminals.py for the finite Jacobian certificates,
literal sources and the exact-target witness. Run python
compile_terminals.py for constructive local lifting and CRT; it uses
modular matrix exponentiation to check sample precision demands
without constructing the enormous resulting tournaments.

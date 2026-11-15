# Exact slide holonomy: a finite graph of integer count changes

A marked-core slide is a useful instance of a closed operation with a
nonzero readout: the surrounding context cannot detect that the local
interface moved, but the Hamilton total detects a definite translation.
This note packages **all tile lengths** at a fixed physical core into a
finite weighted graph. It complements the [matched slide rule](../chanpak_window_descent_20261115/MATCHED_SLIDE.md)
and the all-length [five-core overlap controller](../chanpak_window_descent_20261115/FIVE_CORE_WIDTH.md).
The controller below is for a core returning to the **same marked core**.
A changing-core relation would need a larger boundary state.

## The closed-walk compiler

Fix a nonempty core \(C\), with a literal apex mark set. Let \(h\) be the
ordinary Hamilton count and \(Q\) its two-path-cover count (empty parts
allowed). Let \(s_0,r_0\) be respectively the sums of weighted starts and
ends at **marked** vertices. The weight of an open path is the Hamilton
count of its complement. Let \(t_0=Q-h-r_0\), the unmarked end total.
Use \(z_i,w_i\in\{0,1\}\) for successive marks in two words \(z,w\) of
the same length, read from left to right. At a graph vertex keep
\((s,r,\alpha,\beta)\), with \(s,r\) integers in \([-h,Q]\) and
\(\alpha,\beta\) the preceding two bits. Each pair of bits defines the
possible edge
\[
 s'=2s+h z_i-Qw_i,\qquad
 r'=\frac{r+Qz_i-hw_i}{2},\qquad
 (\alpha',\beta')=(z_i,w_i).
\]
Require an integer \(r'\) and keep the edge only inside the stated bounds.
Give it the integer weight
\[
\omega=(1-z_i)s-w_i(Q-h-r')+
 h\{\alpha(1-z_i)-\beta(1-w_i)\}.
\]
Consider the four base states \((s_0,r_0,\alpha,\beta)\). A
nonempty closed walk based there spells exactly a pair of returning
words with the same unsigned response for \(Cz,wC\), with the
memory bits set to the **last** pair of bits of those words.
Close the literal tournament with guards \(p=\operatorname{last}(w)\),
\(q=\operatorname{first}(z)\). Then the sum of the edge weights is
exactly
\[
 H(pCzq+\mathrm{apex})-H(pwCq+\mathrm{apex}).
\]
There is no reduction modulo a selected integer in this equality.

## Proof of the accounting

The two equations for a match are
\(2^d s_0+h z=Qw+s_0\) and
\(Q\operatorname{rev}(z)+r_0=2^d r_0+h\operatorname{rev}(w)\).
Iterating the equations for \(s',r'\) proves the closed-walk equivalence.
For the bounds, work backward from \(s_d=s_0\) by
\(s_i=(s_{i+1}-h z_{i+1}+Qw_{i+1})/2\), and forward from \(r_0\)
by the displayed update. The interval \([-h,Q]\) is invariant
under each relevant half-step. For an integral matching word,
\(r_i\) is integral as well, using its backward expansion from \(r_d\).

For the count, first append \(z\) to \(C\) and, separately, prepend \(w\)
to \(C\). The first operation contributes the marked start total and
signed start at every 0 in \(z\). The second, by applying converse
and complementing marks, contributes the unmarked end total and
signed end at every 1 in \(w\). In the unsigned part, the totals at
step \(i\) can be written
\[
 S_{i-1}=s_{i-1}
    +Q\sum_{j<i}w_j2^{i-j-1},
 \qquad
 T_i=Q-h-r_i+
    Q\sum_{j>i}(1-z_j)2^{j-i-1}.
\]
These follow by expanding the start recurrence from the left and
using the end match to run the end recurrence back from the right.
In the difference of the accumulated counts, every double term
\(Qw_j(1-z_i)2^{i-j-1}\), \(j<i\), cancels. Only the first two terms
of \(\omega\) remain. The signed contribution without guards uses
the core's signed start at the first \(z\) and the core's signed end
at the last \(w\). The right guard replaces the first boundary
term by \(h z_d(1-z_1)\); the left guard replaces the last by
\(-h w_d(1-w_1)\). In other words all signed terms can be read
cyclically, giving the memory term in \(\omega\). This proves
the integer, guarded count identity.

## The arithmetic of returns

For one base state, restrict to vertices that can be reached from it
and can return to it. Choose one walk from the base to each such
vertex \(v\), of weight \(P(v)\). Then the gcd of *all* genuine slide
displacements through this base equals
\[
 g=\gcd_{u\to v} \bigl(P(u)+\omega(u,v)-P(v)\bigr).
\]
The sums along closed walks telescope, giving one direction.
For the other, append any walk \(v\to\text{base}\) both to the
chosen walk to \(v\) and to the chosen walk to \(u\) followed by
\(u\to v\). Their difference is that edge's reduced weight.
This is a finite certificate, even for arbitrarily long tape words.

* If \(g=0\), the weight is a vertex-potential difference on this
  component: every matched return is count-neutral.
* If both positive and negative closed walks exist, **every** multiple
  of \(g\) is physically realizable as a single guarded return.
  Choose finitely many walks with gcd \(g\), including both signs.
  Integer Bézout coefficients can all be made nonnegative by adding
  zero relations between an opposite-sign pair, then concatenate.
* If only one sign is possible and \(g>0\), the attainable nonzero
  magnitudes contain all sufficiently large multiples of \(g\).
  The additive semigroup of closed-walk weights contains a finite
  subset with that same gcd; normalize by \(g\) and take a least
  representative in each class modulo one positive generator.

In particular, obtaining a unit half-step is not a conjecture about
scanning to sufficient tape depth. In any component with both signs,
it is equivalent to \(g=2\). Concatenating to reach a prescribed
integer displacement can be expensive; the criterion does not bound
the absolute center of its resulting tournament.

## Physical contrast and replay

The [replay](replay_slide_cocycle.py) reconstructs the reachability
restriction and reduced-weight gcd at selected **physical** marked
profiles. Results include:

* At the earlier 13-core profile \((h,Q,s,r)=(13,104,18,18)\),
  a return component has gcd 10 and actual shifts of both signs.
  This fixes the bottleneck for this profile at every word length,
  not just at the old 12-vertex tile.
* At the new 45-core profile \((45,1440,517,818)\),
  the live signed component at memory \((1,1)\) has gcd 2.
  Two explicit 12-bit returns have displacements \(2,-16\).
  Both have left guard 1 and right guard 0. Concatenations give
  any even displacement; in particular the positive unit return
  itself requires no long Bézout combination.

This is a tournament version of an integer cycle invariant: local
changes of potential do not alter the weight of a closed walk.
The finite graph captures which arithmetic translations are real
at one guarded port. Locating such a port in a witness with a
prescribed absolute count remains a different problem.

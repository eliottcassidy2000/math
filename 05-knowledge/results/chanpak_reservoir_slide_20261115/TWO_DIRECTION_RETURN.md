# Two-ended fans: both arms can be ready, but stable menus never refill

A forward move \(F_q\) adds a vertex beating exactly \(q\); a reverse
move \(R_q\) adds a vertex losing exactly to \(q\). Put
\(\delta_F(q)=H(F_qT)-2H(T)\) and define \(\delta_R\) analogously.
A forward/reverse arm is **repeat-stable** if making its own move
preserves its defect, even if the defect is not a unit. This is
different from requiring that the _other_ arm survive the move
with the same defect. Reversing all arcs interchanges \(F,R\) and
preserves \(H\).

## Exact transport across the other direction

By the [fan and cross-port laws](../chanpak_leaf_fans_20261115/CROSS_PORTS.md),
\(q\) is repeat-stable forward exactly when every Hamilton path of
\(T-q\) ends at an in-neighbor of \(q\). Repeat-stable reverse means
all paths of \(T-p\) start at an outneighbor of \(p\).
The set of possible starts of a nonempty tournament equals its initial
strong component. Indeed strong components are linearly ordered and
each nonsingleton strong tournament has a Hamilton cycle; break it at
any desired start and concatenate the component paths.

Assume \(T\) is strong. Given \(F_qT\), a distinct old reverse arm
\(p\ne q\) is repeat-stable exactly when it was repeat-stable
beforehand and \(q\) lies outside the initial strong component of
\(T-p\). To see this, after deleting \(p\), the new vertex reaches
only the suffix starting at \(q\)'s component. If that is not the initial
component, the initial component is unchanged. Otherwise \(q\)
reaches every old vertex of \(T-p\); since there are at least two
old vertices, it also reaches the new vertex via one of the others.
Thus the whole deleted graph becomes strong. A strong original
\(T\) has some in-neighbor of \(p\), so in this latter case reverse
stability fails. If \(p=q\), deleting \(q\) instead makes the leaf
a sink and preserves the old reverse criterion. The newborn could
be reverse-stable only if every Hamilton path of \(T\) started at
\(q\), impossible for strong \(T\) with at least three vertices.
The dual statement transports forward readiness across \(R_p\).

When the performed forward move is itself repeat-stable, it
does not create any repeat-stable forward arm either. In the
known forward transport law, an old arm \(p\ne q\) can survive
exactly when \(N_T^+(q)=\{p\}\), in which case it was already
repeat-stable: \(q\) was the sink of \(T-p\) and \(q\to p\).
The newborn cannot be forward-ready because \(q\) belongs to the
terminal component. Apply the dual argument to a repeat-stable reverse
move. Thus **neither direction's menu ever gains a new vertex in a
repeat-stable two-ended trajectory on a strong seed.**

Each direction can change its used arm at most once, even interleaved
with moves of the other direction. For the forward assertion, at
the last forward move on \(q\) before a switch to \(p\), the
transport law and the noncreation of arms force \(N^+(q)=\{p\}\).
The tournament is strong, hence \(p\) already has an outneighbor;
the new vertex is an additional one. Its outdegree cannot decrease
under adding vertices, so the first forward move at \(p\) eliminates
every other forward arm. Reverse is dual.

Therefore a finite seed bank, even permitting all defects, yields
only \(O((\log X)^4)\) distinct counts through \(X\) in this grammar.
Old roots always belong to the seed and every fixed-root
forward/reverse pair commutes physically: the reverse newcomer
beats the forward newcomer regardless of their construction order.
One can group the two directions for counting outcomes, with at
most two fixed-root runs in each. All midpoints from a strong seed
stay strong, and [Busch's order bound](https://doi.org/10.37236/1141)
makes the combined length \(O(\log X)\).
For a reducible seed, forward repeat-stable moves stay in the
terminal component, and reverse moves stay in the initial component.
Inside each nonsingleton boundary component they still have at most
one switch; at a singleton boundary they are neutral and keep that
boundary a singleton. The two boundaries remain distinct and other
components are unchanged. A singleton seed starts with a neutral
two-vertex tournament. By component multiplicativity the same
upper bound applies to arbitrary finite seed banks.
This does not forbid unstable unit moves, non-leaf rewiring, or using
an infinite seed bank selected from the target.

## Two physical examples, including distinct armed vertices

**Same site.** The strong four-vertex tournament with outgoing
bitmasks \((12,1,2,6)\) has count 5. Site 0 is both repeat-stable
forward \(+1\) and repeat-stable reverse \(-1\), giving children
11 and 9. It also has two distinct forward repeat-stable arms:
\(+1\) at 0 and \(-1\) at 1. This is the smallest strong order
at which either phenomenon occurs: the directed triangle has no
positive defect. A pure forward history starts at label 2 and
appends labels \(0\to2,\ 3\to2,\ 1\to0\), where the right-hand
label denotes each newborn's unique beaten old vertex.

There is a general same-site formula. Treat \(q\) as apex over
\(T-q\), with its outneighbors marked. Let \(S,d\) be the
total and signed marked starts, \(t,f\) the total and signed
unmarked ends, and \(Q\) the ordered two-path count of the spine
(allowing either slot to be empty). Descending into this interface
also keeps the ordinary spine count \(h=H(T-q)\).
With \(a\) forward leaves and \(b\) reverse leaves at \(q\), the
physical marked spine is an all-one transitive prefix of length \(b\),
then the old spine, then an all-zero suffix of length \(a\).
The [ordinal cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md)
therefore gives
\[
H(a,b)=H+(2^b-1)t+f[b>0]
 +(2^a-1)\bigl(Q(2^b-1)+S\bigr)+d[a>0].
\]
Here the all-one prefix has start weight \(2^b-1\) and
signed weight \(1[b>0]\); the all-zero suffix is dual.
Repeat-stability of both original site arms means precisely
\(d=f=0\), by comparing the first and second defect on each axis.
Then the forward and reverse defects at \((a,b)\) are
\[
\delta_F(0,0)+(Q-t)(2^b-1),
\qquad
\delta_R(0,0)+(Q-S)(2^a-1).
\]
The two coefficients \(Q-S,Q-t\) are positive because the
two-path totals include the \(h>0\) covers omitted from each
weighted endpoint sum. They are even under joint repeat-stability,
since \(S,t\) are even when \(d=f=0\), and the two nonidentical
slots can always be exchanged freely on a nonempty spine.
Thus the first opposite move already raises the affected unit
defect by at least two. For the literal count-5 seed, \(Q=8\),
\(S=6\), \(t=4\); making the forward move changes the reverse
\(-1\) to \(+1\), and making the reverse changes the forward \(+1\)
to \(+5\).

**Distinct sites.** The strong count-59 tournament with outgoing
rows \((114,16,3,7,12,30,62)\) has exactly forward repeat-stable
\(-1\) at label 4 and reverse repeat-stable \(+1\) at label 0
among its ready unit arms. It too has a pure forward one-parent
history: start at 4, then append \(0\to4,\ 6\to4,\ 5\to4,\
3\to0,\ 2\to0,\ 1\to4\). With forward/reverse repetitions
\(a,b\) at the two indicated sites its exact count is
\[
H(a,b)=59+58(2^a-1)+60(2^b-1)+56(2^a-1)(2^b-1).
\]
To prove this identity, use the [two-module slot proof](../chanpak_leaf_fans_20261115/OPPOSITE_FANS.md):
each fan's transitive vertices can appear in at most two runs,
with linear coefficients in its own \(2^a\) or \(2^b\) on the
nonempty range. On the axes and the interior the exact polynomial
is fixed by the \(3\times 3\) literal grid
\((59,119,239),(117,233,465),(233,461,917)\).
Both markers remain repeat-stable regardless of mix: at the
origin they satisfy the two cross-direction transport tests;
along an axis each test keeps the old relevant initial/terminal
component, and in the interior deleting the reverse (forward)
root leaves reverse (forward) fan vertices as a source (sink)
module excluding the other root. Nevertheless their defects
are now
\[
-1-4(2^b-1),\qquad +1-2(2^a-1).
\]
In particular both choices exist at the root (117 and 119), but
neither child retains both signs. This is the distinction between
cross-survival of _readiness_ and cross-survival of the _unit
digit_. It holds even in a literal one-parent tree seed.

The accompanying C++ enumerates every strong class at order 7
and covers every strong class at order 8 by extension of all
order-7 classes. It finds no depth-two **raw** unit tree even
when both orientations are allowed and no repeat-stability is
required. Only the stated finite orders are tested; there are
raw one-sign surviving children, and hypotheses about higher-order
unstable arms remain open.

## All-order phase law for a compatible pair of fixed roots

Suppose a strong seed has a forward repeat-stable root \(q\) and
reverse repeat-stable root \(p\), neither necessarily unit.
If \(p=q\) they are automatically compatible. If \(p\ne q\),
call them compatible when \(q\) is outside the initial component
of \(T-p\) and \(p\) is outside the terminal component of \(T-q\).
These conditions are necessary and sufficient for *both*
readiness flags to survive arbitrary usage of just these two
moves. Necessity is the one-step transport test. For
sufficiency, while a root's own leaves are absent, repetitions of
the opposite move preserve the relevant component exactly by that
test. Once the root's own leaves are present, deleting it sees
those leaves as a source module (reverse) or a sink module
(forward); the other old root is outside it. The same-site case
follows directly from the \(p=q\) part of transport.

For a compatible pair put \(N=H(T)\),
\(U=N+\delta_F(q)\), and \(V=N+\delta_R(p)\).
For some integer \(C\) the **entire** opposite-direction grid is
\[
H(a,b)=N+U(2^a-1)+V(2^b-1)+C(2^a-1)(2^b-1).
\]
Indeed, for any fixed \(b\), forward stability makes \(H(a,b)\)
affine in \(2^a\) including \(a=0\) (iterate the unchanging own
defect). Reverse stability similarly gives affinity in \(2^b\)
including \(b=0\). The four values at \(a,b\in\{0,1\}\)
determine the formula. Consequently the two defect drifts are
\[
  \delta_F(a,b)=\delta_F(0,0)+(C-V)(2^b-1),\qquad
  \delta_R(a,b)=\delta_R(0,0)+(C-U)(2^a-1).
\]
If they start at \((-1,+1)\), then \(U<V\). For positive \(a,b\)
the initial opposite pair could persist only with \(C=U=V\);
to become the swapped pair \((+1,-1)\) would require both
\(C>V\) and \(C<U\). Thus **no interior point can ever rearm these
two fixed roots with opposite unit defects, even if the signs are
allowed to swap.** The case \((+1,-1)\) is the converse.
On the pure forward and reverse axes the initial opposite menu
persists exactly when \(C=U\) or \(C=V\), respectively, so it can
persist along at most one of the two axes. This phase law is
stronger than a fixed-sign commutation obstruction and does not
exclude moving to different non-ready roots.

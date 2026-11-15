# Physical routes: moving-modulus calibration and a single-site obstruction

A fixed pair of witnesses can often be compared much more strongly than
their shared residue: one can keep an actual marked square and tune the
two edge responses into an arithmetic relation. This yields infinite
exact physical routes for a specified family of targets. It does not
yet join the two arbitrary witnesses produced by the older bracketing
theorem.

## A marked four-vertex route at a varying modulus

Use the 11-bit center of the [interleaved
construction](../chanpak_mixed_defects_20261115/INTERLEAVED.md):
the base word is 1158 and the two independent xor masks are 32 and
220, denoted A and B. The center counts are 721, 723, 725, 727 in
setting order ∅, A, B, AB. Let a be the integer of the left exterior
and b the reverse-complement integer of the right exterior. These
can be any nonnegative integers. In all four settings the end
bits of the center stay fixed, both seam parity guards vanish,
and the apex marks the first and last spine vertex in opposite
directions, making the tournaments strong. The mark steps and base
are exactly

    A(a,b) = 2 + 32(b−a)
    B(a,b) = 4 + 56a − 44b
    L(a,b) = 721 + 1270a + 1158b + 2048ab + p(a) + p(b),

where p(w)=H(w)−1 on a binary word. The mixed term is zero for
every exterior. The [prior seam formula](../chanpak_mixed_defects_20261115/INTERLEAVED.md)
proves the displayed identities, including the quadratic base term.

Instead of trying to keep both numerical steps fixed, impose the
homogeneous relation B=2A. It gives exactly a=9t,b=10t, t≥0.
Let M=2+32t. The actual route is

    ∅  --A-->  A  --A⊕B-->  B  --A-->  AB,
    L          L+M           L+2M        L+3M.

Each step changes a fixed marked set of one or six apex arcs;
no arbitrary graph replacement is used. Counts stay in one residue
modulo M, and the signed crossings hit both internal targets. At
t=1, the four physical counts are 208063, 208097, 208131, 208165.

There is a second lane: B=−2A gives a=5t+1,b=2t, t≥0. Now
M=30+96t and the sorted route is A, ∅, AB, B, again at steps M.
These are the only four-distinct equally spaced calibrations of
this fixed additive block under the stated exteriors. A four-term
additive rectangle is equally spaced iff the absolute steps have
ratio 1:2 or 2:1. The other two algebraic cases A=±2B are impossible
by reducing modulo 4: A≡2, B≡0.

**Limits.** This is a coupled family of physical witnesses, not a
union of independent approximate hits. But M grows on the order
of the square root of the count and both lanes have v₂(M)=1.
They do not thread the older high-dyadic congruence from arbitrary
prescribed N. Moreover the middle move flips both controls.
If only the two individual switches are allowed, a strictly increasing
three-edge path around a Boolean square with all steps M forces mixed
difference ±2M: the repeated control's first and last changes must
both be +M, and are opposite-oriented parallel edges. Exterior
padding modifies the linear coefficients, not the mixed coefficient,
so such an elementary route cannot have growing M on a fixed square.
The [earlier nonlinear escape](../chanpak_mixed_defects_20261115/NONLINEAR_ESCAPE.md)
does give the elementary A,B,A route 213→215→217→219 at M=2.

## All-order small-component theorem for single-site routes

This obstruction is about routes in the one-apex transitive-spine model
which flip only one apex mark at a time, not about flipping arbitrary
tournament arcs. Allow neutral leading zeros or trailing ones to be
added or deleted at any time. If these are stripped off, the effective
word starts in 1, ends in 0, and has length m. Let each mark flip have
absolute count change at most 2g, even allowing changes of zero.

**Theorem.** At odd m=2k+1 with g<2^(k−2), and at even m=2k with
g<(2^(k−2)−1)/2, every component has at most two effective words.
In particular the small-step graph at any large one-apex count does
not provide long routes when g=o(sqrt(H)). The count can neither
run across a long residue interval nor move to a different effective
word length using such flips.

**Proof.** The [single-site formula](../chanpak_mixed_defects_20261115/SINGLE_SITE.md)
says the gap at pivot U|bit|V is
2(ceil(v/2)−ceil(u/2)), with u=bin(U), v=revcomp(V).
In particular |u−v|≤2g+1 on every allowed edge.

A pivot outside the strict interior of the effective word has one
address zero and the other at least 2^(m−2): at or before the first
1 there is no preceding 1, while the last effective 0 remains at
distance at least m−2; the opposite side is symmetric. Thus no
small flip can change m or touch the neutral padding nontrivially.
Inside the word, at index p starting from 0, u and v have exact
bit lengths p and m−1−p, since the first effective 1 and last
effective 0 are present. At odd m, any non-middle pivot thus gives
|u−v|≥2^k−(2^(k−1)−1)=2^(k−1)+1. The middle pivot is the
only possible control and its outer data do not change across it.

At even m only the middle two pivots can be allowed; any other
pivot gives |u−v|≥2^k−(2^(k−2)−1). Let x,y be the fixed
left ordinary and right reverse-complement addresses outside those
two sites; each is a full (k−1)-bit integer, at least 2^(k−2).
Eligibility of the left middle pivot, regardless of right middle bit,
implies y≤ceil(x/2)+g. Eligibility of the right, regardless of left
middle bit, implies x≤ceil(y/2)+g. If a component ever uses
both, adding gives x+y≤4g+2, contradicting the stated bound.
So at most one pivot can be used on a component.

Finally 2^(m−2)≤H≤2^m for a nonempty effective word: the first
and last sites alone have odd-cycle coefficient 2^(m−2) in
H, while H=1+p(w)≤1+w≤2^m. Thus the thresholds are on the order
of sqrt(H). In particular the polynomial precision 2^k≈N^(1/160)
does not by itself make a one-mark walk rich enough at height N.

This describes a difference between *generating* witnesses and
*transporting* them. To join arbitrary bracketing windows at small
modulus, at least one move must change multiple marks or leave the
one-apex architecture; neutral padding and repeated one-site edits
are insufficient.

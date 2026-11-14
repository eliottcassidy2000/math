# Bounded-feedback tournaments: logarithmic height and an explicit repeat spectrum

The earlier two-apex chart tracks 27 states for each fixed ordering of
the two apices. Two all-order results here hold for arbitrary k:
a logarithmic height bound for an absolute target, and a small
annihilator for repeat-count sequences. In the second result, the
powers of the characteristic roots come from gap occupancy, not
from an experimentally fitted recurrence.

## An absolute height bound by keeping the cycle's ports

Let T be a strong tournament of order n≥3. Suppose a set F of
k vertices can be removed to leave a transitive tournament.
In a Hamiltonian cycle of T, split the non-F vertices into maximal
stretches between vertices of F. There are at most k nonempty
stretches. Retain every member of F and the first and last vertex of
each nonempty stretch, counting a singleton once. Call the number
of other, optional vertices d, so d≥n−3k. The Hamiltonian cycle
exists by Camion's theorem (also from the pancyclicity argument
[here](CONFLICT_GRAPHS.md)).

A stretch is an increasing word in the transitive spine. Any
choice of optional vertices still gives a directed cycle: its
endpoints and the arcs to and from F are unchanged. Every
choice gives a distinct vertex set, and when d≥1 exactly half
of the 2^d choices have odd total size. The retained vertices
number at least three: for k=1 a cycle needs two spine vertices;
for k=2 it needs at least one; for k≥3 this is automatic.
By the [odd-cycle expansion](https://arxiv.org/abs/2307.05569),
the empty packing and these single-cycle packings imply
$$
 H(T)\ge 1+2^{n-3k}.
$$
When d=0, n≤3k and the inequality follows from H(T)≥3
for a nontrivial strong tournament. The proof also works with the
actual number of retained endpoints instead of the upper bound 2k.

Thus for target N>1 the two-apex strong search has the rigorous
cutoff n≤6+floor(log₂(N−1)); arbitrary k gives 3k instead of 6.
The exponential base cannot increase uniformly even for k=1:
reverse the single arc from the first to the last vertex of a
transitive tournament on n vertices. Every cycle uses that arc,
so no two odd cycles are disjoint. The 2^(n−3) choices of a
nonempty odd number of internal vertices give H=1+2^(n−2).

For arbitrary T, discard singleton SCCs and apply the estimate
in each other SCC, using its intersection with F. Path counts
multiply across the ordinal SCC order, so the retained order is
at most 3k+floor(log₂ N). Hence for fixed k, membership in the
whole bounded-feedback spectrum is decidable by searching at
most O_k(N^k) binary contact tapes times a polylogarithmic
factor in N: choose the ≤k named apices and their internal
orientation, then enumerate literal tapes up to this order.
The transfer below evaluates each tape using a constant number
of integer states for fixed k. The case N=1 has a transitive
witness separately.

A fixed template u_1^e_1 ... u_r^e_r of nonempty contact blocks
(with fixed extra positions allowed), restricted to its strong
realizations, consequently produces only O((log X)^r) different
counts ≤X. The same holds for any finite union, with the
largest r. This does not restrict arbitrary words in a fixed
alphabet to finitely many such templates.

## Gap transfers simultaneously triangularize

Fix an order of the k apices in a potential Hamiltonian path
and put K=k+1. A full path is uniquely an assignment of every
spine vertex to one of the K gaps between consecutive apices
(including the exterior gaps). The vertices in each gap appear
in their spine order. An empty internal gap is valid iff its
two consecutive apices form the requisite arc; exterior empty
gaps are always valid.

One gap has three indicator functions e,f,g on its partial
state: e is empty, f is nonempty, and g says the gap is nonempty
and its current last vertex can close the gap. A new spine
vertex has an opening bit a and a closing bit b for each gap;
a missing exterior guard is always 1. When assigning this
vertex to that gap, pull back the functions along the update,
counting an invalid opening as zero. The one-gap map is
$$
 R(e)=0,\qquad R(f)=ae+f,\qquad R(g)=b(ae+f).
$$
Assigning the vertex to exactly one gap means the global letter
operator is the sum of R on each tensor factor and identity on
the others. Its diagonal on a tensor word of length K in e,f,g
is its number j of f's. Every off-diagonal move strictly decreases
some factor of e<f<g, and none increases any factor. The initial
state has all gaps empty, and the final indicator is the product
of g plus the permitted empty indicator for each gap. Summing
over apex permutations gives the exact path count. The direct
sum over permutations has the same diagonal bounds.

## A sharper repeat annihilator

Repeat a fixed nonempty block u of length L, with arbitrary
fixed prefix and suffix. Let M_u be its transfer. It still
preserves componentwise monotonicity, and the diagonal entries
are j^L for j=0,...,K. For the scalar sequence c_t the shift
operator E obeys the uniform annihilator
$$
 \prod_{j=0}^{K}(E-j^L)^{K-j+1}\ c=0.
$$
These are upper bounds on the necessary powers, not claims
of minimality for every block.

To prove the powers, write M_u as diagonal plus strict triangular.
Resolvent expansion of (I-z M_u)^(-1) as a finite sum over strict
chains gives one factor (1-zj^L)^(-1) per visited state carrying
j f's. Along a componentwise decreasing chain each strict move
between states with the same j must lose at least one g:
if the g count stays fixed, no g changes, and any strict
f→e change would decrease j. At the first j-state there are
at most K−j g's, so a chain has at most K−j+1 states of this
same j. For a nonzero eigenvalue this bounds the poles, hence
the Jordan powers. The same argument at the zero eigenvalue
uses the resolvent (xI-M_u)^(-1), whose factor is 1/x.
Therefore the full matrix is annihilated by the stated polynomial;
prefix, suffix, and summing apex orders do not change it.

After the zero-eigenvalue transient (t≥K+1), the exact scalar is
a sum of polynomial multiples of (j^L)^t, j=1,...,K, where
the degree for j is at most K−j. In particular for k=2
the nonzero bases are 1, 2^L, 3^L with powers at most 3, 2, 1.
This shows the difference between a 27-state presentation and
the six modes counting multiplicity that can persist.

Two such repeat sequences at the same L and k either agree for
the whole tail or meet at no more than K(K+1)/2−1 tail indices.
For a nonzero difference, view the tail exponential polynomial
on the real line. Divide by the smallest exponential and
differentiate as many times as that polynomial's degree plus one.
This kills the smallest term, leaves the degrees of all others
unchanged, and loses at most that many distinct zeros by Rolle's
theorem. Induction gives the bound, with the one-term case a
polynomial. For two apices there are at most five tail meetings.
This is a restricted comparison in the spirit of the
[Skolem problem](https://arxiv.org/abs/2607.15510); it neither
allows general complex characteristic roots nor reduces the
unrestricted Chan–Pak question to a fixed recurrence.

The port bound and simultaneous triangularization came from
the background research pass. The individual Jordan bounds,
zero-count consequence and absolute-search application extend
them. See [replay](replay_apex_recurrences.py) for independent
small subset-DP checks and the annihilator controls.

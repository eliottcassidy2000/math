# Matching sections: exactly what small uniform contexts can see

**Proved.** Uniform tournament substitution has a sharp path-cover
observation threshold controlled by a bipartite matching, and a
same-order family of literal strong outputs recovers every observable
coordinate at that threshold. Polynomial recovery does not imply
tournament reconstruction or exact target coverage.

The [weak Lefschetz theorem](https://grothendiecksga.com/read/sga5/en/VII-7.html)
identifies low-degree cohomology under restriction to a smooth
hyperplane section and only gives an injection in the boundary degree.
The analogy here is restricted to the distinction between recovering
specified low-degree data and reconstructing the original object.
This is a combinatorial theorem, not a use of a topological theorem.
The repository's THM-587 is a different Lefschetz connection: the
fixed-point trace of converse, not restriction to a section.
THM-1975 is the correct antecedent for path-cover composition;
Session 329's triangle-cut formula is its two-slot case.

## 1. Matching capacity of every uniform module context

Fix a tournament U on n vertices, and adjoin a set V of b outside
vertices. For each v outside, require that v beats all of U or loses
to all of U. Split V into O (U beats v) and I (v beats U); the induced
outside tournament is otherwise arbitrary. Let ν be the maximum
matching among pairs o→i with o in O and i in I.

Write P_k(U) for the number of ordered covers by k nonempty directed
paths. For every such context,

    H(out)=sum_(1≤k≤n) c_k P_k(U),

where c_k is the number of context words using each outside vertex
once and k copies of one abstract symbol T, with no consecutive T
and with every consecutive direction lawful. A Hamiltonian path
uniquely collapses its k maximal U-runs to these visits. Conversely,
fill the visits with the k paths of an ordered cover: uniformity is
precisely what makes all joins legal. Thus all c_k are nonnegative
integers depending only on the context.

**Exact capacity:** c_k>0 if and only if 1≤k≤ν+1. In any gap between
two visits, the first outside vertex lies in O and the last in I.
Some consecutive pair in this gap crosses from O to I. Distinct
gaps use pairwise disjoint crossing edges and produce a matching, so
k−1≤ν. Conversely take k−1 directed edges forming a matching, and
assemble

    T o_1 i_1 T o_2 i_2 T ... o_(k−1) i_(k−1) T.

Insert every other outside vertex, one at a time. The usual
tournament insertion argument works even with repeated T symbols:
either the new vertex beats the first symbol, loses to the last,
or there is a consecutive pair bracketing a transition from a
symbol beating it to one losing to it. Each outside vertex's
orientation to T is defined, and insertion does not change its
visit count. For k=1 this starts with just T, so no special
connectivity or positivity premise is needed. For a real path k≤n
as well. P_k(U)>0 for 1≤k≤n by cutting any Hamiltonian path.

In particular b outside vertices cannot see above P_(floor(b/2)+1).
At the tight budget b=2(r−1), a nonzero top coefficient requires
|O|=|I|=r−1, then

    c_r=(r−1)! permanent( [o→i]_(o in O,i in I) ).

Every gap is exactly two vertices and every exterior vertex is in
a gap. Choose a perfect matching and order its pairs. Thus the
largest top coefficient is ((r−1)!)², attained by O→I complete.
Even when the number of vertices fits, the matching can fail.
Nonuniform rooted seams do not have this bound; their contact
data must be retained separately.

## 2. Optimal same-order experiments

Let r≥3 and n≥r, and fix b=2r−2. Use cyclic substitution
U→W→X→U with W and X nonempty. All outputs are strong, even
when any individual module is not. For module visit counts
i,j,k, the skeleton multiplicity is 3 if i=j=k, 1 if the
maximum minus the minimum is 1, and 0 otherwise. This follows
because the word must follow periodic cyclic order: an entire
number of cycles permits three start symbols, and one or two
extra symbols specify a unique starting symbol.

For m=2,..,r−1, set W=TT_m and X=TT_(b−m), where TT_s is
transitive of order s. The response has no terms past P_(m+1)
and its coefficient of P_(m+1) is

    m! (P_m(X)+P_(m+1)(X)) > 0.

Indeed the W count must be exactly m and the X count is either
m or m+1; each skeleton occurs once. Here P_j(TT_s)=j! times
the Stirling number of partitions of s objects into j parts.

Take two further probes with W singleton and X of order b−1,
one X transitive, the other C3→TT_(b−4) (empty TT_0 allowed).
Writing h=H(X), q=Q_X(2), their responses are

    (q+h)P_1(U) + (q−h)P_2(U).

Their q/h slopes are respectively 2^(b−1) and 2^(b−2), because
the triangle has (h,q)=(3,12) and ordinal singleton padding
doubles q without changing h. Hence these two rows have a
nonzero determinant and recover P1,P2; the rest recover P3,..,Pr
one at a time by subtraction and exact division. Every reading
uses the same source order and the same auxiliary budget. This
is minimal in exterior vertices for accessing the formal P_r
coordinate; r independent readings are necessary if those
coordinates are treated as free linear unknowns.

**Small case r=2.** To get both coordinates from strong
outputs, b=3 is enough. Take the four-vertex strong quotient
with outgoing masks (2,4,9,3), replacing vertex 0 in the first
experiment and vertex 1 in the second by U. Their readings
are 5P1+2P2 and 5P1+3P2, hence the difference is P2. With
b=2 there is only one cyclic quotient, giving 3P1+P2. If a
nonstrong output is allowed, a transitive quotient at the same
size reads P1 directly and b=2 suffices. This is an information
bound on the linear interface, not a claim that every small
physical family has two independent coordinates.

## 3. Full polynomial, exact cutoff, and non-reconstruction

The optional-slot polynomial is

    Q_U(z)=sum_(i≥1) binom(z,i) P_i(U).

The inherited disjoint-odd-cycle formula makes this a monic
polynomial of degree n, with only powers congruent to n mod 2
and no constant when n>0. Thus it has r0=floor((n−1)/2) lower
coefficients. The values Q(1),...,Q(r0) determine it: subtract
z^n and divide by z (odd n) or z² (even n), then interpolate the
result as a polynomial in z² of degree at most r0−1. These
evaluations are successively invertible triangular combinations
of P1,..,P_(r0). Hence for n≥7, r0 all-strong experiments at
exactly n+2r0−2 vertices recover the full polynomial. Order is
part of the metadata.

For fixed order, equality of the first r path-cover counts is
**exactly** indistinguishability under all uniform module contexts
of at most 2r−2 outside vertices if r≥3, even if tests are
restricted to strong outputs. It implies equality of all such
responses by the matching upper bound, and the triangle sections
are witnesses to the converse. At r=2 use the small-case boundary
above.

Two complementary physical controls check the two types of
limitation:

- Session 326's nonisomorphic strong order-eight sources share
  (H,Q2)=(345,5640) but have Q3=44505 and 45225. Thus no
  uniform context below four outside vertices can distinguish
  them. With two two-vertex transitive reservoirs the top
  coefficient is 4, and their output counts differ by
  4·720=2880. This is a sharp *physical* activation of the
  third layer, not only a coefficient count.
- Session 328's ten-vertex outputs after its Fermat triangle
  have identical entire slot polynomials but are nonisomorphic
  even up to converse. Thus even arbitrary uniform module
  tests cannot distinguish them. A purposefully nonuniform
  observation (e.g. rooted slot data) is another medium.

**Arithmetic jump at the frontier.** If two sources agree in the first
r−1 coordinates, substituting them separately into the same context at
budget 2r−2 gives the exact difference

    (r−1)! permanent(M) ΔP_r,

or zero if the exterior is unbalanced. Since P_r is r! times an
unordered cover count, every difference is a multiple of r!·(r−1)!,
and in a fixed source pair it is a multiple of (r−1)! ΔP_r. The
universal factorial obstruction persists for **every** exterior
size: freely permute the k−1 internal exterior excursions in a
context word with k visits. The prefix and suffix stay put and all
joins remain legal by uniformity; the disjoint nonempty vertex sets
make the action free. Thus c_k is divisible by (k−1)!. Every term
with k≥r then gives a multiple of r!·(r−1)!, even if further slots
have opened. This is an exact absolute-target obstruction: making a
layer observable does not make the repair arbitrarily fine. In the cited eight-vertex pair
and r=3, the only nonzero magnitudes at budget four are 1440 and
2880. Both occur on strong outputs: balanced complete crossing gives
the latter, and reversing one crossing leaves one perfect matching
and gives the former. The two checks use the same number of vertices
and the identical low source profile.

**A binary interval of normalized jumps at one exact budget.**
Put d=r−1. At budget 2d, draw the crossing matrix M with ones on
and above its diagonal, and selected entries (j,0) below it for
j=1,...,d−1. All other entries are zero. This still gives a strongly
connected quotient: the diagonal gives T→o_j→i_j→T for every j.
A nonidentity permutation of the matching matrix decomposes into
fixed points and precisely one nontrivial cycle: every nontrivial
cycle must use a below-diagonal arc, and all those arcs end at 0.
For the chosen return j→0, the remainder of the cycle is an arbitrary
increasing path from 0 to j, with 2^(j−1) choices. Consequently

    permanent(M) = 1 + sum_(j selected) 2^(j−1).

At exactly the same order, the normalized top coefficient c_r/d!
therefore realizes every integer 1,...,2^(d−1), with an explicit
inverse by binary digits. The proof uses the same increasing
path/one back arc mechanism as the near-transitive Fermat block.
For a fixed pair of source profiles whose first difference is
P_r, this controls their output difference through that entire
normalized interval at the very first size where the difference
can be seen. It does not keep the baseline output fixed when the
context changes.

Nothing here makes a prescribed target recover a physically
available source, since a linear inverse need not be realized
by a tournament. The result certifies what a context can
observe and precisely when more information can appear.

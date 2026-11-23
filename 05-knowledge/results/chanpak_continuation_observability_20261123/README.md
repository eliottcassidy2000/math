# What is visible depends on which continuation you allow

The two Berggren height-nine cores are OOOI and OIOI, or canonical
five-vertex masks 8 and 10, in the [preceding packet](../chanpak_berggren_divisor_twins_20261123/README.md).
The continuation problem has different exact answers depending on
how the outside sees the core.

- **Homogeneous contacts:** even all ordered-slot counts agree, so
  no fixed homogeneous shell distinguishes the cores. Conversely,
  cyclic shells with increasing numbers of gaps recover the whole
  ordered-slot vector; it is the exact equivalence invariant for
  homogeneous continuations, not just a sufficient statistic.
- **One individualized apex, then arbitrary ordinal flank words:**
  [five directional reads](TWIN-ATLAS.md) exactly characterize all
  responses. They give 60 profiles on the 64 possible rooted states.
- **One further individualized apex, and one marked edge reversal:**
  a single pair of endpoint heights recovers all 64 states. Two of
  the available flips even work when read from the wrong end.
- **All labelled external continuations:** full adjacency can be
  recovered, even using gap-only data for a marked arc; see
  [arbitrary-interface tomography](INTERFACE.md). This is distinct
  from the smaller sufficient statistics above.

## The shortest decoder in the one-extra-apex probe family

Fix the core's canonical labelling \(0,\ldots,4\). The unknown is
\((A,S)\), where \(A\) is one of the two cores and \(S\subseteq
\{0,\ldots,4\}\) is the set beaten by the designated first apex
\(v=5\). There are 64 such states. Adjoin a second apex \(w=6\)
which beats **only vertex 1** and loses to \(v\) and the other
core vertices. Measure the height. Then reverse \(vw\) and
measure again. In bit masks for the out-neighbors of \(w\), this
means 2 and 34. The resulting ordered pair is injective on all 64
states. In fact the unordered pair, with multiplicity, is also
injective: one can start at either end of this particular flip.

The exact replay scans all 64 individual contact masks for \(w\)
and all \(\binom{64}{2}\) fixed pairs of masks. No individual mask
works (at best 37 distinct heights). There are 74 successful
two-mask choices; exactly three are single-arc reversals.
Thus two outputs are optimal **in this specified seven-vertex
probe family**, without claiming minimality for larger shells or
different observables.

| Two masks on \(w\) | Distinct ordered pairs | Distinct unordered pairs |
|---|---:|---:|
| 2 and 34 | 64 | 64 |
| 10 and 42 | 64 | 62 |
| 23 and 55 | 64 | 64 |

The middle row isolates the danger in reading from the wrong end:
\((A,S)=(8,16)\) gives \((123,101)\), but \((10,26)\) gives
\((101,123)\). They are different cores, not merely two labels
for the same track. For masks 2 and 34, direction can be
discarded, but the difference alone cannot: it takes only 28
values, and 12 occur on **both** cores. Thus the needed data in
this experiment are the two endpoints, not simply a larger
signed derivative. The [replay](short-twin-decoder.py) reports
both the cross-core alias and the complete minimum test.

The [cut equation below](#exact-law-for-every-contact-of-a-single-apex-on-a-general-core)
also identifies what the two probes sample. For the already grafted
six-vertex state, their outputs are \(h+a_1\) and
\(h+a_1+a_v-W_{1v}\). The second probe uses an interaction between
the designated old apex and the one targeted core vertex; it
does not merely measure the old height again.

## Exact law for every contact of a single apex on a general core

Write \(h=H(T)\) and \(q(S)=H(T+x_S)-h\). In the
[inherited cut gauge](../chanpak_mirimanoff_unit_lifts_20261122/CUT-GAUGE.md),
the symmetric edge weight \(W_{ij}\) counts the two orders of
covering the old core by two nonempty paths touching its marked
pair at the open ends. There is an exact directed cut expression
\(q(S)=\sum_{i\to j,\ i\in S,\ j\notin S}W_{ij}\).
Set \(a_i=q(\{i\})\) and \(d_i=\sum_{i\to j}W_{ij}-\sum_{j\to i}W_{ji}\).
Equivalently
\[
 q(S)=\sum_{i\in S}a_i-\sum_{\{i,j\}\subseteq S}W_{ij}.
\]
For an untagged \(i\notin S\), changing its contact to tagged
has height change \(a_i-\sum_{j\in S}W_{ij}\); a tagged
vertex's reverse step has the negative value, evaluated on
the set without that vertex. This is a necessary and
sufficient decision rule for any specified apex contact flip.

The literally opposite view \(S^c\) separates the two registers:
\[
 q(S)+q(S^c)=\sum_{i\in S,j\notin S}W_{ij},\qquad
 q(S)-q(S^c)=\sum_{i\in S}d_i .
\]
The first is an undirected crossing weight, the second a signed
divergence. A target gap tested on both complementary contacts
therefore has one symmetric quadratic test and one signed linear
test. A weighted circulation is invisible to both. For instance,
the two orientations of the triangle cycle on **fixed labels**
are already indistinguishable by all single-apex contacts; their
pair weights agree and the weighted divergence is zero. The
single-apex observer is not universally a graph decoder.

### An optimal fixed linear interface

For \(n\ge4\) the real **linear span of the actual height-response
functions** has dimension \(n+\binom n2\). For \(n\ge5\) their
affine hull has that same dimension. A fixed optimal linear
probe list consists of the empty contact, all \(n\) singletons,
and all but one two-element contact. Indeed
\(W_{ij}=q(\{i\})+q(\{j\})-q(\{i,j\})\); the missing weight
is supplied by \(\sum_i a_i=\sum_{i<j}W_{ij}\), which is
\(q(V)=0\). This is an exact decoder over the integers. Fewer
fixed scalar probes cannot suffice for linear reconstruction
on actual data, or for affine reconstruction when \(n\ge5\).
This is not a lower bound for a nonlinear lookup on a
particular finite collection of cores.

**Proof of the span.** Use coordinates \((h,W,d)\), where
\(d\) is a zero-sum vector. Taking the converse preserves
\((h,W)\) and negates \(d\), so symmetric and signed components
can be studied separately. On a transitive core with labels in
order, \(W_{ij}/2=1\) at distances 1 and 2 and
\(2^{|i-j|-2}\) at larger distances: choose an even-length
path between its endpoints, with complement height 1. Its
divergence is nonzero, and differences of permuted copies
span the zero-sum vertex space (a transposition gives a
nonzero multiple of a coordinate difference).

For the symmetric part, permute the weighted transitive core
averaged with its converse. Applying \( (1-(ij))(1-(k\ell))\)
at four distinct vertices isolates the alternating
quadrilateral on \(ik,i\ell,jk,j\ell\), multiplied by the
weight contrast on that quadrilateral. With four consecutive
labels the contrast of \(W/2\) is \(1-2-1+1=-1\).
Hence the span of orbit differences contains every
alternating quadrilateral. Their span is the kernel of the
symmetric edge-array row-sum map: the orthogonal complement
consists exactly of arrays \(z_{ij}=u_i+u_j\). To see this,
solve for \(u\) on three anchors; quadrilateral equalities
first force the other-to-anchor entries and then all
other-to-other entries to have that form. The row sums of
the transitive profile are not constant (the endpoint and
its neighbor differ), so orbit differences project onto
the full zero-sum row space. This gives every edge-array
direction with total sum zero.

It remains to fill the two scalar directions \(h,\sum W\).
The odd-cycle packing formula implies
\(\sum_{i<j}W_{ij}=Q_2(T)-2h\). Indeed \(W/2\) sums the
even nonempty subset Hamilton counts weighted by their
complements, and even and odd subset contributions are equal:
each fixed odd-cycle packing contains at least one odd-sized
item (a cycle or an uncovered singleton), and sending that
item to the other slot changes the parity. Subtract the
empty subset, of weight \(h\).
For the transitive core and a cyclic triangle padded
ordinally to order \(n\), the pairs \((h,Q_2)\) are
\((1,2^n)\) and \((3,12\cdot 2^{n-3})\), linearly
independent. For the affine hull at \(n\ge5\), the third
pair from a height-nine twin padded to the same order is
\((9,84\cdot2^{n-5})\); the three are noncollinear. This
proves the claimed dimensions and optimality.

## What the exact census says about removing the labels

One-apex labelled indistinguishability is not the same as
nonisomorphic unlabelled ambiguity. The [finite census](apex-atlas-census.cpp)
counts one profile per unlabelled tournament through order eight:
\[
1,1,2,4,12,56,456,6880
\]
distinct profiles, with no nonisomorphic collisions in this range.

For each order, add a vertex with every contact pattern to each
canonically retained tournament of the previous order. This
exhausts isomorphism types by deletion of one vertex. Canonicalize
tournaments by permuting within equal outdegree classes, and cut
profiles by permuting within equal pairs of weighted incoming and
outgoing sums. Each canonical key minimizes over its whole
permitted class. The profile key uses \(H\), the two sums per
vertex, and the symmetric weight array. The endpoint DP counts
even-order paths from each marked start to each end, weighted by
the complementary Hamilton count, giving \(W/2\) by polarization.
No injectivity at larger orders is asserted.

~~~sh
python3 short-twin-decoder.py
g++ -O3 -std=c++17 -Wall -Wextra -Werror apex-atlas-census.cpp -o apex-atlas-census
./apex-atlas-census 8
~~~

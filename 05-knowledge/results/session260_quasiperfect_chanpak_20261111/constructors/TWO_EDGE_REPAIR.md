# Two-edge repair shells, a sharp all-apex gap, and an isolated four-edge repair

**Session260 second wave, 2026-11-11. PROVED: all-apex gap and exact rank-template formula. FINITE-EXACT: stated two/three/four-edge searches. General all-a repair remains OPEN.**

## 1. Actual source and requested target

Use the parent arithmetic lane's source

    T_a=C3(TT_a,TT_(a+1),1),     a≥2,
    n_a=2^a(2^(a+1)−3),
    H(T_a)=n_a+3,               target=n_a+1.

Vertices0,...,a−1 form A; vertices a,...,2a form B; r=2a+1 is the apex.
A→B→r→A, and each block is transitive in increasing vertex order.
This is the excess-two arithmetic family, not a quasiperfect assertion.
The parent/bridge all-a one-edge obstruction is inherited. This wave asks
what a jointly chosen pair or small pattern of arc flips can do that no
single flip can do. Flips never add/delete vertices.

## 2. Complete distance-two shells through a=6

`two_edge_repair.cpp` counts Hamilton paths by exact integer subset DP after
each unordered pair of distinct arcs is reversed. No graph or edge pair is
pruned, and isomorphic copies are retained as labelled presentations.

| a | Order | Source H | Target | Edge pairs tested | Hits | Smallest positive decrease |
|---|---:|---:|---:|---:|---:|---:|
|2|6|23|21|105|0|4|
|3|8|107|105|378|2|2|
|4|10|467|465|990|0|22|
|5|12|1955|1953|2145|0|46|
|6|14|8003|8001|4095|0|94|

The two a=3 repairs are

    {(0,7),(1,6)},       {(1,7),(2,4)}.

They show that the universal one-edge obstruction does not extend to two
edges verbatim. Conversely, these are not all-a repair templates. The
absence at a=4,5,6 is only the **complete finite two-edge shell** statement;
no all-a two-edge exclusion is claimed.

## 3. Stronger all-a result: arbitrary apex reorientation cannot repair

Here allow any number of apex-incident arc flips, but keep the root-deleted
host transitive. Let its length be N=2a+1. Encode the apex cut by a binary
word x_0...x_(N−1), where x_i=1 means r→i. Define

    w(1)=1,      w(d)=2^(d−2) for d≥2.

Then

    H(x)=1+2 sum_(i<j) x_i(1−x_j) w(j−i).                 (1)

**Proof.** Every directed cycle contains r, and all other vertices of a
cycle occur in increasing order. The odd-cycle expansion therefore has
only its constant and single-cycle terms. Fix its first and last old
vertices i,j. There must be an even number of old vertices, so the number
of internal choices is the sum of even binomial coefficients on j−i−1
positions. This is w(j−i), including the adjacent boundary. The apex edges
require x_i=1 and x_j=0. This proves (1).

**Sorting lemma.** Replacing adjacent01 by10 strictly increases H. The new
pair contributes2. Every earlier selected vertex now sees its unselected
endpoint one place farther right, and every later unselected vertex sees
its selected endpoint one place farther left. Since w is nondecreasing,
all other changes are nonnegative.

At fixed k ones, the unique maximum is therefore the sorted word
1^k0^(N−k), with count

    M(k)=(2^k−1)(2^(N−k)−1)+2,

for1≤k≤N−1; k=0,N both give1. Every unsorted word can be sorted by strict
increases. The last step must come from the unique word obtained by
swapping the boundary10 of the sorted word. Thus this boundary-swapped
word is the largest unsorted word at fixed k.

The largest M(k) occurs exactly at k=a,a+1. For these balanced values both
sides have size at least2. Evaluating the one boundary swap in (1) loses

    d_a=2^(a−1)+2^a−2=3·2^(a−1)−2.                       (2)

For completeness its loss is2 from the new inversion, plus
2(2^(k−2)−1) from earlier selected vertices and
2(2^(N−k−2)−1) from later unselected vertices. These are the only changed
summands. The next-largest sorted counts occur at k=a−1,a+2 and lose
3·2^(a−1), which is strictly larger than d_a. More unbalanced counts are
smaller still, because 2^k+2^(N−k) increases away from the center.

Consequently the **largest and second-largest distinct H-values among all
2^(2a+1) apex cuts** are exactly

    n_a+3,           n_a+3−d_a.                          (3)

In particular n_a+1 is absent for every a≥2, even if arbitrarily many apex
arcs are reoriented. This is a genuine all-a theorem, not an inference from
the two-edge census. A repair must alter at least one old host arc or leave
this source construction. `two_edge_repair.py` also checks the two top
levels of the complete cut spectra for a=2,...,8.

An O(N) evaluation of (1) is useful: start h=1, z=p=0; read each bit b;
if b=0 add z+p to h; then replace z by2z+b and p by b. Here z records the
weighted prefix selections. This is the exact algorithm used in the finite
all-cut controls, independently of the graph subset DP.

## 4. A bounded three-edge branch, then its nearest fourth-edge repairs

After the two-edge obstruction, the next search was the explicitly chosen
union of three patterns, not an unbounded expansion of flip radius:

1. three apex-incident arcs;
2. two apex-incident arcs and one old A–B cross arc;
3. all three arcs of any vertex triple.

Duplicate edge triples are removed. The tested counts for a=2,...,6 are
84,331,904,2005,3884. Only a=3 has target hits (seven, listed in the receipt).
For a=4,5,6 the minimum positive decrements are6,4,4 respectively.
This gives a real nearby construction n_a−1 at a=5,6, but not the requested
n_a+1.

We then tested one extra distinct arc at each minimum-decrement triple.
This is a precisely bounded continuation of the near misses. It is not
complete Hamming distance four. The a=4 branch tests210 such continuations
and finds one distinct four-arc repair:

    {(0,9),(2,9),(3,8),(3,9)},     H=465.

The a=5 and a=6 minimum triples have63 and88 fourth-edge continuations;
none reaches the target. Their minimizing triples are respectively

    a=5: {(4,6),(7,11),(9,11)},
    a=6: {(1,13),(3,13),(4,6)}.

The explicit a=4 repaired adjacency rows are in
`two_edge_structure_receipt.json`. All bounded hit and absence scopes are
retained in `three_edge_targeted_receipt.json`.

## 5. Why the a=4 success does not repeat: an exact extra spectral term

The natural rank-preserving continuation of the four-arc repair is

    F_a={(0,r),(a−2,r),(a−1,2a),(a−1,r)},      r=2a+1.

For a≥4 its exact Hamilton-path count is

    H(T_a flipped on F_a)
      =7−45·2^(a−3)−4·3^a+7·2^(2a−3)+2^(a−1)3^a.        (4)

This is proved by a finite compressed path-word calculation, not fitted
from small a. Keep the five fixed vertices

    Afirst, Apenultimate, Alast, Blast, apex.

Let P be the a−3 intermediate A vertices, and S the a B vertices before
Blast. These are transitive modules; the only inter-module direction is
P→S. Between fixed vertices a word has at most one P run and one S run,
so the skeleton language is finite. Expanding its integer Stirling
coefficients into exponentials gives, for module lengths j,k≥0, the terms

    0^j0^k −3·1^j0^k −0^j1^k +7·1^j1^k
    −0^j2^k +2^j0^k −5·1^j2^k −5·2^j1^k
    −4·1^j3^k +7·2^j2^k +4·2^j3^k.

Substitute j=a−3,k=a. For a≥4 all zero-base terms vanish and the result is
(4). The executable explicitly generates every compressed word, uses the
proved acyclic-symbol cutoff, and confirms the graph interface for a=3,...,7.
The exceptional empty-module boundary a=3 has H=9 and is retained in the
receipt rather than passed through (4).

Subtracting n_a+1 from (4) gives

    6−21·2^(a−3)−9·2^(2a−3)+(2^(a−1)−4)3^a.              (5)

It vanishes at a=4. For a≥5, use 2^(a−1)−4≥2^(a−2). Then (5) is at least

    6+2^(a−3)(2·3^a−9·2^a−21)>0.

The bracket is positive at a=5 and stays positive by induction (tripling
the first term versus doubling the second only increases the margin).
Thus this recurring rank pattern hits **only a=4** among a≥3.

The mechanism of failure is informative: changing the host arc creates
three-run channels and hence 3^a terms, while the requested arithmetic
family is purely dyadic. At a=4 these terms happen to cancel the lower
pieces; afterward their positive contribution dominates. A successful
fixed-pattern generalization would have to control these channels, rather
than infer uniform repair from the isolated ten-vertex witness.

## 6. Remaining explicit lemma

The all-apex gap rules out a large and natural class completely. The finite
shells isolate what remains: joint host/apex changes can repair a=3 or4,
but no repeatable rank pattern has been found. The next exact target is a
bounded host-arc template whose compressed-word response has all unwanted
3^a and higher exponential coefficients zero, while its dyadic response
is exactly n_a+1. Alternatively a source replacement could avoid this
near-maximum cut carrier altogether. Neither is supplied here.

## Reproduction

From the repository root:

    g++ -O3 -std=c++17 05-knowledge/results/session260_quasiperfect_chanpak_20261111/constructors/two_edge_repair.cpp -o /tmp/two_edge_repair260
    /tmp/two_edge_repair260 2 6 2
    /tmp/two_edge_repair260 2 6 3
    python3 -B 05-knowledge/results/session260_quasiperfect_chanpak_20261111/constructors/two_edge_repair.py

The C++ loops use unsigned64-bit exact counts; at the largest tested order14
even14! is below2^64, so no Hamilton-path count can overflow. The binary
is a disposable build artifact outside the repository. No external source,
library, numerical approximation, or isomorphism classifier is used.

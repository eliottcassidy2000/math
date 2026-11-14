# Cycle rotors, path–cycle moments, and the closed/open metagraph

Here cycles are directed and unrooted (rotation only); a path is
ordered. Inherited: [THM-4729](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md)
identifies the closed half of the path-tiling cube and proves closed
mass nC(T)/|Aut T| at a class. [The marked-path packet](
../marked_path_tournament_series_20261021/session23-symmetry/README.md)
counts marked paths and rigidly labelled scaffolds, and names a
nonrigid cycle scaffold as a next problem. The three calculations
below resolve that refinement and give new mixed moments. None is
a solution of the Chan–Pak realization problem.

## 1. Exact cycle-marked orbit and symmetry census

Fix n≥3. Complete the cyclic scaffold i→i+1 modulo n, then
quotient its completions by cyclic rotations. These are the
isomorphism classes of pairs (T, Hamiltonian cycle). For each
odd divisor d of n let m=n/d and define

    R(n,d)=n(n−3)/(2d).

A rotation subgroup of order d fixes exactly 2^R(n,d)
completions: every chord orbit has length d (a smaller orbit
would require exchanging two endpoints, forbidden at odd d).
An even rotation subgroup would contain a permutation that
exchanges the endpoints of some pair and fixes no tournament.
Therefore the ordinary unrooted cycle-marked number is

    A(n) = (1/n)∑_{d|n, d odd} φ(d) 2^R(n,d).

The sector at exact rotation stabilizer d has

    A(n,d) = (d/n)∑_{d|e|n, e odd} μ(e/d) 2^R(n,e).

This follows either from divisor inversion of the fixed counts
or by direct Burnside and then keeping the orbit sizes.

For converse merging, reverse cyclic order as you take the
opposite tournament. An anti-reflection is i↦s−i together
with orientation reversal. Among completions fixed by a
rotation subgroup of order d:

* For odd n any reflection fixes 2^D completions, with
  D=(n−3)(m+1)/4.
* For even n the n/2 vertex-axis reflections (even s) fix
  none, while each edge-axis reflection (odd s) fixes 2^D
  completions, where
  D=(n/2−2)m/2+floor(n/4)−1+ceil(m/4).

To prove this, for 2≤l<n/2 let the chord of step l start at i
modulo m. Reflection on its bit sends i to s−l−i. At odd m
there is one fixed position and (m+1)/2 free bits per step.
At even m there are two fixed positions precisely for odd l
and compatible odd s. There are floor(n/4)−1 such odd steps.
For the diameter step at even n (hence even m), the bit at i+m/2
is opposite to the bit at i. An even s forces one pair to equal
its own opposite because 2i=s has a solution. At odd s the
reflection has one fixed antipodal pair iff m/2 is odd; this
pair is fixed without a contradictory sign. All other diameter
pairs are paired by reflection, giving ceil(m/4) free bits.

Put B(n,d) for the self-converse cycle-marked orbit count in
sector d. Invert 2^D over odd multiples e of d as above and
multiply by 1 for odd n, 1/2 for even n. Indeed a completion
with stabilizer d that admits one anti-reflection admits
exactly d, and it has n/d distinct rotated completions. Thus
when summing fixed completions over reflections, every fixed
cycle-marked orbit contributes exactly n. Accordingly the
converse-merged count in sector d is (A(n,d)+B(n,d))/2. In
total B(n) is 2^((n−3)(n+1)/4) for n odd, and
2^(n²/4−n/2−2) for n even.

Changing the cut of a closed Hamiltonian path projects the
closed tilings onto these cycle-marked orbits. The fiber over
an orbit in sector d has m=n/d cuts, not always n. A self-converse
cycle-marked orbit has exactly one blue cut for odd n and exactly
two for even n. To see it, a gap at edge i→i+1 is sent to gap
s−1−i on the quotient of cuts modulo m. The equation
2i=s−1 mod m has respectively one or two solutions. A
non-self-converse cycle-marked orbit has no blue cut. These
facts recover the inherited global closed and blue totals and
distinguish the cycle stabilizer inside their fibers. In
particular no naive division of the closed-tiling total by n
can count *unrooted marked cycle orbits* in general.

## 2. A closed path–cycle joint moment

Let H_k(T) count ordered directed k-vertex paths, so H_n=H.
For 1≤k≤n, the number of joint marked orbit classes is

    ∑_[T] H_k(T)C(T)/|Aut(T)| =
    2^(binom(n,2)−n−k+1)/(n−k)! · J(n,k),

where

    J(n,k) = (n−1)! +
      ∑_{t≥1,2t≤k−1} (n−2t−1)!
        ∑_{r≥1} 2^r binom(t−1,r−1) binom(k−2t,r).

The marks have trivial joint stabilizer: fixing the path
fixes a vertex, and an automorphism preserving the directed
cycle acts by rotation. Fix the marked path on the first k
vertices and count oriented cycle compatibility. For every
spine edge (i,i+1), a candidate cycle has multiplicative
weight 2 if it uses i→i+1, zero if it uses the reverse, one
otherwise: its orientation factor is 1+I_forward−I_reverse.
The free edges provide the powers of two outside J.

Expand that product. At every chosen run of consecutive
spine edges the orientations must be coherent, since
otherwise the cycle would use a vertex with two incoming
or two outgoing edges. A run of length l contributes
1+(−1)^l after summing both orientations. Thus all runs
must be even, and their weights are 2 each. A selected
subset of size 2t split into r runs is specified by a
composition of t into r positive lengths, and the gaps
placing the runs on k−1 available edges. These give the
two binomial factors. Contract each coherent run and the
remaining directed cycles are counted by (n−2t−1)!
(when one block remains its last end joins its first).
The factor (n−k)! removes labellings outside the first
path, justified also by the free joint action.

On the labelled tournament distribution, the resulting
ratio E(H_k C)/(E H_k · E C) is J(n,k)/(n−1)!.
It is strictly larger than 1 for k≥3, and exactly 1 for
k≤2. For full paths its first correction is 2/(n−1),
with higher terms O(n^−2): an even-run subset can be
uniquely divided into adjacent dominoes, and its weight
is at most 2^t. The total relative t-term is therefore
at most 2^t/[t! (n−1)(n−2)…(n−t)], bounded by
(4/(n−1))^t/t! for t≤(n−1)/2.

## 3. Directional closure bias in the tiling-line metagraph

Under the inherited nonpath-arc complement, every path
presentation has one opposite, and a closing arc is
present at precisely one endpoint of that line. Count
edges with multiplicity and compare the mean actual
C(T) on the closed endpoint to that on the open endpoint.
The exact difference is 2^(1−n) times

    2·[n even] +
    ∑_{t≥1, 2t<n} (n−2t−1)!
      ∑_{1≤r≤min(t,n−2t)}
        2^r (2t/r) binom(t−1,r−1) binom(n−2t−1,r−1).

The common stem is a path with its first/last pair free.
Conditional on the stem, pinning the closing edge in
opposite orientations changes its cycle compatibility
weight by 2(I_forward−I_reverse). Expand the remaining
spine weights as above. The selected edge sets are now
circular runs on the n cyclic edges, must contain the
marked closing edge, and again odd runs cancel. For
2t<n, circular subsets with r even runs number
(n/r) binom(t−1,r−1) binom(n−2t−1,r−1), since we choose a
first start and ordered compositions of the positive run
halves and the positive gaps. The fraction containing
the marked edge is 2t/n by rotation symmetry. If all n
edges are selected only the two coherent cycles survive,
with combined sign 1+(−1)^n. This proves the formula.

In particular the average bias is positive for all n≥3
(the two-edge run containing the marked edge already
works). This is an aggregate, not a monotone orientation
of each metagraph line. For example at order six the
closed rows (22,60,8,17,36,13) have two directed Hamilton
cycles, while the open partner (42,4,57,50,35,2) has
three, using zero-based outgoing masks.

The aggregate equality survives merging converse vertices **as
a multigraph retaining all lines and loop incidences**:
(C(T)=C(T^{op})) and the two incidences keep their closure
statuses. It need not survive deduplicating parallel edges or
counting each fundamental black orbit once while leaving every
blue line at weight one. The sum here is over actual original lines.

## 4. All path and cycle lengths still do not determine a merged node

The following order-seven row masks give two rigid tournaments:
$$
 A=(34,68,33,87,71,90,5),qquad
 B=(70,64,98,39,47,3,56).
$$
Both have entire directed ordered path inventory by vertex
length 1 through 7 equal to (7,21,53,107,161,151,65), and
entire directed unrooted cycle inventory by vertex length 3
through 7 equal to (9,15,16,11,3). Sorted score sequences
(2,2,2,2,4,4,5) and (1,2,3,3,3,4,5) distinguish the pair
even up to taking converse (the complement of the first
sequence is (1,2,2,4,4,4,4)). Both nodes are pure black and
have equal ordinary closed/open incidence totals 21 and 44.

Yet even the neighbor H inventory of the metagraph differs:
the minimum H among the other endpoints is 33 for A and 23
for B, reached through closed incidences in both cases.
Literal verification counts all marked paths and reverses all
off-path arcs at each one. Merging converse cannot remove this
difference: H is unchanged by converse, all neighboring marked
classes coming from the converse representative have the same
neighbor H inventory. Thus knowing every simple path and cycle
count and automorphism size loses even a coarse piece of actual
wiring, not merely the exact names of the neighboring nodes.

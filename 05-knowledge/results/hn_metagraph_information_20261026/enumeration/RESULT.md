# A complete family enumerator, with owners and literal metagraph reconstruction

**Status: PROVED elementary, independently audited; declared computations
FINITE-EXACT.** See the [publication audit](../AUDIT.md). The underlying
weighted-core classification, alternating core birth,
and abstract canonical-augmentation method are **INHERITED PROVED**. This packet
materializes their combination without duplicate generation, exports literal
representatives and birth records, and proves a quantitative family-coverage
bound. The declared finite checks are **FINITE-EXACT**, not an infinite atlas.

There is already an all-order structural answer in this repository: exact core
objects, their automorphism actions, positive chain weights, and alternating
attachments generate every tournament class. There is no need to postulate a
new scalar recurrence or to call the known construction open. What remains
expensive is generating and identifying the new core objects at each order.

## 1. Inheritance, corrections, and the actual addition

The closest mechanism is
[weighted_core_transport/core_count.md, sections1--3 and7](../../weighted_core_transport_20261020/core_count.md).
It already proves the weighted-core bijection, its rational response, an
alternating next-core construction, and the symmetry-corrected child count.
Its [growth_check.py](../../weighted_core_transport_20261020/growth_check.py)
already generated all classes through eight. The corresponding
[inheritance.md, sections3--5](../../weighted_core_transport_20261020/inheritance.md)
contains the labelled transition kernel, successful-phase counts, and the
canonical-root option. None is claimed as new here.

The underlying canonical chains are
[THM-4463, homogeneous tournament contacts, section2](../../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md).
The later [canonical-owner forest](../../ports_owners_f1_20261025/anchor.md)
gives a different feedback-rank-preserving construction. Its regular-pentagon
hostile shows why a canonical child need not project to a canonical parent.
Here the parent's weighted-core owner and the child's distinguished vertex
orbit are separate tests.

The canonical hostile is two six-vertex cores with the same automorphism-group
order but different actions. The stronger same-cycle-index/future-child hostile
is already in the inherited growth checker. The corrected historical near miss
is `SC(n)=A000568(n-1)`: see
[the correction](../../../../07-reflections/sc-equals-a000568-identity-corrected.md)
and MISTAKE-049 in `01-canon/MISTAKES.md` (the correction note's older ID is
stale). At order seven the actual figures are
`SC(7)=88`, `A(6)=56`; these are different quantities.

The least-used sidecar for making the structural recurrence executable is the
**entire set of minimizing labelings**. A single chosen labeling is not an
invariant choice of a physical vertex. Its orbit is invariant.

The Anchor is literal complete enumeration; the Niche is recovery of metagraph
line incidence from family data; the Wildcard is a rigorous coverage budget for
old families. The live board is core objects / weight actions / root ownership /
path fibres / high-order count gauges / truncation loss.

No literature-priority claim is made. [OEIS A000568](https://oeis.org/A000568)
was checked on2026-10-26 for the sequence and Davis's classical odd-partition
formula. The formula is not a new result of this session.

## 2. Exact state, map, and complete recurrence

For a tournament `T`, join `u,v` when every outside vertex sees them alike.
THM-4463 says this graph is a disjoint union of consistently directed paths.
Their vertex sets are the unique maximal transitive modules, called chains
below. Contract them to a pair-free tournament `Q`, meaning no homogeneous pair.
Pair-free does **not** mean indecomposable with respect to arbitrary modules.

For a fixed representative `Q` on `q` vertices, let `G=Aut(Q)`. A family member
has positive weights `w=(w_i)` and is the inflation `Q[w]`: vertex `i` becomes
the transitive chain of length `w_i`. The exact class parameter is

    (Q, G-orbit of w),       sum_i w_i=n.                       (E1)

Every class has exactly one such parameter with one representative per core
class. Inside each chain the order is rigid, so

    Aut(Q[w]) = Stab_G(w),                                    (E2)

with its action on the chains and the unique order-preserving lift retained.
This is the inherited class bijection, not only an equality of counts.

Choose the least tuple in each weight orbit. To generate new cores of order
`n+1`, traverse all family members of order `n`, writing their ordered chains
as `L_i=(v_i0,...,v_i,w_i-1)`. A candidate attachment is a phase vector
`b in {0,1}^q` and has arcs

    x -> v_ij  iff b_i xor (j mod2)=1.                         (E3)

Retain one phase representative under `Stab_G(w)`. Reject a candidate exactly
when `x` has a homogeneous neighbor.

**Inherited alternating criterion, recalled with both directions.** A pair of
old vertices cannot become homogeneous after adding `x` unless it already was
homogeneous in the parent. The old homogeneous pairs are precisely adjacent
positions within the chains. They all cease being homogeneous iff the new
incidence bits alternate on every chain. After this, the only remaining
possible homogeneous pairs involve `x`. Thus the stated filter is necessary
and sufficient for the child to be pair-free.

Every pair-free child is reached: delete any vertex `x`, take the canonical
weighted-core description of the parent, and read its alternating phases.
Isomorphisms of rooted children `(U,x)` are exactly parent isomorphisms carrying
these attachment bits. Equations(E1)--(E2) therefore prove that parent family
owners followed by phase orbits generate **each rooted pair-free child once**.
Root forgetting still has nonconstant multiplicity and needs the next step.

## 3. A concrete invariant deletion owner, without child deduplication

Encode a tournament by its upper-triangular arc bits, increasing first endpoint
and then second endpoint, bit1 meaning the earlier label points to the later.
Build an invariant individualization/refinement tree. Start with cells in
nondecreasing outdegree order. Repeatedly split each cell by the vector of
outneighbor counts into the current cells, sorting these vectors and preserving
the old cell order. When stable, choose the first nonsingleton cell and branch
over every choice of a vertex made singleton before its remainder. Refine again.
At a discrete partition read its singleton order as a labeling. Let `c(T)` be
the least integer code among these leaves, and let `P(T)` be every leaf labeling
achieving it. A labeling maps new positions to old vertices. Define

    O(T)={p(0): p in P(T)}.                                   (E4)

For the one-vertex object this is its sole vertex. The empty object is only the
initial parent and has no deletion test.

**Owner lemma.** `O(T)` is exactly one automorphism orbit, is nonempty, and is
equivariant under all tournament isomorphisms.

**Proof.** Any two minimizing labelings differ by an automorphism of the
canonical representative. Conversely composing a minimizing labeling with any
automorphism of `T` gives another minimizing labeling. Consequently their
position-zero images are exactly one orbit. Relabeling `T` transports the whole
set of minimizing labelings, proving equivariance. Each initial cell, signature,
first nonsingleton-cell choice, and branch family is equivariant. Therefore
isomorphisms carry the whole search tree bijectively to the target search tree.
Every automorphism carries leaves to leaves, so no automorphism is lost by
refinement. The code need not be the minimum over all degree-sorted labelings;
the invariant search tree, and retaining all its minimizers, is sufficient. QED.

Keep a rooted candidate only if its new vertex `x` lies in `O(U)`. By section2
every rooted-child orbit appears once; (E4) chooses exactly one of those vertex
orbits for every unrooted child. Therefore:

**Complete owner theorem.** Starting from the empty tournament, repeatedly
generating cores by(E3)--(E4), then generating family members by(E1), produces
every tournament isomorphism class exactly once at every finite order. The
number of output records is `A000568(n)`. Each stage is a finite terminating
procedure; it does not use known sequence values to decide completeness.

The implementation does not deduplicate accepted children: a repeated code
raises an exception. Only canonical phase and weight orbit representatives
are traversed. Codes are retained to verify uniqueness and support queries.
All minimizing labelings also give the exact core automorphism group. The
cheap test that `x` has minimum outdegree may precede full canonicalization,
because every member of(E4) has minimum outdegree; it is only a necessary test.

The new executable artifact is this explicit ownership specialization and
portable representative/birth bank. The abstract canonical-construction-path
principle and inherited core-birth theorem are not re-promoted as discoveries.
No polynomial complexity or efficient general graph-isomorphism claim follows.

## 4. From family members to every literal metagraph line

The exact source is(E1) with the actual arcs of `Q`, not just its cycle index.
For every member `T`, enumerate its directed Hamiltonian paths `p`. Relabel
positions along `p` so that `0->1->...->n-1` is fixed, and record the other
`m=binom(n-1,2)` arcs as `t in X_n={0,1}^m`.

Two paths yield the same `t` iff their position-by-position map is an
automorphism of `T`. The action on ordered paths is free. Hence the fibre of
the ordinary class `[T]` has exactly

    F(T)=H(T)/|Aut(T)|                                       (E5)

distinct tilings. Different classes have disjoint fibres. Every tournament has
a directed Hamiltonian path: insert a new vertex before the first vertex it
beats in an existing directed path, or at the end if it beats none; the
preceding vertex, when present, beats it. Induction proves existence. Thus the
family fibres cover all `X_n`, and

    A(n)=sum_(t in X_n) |Aut(T_t)|/H(T_t).                   (E6)

For `n>=3`, complementing every free coordinate gives the fixed-point-free
involution `C`. Its two-element orbits are literal lines. Reflection permutes
free coordinates `(i,j)->(n-1-j,n-1-i)` and is denoted `R`. It corresponds to
converse followed by reversal of the path, so it preserves converse-merged
class nodes. A line is blue iff `Rt=t`, and black otherwise. Preserve loops
and parallel lines when pushing endpoints to these nodes.

The involutions commute. Moreover `Rt=Ct` is impossible because the outer
chord `(0,n-1)` is a fixed coordinate of reflection and cannot equal its own
complement. Hence reflection freely pairs black lines, even within a fixed
unordered pair of merged endpoints. Every black multiplicity is even. These
are inherited conventions and laws from
[THM-643, blue/black line parity](../../../../01-canon/theorems/THM-643-blue-black-line-parity-structure.md),
[THM-796, three-sorted incidence](../../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md),
and [the paired-black packet](../../merged_metagraph_black_quotient_bridges_20261020/README.md).

Converse is an involution on ordinary classes, with `SC(n)` fixed classes.
The exact merged-node count is

    M(n)=(A(n)+SC(n))/2,       A(n)=2M(n)-SC(n).             (E7)

It is not `A(n)`, and `SC(n)` is not generally `A(n-1)`. The exported lines
include their endpoint codes, colour and multiplicity. This is complete finite
incidence reconstruction, not a proof of any conjectured black connectivity.
Quotienting families to their weight-count response loses the arc and path
data used by(E5)--(E6); no line transport is claimed for that quotient.

## 5. A quantitative obstruction to a fixed complete family bank

Let `B(n,m)` count order-`n` classes whose canonical pair-free core has order at
most `m`, with `1<=m<=n`. Allow inflation of **every** order-`m` tournament,
including those that themselves have clones, by positive chain weights of
total `n`.

**Coverage theorem.** The image is exactly the set counted by `B(n,m)`, and

    B(n,m) <= binom(n-1,m-1) A(m).                          (E8)

**Proof.** Inflation can only concatenate or lengthen the canonical chains of
the parent, so its pair-free core has order at most `m`. Conversely take a
tournament with `q<=m` canonical chains. It has `n-q` intrinsic adjacent gaps.
Cut at any `m-q` of them. The `m` resulting nonempty intervals are transitive
modules; their quotient is an order-`m` tournament and their lengths are a
positive weight vector reconstructing the original object. A fixed labelled
representative of each parent class and all `binom(n-1,m-1)` compositions
therefore give a surjection. Collisions only improve the upper bound. QED.

For a fixed output `T`, the exact inflation fibre consists of `Aut(T)`-orbits
of these `(m-q)`-element cut subsets. This is the fixed-size layer of the
[inherited clone-cut fibre theorem](../../ports_regular_actions_hn_20261025/counting/CLONE_INVOLUTION.md),
not a claim that all fibres have the same cardinality.

In particular, if `c(n)` counts new pair-free cores,

    A(n)-c(n) <= (n-1) A(n-1).                             (E9)

Exact count arithmetic certifies that families with core order below100
cover at most `2^-85` of the order100 classes; below128 they cover at most
`2^-113` of the order128 classes. These are rigorous upper bounds from(E9),
not approximate proportions or literal high-order core lists. The remaining
fractions are new pair-free cores.

For any fixed cutoff `m`, (E8) is polynomial in `n`, whereas
`A(n)>=2^binom(n,2)/n!` (one isomorphism class has at most `n!` labelings).
Thus every fixed finite core bank has asymptotic coverage fraction zero.
This decisively distinguishes complete all-order **generation of new families**
from trying to expand a fixed finite bank forever. It does not rule out an
operation-stable compressed count state or the existing fast arithmetic DP.

## 6. Counting history and current engine, with the types kept separate

The count-only checks call the existing
[lca_ports_reverse/counting/reverse.py](../../lca_ports_reverse_20261025/counting/reverse.py).
It uses the retained defect `D(d)=d*m-sum gcd(c,d)`, an integral odd-factorial
gauge, and the reverse minimum valuation potential. Its normalized exact
integer is `A(N)/2^B`, where `B=floor(N/2)` for `N>=3`. This baseline is not
always the final valuation: `N=50` has one additional factor two. The script
checks the full exact value and a32-bit normalized residue, with every
reachable intermediate gauge replayed.

The historical unit-tail kernel is
[THM-514, A000568 unit-tail core kernel and prime-tail ladder](../../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md).
Its “core” is an odd permutation partition without1-parts, not an actual
pair-free tournament. Its `[mass,cycle count]` outer compression still requires
the inner gcd response. The odd-order valuation theorem is
[THM-305, two-adic valuation of the tournament count](../../../../01-canon/theorems/THM-305-two-adic-valuation-tournament-count.md).
Neither object supplies literal representatives or metagraph incidence.

The signed, clone-excess integer recurrence in
[port_residual_integral/counting](../../port_residual_integral_20261025/counting/README.md)
computes pair-free counts after the cycle-index inverse substitution. Actual
cycle degree and clone excess are different coordinates. Our positive literal
family construction neither replaces that signed engine nor identifies its
intermediate integer states with class sets.

Likewise the old SCC/Lyndon route
[tournament_even_bijection_20261011c.md](../../tournament_even_bijection_20261011c.md)
already connects ordered strong-component assembly and the even-graph count.
It is complementary inherited information, not a new clone-family proof.

## 7. Artifacts, finite universe, and reproduction

`enumerate.py` is standard-library-only and imports no previous tournament
enumerator. The default receipt ends at eight; the completed extended receipt
certifies these exact census rows through nine:

|Order|All classes|New pair-free cores|Accepted rooted phase orbits|
|---:|---:|---:|---:|
|1|1|1|1|
|2|1|0|0|
|3|2|1|1|
|4|4|2|4|
|5|12|5|17|
|6|56|26|132|
|7|456|255|1693|
|8|6880|4602|36348|
|9|191536|147095|1319141|

The last column is **not** the new-core count. It equals the sum of vertex-orbit
counts of the new cores, and the deletion-owner test removes precisely that
root multiplicity. At order eight the code starts with42890 raw phases,
39866 parent-stabilizer phase orbits,36348 accepted rooted objects, and4602
owned unrooted cores. The eight-vertex family contributions by core order are
`1:1, 3:7, 4:26, 5:107, 6:444, 7:1693, 8:4602`.

At order nine the contributions are
`1:1, 3:10, 4:40, 5:210, 6:1168, 7:6664, 8:36348, 9:147095`.
There are 1427130 raw phases and 1393530 stabilizer phase orbits. The completed
final normal replay took about 64 seconds on this host, with every class and
birth retained in `representatives9.json.gz` (3653787 bytes). Gzip uses `mtime=0`.
This extends the inherited literal family census through eight by one order;
the value191536 itself is classical, not a new sequence term.

The earlier direct permutation canonicalizer reached the five-minute resource
budget before completing order nine. No census claim was taken from that
interrupted wrapper. A stale child-process overwrite was caught and repaired
before freezing; [REPLAY_NOTE.md](REPLAY_NOTE.md) records the provenance repair.
The final invariant refinement tree replaces that bottleneck;
the independent unrestricted-permutation check remains as a hostile control.

Independent finite controls are:

- all labelled tournaments and all unrestricted relabelings through order five,
  compared as complete class sets with the family output;
- labelled orbit mass and each labelled core stratum through the output bound,
  independently checked against the signed binomial formula;
- all32768 fixed-path tilings at order seven, decoded directly and compared
  with the Hamilton-path fibres reconstructed from the456 family members;
- exact order-seven metagraph data:88 self-converse classes,272 merged nodes,
  16384 literal lines,512 blue tilings, and every black endpoint multiplicity
  even;
- explicit same-order/same-group-order action hostile, with weight responses
  `1,4,11,26,54,102,180` versus `1,2,7,20,42,84,156` at totals6..12;
- root-forgetting hostile: the transitive triple has three vertex orbits and
  the cyclic triple has one, so a rooted count is not an unrooted count;
- exact count-only comparisons at0..9,50,76,99,100,127,128 with the repository
  b-file,32-bit normalized residues, and118949 intermediate gauge checks;
- all finite(E8) inequalities for the complete census through nine.

Run from the repository root:

```text
python 05-knowledge/results/hn_metagraph_information_20261026/enumeration/enumerate.py
python -O 05-knowledge/results/hn_metagraph_information_20261026/enumeration/enumerate.py --output 05-knowledge/results/hn_metagraph_information_20261026/enumeration/optimized.json
python 05-knowledge/results/hn_metagraph_information_20261026/enumeration/enumerate.py --max-order 9 --output 05-knowledge/results/hn_metagraph_information_20261026/enumeration/order9-normal.json --bank 05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz
python -O 05-knowledge/results/hn_metagraph_information_20261026/enumeration/enumerate.py --max-order 9 --output 05-knowledge/results/hn_metagraph_information_20261026/enumeration/order9-optimized.json --bank 05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz
python 05-knowledge/results/hn_metagraph_information_20261026/enumeration/check_counts.py
python -O 05-knowledge/results/hn_metagraph_information_20261026/enumeration/check_counts.py --output 05-knowledge/results/hn_metagraph_information_20261026/enumeration/counts-optimized.json
```

All gates use explicit exceptions and remain active under `-O`. The JSON
receipts pin the producer bytes; the count receipt also pins the inherited
engine and reference file. `representatives.json` retains canonical codes,
weighted family coordinates, exact automorphism orders, and one parent/phase
birth record for every core through eight. The compressed bank extends the
same schema through nine. Codes decode directly into adjacency rows. Normal
and optimized receipts and decompressed banks agree exactly; no duplicate
optimized bank is retained.

The all-order theorem concerns termination and exact coverage at every chosen
finite bound. The saved atlas and metagraph have the explicit finite bounds
above. Count-only values at100 or128 do not extend those literal atlases.

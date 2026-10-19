# Weak Gold from conditional symmetry and exchange overlap

2026-10-19. **PROVED reductions and sufficient conditions; FINITE-EXACT
selected witnesses.** The general seven-thin
Gold question remains open. The cyclic-triple criterion is explicitly
attributed to Peczarski; no historical-novelty claim is made for any part.

## 1. Convention and exact two-comparison criterion

Let P be a finite non-total poset, let `L(P)` be its linear extensions,
and write `E=e(P)=|L(P)|`. A comparison asks which of two distinct
vertices comes first in the hidden extension. Following the weak Gold
convention of
[packet 94](good_pair_gold_synthesis_20261019.md), the requirement on
every realized two-comparison path is

    E ≥ t₁+t₂.

Determined and repeated comparisons are permitted. Counts below are
uniform extension counts unless a weighted statement is explicit.

Fix a first comparison and orient its larger branch as the event B.
Put `b=|B|≥E/2` and `a=E−b`. For a proposed second comparison C, let
`c=|B∩C|`. The smaller first branch is automatic, because its surviving
count after any second comparison is at most a, and `2a≤E`.
Consequently the following are equivalent:

1. This first comparison can be completed to weak Gold.
2. Some comparison C has both larger-branch cells at most `E−b`.
3. Some comparison C satisfies

       2b−E ≤ c ≤ E−b.                              (1)

4. Some comparison C satisfies

       |2c−b| ≤ 2E−3b.                              (2)

Indeed the two nontrivial path inequalities are `b+c≤E` and
`b+(b−c)≤E`. A zero cell causes no problem. A determined first
comparison cannot give Gold: it leaves E extensions before a positive
second-path cell. Thus any successful first comparison is genuinely
incomparable and both first branches have positive integer counts.

In particular, a possible first comparison must have `b≤2E/3`.
Conversely, that bound alone need not supply the second comparison.
At `b=E/2` every second comparison works. At `b=2E/3`, the larger
branch must split exactly in half.

**Adaptivity gives no extra existence power at depth two.** If an
adaptive strategy succeeds, use its second comparison from the larger
branch on both branches. It still succeeds, since the smaller branch
is automatic. Thus an adaptive weak Gold certificate can always be
replaced by two comparisons fixed before either answer. This observation
aligns (1) with the fixed-comparison formulation rather than modifying
the conjecture. The argument also holds for arbitrary nonnegative
weights, with positive-weight paths as the realized paths.

## 2. Joint-count, covariance, and conditional-balance forms

Let `m=|C|` in the original extension set and put

    Δ=Ec−bm.

This is the determinant of the two-by-two joint count table. Equation
(2) is equivalent to the exact integer condition

    |b(2m−E)+2Δ| ≤ E(2E−3b).                        (3)

With `p=b/E`, `q=m/E`, and
`κ=Cov(1_B,1_C)=c/E−pq`, the same condition is

    |p(2q−1)+2κ| ≤ 2−3p.                            (4)

In particular, the triangle inequality gives the sufficient condition

    p|2q−1|+2|κ| ≤ 2−3p.

This last condition can lose useful cancellation. A conditional half
pair need not be an unconditional half pair: covariance can cancel its
unconditional bias exactly. No universal sign or smallness claim for
comparison covariances is being assumed.

More generally, if a comparison in the larger-branch poset
`R=P+B` has both conditional probabilities at most `θ`, then the
original first comparison gives Gold whenever

    (1+θ)b≤E.                                       (5)

Here `P+B` means adjoining the observed comparison and taking transitive
closure, not an induced subposet. Its uniform extension distribution is
precisely the original distribution conditioned on B. The cases
`θ=1/2` and `θ=2/3` give the thresholds `b≤2E/3` and `b≤3E/5`,
respectively. The conditional balance hypothesis is essential.

## 3. An unconditional exchange-overlap inequality

Let R be any finite poset, with `b=e(R)`. Let σ be an involution of its
vertex labels interchanging distinct vertices u and v. It need not be
an automorphism of R and may move other labels as well.

Relabel R by σ: `σR` has a relation `σx<σy` for each relation `x<y`
of R. Let J be the transitive closure of the union of the relations of
R and σR. If the union is cyclic, define `e(J)=0`. Put

    M=e(J),   d=b−M,   c=e(R+uv).

Then

    |2c−b| ≤ d=b−e(J).                              (6)

**Proof.** Relabeling gives `L(σR)=σL(R)`. Therefore

    L(J)=L(R)∩σL(R),

including the empty-intersection cyclic case. Since `σ²=id`, this
intersection is σ-invariant. On it σ interchanges the orientations
`u<v` and `v<u`, so each orientation has exactly M/2 extensions.
The remaining d extensions of R can contribute at most d to the
absolute orientation imbalance. This proves (6). ∎

The inequality is sharp even on three vertices. Let R have the sole
relation `u<z`, with v isolated, and let σ exchange u and v. Then
`b=3`, `M=2`, and `c=2`; hence `|2c−b|=b−M=1`. No universal smaller
coefficient multiplying d can replace the one in (6).

**Conditional exchange certificate.** For the larger branch R of the
first comparison in Section 1, any such involution gives weak Gold if

    d≤2E−3b,

or equivalently

    e(R∪σR)≥4b−2E.                                 (7)

The union notation in (7) always includes transitive closure and the
cyclic-zero convention. Both larger-branch Gold margins are at least

    E−(3b+d)/2.                                     (8)

On the smaller branch every second comparison works; repeating the
first gives the explicit margin `2b−E` there.

The statement is a sufficient structural certificate, not an equivalence:
invalid exchanges of opposite orientations can cancel, so actual
imbalance may be much smaller than d. It also does not posit that a
suitable involution exists in every seven-thin poset.

## 4. Restored symmetry after one comparison

If σ is an automorphism of R interchanging u and v, then `M=b` and
`d=0`. Thus any first comparison whose larger branch has size
`b≤2E/3` and acquires such a symmetry gives weak Gold. Its two
larger-branch margins are exactly

    E−3b/2.                                         (9)

They are strictly positive if `b<2E/3`. If also `b>E/2`, even a
repeated comparison on the smaller branch has strictly positive margin.
The original poset P can be rigid and have no good pair.

The branch-size restriction cannot be omitted. Take a three-element
chain together with an isolated point, then put a two-element antichain
above all four vertices by ordinal sum. There are eight extensions.
Comparing the chain minimum with the isolated point has a larger branch
of six extensions. That child has a symmetry swapping the two upper
vertices, but its equal split is `3+3`, and `8−6−3=−1`.
Indeed no second comparison can repair that chosen first comparison,
since `6>2·8/3`.

### Structural explanation of the rigid packet-117 witness

Let Q be the ten-point poset from
[packet 117](height_three_seven_thin_no_good_pair_20261019.md), with
strict predecessor masks

    Q=(0,0,0,6,0,19,2,18,86,151).

It is rigid, seven-thin, and has no good pair in either orientation of
the poset. Its known first comparison is `3<5`. Adjoining that relation
changes only the predecessor set of vertex 5, producing

    R=(0,0,0,6,0,31,2,18,86,151).

The permutation

    σ=(2 4)(3 7)(5 9)

fixes vertices `0,1,6,8` and is an automorphism of R. This follows
directly from downsets, without using the cell-count equality:

- `D(3)={1,2}` is carried to `D(7)={1,4}`;
- `D(5)={0,1,2,3,4}` is carried to
  `D(9)={0,1,2,4,7}`;
- `D(6)={1}` and `D(8)={1,2,4,6}` are invariant;
- the four empty downsets are carried among themselves.

The converse mappings hold because σ is an involution. Thus relabeling
R's extensions swaps `2<4` with `4<2`. The exact totals are

    E=11971,   b=6562,   c=b/2=3281.

The larger-branch margin in (9) is `2128>0`. Repeating on the smaller
branch gives margin `1153>0`. Using `2<4` on both branches instead
recovers packet 117's four cells

    (3281,3281,1899,3510)

and four margins `(2128,2128,4663,3052)`.

The second comparison is not unconditionally half-balanced:
`e(Q+2<4)=5180`, whereas `E/2=11971/2`. This is a concrete instance
of the covariance cancellation in (4).

For any ordinal sum containing a Q block, with arbitrary other finite
blocks, every count involving only that block is multiplied by the
product of the other blocks' extension totals. Its conditional symmetry
extends by fixing all other blocks. Consequently this mechanism covers
the entire rigid no-good-pair ordinal-sum family of packet 117, not
merely its ten-point seed.

## 5. Exchange without a symmetric larger branch

There is a closely related seven-thin witness where the chosen larger
branch itself is rigid. Adjoin `6<7` to Q and close transitively. This
also adds `6<9`, and gives

    Q′=(0,0,0,6,0,19,2,82,86,215).

Its incomparability degrees are

    (7,3,6,7,5,6,5,5,5,3).

It has height four and `E=7249` extensions. It has no good pair:
exactly 42 ordered orientations in Q′ and its dual meet the structural
conditions, and their smallest numerator is `3804>7249/2`. The
accompanying receipt lists all 42, and full extension enumeration
independently verifies the counts. Absence of a very good pair is not
being substituted for this probability-sensitive test.

Q′ is rigid because its ten signatures `(|D(v)|,|U(v)|)` are distinct:

    (0,2),(0,6),(0,3),(2,0),(0,4),
    (3,0),(1,3),(3,1),(4,0),(6,0).

Again take `3<5` as the first comparison. Its larger branch has

    R′=(0,0,0,6,0,31,2,82,86,215),   b=4072.

R′ is also rigid. Its signatures are

    (0,2),(0,6),(0,4),(2,1),(0,4),
    (5,0),(1,3),(3,1),(4,0),(6,0).

Only vertices 2 and 4 share a signature. Every other vertex is
therefore fixed; fixed vertex 3 lies above 2 but not 4, so those two
cannot be interchanged either.

Nevertheless the simple transposition `τ=(2 4)`, which is not an
automorphism of R′, gives

    M=e(R′∪τR′)=3192,
    d=b−M=880,
    2E−3b=2282.

Thus (7) certifies Gold despite the absence of any nontrivial symmetry
in this chosen larger branch. The exact heavy-branch cells are
`2221` and `1851`, with imbalance `370≤880`; the defect estimate is
sufficient but not tight. Equation (8) guarantees each heavy margin
at least 701; the actual two margins are 956 and 1326.

Using `2<4` on both first branches gives the complete four cells

    (2221,1851,1224,1953)

and Gold margins

    (956,1326,2848,2119).

Thus this is another no-good-pair poset with strict Gold, and it
illustrates a genuinely non-automorphism use of the overlap theorem.
It is not claimed to exclude every possible restored-symmetry strategy
with some different first comparison.

The same count multiplication applies to ordinal sums of copies of Q′,
with optional inserted chains. These remain seven-thin, rigid, without
good pairs or half pairs; their selected larger branch remains rigid.
For rigidity, the downset-size intervals of successive blocks are
disjoint, so every automorphism preserves each block. For absence of
good pairs, within-block downset inclusions, upset differences, and
comparison probabilities are unchanged. The selected certificate has
minimum normalized fixed-comparison margin `956/7249` for every number
of copies. This is an immediate family corollary, not a minimum-order
claim.

## 6. Classical cyclic triples and the majority-order residual

Peczarski's cyclic-triple Gold criterion is stated explicitly as Fact 3
on page 108 of
[*Comments on the Golden Partition Conjecture* (2017)](https://cdm.ucalgary.ca/article/download/62349/46785/178291),
where it is credited to Lemma 3 of *The Gold Partition Conjecture*,
Order 23 (2006), 89–95. Page 109 also uses the resulting transitive
majority-order argument. These are inherited results, not new claims.

For clarity, the following short proof works for any nonnegative
weights on the six orders of three labels, not just uniform poset
extensions. Normalize their total to 1, and write

    a=w(xyz), b₀=w(yzx), c₀=w(zxy),
    d₀=w(xzy), e₀=w(yxz), f₀=w(zyx).

Suppose all three cyclic comparison probabilities exceed one half.
Cyclically relabel so

    p=Pr(x<y)=a+c₀+d₀
      =min{Pr(x<y),Pr(y<z),Pr(z<x)}.

The sum of those three probabilities is
`1+a+b₀+c₀≤2`, so `p≤2/3`. Since
`b₀+e₀+f₀=1−p`, the other two majority probabilities being at least p
gives

    a≥2p−1+f₀≥2p−1,
    c₀≥2p−1+e₀≥2p−1.

At least one of a and c₀ is at most `1−p`: otherwise
`p≥a+c₀>2−2p`, contradicting `p≤2/3`. In the larger branch `x<y`,
comparison `y<z` isolates a and comparison `z<x` isolates c₀. Choose
one in the interval `[2p−1,1−p]` and apply (1). The smaller branch is
automatic. This proves the weighted cyclic-triple statement directly.

It follows that any weak Gold counterexample has no half pair and no
directed triangle in its strict majority tournament. A tournament
without a directed triangle is transitive, so that tournament is a
linear order extending the poset. The term tournament is literal here:
vertices are poset elements, each edge records the strict comparison
majority, tied pairs have already been excluded, and the retained
numeric probabilities remain necessary for (1) and (7). The bare
orientation discards those margins and cannot alone certify Gold.

## 7. Minimal-counterexample reduction and the exact remaining obstacle

Consider any class of finite posets closed under adjoining an orientation
of an incomparable pair. Seven-thin posets form such a class: adding
relations cannot increase an incomparability degree and strictly reduces
the total number of incomparable pairs.

Suppose the class has a non-total weak Gold counterexample. Choose one
with the least number of incomparable pairs. Then it has no comparison
whose probability lies in `[2/5,3/5]`.

**Proof.** Such a comparison would have larger branch size `b≤3E/5`.
If its larger-branch poset R is non-total, it has fewer incomparable
pairs and hence has Gold by minimality. The first comparison of any
Gold strategy in R is `1/3–2/3` balanced: its larger cell is at most
`2b/3`, by the necessary bound already proved in Section 1. Use that
comparison as the second comparison in the original P. Equation (5)
with `θ=2/3` gives Gold, a contradiction.

If R is total, then `b=1`. Both orientations of the original incomparable
pair occur, so the smaller branch has at least one extension and at most
b. Hence `E=2` and both branches have size 1. Comparing and repeating
gives Gold with equality. No Gold assertion for a total child is used.
This handles the terminal case explicitly. ∎

For the seven-thin class, packet 94 also rules out a good pair in a
minimal counterexample. Combining this with Section 6 gives a precise
necessary residual description:

- no good pair in the poset or its dual;
- a transitive strict majority order extending the poset;
- every majority probability strictly greater than `3/5`;
- for every potential first comparison with `b≤2E/3`, every second
  comparison violates (1), and every candidate exchange involution
  fails the sufficient inequality (7).

For a root probability `p∈(3/5,2/3]`, a conditional `1/3–2/3` pair
is no longer enough in general: the required larger conditional cell
is at most `(1−p)/p`, which decreases from `2/3` to `1/2`. Controlling
this narrower conditional window is the remaining obstacle. Restored
symmetry reaches its endpoint; exchange overlap interpolates toward it.
Neither theorem guarantees that the needed root or exchange exists.

A sharp, falsifiable candidate sufficient for completing the seven-thin
argument is:

> **OPEN candidate dichotomy.** Every non-total seven-thin poset has
> either a good pair (possibly in its dual) or a comparison whose
> probability belongs to `[2/5,3/5]`.

If this statement were proved, the minimal-counterexample argument and
packet 94 would prove weak Gold for every seven-thin poset. It is
stronger than anything established here and may be false. A refutation
would be a seven-thin no-good-pair poset with every comparison outside
that closed interval; it need not be a Gold counterexample.

The good-pair exception is essential: a two-element chain plus an
isolated point is seven-thin and has only `1/3` and `2/3` probabilities
on incomparable pairs, but has a good pair. Q, Q′, and their ordinal
families do not test the candidate's hard gap: their displayed first
comparisons already lie strictly between `2/5` and `3/5`. Cheap checks
of three previously recorded examples (packet 117's twelve-point
witness, a nine-point collar witness, and an eight-point rooted-five
obstruction) also do not falsify it. This small control set is not
evidence for a universal conclusion, and no new broad census was run.

The four-page 2017 primary paper was checked in full. It states the
cyclic-triple and majority-order reductions, but not the `3/5`
minimal-counterexample reduction. A bounded access check did not obtain
the full 2006 paper, so no historical-priority inference is made from
that absence; the `3/5` argument is given above as an elementary proof.

## 8. Exact controls and scope

Run, from the packet directory,

    python3 05-knowledge/experiments/verify_conditional_gold_exchange_20261019.py

The standard-library checker has 1,993 active exact guards. It compares
ideal dynamic programming with full extension enumeration for Q and Q′,
checks every eligible good-pair orientation, verifies overlap both as a
relation-union extension count and as an explicit intersection of word
sets, and verifies the paired orientations. It checks two-block ordinal
count multiplication, the rigidity signatures, the sharp three-point
defect example, the symmetric-but-too-heavy hostile example, and the
`E=2` total-child base case.

Supplementary algebra controls comprise 354 nonzero four-cell integer
tables, 99 cyclic six-order weight laws with entries from 0 through 3,
and 224 conditional-transfer integer cases. They are regression checks
for the uniform proofs, not poset classifications. Normal and optimized
runs produce identical JSON receipts; all guards remain active under
`python -O`.

The good-pair definition is Zaguia's probability-sensitive Definition 1
in the [author's primary paper](https://arxiv.org/pdf/1610.00809v3),
checked directly on 2026-10-19. The inherited good-pair-to-Gold theorem
is packet 94; the Q construction and its original no-good certificate
are packet 117. The research focus of this packet is the
conditional symmetry/overlap mechanism and its explicit rigid witnesses,
without claiming historical novelty or resolving the general conjecture.

The full proof was read and accepted by root on 2026-10-19, including
the exact criterion, adaptive-to-fixed reduction, overlap inequality and
sharpness, structural symmetry and rigidity, and minimal-counterexample
induction. The final normal and optimized control runs each passed all
1,993 exact guards. The candidate dichotomy remains explicitly OPEN.

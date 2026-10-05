# A guarded residue tree at literal 27 mod64

2026-10-05. **PROVED:** the finite residue partition, the all-height
five-step first-descent family, its new half relative to the explicitly
frozen comparison bank, and its power-three exponent phase. **FINITE-EXACT:**
the independent replays and saved-row overlap checks. **INHERITED PROVED:**
the comparison banks and proper Collatz rank. **OPEN:** universal paid
coverage and root certificates for arbitrary resulting children.

The concrete new half is

    219+512t --(1,2,1,1,3)-->209+486t,  t>=0.

For sources 3^k it is exactly k=83 mod128. This contributes another1/16
of the specified exponent domain k=3 mod8 beyond the frozen four-slot
comparison bank. Its48.9706155834282 percent subtotal is an intermediate
subtotal, not the larger concurrent minimum-exit bank's coverage.
In particular, the concurrent bank includes this row: do not add it again.

## 1. Inheritance, types, and the cheapest hostile

The closest mechanism is the switched bank in
[paid portraits](collatz_paid_portrait_controllers_20261004.md), section6,
followed by the first sharper exit and finite sharper-exit probe in
[four-slot compression](collatz_four_slot_compression_20261004.md),
sections6--7. Its q=1,r=2 factor-two threshold was4. The old finite probe
already tested threshold3; the present result gives the direct all-height
proof and exact comparison accounting for this row. There is no claim
that its numerical occurrence was previously unseen.

The canonical hostile is the source27 itself: its first five valuations
are12111, and the fifth endpoint107 still exceeds27. The corrected near
miss is the notation27 mod64. In the preceding four-slot theorem this
congruence applies to the **exponent k** of3^k; it maps to source187 mod256
and word1213. Here the initial congruence applies to the **source n** and
forces1211. Indeed 3^27=59 mod64, not27 mod64.

The least-used sidecar is the original low-bit address together with the
final valuation. Our board is **source / exact word / coarse tail /
carry / payment / retained address**. A property-preserving merge can
discard the exact valuation; retain that sidecar if exact word recovery
is required.

## 2. Four leaves, one exact local partition

Write U(n)=(3n+1)/2^v2(3n+1) on positive odd integers. For
n=27+64t, t>=0, the first four steps are exactly

    n ->41+96t ->31+72t ->47+108t ->71+162t,
         1        2         1         1.

Every displayed checkpoint strictly exceeds n. The next numerator is

    3(71+162t)+1=2(107+243t).

Consequently the next valuation and its source guard have the exact
partition

| Low-bit parameter guard | Source guard | Fifth valuation | Status here |
|---|---|---|---|
| t=0 mod2 | n=27 mod128 | 1 | Still above source at step5 |
| t=1 mod4 | n=91 mod256 | 2 | Still above source at step5 |
| t=3 mod8 | n=219 mod512 | 3 | Newly paid relative to frozen bank |
| t=7 mod8 | n=475 mod512 | at least4 | Already in the old switched bank |

Read low-order bits first. The four leaf addresses are respectively
0,10,110,111. Their parameter densities are1/2,1/4,1/8,1/8; they are
disjoint and exhaust every nonnegative parameter. Merging the final two
leaves gives the coarser guard t=3 mod4, equivalently n=219 mod256.
This is the requested three-branch fifth-reset split:1,2,or at least3.

For the first two rows the fifth multiplier is respectively243/64 or
243/128, greater than1, with positive carry. Hence those checkpoints
still grow for every permitted positive source. This says nothing about
their later steps or other paid dependency rules. For example155 is in
the first row but its inherited checkpoint child111 is already paid.

## 3. The fifth-reset tail pays at every height

The formal affine coefficients for1211 are(81,32,85). Appending a
valuation3 gives

    F_12113(n)=(243n+287)/256.

The exact word cylinder is219 mod512. Substitution gives

    F_12113(219+512t)=209+486t,
    (219+512t)-(209+486t)=10+26t>0.

The earlier checkpoints grow, so this is the first odd checkpoint below
the original source, after exactly five odd steps. The exact word has
eight halvings. Here5 is a word length,8 is its valuation cost, and9/8
is the multiplier of its inherited initial12 block. These are typed
properties of this particular arithmetic word, not interchangeable
geometric counts.

The same proof covers the whole coarse tail. Put n=219+256t. Its first
four valuations remain1211, and the fifth numerator is

    8(209+243t).

Thus the actual endpoint and fifth valuation are

    y=oddpart(209+243t),   a=3+v2(209+243t).

Since y<=209+243t<n with difference at least10+13t, every member
first descends at the fifth odd step. The inequality is strict already
at the least source219; it does not extrapolate from finite samples.
The equivalent coefficient test is(256-243)n>287, valid throughout
this cell. All earlier states exceed n>1, so no earlier root has been
silently padded. If the fifth endpoint is1, a first-hit certificate
stops there.

The inherited proper rank has first component
E(n)=3^b((n-1)/2^b)^2, b=v2(n-1), and rank(1)=(0,0).
Here n=3 mod4, so E(n)=3(n-1)^2/4. Every smaller odd y>1 satisfies
E(y)<E(n), because(3/4)^v2(y-1)<=3/4. Root1 is immediate.
Thus this is original-source payment in both numerical size and the
inherited lexicographic rank. It is a dependency reduction:
a supplied certificate for y can be prepended with the actual word.
The cell alone does not supply such a certificate for every y.

## 4. Exactly what is new in the frozen comparison

The old q=1,r=2 switched rule required a>=4, the cell475 mod512.
The new part of the coarse tail is exactly its sibling219 mod512.
It has natural density1/512 among positive integers, relative density
1/256 among odd integers, and1/8 within literal27 mod64.

Disjointness from the **all-height named dyadic bank**, not just its
saved truncation, follows from the first distinct valuation:

* Rational repeated-anchor rows with initial run longer than one, the
  -17 rows, and permuted1^h2a rows start11. The h=1 rational rows repeat
 12 at least twice and start1212.
* Switched rows(12)^q1^r a have q=1,r=2 only if their departure matches
 1211; that row requires a>=4. All other q,r depart earlier or later.
* The checkpoint155 mod2048 starts121112, so its fifth letter is1.
* The old four-slot cell187 mod256 starts1213. Other paid four-slot
  permutations cannot have first four letters1211: their multiset
  requires an exit of at least3.
* In binary16 the initial one-run here is1, followed by one reset2.
  Its two e=1 continuation tests are16 or6. Our continuation starts11,
  and every other debt depth already disagrees at the next valuation.
  The earlier reset-at-least3 rule and immediate-descent rule also fail
  their respective guards.

The program also checks congruence disjointness against all2217 saved
prior dyadic rows and verifies the binary16 continuation metadata.
These finite checks are controls for the prefix proof, not substitutes
for its all-height scope.

The strengthened original-source ternary bank can meet the new cell
on ordinary integers. Therefore1/256 is a new **dyadic** odd-relative
increment, not an increment over the fused dyadic/ternary union.
That bank excludes every power of3, so it does not subtract from the
following exponent-domain increment. Other repository methods are not
claimed exhausted.

## 5. Source cells and exponent cells have different measures

For b>=3,3 has exact order2^(b-2) modulo2^b, with image the units
congruent to1 or3 modulo8. At modulus512 the exact phase is

    3^k=219 mod512  iff  k=83 mod128.

The half already paid by the old switched bank is k=19 mod128.
Together the coarse source219 mod256 gives k=19 mod64. The preceding
four-slot addition k=27 mod64 is a separate phase.

Within the specified incoming exponent domain k=3 mod8, the new phase
has relative density8/128=1/16. Adding this to the frozen four-slot
interval gives approximately48.9706155834282 percent. The script saves
the exact rational interval and input-file hashes. This is not a
percentage of ordinary sources, not a root-coverage statistic, and not
an increment to add again to the concurrent larger threshold bank.

## 6. Lossless guards versus a paid quotient

The tree's source is a nonnegative integer parameter t; its target is
one labelled leaf and the unused higher bits. Reading one low bit at
each split is a bijection between the integers in a leaf and its
nonnegative quotient parameter. Retaining the addresses and leaf labels
recovers the partition, its original sources, and the corresponding
predicate. Leaf counts or densities alone lose their placement.

The merge of addresses110 and111 preserves the predicate
\"the fifth step pays\", but it loses the exact fifth valuation and
ordinary route length. There is an exact repair: for the merged source
n=219+256t retain v2(209+243t), recovering the fifth exponent and endpoint.
Even without storing it in advance, this sidecar can be computed from
the retained exact source; it cannot be read from the merged leaf label
alone. No arbitrary ties or cosmetic tournament are introduced.

The complementary low-reset leaves remain explicit obligations in this
small tree. Infinite refinement might discover further paid cylinders;
this finite partition gives no theorem that every branch eventually
acquires a paid label.

## 7. Reproduction and boundary controls

Run, from the repository root:

    python 04-computation/experiments/paid_guard_cover_tree_20261005.py
    python -O 04-computation/experiments/paid_guard_cover_tree_20261005.py

The [script](../../04-computation/experiments/paid_guard_cover_tree_20261005.py)
uses explicit checks that survive optimized mode. The
[saved output](paid_guard_cover_tree_20261005.out) records4096 consecutive
local parameters,1024 exact-half parameters, nine literal sources
3^(83+128j), j=0..8, and2217 saved-row disjointness controls.
The exact coefficient derivation and literal valuation replay are
independent paths. The two unpaid fifth-reset branches, source155
already paid elsewhere, the exponent/source27 mismatch, and invalid
source types are explicit hostile controls.

The concurrent original-source credit refinements are developed
separately in [paid guard budgets](paid_guard_budget_20261005.md).
Their smaller common-future children and credit funding are different
receipts from this direct fifth descent. The all-q,r sharper-exit bank
also includes this q=1,r=2 row; its own artifact carries the larger
union accounting.

# Ranking rational charts by added coverage and reusable endpoint certificates

2026-10-04 (America/Denver).

**Status:** PROVED elementary separation, density-tail and endpoint-completion
mechanisms; FINITE-EXACT certificates for the frozen 48-chart atlas and the
explicit universes below; OPEN global positive Collatz coverage. No novelty
claim is made for dyadic cylinders, periodic-word comparison, affine Collatz
composition, or inverse iteration.

The useful outcome is **47 mutually disjoint all-height additions** to the
named baseline. Their combined added natural density among all positive
integers is

    0.00033038280540161839...,

with an exact lower sum through repetition 30 and a proved omitted tail
smaller than `2.087e-59`. The two largest families are the rational charts
`1113` and `1122`. The companion endpoint work retains 432 completed routes
to 1 and a shared, fully checked suffix DAG. A general guarded completion
rule also produces infinitely many sources ending at any specified certified
odd hub not divisible by 3. These are distinct coverage and completion claims.

## Inheritance, comparison universe, and hostile

The input is precisely the 48 negative rational-only cycles retained by
[the arithmetic realization filter](denominator_arithmetic_filter_20261004.md),
at golden denominators `1..40,64,76,81,105`. The marked negative anchor is the
cycle member of smallest absolute value. Its exact positive-prefix growth
was proved there; the compiler is inherited from
[rational-anchor returns](rational_anchor_returns_20261004.md).

The entire comparison baseline is the union of these named sets:

1. All 171 rows / 65 disjoint cylinders of
   [the old swap-lift bank](reset_20260926_swaplift.md).
2. The all-height `-5` family in
   [recursive entry](entry_20260927_recursive.md).
3. The all-height [pure `-17` family](collatz_minus17_return_20261003.md).
4. Both named mixed heads `(1,2)` and `(1)` into that `-17` cycle in
   [the mixed compiler](collatz_mixed_return_compiler_20261003.md).
5. The already-added all-height rational `112` family, anchor `-19/11`,
   in [rational-anchor returns](rational_anchor_returns_20261004.md).

The comparison does not purport to subtract every previously conceivable
head, inverse family, or known convergent source. The script regenerates the
171 old rows and compares them individually with their saved table. Input
file hashes, all 48 rows, all separation certificates, and every new finite
certificate are retained in [the JSON artifact](anchor_coverage_refinement_20261004.json).

Closest mechanism: exact first-descent cylinders with the original source
retained. Canonical hostile: the weakened `112` exit sends `487` to `695`
after six legal shadow steps, still above `487`. The corrected near miss is
the lost terminal oddness bit in
[THM-4512, coefficient descent classes](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md).
The underused sidecar here is **first-descent time as a separation label**.

The concept board is a dyadic source cell, a periodic nominal word, the
original-source descent clock, a rational center, a certified endpoint, and
a retained suffix pointer. The anchor is added coverage; the niche is shared
completed certificates; the wildcard is an explicitly labelled geometric
refinement. Neither a geometric picture nor a compressed word pays a missing
source-relative bound.

## 1. Exact unions and differences, including partial overlap

Write `C(r,K)={n>0:n=r mod2^K}`, with `0<=r<2^K`. Its natural density is
`2^-K`. Two such cylinders are disjoint or one contains the other; they
intersect iff their residues agree modulo `2^min(K,L)`.

To subtract a deeper cylinder from its parent, follow the deeper residue's
binary address and retain the unused sibling at each split. This gives a
finite disjoint difference. For example,

    C(1,1) minus C(1,3) = C(3,2) union C(5,3),
    density = 1/4+1/8 = 3/8.

Discarding the entire parent merely because it intersects a removed cell
would lose this mass. The implementation supports general partial overlaps;
3,969 independent residue-set tests exhaust all ordered pairs of cylinders
of depths at most five. It also compares exact union measures with a separate
finite residue enumeration. No sampled frequencies enter the reported mass.

In the present 48-chart universe every old-bank overlap happens to be whole
containment. This is an observed and exactly checked property of these
inputs, not an assumption of the subtraction engine.

## 2. The descent clock makes infinite-family subtraction finite

For a chart word `w=(a_1,...,a_p)`, let `A=sum(w)`, `P=3^p`, `Q=2^A`, and
write its reduced negative anchor as `-h/d`. The inherited compiler accepts
every `m>=m0`, where `m0` is least with `Q^m>h`. Its cylinder uses the least
`t>=1` for which

    2^t*(Q^m-h)>P^m-h,
    beta=h*(P^m)^(-1) mod2^t,
    C_m=C((beta*Q^m-h)*d^(-1) mod2^(Am+t), Am+t).

Every source in `C_m` has exact first odd descent at time `T=pm`. Its first
`T-1` valuations are precisely the first `T-1` letters of the infinite
periodic word `w^infinity`. Only the final valuation is enlarged by the
actual exit division. The source decoder is `v2(d*n+h)=Am`, so different
repetitions of one chart are already disjoint.

**PROVED periodic separation lemma.** Consider distinct marked primitive
words `w,v`, lengths `p,q`, with compiler cylinders at repetitions `m,n`.
If their cylinders overlap, their exact first-descent times agree:

    pm=qn=T.

Put `L=lcm(p,q)`. The two nominal infinite sequences have common comparison
period L. If their first mismatch has zero-based index `j<L`, overlap
requires `T<=j+1`, because all valuations before the final step are exact.
Since T is a positive multiple of L, the only possible exception is

    j=L-1, T=L, m=L/p, n=L/q.

That isolated exception can be checked against each compiler's lower bound
and, if admissible, by exact dyadic intersection. This argument keeps the
final-valuation exception explicitly; comparing complete nominal words
without it would not be sound.

The saved finite certificate checks **all 1,128 unordered pairs** of the
48 primitive marked words. Their infinite words are distinct, and none has
an admissible exceptional time. Therefore their all-height source families
are pairwise disjoint. This is an all-height consequence of a finite word
comparison, not an extrapolation from repetition 30.

For a mixed baseline with head length H and cycle length q, compare
`w^infinity` with `head · cycle^infinity`. Comparing the first
`H+lcm(p,q)` letters either finds the first mismatch or proves the two
infinite sequences equal. An overlap additionally requires

    pm=H+qk <= j+1.

Only finitely many `(m,k)` can satisfy this after a mismatch. Across all
240 chart/baseline comparisons there is exactly one equality: the `112`
chart is the already-inherited `112` family. Every other comparison has no
admissible clock-compatible exception. Thus the other **47 whole families**
are disjoint from all five named infinite baseline families at every height.

## 3. All-height bank escape and the exact ranking

Every compiler cylinder at repetition m lies in

    d*n+h=0 mod2^(Am).

For a fixed chart, find a depth H such that its center residue
`-h*d^(-1) mod2^H` is disjoint from every old-bank row. Then all repetitions
with `Am>=H` are entirely outside the bank. This is a check against the
finite bank, followed by a proved containment implication.

For all 48 centers such a missed parent exists. The least depths in this
test range from 10 through 15. Every chart is outside the bank from
repetition three onward; most are already outside at their first admissible
repetition. The only early bank-contained cylinders are:

| Word | Anchor | Covered repetition |
|---|---:|---:|
| `1112` | `-65/49` | 2 |
| `112` | `-19/11` | 2 |

The latter whole family is independently removed as an infinite baseline.
There are no partial bank overlaps in this atlas. Consequently the new
union consists of all admissible repetitions of 46 charts plus repetitions
`m>=3` of `1112`: 47 disjoint all-height families.

The finite ranking includes **every admissible m<=30**, totaling 1,394
input cylinders and 1,364 new disjoint cylinders. For any one chart,

    2^t > (P/Q)^m,
    sum_(m>M) 2^(-Am-t_m) < 1/((P-1)*P^M).

The first inequality follows because `(P^m-h)/(Q^m-h)>(P/Q)^m`.
These are strict rational bounds. The JSON records each exact lower sum,
upper bound and uncovered cell. Their total is

    L = 64460752738239543118337511721732842656105763433642785083363134028344828679904765711412997803737476414057852160521106062898244500922191380481
        / 195109284394749514461349826862072894109287383916560696928697309976585733676235351257519131441468248197489183195087913930965498479955517831643136,

    L < delta_new < L + E,
    E = sum_(47 new charts) 1/((P_chart-1)*P_chart^30)
      < 2.087*10^(-59).

Existence of the natural density is also proved: the tail of each family
after m=M is contained in its single center cylinder of depth `A(M+1)`.
The union of finitely many such tail cylinders has upper density tending
to zero. Finite disjoint unions therefore approximate the countable union.
No general countable-additivity rule for natural density is assumed.

The strongest candidates are:

| Word | Anchor | First added cylinder | First-descent endpoint of least source | Added density, approximate |
|---|---:|---|---:|---:|
| `1113` | `-65/17` | `719 mod8192`, m=2 | 577 | 0.00012303912265610248 |
| `1122` | `-73/17` | `6983 mod8192`, m=2 | 2797 | 0.00012303912265610248 |
| `11113` | `-211/115` | retained in JSON | retained in JSON | 0.000015318627450980392 |
| `11122` | `-227/115` | retained in JSON | retained in JSON | 0.000015318627450980392 |
| `11212` | `-251/115` | retained in JSON | retained in JSON | 0.000015318627450980392 |

Both leading lower bounds exceed every other chart's all-height upper
bound, so their superiority is certified beyond the cutoff. Their finite
lower sums through 30 coincide. An all-height equality of their densities
is **not** asserted. Since they have the same P,Q and `65<73`, monotonicity
of `(P^m-h)/(Q^m-h)` in h gives `t_m(65)<=t_m(73)`, hence the `1113` density
is at least the `1122` density. Their first cylinders alone add exactly
`1/4096` of all positive integers. The respective bank-missed parents are
`719 mod2048` and `839 mod2048`.

The full table ranks exact finite added mass. Ties or overlapping infinite
intervals elsewhere are not promoted to an unproved total ordering of
infinite densities. Densities concern starting integers, not visit rates
along a Collatz orbit.

## 4. Completed finite routes and shared suffixes

The endpoint universe is fixed in advance: for each of the 48 charts take
the first three admissible repetition counts and the three source lifts
`r, r+2^K, r+2*2^K`. This is **432 distinct certificate instances**.
The inherited independent compiler verifier checks each first descent,
every exact valuation and the closed-form endpoint. A second literal
source replay stops at its first visit to 1.

All 432 endpoints are distinct and all 432 sources reach 1 in this universe.
The largest source has 53 binary digits; the longest full first-hit route
uses 232 odd steps. The declared resource caps are 10,000 odd steps and
10,000 bits; no source reached either cap. These are finite certificates,
not a statistical estimate or a coverage theorem for the source cylinders.

Each suffix DAG node stores the actual integer, its next odd integer, the
exact valuation, and the remaining first-hit distance. Edges lower that
distance by one. The only terminal is 1, with distance zero. The graph
retains **20,774 distinct suffix edges**, versus **32,027** edges when the
432 suffixes are separately expanded. Every saved edge is independently
checked against literal odd iteration. For example, 402 routes share node
5, 197 share node 17, and 176 share node 577. Sharing means identical actual
states and future routes, not equal colors, word lengths or residues.

The per-chart table also retains the sum and maximum suffix length over
its nine specified sources. This is a resource comparison on a declared
finite sample, not an expected running time for random members of a chart.

## 5. A source family can carry a reusable completed suffix

**PROVED completion rule.** Fix any compiler chart, admissible m, and
positive odd hub rho with `3` not dividing rho. Suppose rho already has
a retained finite first-hit route to 1. Select an exponent tau with

    tau>=t_m,
    d*rho*2^tau = -h modP^m.

Then define

    b=(h+d*rho*2^tau)/P^m,
    n=(b*Q^m-h)/d.

These are integers: the congruence gives b, and `P=Q mod d` gives
`b*Q^m=h mod d`. Also b is positive and odd, and `Q^m>h` makes n positive.
The dyadic guard is exactly the compiler guard, and

    b*P^m-h=d*rho*2^tau

has valuation exactly tau. Thus n has its exact first descent at pm,
directly to rho. Since every earlier iterate exceeds the retained source
and rho<n, none visits 1 early. Appending rho's certified suffix gives an
actual first-hit certificate to 1.

There is exactly one exponent class modulo `2*3^(pm-1)`: 2 generates the
units modulo `3^(pm)`. The elementary order proof is inherited from the
ternary lifting argument in rational-anchor returns: order two modulo 3,
with `v3(2^(2*3^j)-1)=j+1`. Choose the first exponent in that class at
least `t_m`; every later exponent in the class supplies another source.
The script finds the class by testing three lifts at each ternary digit,
without expanding the resulting enormous source.

For all 48 charts at their first admissible m, the script compiles the four
already-certified hubs `1,5,17,89`: **192 symbolic infinite families**.
All modular guards and exact exponent periods are checked. The 21 instances
whose least terminal exponent is at most 4,096 are additionally expanded
and replayed literally. Examples illustrate why endpoint choice matters:

| Chart, m | Hub | Least tau | Exponent period | Hub suffix length |
|---|---:|---:|---:|---:|
| `1113`, 2 | 1 | 79 | 4,374 | 0 |
| `1122`, 2 | 1 | 3,102 | 4,374 | 0 |
| `1122`, 2 | 17 | 153 | 4,374 | 3 |

For `1122`, sharing the three-step suffix `17 -> 13 -> 5 -> 1` reduces
the least tested terminal exponent from 3,102 to 153. The actual odd suffix
is checked in the DAG; neither scalar exponent size nor first-descent time
alone measures the entire certificate cost.

For fixed chart, m and hub, these completed sources grow exponentially in
the exponent-class index. Each such family has natural density zero.
The 192 sparse completions therefore must not be added to the positive
density of already-counted descent cylinders. Their benefit is a sealed
endpoint, not additional source mass.

## 6. A precise, limited refinement picture

The arithmetic refinement is

    C(r,K) = C(r,K+1) disjoint-union C(r+2^K,K+1).

An explicit geometric display can group two consecutive least-significant
address bits and assign `00,01,10,11` to the three corner triangles and
the central triangle of a midpoint subdivision, in a fixed labelled order.
Repeat the same labelled rule inside each child. Depth `2k` cylinders then
correspond bijectively to the `4^k` labelled depth-k triangles. Ancestor
containment, disjoint interiors, and normalized area `4^-k` match source
cylinder containment, disjointness, and density. Odd binary depths correspond
to a union of two children, so this is not a one-bit/one-triangle claim.

This finite address map preserves refinement and measure. It discards the
actual integer, ordered valuation word, original-source threshold and
endpoint unless those remain attached labels. At infinite depth, boundary
points can have multiple addresses and need that address sidecar. No
intrinsic geometric Collatz map, tournament orientation, or noble-polyhedron
claim follows from the display. The productive transfer is an exact way to
visualize partial cylinder subtraction while retaining the arithmetic data.

## Reproduction

    python -X utf8 -B 04-computation/experiments/anchor_coverage_refinement_20261004.py --json 05-knowledge/results/anchor_coverage_refinement_20261004.json
    python -X utf8 -B -O 04-computation/experiments/anchor_coverage_refinement_20261004.py

The [saved output](anchor_coverage_refinement_20261004.out) and JSON retain
the exact finite ranking, all separation and parent-cylinder certificates,
the partial-overlap hostile, every completed source, the suffix DAG, and
the symbolic hub families. All mathematical checks remain active under
optimized Python. The outstanding obligation is to assign every arbitrary
positive source a successful certificate; neither this ranked atlas nor its
shared suffix graph establishes that totality.

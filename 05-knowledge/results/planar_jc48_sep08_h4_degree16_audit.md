# Independent audit: degree at least sixteen for the marked mixed-cusp covering

**Status: PROVED / INDEPENDENT ANALYTIC, ENUMERATION, AND SOURCE AUDIT PASS.**
The geometry sibling independently read
[the full primary proof](planar_jc48_sep08_h4_degree16.md),
its [standalone source](../../04-computation/planar_jc48_sep08_h4_degree16.py),
and the load-bearing marked-retention suppliers. The primary was
`RESERVED` at acceptance; root owns its promotion. No producer file
was edited. No mathematical correction is required.

## 1. Domain and dependency audit

The accepted conclusion is **mapping degree `D>=16`**, under the actual
geometric hypotheses inherited in the primary from
[odd retention, Section 1](planar_jc48_sep08_odd_retention.md): a
nonautomorphic polynomial Keller map whose whole irreducible
nonproperness support has normalization `A1`, exactly two ordinary
`(2,3)` cusps, one `(2,5)` cusp, at least three ordinary nodes and
no other affine singularities. The simultaneous actual local marking
and generation hypotheses are required. This is not a theorem about
arbitrary unmarked H4 actions or about supports with additional
unlisted singularities.

The four positive meridians generate the actual transitive action and
have equal cycle type. In the primary's letters, the cusp pairs are
`(a,b),(b,c),(c,d)` with exponents `3,3,5`; three distinct actual
node pairs are `(a,c),(a,d),(b,d)`. There is a common positive actual
retained-sheet count `k`, and common moved count `t`, with
`1<=k<=D-t`. The nonautomorphic setting pays `D>=2`.

I checked the complete list of nonidentity cycle partitions with moved
size at most six. There are exactly ten:

    (2), (3), (4), (2,2), (5), (3,2),
    (6), (4,2), (3,3), (2,2,2).

The primary routes these respectively to the proved single-cycle,
arbitrary-retention involution, and mixed `(3)(2)`, `(3)(3)`,
`(4)(2)` results. The mixed `(4)(2)` primary is currently promoted
and independently audited; it is not a reserved dependency. The
single-cycle abstract equality theorem is used with the separately
paid actual generation and positive retained count. A single common
cycle moving fewer than `D` letters is not transitive. The mixed
equality results are likewise valid in arbitrary ambient degree,
rather than only the ambient size of their finite pair certificates.
For involutions the actual retained-subset consumer is the proved
one, not a full-fixed-set specialization.

An identity meridian action cannot be transitive for `D>=2`, and a
permutation cannot move one letter. These facts justify `t>=7` before
the finite degree head begins. No new small-cycle theorem is silently
assumed in the final three cases.

## 2. Actual marking and the inequalities it permits

I reread [mixed-cusp braid, Sections 4–5](planar_jc48_sep08_mixed_cusp_braid.md).
Its positive free basis `(z,y,w,x)` is the present `(a,b,c,d)`.
The cusp and node identifications use exact simultaneous conjugations
and the certified inside-pair inverse Hurwitz corrections. Those
corrections send

    (sigma,tau) -> (tau,tau^-1 sigma tau),
    (A,B) -> (B,tau^-1 A).

The corrected retained subsets are still subsets of their respective
fixed sets and still have size `k`. The actual cusp intersection,
deleted-overlap cardinality, joint fixed count, and moved-support
intersection are preserved in the manner proved by the supplier.
An ordinary-cusp injection is applied to the directly accessed pair
before transporting its count; no fresh re-access equation is
assumed from a formal basis change alone.

Thus the local actual subsets at different singularities need not be
one global choice of retained subsets for all four letters. The proof
only uses their common sizes and the paid pairwise count maps. This
distinction is particularly important in the partial-retention case.

The inherited Euler ledger and odd-retention bounds give

    1=-2k+n_ab+n_bc+n_cd+W,
    n_ab,n_bc >= max(0,k-floor(t/2),ceil((3k-D)/2)),
    n_cd >= max(0,k-floor(2t/3),ceil((5k-2D)/3)),
    W>=3 max(0,D-2k).

Here `W` includes every actual node. Extra nodes add nonnegative
terms. The exponent-five cusp uses the proved weaker exponent-five
inequality, not the false ordinary half-overlap extrapolation.

For each marked node, the actual deleted sets contain the moved
supports, so its overlap is at least the moved-support intersection.
If its meridians commute, that intersection is invariant under each
meridian and is a union of nontrivial cycles of either one. A nonempty
intersection therefore has size at least two. This assertion is
analytic; the source's elementary minimum-length control is not being
used as a substitute for it.

The two ordinary support bounds are

    |S_a intersect S_b|, |S_b intersect S_c|>=ceil(t/2).

When retention is full, `k=D-t`, every retained subset is its entire
fixed set. Cusp counts are consequently `D-2t+j` and the three
marked node counts are `j`, where `j` is the corresponding support
intersection. The six marked pairs are all six pairs of the same
four supports. This pays the full-retention incidence calculation,
including after the allowed local basis corrections.

Explicitly,

    1>=D-4t+sum_(i<j)|S_i intersect S_j|>=4t-2D.

For the second inequality, each of the `D` letters contributes
`binom(m,2)>=2m-3` when it occurs in `m` of the four supports;
their total incidence is `4t`. This holds even for `m=0`, so unused
letters need not be discarded. Integer parity now gives `D>=2t`.
The inequality direction and the rounding step are correct.

For partial retention, if `delta=(D-t)-k`, the two local retained
subsets at a given cusp can together omit at most `2delta` letters
of the joint fixed set. Hence

    n_ab>=D-2t+|S_a intersect S_b|-2delta.

This is a conservative ordinary subset bound. It is valid without
full retention or a common globally chosen retained set, and does
not depend on transporting an unproved re-access identity.

## 3. Independent complete scalar enumeration

I did not import the producer or repeat only its floor formula.
My independent enumeration reverses the outer loop order to run
over `D`, then `k`, then `t`, and explicitly checks existence of
all three individual cusp counts. It uses the equivalent integer
conditions

    2n3>=3k-D,       2(k-n3)<=t,
    3n5>=5k-2D,      3(k-n5)<=2t,

with each count in `0..k`. For every possible triple of such counts
it sets `W=1+2k-n_ab-n_bc-n_cd` and retains it exactly when
`W>=0` and `W>=3(D-2k)`.

The exact universe `D=2..16`, `k=1..D-1`, `t=7..D-k` contains
**165 integer triples**. It has **234 feasible cusp-count tuples**
before the full-retention deletion. This is a scalar universe only;
it does not assert those tuples are realized by sets or permutations.

It independently reproduces these complete raw rows below sixteen:

| D | t | k | minimum ordinary count | minimum fifth count | maximum W |
|---:|---:|---:|---:|---:|---:|
|12|7|5|2|1|6|
|13|7|6|3|2|5|
|13|8|5|1|0|9|
|14|7|7|4|3|4|
|14|8|6|2|1|8|
|15|7|7|4|3|4|
|15|7|8|5|4|3|
|15|8|7|3|2|7|
|15|9|6|2|0|9|

There are no raw rows at `D<=11`. Full retention is present in
every displayed row except `(15,7,7)`. Applying `D>=2t` leaves
exactly `(14,7,7)`, `(15,7,7)`, and `(15,7,8)`.

At degree sixteen, the raw enumeration additionally has the row
`(t,k,n3,n5,Wmax)=(9,7,3,1,8)`, which is correctly removed by
full retention. The four surviving scalar rows are precisely the
four reported by the primary. The formal Euler control is honestly
typed as a scalar survivor, not a set system or covering realization.

For reproducibility, the independent core was:

```python
from itertools import product
raw, post = {}, {}
for D in range(2, 17):
    raw[D], post[D] = [], []
    for k in range(1, D):
        for t in range(7, D-k+1):
            ordinary = [n for n in range(k+1)
                        if 2*n >= 3*k-D and 2*(k-n) <= t]
            fifth = [n for n in range(k+1)
                     if 3*n >= 5*k-2*D and 3*(k-n) <= 2*t]
            feasible = []
            for na, nb, nc in product(ordinary, ordinary, fifth):
                W = 1+2*k-na-nb-nc
                if W >= 0 and W >= 3*(D-2*k):
                    feasible.append((na, nb, nc, W))
            if not feasible:
                continue
            row = (t, k, min(a for a,b,c,w in feasible),
                   min(c for a,b,c,w in feasible),
                   max(w for a,b,c,w in feasible))
            raw[D].append(row)
            if k < D-t or 2*t <= D:
                post[D].append(row)
    raw[D].sort()
    post[D].sort()
```

The independent record `{'raw':raw,'post':post}`, serialized with
sorted keys and compact JSON separators, has SHA256
`d5f205ccd596b6a2b94c59c78e837aa3d2b9342eeece838c8117f2da19700627`.

## 4. Independent support proof for the three surviving rows

For all three rows `t=7`. The ordinary intersections each have at
least four letters in the seven-set `S_b`, so `S_a intersect S_c`
is nonempty. Commutation strengthens its size to at least two.

### Full retention at degree fourteen

Here `D=14`, `k=t=7`, and `W<=4`. If both `ad` and `bd`
support intersections were positive, they and `ac` would already
contribute at least six. Thus one of `ad,bd` is zero.

If `ad=0`, the seven-sets `S_a,S_d` are complementary in all
fourteen letters. Consequently `bd=7-ab` and `cd=7-ac`.
The Euler expression from the three cusps and marked nodes is

    -14+ab+bc+cd+ac+ad+bd=bc>=4.

If `bd=0`, the complementary pair is `S_b,S_d`; then
`ad=7-ab`, `cd=7-bc`, and the same expression is `ac>=2`.
Any extra nodes only increase these values. Both alternatives
contradict Euler value one. The complement identities were also
checked over their complete relevant integer intervals by the
producer; the proof itself uses actual complementary supports.

### Partial retention at degree fifteen

Here `D=15`, `t=k=7`. Each fixed set has eight letters, so each
local actual retained set has one omitted fixed letter. Since
`omega_ac>=2`, every node has `omega>=D-2k=1`, and `W<=4`,
the three marked overlaps must be exactly `(2,1,1)`. Additional
nodes would make this case still more restrictive.

The moved intersections `ad,bd` are at most one. Commutation
forbids positive singletons, so both are zero. Hence the seven-set
`S_d` is disjoint from `S_a union S_b`, which has size at most
eight; thus `ab>=6`. The common fixed set of `a,b` has size
`D-2t+ab>=7`. Removing at most two letters, one from each actual
retained subset, leaves `n_ab>=5`. This is the required explicit
partial-retention calculation, not the stronger false assertion
that both actual subsets contain all eight fixed letters.

Together with `n_bc>=4`, `n_cd>=3`, and the marked node total
four, the Euler expression is at least

    -14+5+4+3+4=2.

This eliminates the partially retained row.

### Full retention at degree fifteen

Here `D=15`, `t=7`, `k=8`, and `W<=3`. The marked node
counts are their actual moved intersections. Since `ac>=2`,
neither `ad` nor `bd` can be positive. Again `ab>=6` follows
from the seven-set `S_d` avoiding `S_a union S_b`.

Full retention yields `n_ab=D-2t+ab>=7`. The other lower bounds
are `n_bc>=5`, `n_cd>=4`, and `W>=2`. Thus Euler is at least

    -16+7+5+4+2=2.

All surviving rows are eliminated. This proves the necessary
mapping-degree bound, with no claim about the next degree's
realizability or about closure of all H4 cycle types.

## 5. Exact-source read and independent replay

I read the complete source, including its ten small cycle partitions,
negative-numerator-safe ceiling formulas, full declared integer
ranges, retention filter, complement identities, exact node-count
allocation checks, and the next-degree stopping controls. All checks
use explicit exceptions and remain active in optimized Python.

Both independent replays completed successfully:

```sh
python3 -B 04-computation/planar_jc48_sep08_h4_degree16.py
python3 -B -O 04-computation/planar_jc48_sep08_h4_degree16.py
```

Each passes **561 gates** and is byte-identical to the frozen output.
This agrees with, but does not replace, the independent actual-count
enumeration and analytic support arguments above.

| File at acceptance | Bytes | SHA256 |
|---|---:|---|
| Producer source | 5,256 | `40e6c9ac619658ba714aab9987e19357fd3de7450cca3d4b4e688cfe5e81edb8` |
| Frozen output and both independent replays | 696 | `b83dfcb3a3c94816840540c5e35627754d5be3f61629a944d686198fb8383623` |
| Primary before promotion | 10,878 | `15518f9346ad5e4b20e6fca68a55ea85100c93e859d715026b9db20abe7cab41` |

The producer's semantic digest is
`590353fa49ea3ff283a45804ff0d6960c767cc259e95a0883bebc6dc689d7039`.
The independent record digest above has a different specified payload
and is not being compared with this producer digest.

**Final audit: PASS.** Actual partial-retention marking, finite-head
completeness, all three support contradictions, and the retained
degree-sixteen stopping scope are accepted. Source/output were left
unchanged; root owns status promotion and Git integration.

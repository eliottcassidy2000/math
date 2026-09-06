# Leaf cofactor bounds remove the unit assumption from five decoder-entry subclasses

**Status: PROVED in the stated scope; independently audited.** The elementary tree
inequality has no LRC dependency. Entry conclusions are relative to the
explicit inherited actual-decoder, signed-box, hereditary-gcd and proper-
component phase suppliers below. LRC(14) remains **OPEN**.

## Inheritance and connection

The closest proved result is the [larger-unit entry theorem](overnight15_20260906_lrc_larger_unit.md),
which closes five unbalanced split types when the larger primitive component
contains one. The [signed-box theorem](overnight12_20260906_lrc_gcd_semigroup.md)
already replaces a unit by a maximum-endpoint gcd in the `2+11` orientation.
We recover its exact relation mechanism, then use a less-used sidecar:
**which tree cofactors correspond to leaves**. This sharpens the bound on
the smaller component's minimum and transfers the gcd replacement to all
split sizes. No claim is made that a primitive component always has a
qualifying pair.

The inherited decoder source is [THM-3818 / inert pair packet](../../01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md),
particularly its actual crossing prohibition, internal pair height and
tree cofactors. Its box supplier is [THM-2052 / finite-height relation code](../../01-canon/theorems/THM-2052-finite-height-forces-high-rank-bounded-relation-code.md).
The completed [joint-shadow sieve](lrc14_joint_shadow_empty_core_next_sep06.md)
supplies the large-subset gcd bounds, not its reserved namespace. Proper-
component LRC is a **CITED** input exactly as in the larger-unit theorem;
no new literature claim or external priority claim is made.

The canonical hostile is a primitive connected cofactor shape with large
minimum; gcd one is not a unit. The corrected near miss is treating the
largest pair label as the maximum of the full protected body. The full
maximum stays in every phase estimate. Another boundary is that an attained
real-weight cofactor bound need not be sharp after integral common-content
division or a distinct-label restriction.

The concept board is: rooted cofactors, leaf-cut proportions, primitive
minimum, endpoint gcd, signed crossing radius, and a protected physical
phase interval. The maps and losses are explicit in the proofs below.

## 1. A sharp minimum-cofactor inequality

Let T be a tree on n>=3 vertices. Each edge ij carries a relation
`a_ij w_i=b_ij w_j`, with positive real coefficients and
`a_ij+b_ij<=C`. For each vertex v let c_v be the positive rooted cofactor:
orient the tree toward v and multiply, for each edge, the coefficient of
the vertex farther from v. These are the absolute maximal minors of the
weighted incidence matrix, and form a positive kernel vector.

**Theorem 1.**

    min_v c_v <= C^(n-1) (n-2)^(n-2)/(n-1)^(n-1).       (1)

The constant is sharp for this positive-real cofactor class. If T has ell
leaves, there is the stronger, topology-dependent bound

    min_v c_v <= C^(n-1)
       *[(ell-1)^(ell-1)/ell^ell]^((n-1)/ell).           (2)

For n=2, the bound is C/2. If the edge coefficients are positive integers
and w is the primitive positive integral kernel vector, the same bounds
hold for `min w`, with the right side rounded down if desired.

**Proof.** Every edge cut contains k original leaves on one side and ell-k
on the other, with `1<=k<=ell-1`. Take the geometric mean of the cofactors
rooted at the ell leaves. The edge contributes, up to interchange,
`a^(k/ell) b^(1-k/ell)`. For `alpha=k/ell`, weighted AM-GM gives

    a^alpha b^(1-alpha)
      <= (a+b) alpha^alpha(1-alpha)^(1-alpha).

The symmetric function `phi(alpha)=alpha^alpha(1-alpha)^(1-alpha)` is
largest on `[1/ell,1-1/ell]` at its endpoints. Multiplying over the n-1
edges gives (2), since the minimum over all vertices is at most the
geometric mean over the leaves. For n>=3, `2<=ell<=n-1`; `phi(1/ell)`
increases with ell, giving (1).

For sharpness take a star with all center coefficients `C/(n-1)` and
leaf coefficients `C(n-2)/(n-1)`. The leaf cofactors are precisely the
right side of (1), while the center cofactor is n-2 times as large.
For integer coefficients, the primitive vector is the cofactor vector
divided by its positive integer content. This can only lower the minimum.
Thus cofactor sharpness is not asserted to survive primitive normalization
or distinct positive labels. QED.

The proof also retains every individual cut proportion before applying
the uniform leaf bound: their product can improve (2) for a specific tree.
This is a bound on a minimum, not on the component maximum.

## 2. Unit-free native entry

Set `Q=91^6=567869252041`. Assume an actual primitive thirteen-speed row
in the physical box `sum speeds<=Q^2` has exactly two THM-3818 decoder
components and satisfies the full relation-span equality `W_(Q,3)=V_dec`.
Write its actual components as

    tV union gU, gcd(V)=gcd(U)=gcd(t,g)=1,
    |V|=a, |U|=b, a+b=13, 1<=a<=6,
    v=min V, L=max U.

Neither shape is required to contain one. For a pair `u<L` in U put

    D=gcd(u,L), A=u/D, B=L/D,
    delta=gcd(gD,tv), c=gD/delta,
    R=Q(A+B)-(A-1)(B-1).

Actual entry supplies `B<=Q`. If `c<=Q`, the signed-box theorem implies

    tv/delta > R > QL/D.                               (3)

Indeed a point in the central box would give a nonzero mixed relation
`c(tv)=r(gu)+s(gL)` of height at most Q, outside the decoder span. This is
excluded by the actual equality, not by a guessed partition.

**Theorem 2 (native criterion).** If

    c<=Q,             a delta R >= 7(b+1)Lv,            (4)

then the physical row has a common phase of clearance strictly greater
than `1/14`.

**Proof.** Equation (3) and (4) give `t>7(b+1)L/a`. The proper-component
phase supplier gives U a phase of clearance `1/(b+1)`. Its L-Lipschitz
minimum-clearance function is strictly above `1/14` in an open interval
of full length `a/[7(b+1)L]`. Lifting any V phase of clearance `1/(a+1)`
to the t phases `(eta+j)/t` preserves every V clearance. Since gcd(t,g)=1,
their images in the U clock form a complete translated t-grid. The
protected interval is longer than its spacing, so one image is strictly
inside. This proves the full thirteen-speed consequence. QED.

The weak-boundary refinement is also exact. Put `z=gcd(delta,v)`,
`d=delta/z`, `e=v/z`. Since delta divides tv, d divides t. Therefore (3)
forces

    t >= d*(floor(R/e)+1).

If this integer lower bound is at least `ceil(7(b+1)L/a)`, the closed
protected interval supplies weak clearance `1/14`. The strict assertion
in (4) does not require this refinement.

## 3. Explicit endpoint-gcd subclasses at every unbalanced split

For hypothetical strict failure the inherited large-subset table implies

    g<=C_b,  (C_12,C_11,C_10,C_9,C_8,C_7)=(1,2,4,9,30,90).

**Corollary 3.** If some maximum-endpoint pair of U has

    D<=Q/C_b,              D min V <= aQ/[7(b+1)],       (5)

then the actual row has a weakly safe phase. Under `g<=C_b` the argument
constructs a strictly safe phase. For `g>C_b`, only the inherited weak
safety consequence is claimed.

To prove this, (5) gives `c<=gD<=Q` in the remaining branch. Equation (3)
then yields `t>QL/(Dv)>=7(b+1)L/a`, so Theorem 2's phase argument applies.

Each edge of a connected primitive V has positive primitive coefficients
whose sum is at most 356. A spanning tree and Theorem 1 give the following
uniform maximum-endpoint gcd cutoffs; they imply (5):

| Split a+b | Bound on min V | Sufficient maximum-endpoint gcd D |
|---|---:|---:|
| 1+12 | 1 | `6,240,321,451` |
| 2+11 | 177 | `76,388,115` |
| 3+10 | 31,684 | `698,294` |
| 4+9 | 6,684,150 | `4,854` |
| 5+8 | 1,694,040,507 | `26` |

The two-vertex bound uses the actual primitive atlas `p+q<=356`, which
gives `min(p,q)<=177`. For three through five vertices, the exact real
bounds before integer rounding are

    31684, 180472064/27, 1694040507.

For six vertices the cofactor bound is `1463827680198656/3125`; its floor
is `468424857663`, still above `3Q/28=60843134147`. No automatic balanced
closure follows. The exact native criterion remains available there.

The `2+11` row is the inherited endpoint-gcd theorem, recovered with its
original constant; it is not a new discovery. The other rows remove the
larger-unit requirement on the stated gcd subclasses. Merely substituting
the old `355^(a-1)` bound would give cutoffs `175558,725,2` in the last
three rows. The leaf calculation strengthens them to `698294,4854,26`.

## 4. Actual controls and boundaries

### Incoming cross-divisor strengthening

Concurrent commit `8e560f2142` supplies the audited
[cross-component divisor criterion](open_frontier_sep06_decoder.md), Sections
2--3. Its native formula overlaps Theorem 2 and also permits a selected pair
away from the maximum, while retaining the full maximum in the phase bound.
That native mechanism is therefore credited as concurrent overlapping work.
The leaf-cofactor inequality and its sharper numerical bounds are the
additional input here.

**Corollary 4.** Every numerical gcd cutoff in the table remains sufficient
after replacing D by

    H=D/gcd(D,min V).

Indeed put `v=min V` and `d=gcd(D,v)`. The divisor d divides delta, so
`c<=C_b D/d=C_b H` in the remaining gcd branch, and (3) gives
`t>Q L/(H v)`. Equivalently, the incoming conditions are
`C_b H<=Q` and `lcm(D,v)=Hv<=aQ/[7(b+1)]`. Each table cutoff pays both by
the same proved bound on min V. Thus the original D-cutoffs are valid
weaker subclasses, and the sharper cancellation form is an immediate
combination of the two independently audited mechanisms. The
[independent entry audit](second_20260906_entry_audit.md) accepts this combination.

The [balanced actual hostile](second_20260906_decoder.md) also settles the
incoming question of whether actual equality, height and all retained
profiles always force a qualifying cross-divisor score: **they do not**.
For every smaller label v and every maximum-endpoint gcd D, even allowing
all shared-factor cancellation,

    lcm(D,v) >= v >= min V > 3Q/28.

This fails the balanced score's second condition for every choice. It is
a refutation of that sufficient-score implication, with an explicit safe
physical row, not a refutation of safety or of the criterion.

### Actual controls

The standalone exact certificate gives all six actual split types
with a unit-free larger shape. The common nested larger bank is

    (6,10,15,20,24,30,60,18,2,3,40,45),

using its first b labels. Its maximum is 60, every maximum-endpoint gcd
exceeds one, and every shape is primitive and connected in the actual
atlas. The smaller shapes are the inherited connected controls of sizes
one through six. Taking `g=1`, `t=120Q+1` fits the physical box and
excludes every mixed relation by dominance. Literal checks recover the
two actual graph components, verify both mixed-support orientations, and
construct an exact full-row safe phase. Thus the new nonunit condition
has actual equality entries; no virtual partition is used as evidence.

The map from a signed coefficient box to entry discards coefficient signs
only after obtaining a forbidden relation; it then retains the two scales,
the full larger maximum and a complete translated grid. The cofactor map
forgets integer content only in the safe direction of an upper bound.
The current balanced hostile is handled separately, retaining its actual
graph, box, mixed relations and all large-subset filters.

No universal implication from primitive gcd one to a small endpoint gcd
has been proved. These are proved sufficient subclasses, not
a classification of every two-component row or a solution of LRC(14).

## 5. Reproduction and independent audit

Run from the repository root:

```bash
python3 -B 04-computation/second_20260906_entry.py
python3 -B -O 04-computation/second_20260906_entry.py
```

The [source](../../04-computation/second_20260906_entry.py) and
[frozen output](second_20260906_entry.out) give **134,477 weighted tree
cases and 884,724 always-active exact gates**. The universe is all labelled
trees on three vertices with coefficient sums at most eight, four vertices
with sums at most six, and five vertices with the five directed weight
pairs `(1,1),(1,2),(2,1),(1,3),(3,1)`. Exact star controls cover n=3,...,9.
The six actual controls retain all **1,001 mixed supports** and test the
constructed phase against all thirteen physical speeds. Normal and
optimized output bytes agree.

The [independent audit](second_20260906_entry_audit.md) accepts the entire
analytic argument and recomputes **26,800** additional weighted trees by
connected edge subsets and literal Bareiss minors, using no producer import.
It independently reconstructs the actual entries and phases, with **207,160**
gates. These finite controls validate implementations and adversarial
boundaries; Theorems 1 and 2 supply the universal statements. Raw source and
output hashes are fixed in the [session manifest](second_20260906_manifest.json).

The [balanced actual hostile](second_20260906_decoder.md) fails (4) for
every maximum-endpoint pair in its larger component, while satisfying the
actual entry, box and retained joint-shadow filters. Its explicit safe phase
shows that failure of these sufficient tests is not an unsafe row. This is
the precise stopping boundary of the current entry argument.

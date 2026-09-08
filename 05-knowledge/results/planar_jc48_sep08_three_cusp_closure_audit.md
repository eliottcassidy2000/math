# Independent audit of the all-degree marked three-cusp closure

**Status: full independent analytic/source audit PASS; normal, optimized
and frozen output agree at 43,138 always-active exact gates.** This
accepts [the complete three-cusp closure](planar_jc48_sep08_three_cusp_closure.md)
for its specified whole nonproperness support and actual local marking.
It also accepts the corollary for the declared birational polynomial
normalization class of degrees `(4,6)`. No general solution of JC(2),
exclusion of an individual component in a larger nonproperness set, or
unmarked Artin-group statement follows.

## 1. Precise entry and status of the suppliers

The hypothetical map is a nonautomorphic complex planar polynomial
Keller map. Its entire nonproperness set is one irreducible curve with
normalization `A1`, exactly three ordinary cusps and `N>=3` ordinary
nodes, with no other affine singularities. Four positive meridians
`a,b,c,d` generate the affine complement group and satisfy the marked
Artin D4 relations. With `e=bcb^-1`, `f=e^-1ae`, the actual local
pairs are, up to the paid count transport,

    cusps: (b,c), (b,d), (a,e);
    three distinct nodes: (f,b), (a,d), (c,d).

The distinction between a marked quotient and an unmarked presentation
is essential. I read the final primary against the proved
[passport](planar_jc48_sep08_three_cusp_passport.md), the
[representative braid and local-cluster supplier](planar_jc48_sep08_three_cusp_braid.md),
the [all-retention involution exclusion](planar_jc48_sep08_d4_retention.md)
and the [single-cycle theorem](planar_jc48_sep08_cycle_support.md).
The latter two currently have PROVED / FINITE-EXACT / INDEPENDENTLY
AUDITED status. Only their small moved-support cases are used here.
They retain the classical geometric mapping-degree-four exclusion;
its indispensable role is exposed by the S4 control below.

For the geometric corollary, I also recovered the precise statement and
marked transport in
[the complete three-cusp family](planar_jc48_sep08_three_cusp_family.md),
particularly its Sections 1 and 6 and its independent audit. That proved
supplier covers both the infinity-`(2,7)` stratum with four nodes and
the infinity-`(2,9)` stratum with three nodes. It transports four actual
positive meridians and the cusp/node neighborhoods and common accesses,
not a hypothetical Keller source. The new theorem's `N>=3` form
therefore applies to both strata. The older degree floor in that frozen
family theorem is not itself being used as an all-degree conclusion;
the new incidence argument supplies that additional consequence.

## 2. Actual retained sheets and the exact Euler identity

Let `D>1` be mapping degree and `k` the actual affine fibre count
on the whole smooth stratum. The inherited suppliers give `1<=k<D`.
The primary's direct positivity explanation is valid: the nonconstant
pullback of an equation of the curve has a nonempty divisor, and
quasi-finiteness prevents any divisor component from contracting to a
point. Thus a dense open subset of the irreducible curve is in the
image. Constancy on the smooth stratum is the separate inherited page
statement.

The off-curve cover is transitive. Conjugacy of positive meridians gives
a common moved-support size `t`. Since the four meridians generate,
`t=0` would make the transitive cover trivial; a permutation cannot move
one label. Therefore `t>=2`. Define

    a0=D-t, delta=a0-k>=0, D=t+k+delta.

Here `a0` is the full fixed count, while `k` is the actual retained
count. The proof preserves their difference throughout.

At a cusp the actual fibre count is `n_i=|A intersect B|` for its
correctly accessed retained sets. At a node the two deleted sets overlap
in `omega`, with actual fibre count `2k-D+omega`. Each `omega`
is nonnegative and at least `D-2k`. With `W` the sum over all `N`
nodes, the exact Euler calculation is

    chi(C)=1-N, chi(C_smooth)=-2N-2,
    chi(A2 minus C)=N,
    1=ND+(-2N-2)k+sum n_i+N(2k-D)+W
     =-2k+sum n_i+W.

This uses actual singular-fibre specialization and Euler integration,
not an arbitrary permutation count in place of the fibres. All node
contributions remain present. No properness of the original Keller map
is assumed.

## 3. Original pairs and the support incidence coordinates

The exact braid identities are

    e b e^-1=c, c e c^-1=b.

Consequently simultaneous conjugation of node0 `(f,b)` by `e` gives
`(a,c)`, and simultaneous conjugation of the third cusp `(a,e)` by
`c` gives `(a,b)` because `a` and `c` commute. The other formal
node pairs are `(a,d),(c,d)`. Thus the three node support intersections
are exactly the three intersections among the leaf supports of `a,c,d`,
and the three cusp joint fixed counts correspond to the central-leaf
pairs involving `b`.

These are identities in the same global marked group. Separate
simultaneous conjugations preserve the relevant scalar counts and
produce the stated common leaf supports; they are not arbitrary choices
of differently rebased leaf permutations. The directly certified node2
pair is `(e,bdb^-1)` and its already paid conjugation by `b^-1`
gives `(c,d)`.

The two inside-pair inverse Hurwitz changes require the additional
local-cluster proof already audited in
[the braid audit](planar_jc48_sep08_three_cusp_braid_audit.md). They act
on actual retained data as

    (sigma,tau) -> (tau,tau^-1 sigma tau),
    (A,B) -> (B,tau^-1 A), B subset Fix(tau).

They preserve retained intersection size, deleted-overlap size and the
joint generated subgroup, hence its full fixed set. The cusp inequality
is derived at the actual pair before transporting these counts. No
unsupported reaccess equation is inferred after an abstract Hurwitz
change. At every node a moved label is necessarily deleted, so the
actual `omega` dominates the corresponding support intersection.
Full retention is not needed for any of these comparisons.

## 4. The arbitrary-cycle cusp lemma

For braided permutations `sigma,tau`, put `g=sigma*tau`. Their
braid relation gives `g sigma g^-1=tau`, hence
`g Fix(sigma)=Fix(tau)`.

More generally let `A` be pointwise sigma-fixed and let `B=gA` be
pointwise tau-fixed, with both sets of size `k` and intersection size
`n`. The primary proves `tau(A minus B)` lies outside `A union B`.
Its argument is correct. An independent short proof uses the opposite
injection: `sigma(B minus A)` lies outside that union. It cannot meet
`A`, because sigma fixes `A`; and for `x in B minus A`,
`sigma(x)=g(x)`, so if `sigma(x)` belonged to `B=gA`, injectivity of
`g` would give `x in A`. Thus

    k-n <= D-2k+n, or 2n>=3k-D.

This proof is valid for empty sets as well and does not assume a
particular cycle type. Applying it to the full fixed sets, whose
cardinalities are `D-t` and whose intersection has size `D-2t+j`,
gives `j>=ceil(t/2)` for the two moved supports. This purely group
theoretic use of full fixed sets does not replace actual retention.

For an actual cusp, every retained point in the full joint fixed set
is fixed by `g`, so is also in the other retained set. Hence

    n>=k-(D-t)+(D-2t+j)=k-t+j
      >=k-floor(t/2).

The actual-pair transport of Section 3 makes this available at all
three cusps. A negative right side is allowed: it is a valid weaker
lower bound, and never treated as an actual negative cardinality.

Write `F=floor(t/2)`, `C=ceil(t/2)` for the two scalar halves, so
`t=F+C`. Their letters in these formulas denote numbers, not the
geometric curve or the meridian named c. The Euler identity gives

    W<=3F+1-k.

## 5. All three lower bounds are needed and correctly directed

Let `L1,L2,L3` be the moved supports of `a,c,d`, and let `B0`
be that of `b`; all have size `t`. Define
`P=sum_(i<j)|Li intersect Lj|`. The three marked nodes give `W>=P`.

Exact inclusion-exclusion gives

    P=3t-|L1 union L2 union L3|+|L1 intersect L2 intersect L3|
     >=3t-D.

Inside the central support put `Ti=Li intersect B0`. The cusp lemma
gives `|Ti|>=C`, while their union has at most `t` elements. Applying
the same identity there gives

    P>=sum_(i<j)|Ti intersect Tj|>=3C-t.

Finally the actual node bound gives `W>=3(D-2k)`. For `D-2k>=0`
this follows from at least three nodes; for `D-2k<0`, nonnegativity of
`W` is already stronger. Thus the proof does not multiply a negative
bound by an inequality for `N` in the wrong direction.

Comparing these three inequalities separately with the upper bound on
`W` yields exactly

    delta>=2t-3F-1,
    2k>=3C+3delta-1,
    k<=1+4F-2C.

The ambient leaf union, central incidence and actual retained count are
separate inputs. Omitting one of them would not prove the same reduction.

## 6. Unbounded parity proof and the last equality case

For odd `t=2m+1`, the inequalities give
`delta>=m+1`, `k>=3m+3`, `k<=2m-1`, impossible. For even
`t=2m`, they give `delta>=m-1`, `k>=3m-2`, `k<=2m+1`,
so `m<=3`. This is an unrestricted integer argument, with no ambient
mapping-degree cutoff or experimental extrapolation.

At moved count two the meridian is a transposition. At moved count four
its nontrivial cycle decomposition is precisely a four-cycle or two
transpositions. The currently proved all-retention involution and
single-cycle theorems have the same actual marking and `N>=3`
hypotheses, and exclude these cases. Their classical degree-four
geometric input remains necessary; the existence of the S4 group
control prevents omitting it.

For `t=6`, the inequalities force `delta=2,k=7,D=15`:
`delta>=2`, `k<=7` and `2k>=8+3delta` leave no freedom.
The upper bound is `W<=3`; the actual node bound is
`W>=N>=3`. Thus `N=3`, `W=3`, and every node has `omega=1`.
The cusp bounds and Euler identity also force all three `n_i=4`.

For commuting permutations, their moved-support intersection is
invariant under either permutation. It contains no fixed point of that
permutation, so it is a union of nontrivial cycles. Its cardinality
cannot equal one. Each marked support intersection is at most its
actual node overlap one; hence every such intersection is zero. Through
the exact pair map these are all three leaf intersections, giving
`P=0`. But central incidence requires `P>=3C-t=3`. This is the
final contradiction. It applies to every cycle decomposition with
moved count six; no classification of those decompositions is needed.

The complete post-reduction integer relaxation is independently
recovered as

    (t,delta,k,D) =
    (2,0,1,3), (2,0,2,4), (2,0,3,5), (2,1,3,6),
    (4,1,4,9), (4,1,5,10), (6,2,7,15).

These seven rows are not alleged polynomial passports. The analytic
parity argument proves why their finite enumeration is complete.

## 7. Independent finite-control audit

I read the full standard-library source. Its explicit checks remain
active under optimized Python. Its formal affine coefficient identities
match the comparisons above. Its bounds for the small integer universe
are derived from the analytic reduction and recover exactly the seven
rows, including their equal lower and upper W values.

The set-control universe is all triples, unordered with repetitions,
of t-element leaves in D labels, for `D=2,...,7`, `2<=t<D`, each
meeting a fixed central t-set in at least `ceil(t/2)` points.
These inequalities are symmetric in the leaves, so this representative
choice loses no tested configuration. Repetitions are expressly kept.
I independently counted the eligible leaves by

    M(D,t)=sum_(j>=ceil(t/2)) binom(t,j)*binom(D-t,t-j).

There are `binom(M(D,t)+2,3)` triples for each `(D,t)`; the total
is exactly **9,441**, matching the producer. Both inclusion-exclusion
identities are proved analytically, so this bank is a check of their
implementation and equality boundaries, not a proof outside its range.

The permutation controls cover `D=2,...,6`, one canonical sigma for
each cycle partition and every tau of the same type. Every same-type
ordered pair is represented up to simultaneous relabeling; no claim of
a free orbit action or unique representative is needed. Identity,
empty-fixed-set and mixed-cycle cases are retained. In a separate
implementation, I recovered all counts:

| Control | Count |
|---|---:|
| Same-type permutation pairs | 872 |
| Braided pairs | 152 |
| Subsets of their full fixed sets | 673 |
| Commuting pairs | 74 |

The source checks the full-fixed conjugator, arbitrary retained-set
injection, deficit bound, half-support bound, valid inside-pair count
transport and absence of singleton commuting intersections. The
same-type restriction is appropriate for actual positive meridians.
Separately, I also checked all eleven canonical sigma types on six
letters against all 720 tau permutations, without a same-type filter:
901 commuting pairs, 97 braided pairs and 460 retained-subset cases.
Those independent checks agree with the general injection and commuting
lemmas; they are supplementary finite controls only.

Both principal hostiles are correctly scoped. The degree-fifteen four
support sets give half-overlaps three, singleton leaf intersections,
full ambient union and a formal Euler value one. They do not define
commuting permutations or actual reaccess data. The missing commuting
cycle invariant is exactly why that numerical survivor fails. Conversely
the marked transitive S4 tuple is an actual finite group action with
full retained/fixed sets, cusp counts `(1,1,1)`, node overlaps
`(0,2,0)` and Euler one. Its failure to be a Keller realization uses
the separately cited geometric degree-four theorem, not the numerical
bounds alone.

## 8. Final replay pins and accepted consequence

I independently ran

    python3 -B 04-computation/planar_jc48_sep08_three_cusp_closure.py
    python3 -B -O 04-computation/planar_jc48_sep08_three_cusp_closure.py

Both completed with **43,138 always-active gates** and matched all
764 bytes of the frozen output. Fresh source and output hashes match:

| Artifact | Bytes | SHA256 |
|---|---:|---|
| [Source](../../04-computation/planar_jc48_sep08_three_cusp_closure.py) | 11,027 | `3b7a0978e0ecfa11cdf5c6d44cb66fb66109ffff965361bf0b5177061a036423` |
| [Output](planar_jc48_sep08_three_cusp_closure.out) | 764 | `bca03a5fd8697bb66984957caea989f299bdfd1b6eb4e4b36c227dd25c849377` |

Semantic digest:
`701ae302f997498f7168d77fda33452e5d73fc35b3f555b9ca5f946e4d392710`.

No mathematical, source or witness correction was needed. The complete
primary proof is accepted for owner-controlled promotion. It rules out
all mapping degrees for the stated whole marked-D4 three-cusp support;
the proved family and braid suppliers pay those hypotheses for the
specified ordinary three-cusp `(4,6)` normalization class. The argument
does not promote that target-complement condition to all Keller curves,
transport an actual source across a family, or exclude one component
inside a larger nonproperness locus.

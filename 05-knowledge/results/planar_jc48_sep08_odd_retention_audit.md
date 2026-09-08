# Independent audit: odd-cusp retention and the marked H4 degree floor

**Audit status: PASS — complete analytic/source review and independent
normal, optimized and frozen-output comparison.** September 8, 2026.
This audit accepts the statements and scopes below; root owns primary
status promotion. No primary source, output or inherited artifact was
modified during this audit.

Primary: [odd retention](planar_jc48_sep08_odd_retention.md).
Source: [exact controls](../../04-computation/planar_jc48_sep08_odd_retention.py).
Output: [frozen replay](planar_jc48_sep08_odd_retention.out).

## 1. Accepted statements and their boundaries

For every integer r>=1 and every finite permutation pair satisfying the
alternating relation of length 2r+1, with the rightmost-first convention,
put g=(sigma tau)^r. For **any** A contained in Fix(sigma), including the
empty and full sets, put B=gA. If D is the ambient cardinality, k=|A| and
n=|A intersect B|, the accepted conclusions are

    (r+1)n >= (2r+1)k-rD,
    n >= k-floor(rt/(r+1)),

where t is the common number of moved letters. The first inequality is
sharp for every r. The stronger bound retains only one deficit from the
full fixed set. Neither statement assumes transitivity or a geometric
realization.

The accepted geometric consequence is a **necessary degree floor D>=12**
for a nonautomorphic complex polynomial Keller map with the complete
whole-support package specified in the primary: irreducible nonproperness
curve, normalization A1, exactly two ordinary cusps, one (2,5) cusp and
N>=3 ordinary nodes, no other affine singularities, and the actual common
positive-meridian marking of the Artin H4 chain. Its three marked node
pairs must be distinct actual nodes. Additional nodes contribute to the
same Euler ledger with nonnegative overlaps.

The literal mixed-cusp curve and its separately paid intrinsic good locus
satisfy this marking. The conclusion does not apply to an arbitrary curve
merely because it has the same singularity list. It does not exclude all
H4 monodromy, assert realizability of the residual numerical rows, or solve
JC(2). In particular the ordinary half-support bound is not used at the
five-cusp.

## 2. Independent check of the local proof

Write U=A\B and O=Omega\(A union B). The odd relation says
g sigma=tau g, so B is fixed pointwise by tau. A point of A fixed by tau
is jointly fixed and hence fixed by g; it is therefore in B. Thus U has
no tau-fixed point. Every moved tau-cycle that meets U avoids B and is
contained in U union O.

Suppose x,tau x,...,tau^r x all belong to U. Each is sigma-fixed.
Evaluating g sigma=tau g at x, in the declared composition convention,
gives tau^r x=tau^(r+1)x, contradicting that the cycle is moved. Repeated
positions are permitted in this argument. In particular it also excludes
a short moved cycle lying entirely in U; a cycle-length assumption has
not been smuggled into the separator argument.

Every such cycle consequently has an O label after every run of at most
r U labels. Counting cyclic runs gives |U|<=r|O|, with O labels on unused
cycles only enlarging the right side. Substituting k-n and D-2k+n proves
the first inequality. This is an unbounded proof, not an inference from
the finite permutation controls.

For the stronger bound put J=Fix(sigma) intersect Fix(tau). Both directions
of the identity

    A intersect B = A intersect J

hold: an intersection point is fixed by both meridians, while a jointly
fixed point of A is fixed by g and hence retained in B. Conjugacy gives
equal fixed-set sizes a0=D-t. Applying the already proved inequality to
the full fixed sets gives the moved-support overlap bound
j>=ceil(t/(r+1)). Since |J|=D-2t+j and A misses only delta=a0-k labels
from Fix(sigma),

    n >= |J|-delta = k-t+j
      >= k-floor(rt/(r+1)).

The intermediate inequality with the actual j is preserved and is
potentially stronger. Subtracting two independent deficits would discard
the same-point re-access information unnecessarily. The displayed floor
bound is at least as strong as the first inequality because t<=D-k.

The sharpness family uses two (r+1)-cycles with a single common label,
sigma=(o,b1,...,br), tau=(o,a1,...,ar). Their product has length 2r+1;
its r-th power maps the a-block to the b-block and conjugates sigma to
tau. With D=2r+1,k=r,n=0 the first inequality is equality. Disjoint copies
and jointly retained fixed labels preserve equality, with all quantifiers
and empty/full-set boundaries consistent.

The named five-letter control sigma=(1 2 3), tau=(3 4 5), A={4,5},
B={1,2} has g=(sigma tau)^2 and n=0. Its support intersection is one.
It attains 3n=5k-2D and disproves both the proposed 4n>=5k-D and an
ordinary half-overlap assumption at exponent five. It is not a global
H4 action or a Keller map.

## 3. Actual retained pages, pair changes and the Euler ledger

The primary inherits the actual left-action access convention from
[odd-cusp passports](planar_jc48_sep06_odd_cusp.md), using the normal
finite extension and pure deleted boundary in the
[ordinary cusp supplier](planar_jc48_sep06_cusp_passport.md) and
[THM-3578, Zariski-main boundary rank and sheet debt](../../01-canon/theorems/THM-3578-zariski-main-boundary-rank-and-sheet-debt.md).
Those suppliers identify the actual re-accessed sets as B=gA and the
number of actual inverses over the cusp as |A intersect B|. The primary
does not replace these sets by all inertia-fixed labels.

The local Shimada marked-meridian presentation is inherited through those
already audited suppliers. This audit makes no new literature-priority
claim and does not treat the source's finite bank as a replacement for
that topological input.

I checked the marking against the previously independently audited
[mixed-cusp braid supplier](planar_jc48_sep08_mixed_cusp_braid.md) and its
[audit](planar_jc48_sep08_mixed_cusp_braid_audit.md). The letters
(z,y,w,x) become (a,b,c,d); the resulting cusp pairs are (a,b),(b,c),(c,d)
with exponents 3,3,5, and the three marked node pairs are
(a,c),(a,d),(b,d). The actual five-cusp pair is directly accessed. For a
permitted inside-pair inverse Hurwitz change, actual sets transform as
(A,B) to (B,tau^-1 A), with B fixed by tau. Their intersection cardinality
is unchanged. The generated local subgroup, joint fixed set, and moved
support intersection are unchanged as well. Therefore one applies the
local theorem to the paid direct access, then transports these counts.
No unproved re-access formula is inferred after an arbitrary free-group
basis change.

Whole irreducibility gives a common actual smooth-support retained size
k>=1 and a common positive-meridian cycle type, hence a common moved
count t with k<=D-t. Euler integration over the normalization-A1 curve
gives

    1=-2k+n_ab+n_bc+n_cd+W,
    omega >= max(0,D-2k)

at every node. The (2,5) cusp is unibranch and therefore has the same
Euler-stratum coefficient as an ordinary cusp; only its local group
inequality changes. The source applies the correct r=2 inequality there.

At each marked node the moved support is contained in the corresponding
deleted sheet set, so the actual overlap omega dominates the moved
support intersection. Commutation makes this intersection an invariant
union of nontrivial cycles; it cannot have cardinality one. For t=5 the
only remaining cycle type after the proved small-type exclusions is
(3)(2), so every positive marked node support intersection is at least
two. Two ordinary support overlaps of at least three inside the five
letters of S_b force S_a intersect S_c to be nonempty, and hence to have
size at least two.

With full retention k=D-t, cusp counts equal D-2t+j and marked node
overlaps equal j. All six pairs among the four supports occur exactly
once among the three marked cusps and three marked nodes. Extra nodes
are nonnegative, so

    1 >= D-4t + sum_{i<j}|S_i intersect S_j|.

For every possible incidence multiplicity m=0,...,4,
binom(m,2)>=2m-3. Summing over all D labels yields
sum intersections>=8t-3D and therefore 1>=4t-2D.
Because 4t-2D is even, this implies D>=2t. This argument does not
require transitivity or exclude unused labels by assumption. It is
valid only with the declared full-retention hypothesis; it is not
applied to the later partial-retention stopping example.

## 4. Complete small-degree reduction

The current primary statements of
[H4 involutions](planar_jc48_sep08_h4_involution.md) and
[H4 single cycles](planar_jc48_sep08_h4_single_cycle.md) are PROVED and
independently audited. They apply to the paid actual generating tuple
with arbitrary retention. Involutions and single nontrivial cycles exhaust
all nonidentity cycle types moving at most four letters. The identity
action cannot be a nonautomorphic transitive cover. Thus t>=5, and at
t=5 the sole remaining type is (3)(2). This reduction has not been
deduced from a bounded symmetric-group search.

An independent direct integer reconstruction, described in Section 6,
recovers exactly these scalar survivors, before support refinements:

| D | Pairs (t,k) surviving cusp and node cardinalities |
|---|---|
| 2 through 8 | none |
| 9 | (5,4) |
| 10 | (5,5), (6,4) |
| 11 | (5,5), (5,6), (6,5) |
| 12 | (5,6), (5,7), (6,5), (6,6), (7,5) |

I also checked the primary's analytic derivation of this finite head.
For D<=7, 1>=3D-8k>=40-5D>=5 is impossible. At D=8 the only value
not already excluded by the node term alone is t=5,k=3; two ordinary
counts at least one and W>=6 make Euler at least two. The remaining
support eliminations are as follows.

* D=9, (t,k)=(5,4): full retention violates D>=2t directly. Separately,
  the three commuting marked support pairs all overlap, hence W>=6,
  while the cusp bounds give W<=4.
* D=10, (6,4): full retention gives 1>=4. For (5,5), W<=3 and
  omega_ac>=2 force the other marked support intersections to be zero.
  A five-element S_d disjoint from S_a union S_b in ten labels forces
  S_a=S_b. Then n_ab=5,n_bc>=3,n_cd>=2 and omega_ac>=3 make Euler
  at least three.
* D=11, (6,5): full retention violates D>=2t. For (5,5), every node
  has overlap at least one and omega_ac>=2, giving W>=4 against W<=3.
  For (5,6), full retention gives W<=2. The ac node uses at least two,
  so S_d is disjoint from S_a union S_b. Their sizes force
  |S_a intersect S_b|>=4. Consequently n_ab>=5,n_bc>=4,n_cd>=3 and
  W>=2 make Euler at least two.

These exhaust all D<=11. For D=12, the additional rows (5,7) and (7,5)
are removed respectively by omega_ac>=2 versus W<=1, and full-retention
D>=2t. The three remaining numerical rows in the primary are correctly
labelled necessary data only. No existence or sharpness assertion for the
geometric degree floor is made.

The multiplication by the number N of nodes is used safely: every node
has lower bound max(0,D-2k), and N>=3. Neither the proof nor the source
multiplies a negative lower bound by N>=3 to obtain a stronger estimate.

## 5. The unbounded stopping object retains exactly the claimed data

I independently summed all seven membership cells. Their total is 33q;
each of the four support sizes is 13q; edge intersections are
(7q,7q,5q), and each marked node support intersection is q. With
k=16q the retained deficit is 4q>0. The declared cusp counts
(10q,10q,8q) equal k-t+j on their respective edges and satisfy both
odd-retention inequalities for every integer q>=2, including the floor
rounding. The node counts (2q+1,q,q) dominate both their support
intersections and D-2k=q. None of the node support intersections is a
singleton, and

    -2(16q)+(10q+10q+8q)+(2q+1+q+q)=1.

Thus the formal family is unbounded and satisfies precisely the necessary
cardinality inequalities claimed. Full retention is false, so the
conditional incidence inequality is unavailable. The example supplies
neither permutations, simultaneous retained subsets, actual re-access
maps, transitivity, nor a Keller map. The primary explicitly identifies
the missing common cycle action. It is a valid stopping witness for the
stated relaxation and not a counterexample to a geometric theorem.

## 6. Independent source and reproduction audit

The producer is standalone Python standard library code. Its `check`
function raises an exception and remains active under `python3 -O`;
there are no removable `assert` gates or imported mathematical engines.
I read the complete source and checked its rightmost-first composition,
odd word construction, exact conjugator, subset masks, cycle-run test,
rounding for negative ceiling numerators, and explicit source universes.

The main finite universe consists of every ordered pair in S_D for
D=1,...,5 and r=1,...,4, filtered only by the literal odd-braid relation,
then every subset of Fix(sigma). Empty, partial and full retention are
included. The 60,068 pair trials equal
4 sum_{D=1}^5 (D!)^2. The producer records 2,328 braided pairs, 6,096
subsets, 2,192 proper nonempty retained subsets and 230 direct-bound
equalities. Every relevant conjugator and inequality check is active.

I reconstructed that entire universe by a separate conjugacy-weighted
path, without importing or executing producer definitions. For each
integer partition of D I built one sigma with that cycle partition and
enumerated every literal tau. A sigma class has weight

    D! / product_l (l^(m_l) m_l!),

where m_l is the multiplicity of length l. Simultaneous conjugation
preserves the braid relation and bijects all the retained subsets, so
this weighting recovers the complete ordered-pair universe exactly.
Words were evaluated directly on each label, rather than through the
producer's permutation multiplication routine. Independent cycle walks
checked actual separators and maximum cyclic U-run lengths.

The independent route used 72 partition/r rows, 3,932 literal pair checks
and 673 literal subset checks, and reproduced all five weighted totals
above. It independently verified B subset Fix(tau),
A intersect B=A intersect J, the cyclic-run bound and both retained
inequalities on that exact universe.

For the scalar head I independently enumerated D=2,...,12, every
5<=t<D, every 1<=k<=D-t, and every triple of actual integers
(n_ab,n_bc,n_cd) in {0,...,k}^3. I imposed the three direct odd inequalities,
the one-deficit bounds and
W=1+2k-n_ab-n_bc-n_cd>=3 max(0,D-2k), rather than copying the producer's
precomputed minimum/ceiling predicate. This recovered the full table in
Section 4. The temporary independent scripts and outputs are diagnostic;
no statement depends on their persistence instead of the analytic proof.

The remaining producer controls are correctly scoped: 36 sharp-family
examples for r=1,...,12; the literal five-letter hostile; all 16 support
membership masks; the complete scalar head; the listed exact strict
small-degree comparisons; and five values of the analytic unbounded
set family. These validate the implementation without substituting a
finite census for the all-r theorem.

Independent commands, from the worktree root:

```bash
python3 04-computation/planar_jc48_sep08_odd_retention.py
python3 -O 04-computation/planar_jc48_sep08_odd_retention.py
```

Both completed successfully. Their **34,983-gate**, **1,358-byte** outputs
are byte-identical to the frozen primary output. The semantic control
hash is `ada3df141e8d5fa747aa178a14ff36ddbd8f88d42aa43bbd93cfd8e9002cbaf0`.

Accepted pins:

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Source | 8,207 | `c5d2865731c5b4866a78298ef2760fe9ff0acf85a58442386a78fc87ea41cfff` |
| Frozen output and each independent replay | 1,358 | `6a8d3d36ca5d2bf86874b78a0df93fe33a3a97c9b082ea98b3e3cfb2cb9dea0a` |
| Primary at final audit, before status promotion | 17,943 | `7066b810e1dc0b3a5eaaa90d76a17a99d0047acc319d33a132b4d4f5e80a413e` |

No mathematical correction remains. This audit is frozen for root's
promotion and checkpoint, with the geometric and nonrealization limits
in Section 1 preserved.

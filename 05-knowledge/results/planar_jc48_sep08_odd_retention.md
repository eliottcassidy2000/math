# Sharp odd-cusp retention bounds and a marked H4 degree floor

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** September 8, 2026.

## 1. Statements, inheritance and scope

Let sigma and tau be permutations of a finite set Omega of cardinality D,
with rightmost-first composition, satisfying the alternating braid relation
of length `2r+1`, where r>=1. Put `g=(sigma tau)^r`. Suppose

    A subset Fix(sigma),   B=gA subset Fix(tau),
    |A|=|B|=k,             n=|A intersect B|.

Empty and full retained subsets are allowed in this abstract statement.
Then

    (r+1)n >= (2r+1)k-rD.                              (1)

The coefficients in (1) are attained for every r. If each meridian moves t
letters, there is also the useful stronger bound

    n >= k-floor(rt/(r+1)).                            (2)

In particular the correct five-cusp bound is `3n>=5k-2D`, and (2) gives
`n>=k-floor(2t/3)`. The former proposed `4n>=5k-D` is false.

The geometric consequence here is deliberately narrower. Suppose a
nonautomorphic complex polynomial Keller map has whole irreducible
nonproperness set with normalization A1, exactly two ordinary (2,3) cusps,
one (2,5) cusp, N>=3 ordinary nodes, and no other affine singularities.
Assume its actual marked meridians generate a quotient of the Artin H4
presentation, with the three actual cusp pairs and three distinct actual
node pairs simultaneously marked as

    ordinary cusps (a,b),(b,c);  five-cusp (c,d);
    nodes (a,c),(a,d),(b,d),
    aba=bab, bcb=cbc, cdcdc=dcdcd,
    [a,c]=[a,d]=[b,d]=1.                               (3)

Then its mapping degree satisfies

    D>=12.                                            (4)

The audited [actual mixed-cusp braid supplier](planar_jc48_sep08_mixed_cusp_braid.md)
pays precisely this marking for the literal curve

    U=t^4-(8/3)t^3-2t^2+8t,
    V=t^6-(3/2)t^4-16t^3+48t.

Its letters `(z,y,w,x)` become `(a,b,c,d)` in (3). The separately proved
[boundary good-locus theorem](planar_jc48_sep08_three_cusp_boundary.md)
transports its actual local neighborhoods and marked positive meridians
through the connected intrinsic good (5,3,3) boundary family. Thus (4)
applies throughout that paid family. None of these statements asserts that
an abstract permutation action realizes a Keller map.

The closest mechanism is the proved general odd-braid run lemma in
[H4 single-cycle support, §2](planar_jc48_sep08_h4_single_cycle.md).
The actual access convention is inherited from
[the odd-cusp passport, §2](planar_jc48_sep06_odd_cusp.md), which in turn uses
[the ordinary actual-page/purity supplier](planar_jc48_sep06_cusp_passport.md)
and THM-3578,
[Zariski-main boundary rank and sheet debt](../../01-canon/theorems/THM-3578-zariski-main-boundary-rank-and-sheet-debt.md).
Its local marked-meridian source is Shimada's
[*Lectures on Zariski van-Kampen theorem*](https://www.math.sci.hiroshima-u.ac.jp/shimada/LectureNotes/LNZV.pdf),
§6, pp.24–25, already read and audited in the inherited supplier. We do
not claim a new literature classification.

The hostile is the five-letter braid-five pair below. The repaired near
miss is the false ordinary half-overlap bound at exponent five. The least
used sidecar is that a jointly fixed letter is either retained in both
re-accessed subsets or in neither: there is only one retained deficit.
The concept board is cycle runs, actual joint retention, marked H4 support
incidence, the exact Euler budget, and the lost cycle-order obstruction.
Targeted exact/synonym searches of current canon, result notes and the
mistakes ledger found the new bound as a current board proposal, not an
already proved general-retention theorem. This is a local antecedent
check, not an external priority assertion. Root proposed (1); this lane
independently proves, tests and develops its consequences.

## 2. The arbitrary-retention run proof

Write `I=A intersect B`, `U=A\B`, and
`O=Omega\(A union B)`. A point of A fixed by tau is fixed by both generators,
hence by g. Such a point lies in gA=B. Therefore U has no tau-fixed point.
Every nontrivial tau-cycle meeting U avoids B, because B is pointwise
fixed by tau; that entire cycle lies in U union O.

No r+1 consecutive positions on such a cycle can all belong to U. If
`x,tau x,...,tau^r x` were all in U, sigma would fix all of them. Apply

    (sigma tau)^r sigma = tau(sigma tau)^r

to x. Its two sides become `tau^r x` and `tau^(r+1) x`, respectively,
forcing tau x=x, a contradiction. This argument allows repetitions, so
it also excludes a whole short moved cycle contained in U.

Each of these cyclic words therefore has an O separator between successive
runs of at most r U-labels. Count each run against its following separator.
Summing over moved cycles yields

    |U| <= r|O|.

Since `|U|=k-n` and `|O|=D-2k+n`, this is exactly (1). O-labels on other
cycles only increase the right side. No transitivity, full retention,
involutivity, bound on D, or bound on n was used.

For the actual higher cusp, the inherited marked relation is precisely
`g sigma g^-1=tau`, with `g=(sigma tau)^r`. Re-accessing the same smooth
cusp point with the inherited left-action convention sends its actual
subset A to B=gA. A common retained label is a singleton orbit of the
entire local group; its normal finite extension and pure deleted boundary
identify it with an actual cusp inverse. Thus the geometric cusp fibre
count is n. This is the paid actual-page argument, not a replacement of
retained subsets by all inertia-fixed letters.

The two ordinary pairs in the actual H4 supplier may be changed by a
certified inside-pair inverse Hurwitz move. Apply (1) to the directly
accessed pair, then transport its actual intersection count. The joint
fixed set and support intersection are also invariant under that local
basis change. No new re-access equation is inferred solely from a free
algebraic conjugation. The actual five-cusp pair requires no such repair.

### The one-deficit strengthening

Put `F_sigma=Fix(sigma)`, `F_tau=Fix(tau)`, and `J=F_sigma intersect F_tau`.
The conjugacy gives `g F_sigma=F_tau`, while g fixes J pointwise. Hence

    A intersect B = A intersect J.                    (5)

Let `a0=D-t` and `delta=a0-k>=0`. Applying (1) to the full fixed sets,
or applying the inherited run lemma to moved supports, gives

    j=|supp(sigma) intersect supp(tau)|>=ceil(t/(r+1)),
    |J|=D-2t+j.

At most delta letters of J are missing from A, not twice delta. From (5),

    n>=|J|-delta=k-t+j>=k-floor(rt/(r+1)),

which proves (2). The intermediate inequality `n>=k-t+j` retains the
specific marked support intersection and is often stronger than its
floor relaxation.

### Sharpness and the first failed implication

For each r take disjoint labels `o,a_1,...,a_r,b_1,...,b_r` and put

    sigma=(o,b_1,...,b_r),  tau=(o,a_1,...,a_r),
    A={a_1,...,a_r},       B={b_1,...,b_r}.

Then g=(sigma tau)^r maps A to B and conjugates sigma to tau. Here
`D=2r+1,k=r,n=0`, so (1) is equality. Disjoint unions of these blocks
and jointly retained fixed letters remain equality controls.

In particular `sigma=(1 2 3), tau=(3 4 5)` on five letters, with
`A={4,5}`, `B={1,2}`, has r=2,k=2,n=0. The alternating five-letter words
agree and `gA=B`. Its moved supports meet in one letter. It attains
`3n=5k-2D=0`, but violates `4n>=5k-D` and the ordinary half-support bound.
The first failed implication in those stronger guesses is limiting an
exponent-five cycle run to one, rather than two, outside labels. This
local marked permutation example is not a Keller map or a representation
of the actual global curve complement.

## 3. What the H4 Euler and support counts retain

For the geometric hypotheses of §1, all four positive meridians are
conjugate, so share moved count t. Their actual retained count along the
whole smooth support is the same k, with `1<=k<=D-t`. Let the cusp counts
be n_ab,n_bc,n_cd and write W for the sum of all N actual node overlaps.
The normalization-A1 Euler integration gives exactly

    1=-2k+n_ab+n_bc+n_cd+W,
    omega_p>=max(0,D-2k).                              (6)

The cusp being (2,5) changes its local group, but not the unibranch Euler
stratum. No ordinary-cusp inequality is used there. With
`f=floor(t/2)` and `h=floor(2t/3)`, (1)--(2) imply

    n_ab,n_bc >= max(0,k-f,ceil((3k-D)/2)),
    n_cd >= max(0,k-h,ceil((5k-2D)/3)),
    W <= 1+2k-n_ab-n_bc-n_cd.                          (7)

Every marked node overlap is at least the intersection of its two moved
supports, because those supports are contained in its deleted-sheet sets.
At a commuting node, the support intersection is invariant under either
meridian, hence a union of its nontrivial cycles. In particular it cannot
have exactly one letter. For the cycle type (3)(2), every positive such
intersection has size at least two.

Write S_a,S_b,S_c,S_d for the four supports. The two ordinary edges give

    |S_a intersect S_b|, |S_b intersect S_c|>=ceil(t/2).

For t=5 this forces S_a and S_c to overlap: their intersections with S_b
have total cardinality at least six inside its five letters. Their marked
node therefore has support intersection at least two.

### An additional full-retention incidence bound

If `k=D-t`, each actual retained set is its meridian's entire fixed set.
Consequently every marked cusp count is
`D-2t+|S_i intersect S_j|`, and every marked node overlap equals the
corresponding support intersection. The six marked pairs are exactly all
six pairs of four supports. Additional nodes contribute nonnegatively.
Thus

    1 >= D-4t + sum_{i<j}|S_i intersect S_j|
      >= 4t-2D.                                      (8)

For the last step, if a letter belongs to m of the four supports, it
contributes `binom(m,2)` to the sum, and
`binom(m,2)>=2m-3` for m=0,1,2,3,4. Summing and using total incidence 4t
proves (8). Transitivity in fact removes m=0, but is not required for this
inequality. Since D and t are integers, full retention therefore forces

    D>=2t.                                            (9)

This argument does not assume that the six actual pairs were initially
presented in the same elementary fibre basis. The marked supplier pays
the simultaneous maps, and with full retention all these counts are joint
fixed-set/support cardinalities, unchanged under its permitted local moves.

## 4. Complete exclusion below degree twelve

The PROVED [H4 involution theorem](planar_jc48_sep08_h4_involution.md)
excludes every involutive meridian action, with arbitrary retained subsets.
The PROVED [H4 single-cycle theorem](planar_jc48_sep08_h4_single_cycle.md)
excludes every meridian with just one nontrivial cycle. Thus a geometric
survivor has t>=5; for t=5 its only possible type is (3)(2). These are
proved dependencies, not inferred from a numerical finite-group bank.

For D<=7, use k<=D-5 and nonnegative cusp counts in (6):
`1>=3D-8k>=40-5D>=5`, impossible. For D=8, k<=3; k<=2 is already
impossible by `3D-8k>1`. The only remaining k=3,t=5 has each ordinary
cusp count at least one and W>=6, so the right side of (6) is at least two.

For D=9, k<=4; k<=3 again gives `3D-8k>1`. The remaining t=5,k=4 is
full retention. Its ordinary counts are at least two, the five-cusp count
at least one, and (7) gives W<=4. Every pair of five-element supports in
nine letters overlaps. Each of the three marked commuting nodes therefore
has overlap at least two, giving W>=6, a contradiction. Equation (9) also
excludes this case independently.

For D=10, k<=5 and k<=3 is impossible by the node bound. If k=4,t=5,
(7) gives cusp total at least five, W>=6, and Euler at least three.
If k=4,t=6, retention is full and (8) gives `1>=4`, impossible.
Thus only t=5,k=5 remains. Its two ordinary counts are at least three,
its five-cusp count at least two, so W<=3. The node (a,c) already has
intersection at least two. The two other marked node intersections must
be zero, since any positive one has at least two letters. Hence S_d is
disjoint from S_a union S_b. These three supports all have size five in
ten letters, forcing S_a=S_b. Now the full-retention cusp counts satisfy
`n_ab=5,n_bc>=3,n_cd>=2`, and the node (a,c) has overlap at least three,
since S_a=S_b. Equation (6) is at least `-10+5+3+2+3=3`, impossible.

For D=11, k<=6. The inequalities (7) and the node bound leave only
`(t,k)=(5,5),(5,6),(6,5)`. Explicitly k<=3 is impossible directly;
k=4 gives W>=9 while the cusp lower bounds at t=5,6,7 already make
Euler exceed one; at k=5 only t=5,6 are allowed, and at k=6 only t=5.
The complete declared integer head is reproduced in the source.
The (6,5) case is full retention and contradicts (9). For (5,5), (7)
gives W<=3, whereas every node has omega>=1 and (a,c) has omega>=2,
so W>=4. Finally (5,6) is full retention, with ordinary counts at least
four, five-cusp count at least three and W<=2. As above the node (a,c)
uses at least two, so S_d is disjoint from S_a union S_b. In eleven
letters this forces `|S_a intersect S_b|>=15-11=4`. Full retention now
gives `n_ab>=1+4=5,n_bc>=4,n_cd>=3`, and W>=2. Equation (6) is at least
`-12+5+4+3+2=2`, impossible. This exhausts D<=11 and proves (4).

At D=12 the same scalar and incidence reductions still leave, as necessary
numerical data only, `(t,k)=(5,6),(6,5),(6,6)`. The pre-incidence scalar
head also has (5,7) and (7,5), removed respectively by the (a,c) node bound
and (9). No surviving row is asserted to be a marked permutation passport,
a complement representation, or a Keller map. Degree twelve and higher
remain outside the present exclusion.

## 5. An unbounded stopping control for the cardinality relaxation

The new local inequality does not by itself close H4. Here is an explicit
unbounded family satisfying the retained-count, support-incidence and node
capacity bounds, even with all three node support intersections positive.
It deliberately supplies sets and integers, not permutations or re-access.

For any integer q>=2 use four support labels a,b,c,d. Partition an ambient
set into cells according to exactly which supports contain a point:

| Cell | only a | a,b | only c | b,c | only d | c,d | a,b,c,d |
|---|---:|---:|---:|---:|---:|---:|---:|
| Cardinality |6q|6q|2q|6q|8q|4q|q|

Then D=33q, every support has t=13q, and put

    k=16q, delta=D-t-k=4q,
    (n_ab,n_bc,n_cd)=(10q,10q,8q),
    (omega_ac,omega_ad,omega_bd)=(2q+1,q,q).             (10)

The three edge intersections are (7q,7q,5q), and all three node
intersections equal q. Both ordinary half-support bounds and the
five-cusp third-support bound hold. Each n equals `k-t+j` for its edge,
and satisfies both (1) and (2). Every node omega is at least its support
intersection and `D-2k=q`; no node support intersection is a singleton.
Finally `-2k+sum n+sum omega=1` exactly. The full-retention condition of
(8) is false because delta>0, so it cannot exclude this control.

Thus these cardinality inequalities alone admit arbitrarily large data.
They do not supply the permutations, exact cycle runs, simultaneous actual
retained subsets, or global transitivity. The first unpaid predicate is a
common cycle action realizing the six labelled local relations. This is
the strongest survivor of the current relaxation, not a counterexample to
any geometric theorem or a proof that a stronger H4 obstruction is absent.

The connection map sends the actual finite cover to its marked permutations,
then to support cells and actual cardinalities. It preserves the necessary
Euler/count inequalities, but the second map destroys cycle order and the
actual re-access maps. The cheapest next decisive test is a cycle-compatible
realization or obstruction for a small residual profile, not another generic
increase in the scalar cutoff. The five-letter hostile remains the control
against importing an ordinary braid bound at the fifth cusp.

## 6. Exact reproduction and status

Source: [planar_jc48_sep08_odd_retention.py](../../04-computation/planar_jc48_sep08_odd_retention.py).
Frozen output: [planar_jc48_sep08_odd_retention.out](planar_jc48_sep08_odd_retention.out).
Run from the repository worktree:

```bash
python3 04-computation/planar_jc48_sep08_odd_retention.py
python3 -O 04-computation/planar_jc48_sep08_odd_retention.py
```

The declared complete small universe is every ordered pair in S_D for
D=1,...,5, r=1,...,4, retaining only the literal odd-braid equation, then
every subset A of Fix(sigma), including empty/full subsets. There are
60,068 pair trials, 2,328 braided pairs, 6,096 subsets and 2,192 partial
retention cases. All conjugator, cyclic-run, direct and one-deficit bounds
are always-active checks. This controls the proof; it does not infer an
unbounded theorem from a finite permutation search.

Further controls include 36 sharp equality examples for r=1,...,12,
the named false-bound hostile, all sixteen four-support membership cells,
the complete scalar head through D=12 with t>=5 and all allowed k,
each analytic small-degree contradiction, and five evaluations of the
unbounded formal set family. The last family is explicitly not a group
bank. The proof of (1)--(2), the low-degree case exhaustion and the
unbounded formal formulas are analytic.

The [independent complete analytic/source audit](planar_jc48_sep08_odd_retention_audit.md)
accepts the local run proof, one-deficit strengthening, every D<=11 case
and the nonrealized unbounded set control. Independent conjugacy-weighted
enumeration recovers the entire declared universe. All earlier frozen
supplier artifacts are unchanged.

Frozen normal and optimized replays are byte-identical: **34,983
always-active gates**, 1,358 output bytes. Source: 8,207 bytes, SHA256
`c5d2865731c5b4866a78298ef2760fe9ff0acf85a58442386a78fc87ea41cfff`.
Output SHA256
`6a8d3d36ca5d2bf86874b78a0df93fe33a3a97c9b082ea98b3e3cfb2cb9dea0a`.
Semantic control SHA256
`ada3df141e8d5fa747aa178a14ff36ddbd8f88d42aa43bbd93cfd8e9002cbaf0`.
The source and output are frozen and independently accepted.

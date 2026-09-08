# Independent audit of single-cycle monodromy and the degree-ten floor

**Status: full independent analytic/source audit PASS.**
Root read the complete [primary proof](planar_jc48_sep08_cycle_support.md),
the entire [standalone source](../../04-computation/planar_jc48_sep08_cycle_support.py),
and the inherited actual cusp/node passport. The arbitrary-retention
[involution theorem](planar_jc48_sep08_d4_retention.md) and its separate
full53-subgroup audit are now accepted and promoted. Their Coxeter
supplier is explicitly cited and paid. The degree-ten corollary therefore
has no remaining RESERVED dependency.

The primary uniformly excludes a positive meridian whose permutation
is exactly one nontrivial cycle, under its actual marked three-cusp
geometry. Combining it with the involution theorem excludes mapping
degrees8 and9. This is not a full three-cusp Keller exclusion.

## 1. The support argument is uniform, with the right direction

For braided permutations sigma,tau, g=sigma*tau satisfies
g*sigma*g^-1=tau. It therefore carries the entire fixed set of sigma
onto that of tau, even when that set is empty. The injection proof is
valid for any retained subset A with B=gA fixed by tau. To make its
short set argument explicit: if x is in A minus B, then tau*x cannot
belong to B, since tau fixes B pointwise. If tau*x also belongs to
A minus B, evaluate sigma*tau*sigma=tau*sigma*tau at x. Both x and
tau*x are sigma-fixed, so tau*x=tau²*x. Thus tau*x=x, giving gx=x
and hence x in gA=B, a contradiction. Consequently tau maps A minus B
injectively into the complement of A union B. This gives2n>=3k-D.

Applying that proved inequality to full fixed sets of single m-cycles
gives support overlap at least ceil(m/2). The conclusion is necessary,
not an iff description of braiding. The proof does not extrapolate the
observed pair-bank overlap list to arbitrary lengths.

If two single cycles commute, each preserves the other's support.
If their supports intersect, the one nontrivial orbit of the second
cycle lies entirely in the invariant support of the first. Equal sizes
then force equal supports. This justifies the identical-or-disjoint
leaf-block structure; the statement would fail without the single-cycle
hypothesis and is never used for arbitrary permutations.

For odd m, two different leaf blocks would require more than m letters
of the central support, so only one block exists and the total support
is at most floor(3m/2). For even m there are at most two blocks. In the
two-block case the central support is exhausted by equal m/2 portions
in them, and the total is2m. Equality in the transitive degree bound
D<=2m therefore forces even m and exactly two leaf blocks. Three leaves
distributed over those blocks have multiplicities2 and1. These are
complete support alternatives, not selected orbit representatives.

## 2. Every use of actual retained information is paid

At a genuine cusp, each point fixed by both local meridians is fixed by
their product. If A has k letters inside an a0-letter fixed set, it
contains at least k-a0+f letters of the joint fixed set. Those letters
also lie in the reaccessed retained set B. Thus n>=k-a0+f, without
assuming k=a0. This is the same elementary deficit bound independently
audited in the arbitrary-retention theorem.

At a node, both permutation supports lie in the actual deleted sets.
Their intersection is therefore a valid lower bound on the deleted
overlap. The proof also uses the separate universal bound omega>=D-2k.
No support intersection is asserted to equal the deleted overlap until
full fixed retention is actually proved in the degree-nine case.

The original third cusp has joint fixed count equal to that of (a,b):
simultaneous conjugation by c carries (a,e) to (a,b), using [a,c]=1
and c*e*c^-1=b. The first node is simultaneously conjugated by e from
(e^-1ae,b) to (a,c); the other two node pairs give (a,d),(c,d).
These are exact marked relations, not permutations of unlabeled cusp
counts. At the two actual inside-pair changes, the braid supplier
separately proves preservation of actual intersections and deleted
overlaps. Joint fixed sets are unchanged because the local subgroup
is unchanged. Applying the deficit bound to the actual local pair and
then using that preserved numerical data is legitimate.

## 3. Both global support regimes force the claimed Euler bounds

If D<2m, two m-element supports in D letters cannot be disjoint.
At every node the local meridians commute and are single m-cycles,
so their supports coincide. Hence each actual node overlap is at least
m. With k<=D-m and at least three nodes, the actual Euler identity gives
1>=5m-2D>1. Extra nodes only strengthen this obstruction.

If D=2m, write q=m-k. Every central-leaf cusp joint fixed count is
m/2; the simultaneous conjugation above gives the same count for the
original third cusp. The three cusp counts are each at least m/2-q.
Among the three marked nodes, precisely one leaf pair has identical
support, giving overlap at least m. The other two give at least2q
each from D-2k. Thus W>=m+4q, and

    -2k+sum n_i+W >= m/2+3q.

For even m>=4 this is at least two. Odd m cannot attain D=2m. For
m=2 it forces q=0 and leaves precisely D=4,k=2. The cited classical
geometric degree-four exclusion is essential, since the S4 control has
the required Euler value one. This establishes the uniform conclusion
for every m>=2 without a finite-group classification.

The lower cusp bound is allowed to be negative in an intermediate
estimate; it remains a valid lower bound. The proof never infers a
negative actual count or uses that bound as an equality. Likewise all
unused nodes contribute nonnegative overlap rather than being omitted
from the exact Euler identity.

## 4. Mapping degrees eight and nine

For D=8 the universal node bound excludes k<=2. For k=3, the cusp
injection gives three counts at least one, hence W<=4, whereas the
three node lower bounds give W>=6. Thus k>=4. A nonidentity meridian
then moves at most four labels: its type is involutive or one single
three- or four-cycle. The two accepted theorems exclude all possibilities.

For D=9, k<=3 contradicts the universal node bound, and k>=5 again
leaves moved support at most four. At k=4 the only additional type
after those two exclusions is (3)(2). It fixes exactly four labels,
so actual retention is now the complete fixed set. The cusp injection
gives n_i>=2, hence W<=3. Each actual node overlap is positive because
D-2k=1. Its full-support interpretation is now justified. Commuting
permutations preserve the support intersection, which is a union of
the first permutation's cycles of lengths2 and3. Every positive
overlap is therefore at least two, giving W>=6, a contradiction.

The auditor independently enumerated all integer retained counts and
all partitions of the possible moved-support size for D=8,9. After
the universal Euler bound, single-cycle exclusion and involutive
exclusion, its only remaining tuple was exactly

    (D,k,cycle lengths,full fixed count)=(9,4,(2,3),4).

This independent path checks completeness of the short cycle-type
list, not the existence of any geometric cover. The final overlap
contradiction is analytic and is also controlled by the producer's
complete local mixed-cycle bank.

Combining these two new degrees with the inherited marked passport
and its cited degrees2--4 gives mapping degree at least ten. The actual
braid and whole ordinary-three-cusp family have already passed separate
audits, so the corollary applies to that entire stated family as whole
irreducible Keller support. It says nothing about a component inside a
larger nonproperness locus or about realizing a permutation action by
a polynomial source.

## 5. Finite controls, source, and frozen output

The producer enumerates every single m-cycle against a fixed standard
one on2m letters for m=2,3,4,5. Choosing each support once and fixing
its smallest element as the cycle's first entry avoids duplicate cyclic
notations. Every possible pair support union embeds in that ambient
set, so these controls are complete for the named lengths. The analytic
proof supplies the all-length quantifier.

The three exact transitive controls replay all six original words,
all cusp and node pairs and all retained counts. Their group orders
are24,12,192 and Euler values1,7,2. The last image has generators of
order four; its order192 does not identify its marked action with the
involutive reflection action. It attains the even support/Euler bound
and is the appropriate hostile to a transitivity-only argument.

The local degree-nine bank enumerates all2520 permutations of type(3)(2)
against one fixed such permutation. The commuting overlap counts have
sizes2,3,5, never one. Its complete analytic cycle-union explanation
is recorded above; the bank is not a census of whole representations.

The complete source uses exact integer permutations and always-active
checks. Normal and optimized replays each pass373 gates and match the
frozen output byte for byte. Source:4,836 bytes, SHA256
`762ae98683bcda5116d66c9354e5c8c4ed37af9fd8227c45cebabd8af300d7aa`.
Output:1,872 bytes, SHA256
`316adec4bbd9767c739dc6f59cfdc3201c042f439a17d8c425e729d1108e662a`.

No mathematical or source correction was needed. Root accepts primary
promotion and updating the conditional status paragraph to reflect the
already accepted involution dependency. Higher-degree extensions need
their own proofs and audits and are not part of this frozen result.

# Matching boundaries and two levels of actual ternary returns

**PROVED scoped:** BC1--BC4, BC6--BC7 below. **FINITE-EXACT:** the exhaustive minimum
orders in BC5, backed by two different counting algorithms. **OPEN:** an
unbounded return grammar and the full odd Hamilton-path spectrum. These are
raw directed Hamilton-path counts H, never H/Aut. Native means complementing
every arc except those of an actual directed Hamilton path.

## BC1. The exact boundary controlling clone growth

Let T be a finite tournament, v a named vertex, O=N+(v), I=N-(v).
Form the bipartite graph B_v on O and I with edge oi exactly when o->i in T.
Let nu be its maximum matching size. Replace v by a named ascending
transitive block of m vertices, obtaining T_m, m>=1.

Then

    H(T_m) = sum_{r=1}^{nu+1} c_r r! S(m,r),
    c_r >= 0 integral, c_(nu+1)>0.

Consequently lim H(T_m)^(1/m)=nu+1. In particular, a family with exact
formula A*3^m+B*2^m+C, A>0, has nu=2: the root's directed triangles
have a two-vertex transversal, and no one-vertex transversal.

**Proof.** The ordered-run decomposition is the one-factor specialization of
[THM-3121, path-cover walk-content substitution kernel](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).
Here c_r counts directed skeleton words containing every exterior vertex
once and v exactly r times, with no consecutive copies of v. Assign the
m increasing clone vertices to r ordered nonempty sets: there are r!S(m,r)
choices, and each assignment gives one actual path. This is a bijection.

Between consecutive v occurrences, a skeleton passes from O to I at least
once. Choose one O->I arc in each intervening segment. The chosen arcs have
disjoint exterior vertices, giving a matching of size r-1. Hence r<=nu+1.

Conversely, a matching o_i->i_i gives the directed skeleton

    v,o_1,i_1,v,o_2,i_2,v,...,o_nu,i_nu,v.

Insert every remaining exterior vertex into this directed word. The usual
insertion argument needs a directed word, not distinct occurrences of v:
prepend if the new vertex beats the first term, append if the last term
beats it, otherwise use the first adjacent transition from a vertex beating
it to a vertex it beats. Existing vertices keep their order and the number
of v occurrences stays nu+1. Thus c_(nu+1)>0. For nu=0 start with v.

Finally r!S(m,r)=sum_j (-1)^(r-j) binom(r,j) j^m. The largest base has
positive coefficient c_(nu+1), so it cannot cancel. Bipartite matching/cover
duality gives the transversal statement. For this size-two case one can
also obtain the cover by the standard alternating-path construction from
a maximum matching: take unreached O vertices and reached I vertices.

**Boundary.** This determines the growth base, not the coefficients or legal
return words. A degree bound min(|O|,|I|)+1 was sufficient for the earlier
experiments; the matching is the exact invariant. A two-cover can lie on
one side or have one endpoint on each side. Neither kind is a cosmetic tie
or an arbitrary tournament orientation of another object.

For a cover by k directed paths, the same argument gives the useful upper
bound r<=nu+k on total clone runs. Keep the ordered endpoints and factorial
factors when using that extension; it is not a scalar-H closure theorem.

## BC2. A constant-cost compiler preserving the clone interface

Suppose named tournaments T,U have the same vertex set and the same arcs
incident to v. If d exterior pairs differ, there is a single uniform recipe
for all m that takes T_m to U_m, keeps every old clone address and every
designated exterior pin, uses 2d legal native complements, 3d temporary
births and 3d temporary deletions, and has at most three temporary vertices
present at once. H and strongness are not promised at intermediate steps.

**Proof for one arc a->b, with a,b different from v.** Choose a Hamilton
path P0 of the exterior restriction T-{a,b}; v is included. Add three
temporary vertices in order:

* z loses only to b among the old vertices;
* s is a sink relative to all vertices then present;
* w beats only a among the vertices then present.

The first legal native path is (a,b,z,P0,s,w). After complementing it,
(a,w,b,z,P0,s) is legal: only the two arcs into and out of w need reversal,
and all other consecutive arcs were retained. Restricted to the previous
vertex set, the symmetric difference of the two path spines is exactly ab.
Two complements therefore reverse ab alone. Delete w,s,z by their names.

The ordered predecessor and successor of v agree in the two paths, including
when v is an endpoint of P0. Expand v to its full ascending clone word in
both paths. Both guards remain legal, and the complete two-step endpoint is
the desired cloned tournament. This is the common-frame mechanism of
[LIFTING L1](../adaptive_four_returns_20261103/LIFTING.md), explicitly replayed
here. Repeat for each differing exterior arc.

The inherited [one-helper compiler, METAGRAPH §5](../adaptive_metagraph_reformulation_20261103/METAGRAPH.md)
already reverses arbitrary named arcs more economically at a single size.
BC2 is a uniform clone-compatible version with an explicit constant frame;
it is not a claim that three helpers are necessary. No old relabelling,
converse operation, disappearing pin, or uncharged path choice occurs.

The integer clock here counts native complements separately from births and
deletions. On order L, a native complement reverses exactly binom(L-1,2)
arcs. Thus BC2's d exterior repairs, at peak order n+m+2, cost exactly
2d*binom(n+m+1,2) arc reversals if expanded to that finer clock. A constant
number of native events does not mean a constant number of arc reversals
as the clone grows.

**Birth extension.** If a larger target's root star extends the source star,
first birth its permanent new vertices with their final incident cuts.
Then apply BC2 to the differing old exterior pairs. New target names may be
chosen when designing the construction; they do not authorize a physical
relabel of the source. `returns.json` retains this distinction.

## BC3. A uniform two-level ternary tree

Let A have outgoing rows (10,36,1,22,7,29). Retain named pins 5 and 1,
and replace vertex 3 by TT_m with addresses (3,0),...,(3,m-1). Then

    h_m=H(A_m)=8*3^m+6*2^m-5,                 m>=1.

Its skeleton coefficients are (31,30,8), since

    R1=1, R2=2^m-2, R3=3^m-3*2^m+3.

There are actual uniform programs from this source to each of

    9h_m+8, 9h_m+10, 9h_m+14, 9h_m+16,

through the indicated intermediate 3h_m+2 or 3h_m+4. All endpoint
tournaments are strong. The same old names, pins and ordered clone block
survive. Their boundary responses may change; we do not assert unchanged
NR readiness or unchanged rooted path-cover profiles.

The first two exact programs are rows 1 and 0 of `scout31.json`:

| Offset | P | Q | Final named rows |
|---|---|---|---|
| +2 | (0,3,1,5,4,2) | (0,3,1,4,2,5) | (10,20,33,22,37,75,31) |
| +4 | (0,1,5,3,4,2) | (0,2,1,5,3,4) | (12,97,82,86,3,93,17) |

Each uses two legal native complements and one permanent birth. The clone
frames agree in P,Q. Source and target have at most four clone runs;
the exact skeleton coefficients are respectively
(31,30,8), (95,90,24), and (97,90,24), with higher coefficients zero.
Thus the displayed ternary identities hold for every m, not merely the
sampled sizes.

The four second steps are:

| Input at m=1 | Output at m=1 | Formula | Final base order | Second-step native complements | Permanent births | Temporary births/deletions | Peak base order |
|---:|---:|---|---:|---:|---:|---:|---:|
| 95 | 287 | 3H+2 | 9 | 4 | 2 | 6 / 6 | 12 |
| 95 | 289 | 3H+4 | 9 | 4 | 2 | 6 / 6 | 12 |
| 97 | 293 | 3H+2 | 8 | 2 | 1 | 0 / 0 | 8 |
| 97 | 295 | 3H+4 | 10 | 6 | 3 | 9 / 9 | 13 |

At general m add m-1 to each order; operation counts stay unchanged.
Add the first-step two complements and one birth for whole-branch totals.
The 97->293 step is program 0 of `next_0.json`. The other three are literal
birth-and-repair programs in `returns.json`. Their differing old arcs are,
in original source orientation:

    95->287: (5,1), (4,5);
    95->289: (2,5), (4,5);
    97->295: (1,0), (4,0), (6,0).

None touches the clone root 3. The final named rows are

    287: (10,52,33,22,5,89,159,319,127),
    289: (138,148,1,22,133,95,415,44,191),
    293: (152,145,195,86,4,95,19,120),
    295: (478,224,210,470,130,221,16,64,758,255).

Their skeleton coefficients are (H,270,72), higher coefficients zero.
BC2 proves the actual source-to-target transport, and BC1/run splitting
proves the all-m count identity. Base strongness plus nonempty transitive
substitution proves endpoint strongness for every m. The finite seed and
program checks are reproduced independently in `audit.py`.

**Unresolved step.** The resulting four targets are not proved to belong
to a source grammar closed under both future offsets. Two levels, even
uniformly for infinitely many sources, do not establish the full spectrum.

## BC4. Why the 31/33 comparison needs its own boundary coordinate

Reverse the exterior arc between 0 and 4 in A to get

    B=(26,36,1,22,6,29), H(B)=33.

The two bases have the same score multiset (1,2,2,3,3,4), an explicit instance
of [THM-1865, H not score-determined and inflation neutrality](../../../01-canon/theorems/THM-1865-h-not-score-determined-and-inflation-neutrality.md).
Only that theorem's exact score and source/sink statements are inherited.

At the same clone root, B has coefficients (33,26,6), so

    H(B_m)=6*3^m+8*2^m-1,
    H(B_m)-H(A_m)=-2*3^m+2*2^m+4.

The difference is 2,-6,-34,-126 at m=1,2,3,4. The midpoint and half-difference
are 7(3^m+2^m)-3 and 3^m-2^m-2 respectively (the latter is A minus midpoint).
Thus 32+-1 at the seed is not a persistent one-bit offset. BC2 can implement
the exterior flip at every m, but it cannot make this changing count defect
constant. The missing data are the separate run coefficients, not an
unrecorded sign convention.

## BC5. Certified minimum sizes in this clone representation

For each of the four functions 9h_m+d, d=8,10,14,16, the smallest possible
base order of any marked tournament T with H(T[TT_m at v])=9h_m+d for all m
is respectively

    9, 9, 8, 10.                         [FINITE-EXACT]

This is not a bound on the minimum order realizing the single integer at
m=1. Many smaller bases realize that scalar without the same refinement law.

**Complete universe.** For n<=6, independent permutation-orbit enumeration
covers every unlabelled tournament and every root. For each n=7,8,9 and each
root out/in split a,b with a<=b and min(a,b)>=2, choose an unlabelled
tournament on each neighbourhood and every a*b cross orientation. These
presentations cover every rooted tournament; they need not be pairwise
nonisomorphic. Converse covers a>b and preserves clone counts after reversing
the clone order for the count bijection. BC1 rules out min(a,b)<2.

Require strongness and coefficients (H,270,72,0,...). The independent
permutation-skeleton enumeration agrees with subset-DP clone evaluation:

| a,b | Presentations | Uniform hits at H=287,289,293,295 |
|---|---:|---|
| 2,4 | 1,024 | 0,0,0,0 |
| 3,3 | 2,048 | 0,0,0,0 |
| 2,5 | 12,288 | 0,0,0,0 |
| 3,4 | 32,768 | 0,0,3,0 |
| 2,6 | 229,376 | 0,0,0,0 |
| 3,5 | 786,432 | 1,1,0,0 |
| 4,4 | 1,048,576 | 0,0,0,0 |

All four target sequences have gcd 1 among their first three values. If a
weak tournament realized one, every strong component outside the root's
component would have constant H dividing that gcd. Such components have H=1
and are singletons. Removing them leaves a smaller strong rooted realization
with the same function. Thus the complete lower-order strong census also
excludes weak lower-order representations.

The upper witnesses are in BC3. Only a matching-cover-filtered portion of
order ten was searched: a,b=(3,6) and (4,5), 1,502,368 and 454,368
presentations. This is sufficient for the 295 existence witness, not a
complete order-ten census. The two-cover filter is justified by BC1 and is
implemented by taking the union over all two-vertex covers of cross masks.

## BC6. The next boundary is recorded, not assumed repaired

For an exterior pin q, let F_q count ordered two-path covers, with the first
path nonempty ending at q and the second possibly empty. Keep the defect
D_q=F_q-H used by the inherited [ready-interface rule](../two_entry_funnels_20261103/READY.md).
For each new target family, BC1's two-cover bound implies at most four
clone runs in these covers. Hence D_q is an exact linear combination of
R1,R2,R3,R4; four samples determine it, and a fifth is an independent
degree-bound check. `boundary.json` stores all exterior-pin coefficients.

Among the exterior pins, the complete constant-defect list is:

    287 family: q=7, D_q=0;
    289 family: q=6, D_q=0;
    293 family: none;
    295 family: q=9, D_q=17.

These identities hold for every m>=1. Every other exterior defect is
nonconstant, so neither D_q=1 nor D_q=2 holds uniformly. At m=1 no vertex,
including the clone root, has defect 1 or 2. Consequently even a re-marking
of these particular base outputs does not enable the usual ready +2/+4
return. The inherited NR formula h+2F would use offsets 0,0,34 on the three
constant ports listed above. This is an explicit next obstruction, not a
claim that those H values have no other ready representative.

**Proof of the degree assertion.** In an ordered two-path cover there are
at most two paths meeting the clone. Between successive clone runs within
each path, choose a crossing arc as in BC1. These form a matching of size
at least r-2, so r<=4. The ordered skeleton/partition bijection retains
the pin-ending constraint and the empty-second-path case. The computation
uses the exact cut identity

    F_q = sum_{S containing q} end_q(T[S]) H(T[V-S]), H(empty)=1.

Literal enumeration of every directed prefix independently audits all four
determining samples in `audit.py`. This bound, not a constant-looking finite
plot, proves the identities. It also explains why scalar branch success has
not closed the return grammar.

## BC7. Exact linear moment budget for a matching cap

Fix integers d>=0 and K>=1. Consider, over all finite rooted tournaments
(T,v) with crossing matching number at most d, the functions of m>=1 that
count ordered covers of T_m by k nonempty directed paths, 1<=k<=K.
Allowing additional fixed exterior endpoint constraints does not enlarge
the following space. Its rational linear span is exactly

    span_Q {1^m,2^m,...,(d+K)^m},

of dimension d+K. Thus at most d+K run coefficients or exact samples specify
each such function, and no smaller linear function space contains this
whole family. For our d=2, Hamilton queries use three modes; two-path queries
use four. Letting K grow makes this budget unbounded.

**Upper bound.** With r total clone runs, at most k paths meet the clone.
The inter-run crossing arcs give a matching of size at least r-k, so
r<=d+k. Ordered skeletons, with any fixed exterior endpoint conditions,
give a nonnegative integer combination of R_r(m)=r!S(m,r). These have exactly
the indicated exponential bases. Empty auxiliary paths, if allowed as in
BC6, merely add terms with fewer components and obey the same bound.

**Realizable lower bound.** For each j>=0 use a rooted tournament W_j on
1+2j vertices: the root beats o_1,...,o_j and loses to i_1,...,i_j;
the crossing arcs are exactly o_a->i_a, with every other cross arc directed
from I to O. Orient each side as an increasing transitive tournament.
Its crossing graph is a matching of size j. For j=0 it is the singleton.

The matching skeleton from BC1, together with k-1 separate singleton-root
paths, gives an ordered k-path skeleton with j+k clone runs. Insert all
remaining exterior vertices into its first path as before. Hence the
largest base in the k-cover function of W_j is exactly j+k, with positive
coefficient. For every b=1,...,d+K choose

    j=min(d,b-1),  k=b-j.

Then 1<=k<=K and the resulting function has positive top base b and only
smaller bases below it. The d+K functions are triangularly independent.
Vandermonde independence of the sequences b^m proves the claim.

`moment_budget.py` supplies literal witness tournaments and exact coefficient
matrices for (d,K)=(0,3),(1,3),(2,3), with determinants 1,6,144.
The independent audit uses vertex permutations and path cuts instead of
the subset-cover convolution. These are controls of the general proof.

**Exact scope.** This is a dimension theorem for functions of the clone
parameter over the stated collection of sources and queries. It is not
a lower bound on arbitrary nonlinear encodings, a minimum bit count, or a
closed finite-state native controller. In particular it does not say that
these coefficients retain the legal Hamilton path needed for the next
complement. Source realization and path ownership remain separate duties.

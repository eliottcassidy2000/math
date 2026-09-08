# A degree-uniform obstruction for actual marked three-cusp D4 monodromy

**Status: PROVED / FINITE-EXACT / INDEPENDENTLY AUDITED.**
The complete analytic argument, full source and exact controls have passed independent audit.
The geometric conclusion concerns the declared whole nonproperness support;
it is not a solution of the planar Jacobian conjecture.

## 1. Statement, inheritance, and scope

Suppose a nonautomorphic complex planar polynomial Keller map has whole
nonproperness set an irreducible curve C whose normalization is A1 and
whose affine singularities are exactly three ordinary cusps and N>=3
ordinary nodes. Assume four positive curve meridians a,b,c,d generate
the affine complement group and satisfy the Artin D4 relations

    braid(b,a), braid(b,c), braid(b,d),
    [a,c]=[a,d]=[c,d]=1,                               (1)

where braid(x,y) means xyx=yxy. Retain the following actual local
marking, with e=bcb^-1 and f=e^-1ae:

    cusp pairs: (b,c), (b,d), (a,e);
    three distinct node pairs: (f,b), (a,d), (c,d).      (2)

The marking may include the certified inside-pair Hurwitz changes
described in §2, with their actual retained-intersection and deleted-overlap
counts transported. It must not be inferred from an unmarked abstract
presentation alone. Extra global relations and additional ordinary nodes
are allowed.

**Theorem.** No Keller map satisfying these hypotheses exists, in any
mapping degree.

**Geometric corollary.** No irreducible affine curve admitting a birational
polynomial normalization (U(t),V(t)) with deg U=4 and deg V=6, and having
exactly three ordinary cusps and otherwise only ordinary nodes, can be
the whole nonproperness set of a nonautomorphic polynomial Keller map.
The [proved family supplier](planar_jc48_sep08_three_cusp_family.md),
§§1–4, supplies (1)–(2), normalization A1, and respectively N=4 or N=3
on its infinity-(2,7) and infinity-(2,9) good strata. Its marked
representative is supplied by the [exact six-loop braid proof](planar_jc48_sep08_three_cusp_braid.md).

This is a theorem about a specified whole-support class. It does not
exclude that curve as one component of a larger nonproperness set, does
not assert that all Keller curves admit this polynomial-normalization
chart, and does not assert equality of the complement group with Artin D4.
Here D4 names the rank-four Artin/Coxeter diagram, not the dihedral group.

The closest proved mechanisms are the [marked passport and actual Euler
ledger](planar_jc48_sep08_three_cusp_passport.md), the [single-cycle support
theorem](planar_jc48_sep08_cycle_support.md), and the [arbitrary-retention
involution exclusion](planar_jc48_sep08_d4_retention.md). The latter two
suppliers are PROVED and independently audited; only their moved-support
two/four consequences are needed in §5. Their geometric degree-four
input has the original primary source and correct generic-fibre degree
recorded in the [ordinary-cusp passport](planar_jc48_sep06_cusp_passport.md).
No new literature priority claim is made here.

The named hostile is the transitive S4 action satisfying the whole local
Euler ledger. It shows why the small-cycle suppliers cannot simply be
discarded. A second, purely formal degree-fifteen set configuration
attains every numerical bound before commuting is imposed. The corrected
near miss is replacing actual retained sets by all inertia-fixed labels.
The least-used sidecar is incidence *inside the central moved support*,
combined with the fact that commuting supports cannot intersect in one
label.

The five live concepts are marked local pairs, actual retention deficit,
central support incidence, the Euler budget, and the equality boundary.
Root proposed the central-incidence estimate and the resulting parity
elimination; this lane independently derived the complete inequalities
and all marking maps. Certificate_audit independently checked the
degree-uniform argument before the final artifact audit. The planned
intermediate degree-ten through degree-thirteen tables were left unwritten:
the argument below has no finite degree cutoff.

## 2. Actual sheets, Euler integration, and the pair map

Let D>1 be the mapping degree, and let k be the actual affine fibre
count on the entire smooth stratum of C. The inherited page-constancy
and actual-cover suppliers give

    1<=k<D.                                            (3)

For clarity, positivity of k also follows from pulling back a nonconstant
irreducible equation of C: its nonconstant pullback has a divisor, and
quasi-finiteness prevents that divisor from contracting to a point. The
monodromy away from the whole nonproperness curve is transitive. All
positive meridians of the irreducible curve are conjugate, so their
permutations have a common moved-support size t. Since the four marked
meridians generate, t=0 would give a trivial transitive action, contrary
to D>1. A permutation never moves exactly one label. Therefore t>=2.

Write

    a0=D-t,       delta=a0-k>=0,
    D=t+k+delta.                                       (4)

The number a0 is the full inertia-fixed count; k is the actual retained
count. They need not agree. Let n_i be the three actual cusp fibre
counts, and let omega_p be the overlap of the two actual deleted subsets
at a node p. Put W=sum_p omega_p, over all N nodes. The retained subsets
have size k, so

    omega_p>=max(0,D-2k).                               (5)

The actual node fibre count is 2k-D+omega_p. The normalization A1 with
three unibranch cusps and N two-branch nodes has

    chi(C_smooth)=-2N-2,
    chi(C)=1-N,       chi(A2\C)=N.

Euler integration of the finite actual fibres, using the inherited
local sheet suppliers at the singularities, thus gives

    1=ND+(-2N-2)k+sum_i n_i+N(2k-D)+W
     =-2k+sum_i n_i+W.                                 (6)

Euler characteristic may be read as compactly supported Euler
characteristic; no properness of the original Keller map is being
assumed. The off-C covering and the actual singular-fibre formula are
the proved suppliers, not an abstract permutation substitute.

Here is the precise map from the original local data to support sets.
The braid relation between b,c gives

    e b e^-1=c,        c e c^-1=b.                      (7)

Conjugating the original node0 pair (f,b) simultaneously by e gives
(a,c). The other two node pairs in (2) are already (a,d),(c,d).
Consequently their support intersections are exactly the three pair
intersections of the leaf supports supp(a),supp(c),supp(d).
The original third cusp pair (a,e), conjugated simultaneously by c,
becomes (a,b), because a commutes with c. Thus the three original
cusp joint fixed counts correspond to the three central-leaf pairs.

At the directly certified path level, node2 initially has pair
(e,bdb^-1), jointly conjugated by b^-1 to (c,d). Two further local
pairs in the braid supplier differ from (2) by one inside-pair inverse
Hurwitz move. In the convention

    (sigma,tau) -> (tau,tau^-1 sigma tau),
    (A,B)       -> (B,tau^-1 A),                        (8)

where B is pointwise fixed by tau, the actual retained intersection
size and deleted-overlap size are unchanged. The generated local
subgroup, and hence its full joint fixed set, is unchanged as well.
These are proved local distinguished-basis changes in the certified
cluster, not assumed base-coordinate loops. The cusp re-access
inequality below is applied first to the directly certified actual
pair, then its count is transported by (8). No new re-access equation
is inferred from an abstract Hurwitz change alone.

All simultaneous conjugations preserve support intersection sizes and
transport both actual subsets. The node estimates below use only that
a moved label cannot belong to the actual retained set. Thus each
actual deleted overlap dominates the relevant support intersection,
without requiring full retention or identifying the retained sets at
different singularities.

## 3. Arbitrary-cycle braid overlap and the cusp deficit bound

We need a local lemma that does not assume a single cycle or a particular
ambient degree. Let sigma,tau be braided permutations of a D-element
set Omega, and put g=sigma*tau. The braid relation gives
g sigma g^-1=tau. Suppose A is pointwise sigma-fixed and B=gA is
pointwise tau-fixed, with |A|=|B|=k and n=|A intersect B|. Then

    2n>=3k-D.                                          (9)

To prove this, put U=A\B and O=Omega\(A union B). Since tau fixes B,
tau(U) cannot meet B. If x and tau(x) both lie in U, then sigma fixes
both. Applying sigma*tau*sigma=tau*sigma*tau to x gives
tau(x)=tau^2(x), so tau(x)=x. Hence g(x)=x, which contradicts x not
in B=gA. Therefore tau(U) is contained in O, so
k-n<=D-2k+n. This proves (9), even when the chosen fixed set is empty.

Apply this abstract lemma to the *full* fixed sets A=Fix(sigma) and
B=Fix(tau), which do satisfy B=gA. If their common moved count is t
and their support intersection has size j, then

    k_full=D-t,       n_full=D-2t+j.

Substitution into (9) gives the arbitrary-cycle half-support bound

    j>=ceil(t/2).                                      (10)

This application to full fixed sets is purely group theoretic. It
does not change the actual retained count k in (3)–(6).

For a directly accessed actual ordinary cusp, the retained-page
supplier gives A subset Fix(sigma), B=gA subset Fix(tau), and its
actual cusp count n=|A intersect B|. Let J0 be the full joint fixed
set. Every point of A intersect J0 is fixed by g, hence belongs to B.
With |J0|=D-2t+j and |Fix(sigma)|=a0=D-t, this gives

    n>=k-a0+|J0|=k-t+j>=k-floor(t/2).                  (11)

The bound is valid if its right side is negative; n>=0 is still true
but is not needed for the next argument. The maps of §2 preserve the
joint fixed cardinality and the actual intersection count, so (11)
applies to all three original cusp counts with the common t.

Set

    F=floor(t/2),        C=ceil(t/2),        t=F+C.

Equations (6) and (11) yield the first central estimate:

    W<=3F+1-k.                                         (12)

## 4. Three independent lower bounds on node overlap

Let L1,L2,L3 be the moved supports of a,c,d, and let B0=supp(b).
Each has size t. Write

    P=sum_{i<j}|Li intersect Lj|.

The three distinct marked nodes of §2 give W>=P. Additional node
overlaps are nonnegative and remain in W.

First, exact three-set inclusion-exclusion gives

    P=3t-|L1 union L2 union L3|+|L1 intersect L2 intersect L3|
     >=3t-D.

Second, put Ti=Li intersect B0. By (10), |Ti|>=C. These three sets
lie in the t-element central support, so

    P>=sum_{i<j}|Ti intersect Tj|
      =sum_i|Ti|-|T1 union T2 union T3|+|T1 intersect T2 intersect T3|
      >=3C-t.

Third, (5) and N>=3 give W>=3(D-2k). If D-2k<0, this last inequality
follows from W>=0; if D-2k>=0, it follows by summing at least three
nonnegative lower bounds. Thus no multiplication by N has reversed
an inequality in the negative case. Together:

    W>=3t-D,        W>=3C-t,        W>=3(D-2k).          (13)

These bounds retain different data: the ambient union of the leaves,
their incidence inside the central support, and the actual retained
cardinality at every node. None alone replaces the other two.

## 5. Unbounded parity reduction and the final equality obstruction

Compare each lower bound in (13) with (12), and use D=t+k+delta.
The exact resulting inequalities are

    delta>=2t-3F-1,
    2k>=3C+3delta-1,
    k<=1+4F-2C.                                        (14)

The first comes from the ambient union, the second from the actual
node bound, and the third from central incidence.

If t=2m+1 is odd, (14) gives

    delta>=m+1,       k>=3m+3,       k<=2m-1,

which is impossible. If t=2m is even, (14) gives

    delta>=m-1,       k>=3m-2,       k<=2m+1.

Consequently m<=3. This proves the all-degree reduction to moved
counts t=2,4,6 without any finite cutoff on D, delta or k.

For t=2, the positive meridian is a transposition. For t=4, its
nontrivial cycle type is either one four-cycle or two transpositions.
The PROVED [arbitrary-retention involution theorem](planar_jc48_sep08_d4_retention.md),
§4, eliminates the transposition and two-transposition cases with these
same actual markings and all N>=3. The PROVED [single-cycle theorem](planar_jc48_sep08_cycle_support.md),
§4, eliminates the single four-cycle case. Their classical
geometric-degree-four input is essential, as the S4 hostile below shows.

It remains to eliminate t=6. Now (14) forces

    delta=2,       k=7,       D=15.                     (15)

Indeed k<=7 and 2k>=8+3delta, together with delta>=2, give equality.
Equation (12) gives W<=3, while D-2k=1 gives W>=N>=3. Therefore

    N=3,       W=3,       omega_p=1 at every node,
    n_1=n_2=n_3=4.                                     (16)

The cusp equalities follow from (11) and (6); the contradiction already
follows from the node equalities. For any commuting permutations x,y,
their moved-support intersection is invariant under x. The intersection
contains no x-fixed label by definition. It is therefore a disjoint
union of nontrivial x-cycles, so its size is zero or at least two,
never one.

At each marked node the support intersection is at most its actual
deleted overlap omega_p=1. It must consequently be zero. By the paid
pair map of §2, all three leaf intersections vanish, giving P=0.
But the central-incidence estimate already proves

    P>=3C-t=9-6=3.

This contradiction removes the last case and proves the theorem.

## 6. Sharp formal controls and the preserved information

The inequalities (14), together with delta>=0 and k>=1, have exactly
seven post-reduction integer solutions (t,delta,k,D):

    (2,0,1,3), (2,0,2,4), (2,0,3,5), (2,1,3,6),
    (4,1,4,9), (4,1,5,10), (6,2,7,15).                (17)

This is a description of a numerical relaxation, not a list of actual
Keller passports. The analytic reduction precedes this finite list.

The degree-fifteen survivor has a concrete support-set hostile. On
labels 0,...,14 take

    B0={0,1,2,3,4,5},
    L1={0,1,2,6,7,8},
    L2={2,3,4,9,10,11},
    L3={0,4,5,12,13,14}.

All four supports have size six, each |B0 intersect Li|=3, and each
leaf-pair intersection has size one. The leaf union is all fifteen
labels, so both support inequalities are equalities. The formal values
k=7, delta=2, n_i=4, omega_p=1 satisfy the Euler ledger. These are
only sets and numerical counts: they cannot be the supports of the
required commuting leaf permutations. This is the first failed
implication when the commuting-cycle sidecar is dropped.

For the small-cycle boundary use, on four labels, a=d=(34), b=(23),
c=(12). This transitive S4 action satisfies the full marked D4
relations. With all k=2 fixed letters retained, the original cusp
counts are (1,1,1), the original node overlaps are (0,2,0), and the
Euler value is one. Thus group relations, transitivity, actual fixed
sets and the Euler identity alone have this honest abstract survivor.
The inherited geometric-degree-four theorem is the precise final
supplier excluding it from Keller geometry.

The successful connection is:

| Coordinate | Role |
|---|---|
| Source | Actual marked cusp/node pages of a whole irreducible nonproperness curve |
| Target | Four permutation supports with a common moved count and a retained deficit |
| Map | Apply monodromy, retain k, and transport the original pairs only by the paid simultaneous maps |
| Preserved predicate | Actual Euler identity and lower bounds for actual cusp fibres and deleted node overlaps |
| Lost information | Full permutations, access-path words, individual deleted fixed labels, and most local compatibility |
| Required sidecars | Braid half-overlap, central incidence, commuting nontrivial-cycle sizes, and small-cycle geometric exclusions |
| Cheapest decisive test | The formal t6/D15 support configuration has singleton leaf intersections, forbidden by commutation |

The proof applies to any whole-support curve satisfying §1, even beyond
the supplied (4,6) class. It does not establish those hypotheses for
other curves. The three marked nodes and three cusps are used
quantitatively; changing this inventory requires a new Euler and
incidence analysis.

## 7. Exact controls, reproduction, and audit boundary

The source is
[planar_jc48_sep08_three_cusp_closure.py](../../04-computation/planar_jc48_sep08_three_cusp_closure.py),
with frozen output [planar_jc48_sep08_three_cusp_closure.out](planar_jc48_sep08_three_cusp_closure.out).
It uses only the Python standard library, integer permutations and
finite sets, and always-active checks.

The declared control universes are:

1. Formal affine coefficient identities for the three comparisons and
   both unbounded parity formulas. After that analytic reduction, the
   complete integer universe m=1,2,3 with its derived exact delta and k
   bounds gives precisely (17). There is no experimental degree cutoff.
2. All unordered-with-repetition triples of t-subsets on D=2,...,7,
   for 2<=t<D, whose intersection with the fixed central t-subset has
   size at least ceil(t/2). All 9,441 configurations satisfy both
   inclusion-exclusion bounds.
3. For each D=2,...,6 and every partition of D, the canonical permutation
   sigma of that type against every tau of the same type. This is the
   complete same-type pair universe up to simultaneous relabeling,
   with 872 pairs, including identity, mixed-cycle and empty-fixed-set
   controls. All 152 braided pairs, all 673 subsets of their full fixed
   sets, and all 74 commuting pairs are checked. The source explicitly
   tests the injection, retention deficit, half-overlap, inside-pair
   count transport and no-singleton commuting intersection.
4. The formal degree-fifteen support hostile and the marked transitive
   S4 Euler-one control, including the actual simultaneous pair maps.

These finite banks check the algebra and failure boundaries. They are
not used to infer a theorem outside their universe; §§2–5 provide the
degree-uniform proof.

Reproduce from the worktree root:

    python3 04-computation/planar_jc48_sep08_three_cusp_closure.py
    python3 -O 04-computation/planar_jc48_sep08_three_cusp_closure.py

Both runs pass 43,138 always-active gates. The semantic report digest is
`701ae302f997498f7168d77fda33452e5d73fc35b3f555b9ca5f946e4d392710`.
Final byte hashes and independent replay status are recorded below.
The source does not import any reserved
proof, numerical braid scout, symbolic degree assumption or unpublished
finite permutation classification.

The source and output are frozen, with normal and optimized output
byte-identical:

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Source `.py` | 11,027 | `3b7a0978e0ecfa11cdf5c6d44cb66fb66109ffff965361bf0b5177061a036423` |
| Output `.out` | 764 | `bca03a5fd8697bb66984957caea989f299bdfd1b6eb4e4b36c227dd25c849377` |

The [independent final audit](planar_jc48_sep08_three_cusp_closure_audit.md) accepts the entire proof, source and both43,138-gate replays. Its separate binomial and permutation counts recover all control universes and the seven exact scalar rows. The actual marked family supplier and all small-cycle dependencies are proved.

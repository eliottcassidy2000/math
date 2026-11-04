# Integer census and the boundary of uniform ternary returns

2026-11-03. Local labels IC/SF/FG are scoped to this file. **PROVED** means
the argument below; finite exhaustive conclusions also use the stated
**FINITE-EXACT** receipts and the **CITED** strong-order bound. Same-agent
review with independent algorithms, not an external proof audit.

All sources are Boolean tournaments. H counts directed Hamilton paths,
not paths divided by automorphism order. F_r is the ordered two-path-cover
count whose first, nonempty bin ends at r. A root is b-ready if F_r-H=b.
Cloning q means replacement by an ascending transitive module TT_m.

## IC1. An integral version of the projective Burnside recurrence

**Recovered mechanism, not a new theorem.** The final inheritance pass
located this exact full-factorial gauge, formula, divisibility proof and
unit completion in [Session105 section1](../session105_metagraph_recursion_20261029/enumeration/INTEGER_PROJECTIVE.md).
Sections2--3 there also contain the odd gauge and packed degree buckets
used in IC2. The proof below is retained as a self-contained derivation
and computational reference; this session contributes an independent
implementation/comparison, not a new counting method or compression.

The inherited recurrence is
[Session98, recurrence](../session98_tournament_growth_20261029/recurrence/verify.py),
with its unit-tail mechanism in
[THM-514, unit-tail core kernel and prime-tail ladder](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md).

For a partial multiset mu of odd permutation-cycle lengths, let

    r=N-|mu|,  z_mu=product_p p^(k_p) k_p!,
    h_j=sum_(l in mu) gcd(l,j),
    e=sum_(l in mu) (l-1)/2 + sum_(i<j) gcd(l_i,l_j).

Even permutation cycles cannot fix a tournament: a half-turn reverses an
opposite pair. For odd cycles, e is the number of unordered-pair orbits
within the already placed cycles. The usual Burnside summand is 2^e/z_mu.

The inherited projective state retains r and s_j=h_j-j*h_1 for future
odd lengths j<=r, with rational weight W=(2^e/z_mu) 2^(r*h_1). Replace W by

    I=(N!/r!) W = N!/(r! z_mu) * 2^(e+r*h_1).              (IC1)

This is an integer for every individual history: N!/(r! z_mu) counts the
ways to place the partial cycle type on the N labelled vertices. Therefore
the sum over histories merged into a state is also integral.

At future length p, insert k cycles and put r'=r-kp. Then

    s'_j=s_j+k(gcd(p,j)-j),
    E=k(p-1)/2 + p*k(k-1)/2 + k*s_p + r'*k,
    I'=I * [r!/(r'! p^k k!)] * 2^E.                       (IC2)

The bracket is an integer cycle-placement count. E can be negative; the
division by 2^(-E) is nevertheless exact. To see this, apply IC2 before
merging histories: its result is precisely IC1 for the new partial type,
whose exponent e'+r'*h'_1 is nonnegative. All histories in the merged state
have the same E and placement multiplier, so the divisibility survives
aggregation. The implementation checks each division with a remainder.

At p=1 all r remaining vertices become unit cycles. They contribute
2^binom(r,2), with placement multiplier one. Thus

    a_N = (sum_(terminal states) I * 2^binom(r,2))/N!.      (IC3)

For N=0 return a_0=1. This is exactly A000568. The state quotient is
unchanged; the improvement removes Fraction normalization and completes
the unit tail directly. No new asymptotic complexity or enumeration of
individual isomorphism classes is asserted.

`integer_census.py` compares all N<=30 against independent odd-partition
Burnside, and N=60,100,150,200 against the same-state rational algorithm.
It also matches the original projective implementation at60 and100.
At200 it checks983315 exact negative shifts and2303050 transitions, with
95709 peak states. Observed times were6.83s integral versus19.72s rational;
timings are measurements, not part of the mathematical certificate.

### IC2. Integration with the incoming integral production engines

Incoming [Session175 projective arithmetic](../session175_tournament_spectrum_20261104/a000568/projective.py)
already combines the same projective quotient with a degree-local odd
factorial denominator. The older [Session87 degree buckets](../session87_four_vertices_enumeration_20261028/FAST_COUNTS.md)
already share multiplicity work across response states. Integer arithmetic,
the projective quotient and bucket scheduling are therefore inherited.

For used degree d=N-r, let O(d)=odd(d!) and let
W_odd=O(d)*2^(e+r*h_1)/z_mu be the incoming weight. The exact change to IC1 is

    I=binom(N,d)*2^v2(d!)*W_odd.                          (IC4)

Both sides are integer gauges of the same histories. `bucket_census.py`
reimplements the inherited projective packing and degree-bucket combination,
and implements both gauges. Store states[d][response], then loop over
d,k,response, sharing each multiplicity factor and expired-coordinate
mask. Packing is exact: each digit j*h_1-h_j is nonnegative and at most d,
because every processed length l>j and j-gcd(l,j)<=l. Width bit_length(N)
therefore prevents carries. The k=0 branch still expires coordinates.
Every other transition is IC2 or its algebraically equivalent odd gauge.

The four-engine comparison matches every value tested, including0..30
for both new gauges and100,150,200 against incoming packed arithmetic and
the older full-prefix counter. At200 the first observation was2.361s
incoming packed,1.994s bucket odd,1.850s bucket placement, and1.634s for
the older full-prefix counter. Thus the bucket combination improves that
packed scalar implementation in this test, but **does not beat the older
full-prefix engine**. The2.9x number above is only the comparison with
Fractions. No general speed or memory improvement is claimed.

## SF1. Dyadic growth forces a matching-one boundary

Use the proved ordered-run law
[BC1](../boundary_carry_returns_20261103/PROOFS.md#bc1-the-exact-boundary-controlling-clone-growth).
For O=N+(q), I=N-(q), let nu be the matching number of the bipartite arcs
O->I. Then

    H(T_m)=sum_(j=1)^(nu+1) c_j j! S(m,j),
    c_j>=0, c_(nu+1)>0.                                  (SF1)

Consequently, if H(T_m)=A*2^m+B with A>0 for all sufficiently large m,
then nu=1. This is a structural conclusion, not a fit to two samples.
It also makes the same formula valid for every m>=1. Endpoint-restricted
Hamilton counts, including counts starting at a fixed exterior vertex,
are obtained by restricting the same skeletons. They lie in span{1,2^m}
when nu<=1; their leading coefficient may vanish.

## SF2. Two samples certify exterior one-readiness at every clone size

Let T be strong and let r be exterior to the cloned vertex q. Suppose
r has unique in-neighbour w, q is neither r nor w, and nu_T(q)<=1.
For every m the root star still has unique in-neighbour w. The exact
single-in-port identity gives

    F_r(T_m)-H(T_m)=start_w((T-r)_m).                      (SF2)

Deleting r cannot increase the crossing matching at q. By SF1's
endpoint-restricted form, the right side is alpha*2^m+beta. If its values
at m=1,2 are both one, alpha=0,beta=1. Hence r is one-ready for all m>=1.
No extra empirical samples are a premise of this conclusion.

More generally, suppose a fixed strong cloned family has a one-ready mark
for infinitely many m; the mark is allowed to change. There are finitely
many exterior vertices. If an exterior vertex is used infinitely often,
it cannot receive arcs from the clone module, since its indegree would
grow with m while a one-ready root has indegree one. Thus its unique
in-port is a fixed exterior w. SF2 implies that this root was already
one-ready at m=1 and stays so at every m. If clone marks occur infinitely
often instead, apply the complete clone-mark classification
[RC1](../return18979_ports_20261103/PROOFS.md#rc1-the-complete-ready-locus-of-transitive-root-growth).
Only the first clone can work; its exterior grandparent must be transitive.

This is a compactness statement for a *fixed finite exterior*, proved by
finite pigeonhole and two exponential modes. It gives no license to
identify different sources having equal H or equal phases.

The same obstruction extends to a bounded number of exterior vertices
whose orientation may vary with m: there are finitely many base isomorphism
types, so one recurs infinitely often. SF1's positive leading exponential
forces matching one even if its H values equal the target only on that
infinite subsequence; two distinct dyadic affine laws cannot agree there.
Apply SF2 to that recurring base. Thus a successful one-ready repair of
SF3's target with a single uniform TT_m module needs unbounded exterior
order, or must change the module hypothesis itself.

## SF3. The1182 debt cannot be removed by any fixed exterior with one-ready closure

The previous repair gives h_m=1243*2^m-1061. Its desired +2 family is

    g_m=3h_m+2=3729*2^m-3181.                              (SF3)

**Claim, with finite certificate:** there is no fixed finite strong T and
vertex q for which H(T_m)=g_m for all sufficiently large m and T_m has a
one-ready vertex for infinitely many m. The exterior may have any fixed
number of feedback vertices, and the designated ready vertex may vary.

SF1 forces nu=1. If infinitely many ready marks are clones, RC1 puts the
family in the one-apex carrier
[RC3--RC4](../return18979_ports_20261103/PROOFS.md#rc3-a-complete-binary-carrier-for-this-ready-growth-locus).
Its constant is the binary recurrence B(0)=1, B(2a)=B(a),
B(2a+1)=B(a)-a-(a mod2). At coefficient3729, B=-1999, not-3181.

Otherwise SF2 gives a one-ready exterior root already at m=1, with H=4277.
Every strong one-ready source is an indegree-one newborn of a strong
unique-start parent, by
[READY sections1--3](../two_entry_funnels_20261103/READY.md).
The parent satisfies H_parent<=(4277+1)/2=2139. Busch's sharp strong-order
bound has f(16)=3125, so the parent order is at most15.
[Busch2006, Corollary1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).

The complete unique-start grammar through15 contains196417 classes.
`exterior_search.json` has exactly two positive parent classes:

|Parent order/rank|Child order|Child profile (H,start,deletion,F)|
|---|---:|---|
|10/685|11|(4277,206,207,4278)|
|13/1047|14|(4277,920,921,4278)|

The first has no eligible matching-one exterior clone. The second has
nine; their H(T_2) values are8477 (three),8459,8651,8189 (three),12095.
None is11735=g_2. Cloning the unique in-port is excluded by its growing
indegree at the ready root. This exhausts the exterior alternatives.

`offset_fibres.json` further checks offsets -4,-2,0,2,4 in 3h_m+d:
seven ready base classes and34 eligible clone queries, with no coefficient
match. The clone-mark constant B(3729) also differs from all five required
constants. The same proof therefore excludes one-ready closure in each
of these five fixed-uniform-clone lanes.

**Boundary.** This is not a nonexistence theorem for the H values themselves,
for individual18979 returns, or for two-ready/growing/split-module carriers.
The searched two-ready core bank through9 is only a bounded test. No
complete two-ready exclusion is claimed.

## SF4. Positive exterior families and a new ready-value gap

SF2 certifies five actual strong, one-ready families sharing H(T_1)=4277:

|Growth vertex in saved14-vertex seed|H(T_m), m>=1|
|---:|---|
|4|2100*2^m+77|
|7|2091*2^m+95|
|8|2187*2^m-97|
|9|1956*2^m+365|
|12|3909*2^m-3541|

The seed rows, clone addresses, root13, five Q4 phases at both queried
ports, and a full eighteen-polynomial seed observer at ports(13,8) are
in `family_certificate.json`. The third clone sample is an independent
literal check; the all-size conclusion is SF2. These are explicit new
sources, not source-preserving replacements for the saved18979 tournament.

Two literal births from the seed give the strong count12833=3*4277+2.
However, **no tournament at any order with H=12833 has a one-ready root**.
For the strong case the parent bound is6417, and f(18)=9375, so parents
through17 suffice. `next_fibres.json` exhausts1346268 parent classes and
finds none. Weak sources cannot help: the proved SCC identity
[O1](../ready_return_pentagon_20261103/OBSTRUCTION.md#o1-weak-sources-do-not-repair-a-positive-h-ready-obstruction)
reduces any one-ready source of H>1 to a strong source of the same H.

The adjacent12835=3*4277+4 has exactly one strong one-ready rooted class,
with parent order12/rank2611. Its two matching-one exterior clones have
coefficients (4602,3631) and (9105,-5375), neither matching any of the five
desired +4 continuations. This is an interface gap, not a spectrum gap.

Moving the live mark back into the clone module does not repair these
continuations either. The same binary B recurrence gives:

|Target coefficient|B(coefficient)|Required +2 constant|Required +4 constant|
|---:|---:|---:|---:|
|6300|-1403|233|235|
|6273|-3161|287|289|
|6561|-3421|-289|-287|
|5868|-1267|1097|1099|
|11727|-11161|-10621|-10619|

Every row disagrees. Combine this table with SF1--SF2 and the complete
12833/12835 exterior-root fibres: neither ternary continuation of any
of the five families can return to one-ready sources infinitely often
within a fixed uniform-clone strong family, even with changing marks.
The eventual-readiness quantifier matters: a clone mark could fail at
m=1 and become ready at m>=2, so the empty12833 base fibre alone would
not have justified this stronger conclusion.

## SF5. An actual return at the first stable clone size

The uniform obstruction does not prevent size-dependent returns. At m=2
the named repaired source has H=3911 and F_root=3912. Its literal two-birth
reset gives a strong16-vertex source with H=11735. The complete one-ready
fibre at11735 has exactly two strong rooted classes, of orders12 and14:
their parents have (order,rank)=(11,2277),(13,2056). For11737 the five
parent addresses are (12,2616),(12,2668),(13,2083),(13,2231),(13,2965).
Parent orders through17 suffice, since5869<f(18)=9375. `first_ready.py`
uses the complete graded grammar and checks all seven positives literally.

The same script constructs a return of the actual11735 child to the
order14 one-ready target. It deletes exactly old names2 and5, retaining
the two clone addresses31/31.1, old ports6/21/31, and both reset vertices.
It then reverses11 exterior arcs, using BC2's three-helper common-frame
compiler for each:22 legal native complements,33 helper births and33
helper deletions. Every completed arc repair is strong; peak repair order
is17. The live root's whole star and common guard frame are retained.
Helper-expanded intermediate strongness is not asserted.

This uses no physical relabelling: a target-to-retained-address map designs
the final Boolean relation on the existing names, and every change is
implemented by its saved legal guards. One literal birth realizes the
initial transitive clone, with its address map recorded. Actual starting,
returned and next source rows, root/spare phases and native clocks remain
in `first_ready.json`. `audit_first_ready.py` independently reconstructs
the births, matrix complements, deletions and Hamilton counts.

Two old deletions are **necessary** for a strong one-ready return from this
16-vertex source: every strong one-ready endpoint of H11735 has at most14
vertices, so retaining all16 old names as distinct vertices is impossible,
regardless of extra helpers. The construction attains this loss bound;
22 is not claimed to be the minimum native-event cost.

The next literal two-birth step gives35207=3*11735+2. Thus we obtain

    3911 -> 11735 --(actual neutral ready return)--> 11735 -> 35207.

This is the first stable clone size on the previously repaired source
branch. It is not a return of the separate18979 source, a uniform family
identity, or a ready return at35207. Those distinctions prevent a single
positive program from being promoted to closed ternary induction.

## FG1. A source-word recurrence and the two-feedback test

Fix feedback vertices w,p and a transitive chain. A chain letter
l in{0,1,2,3} records the Boolean arcs w->v and p->v as its two bits.
Retain the feedback edge as a separate bit. This is the older
[feedback-gap construction](../feedback_word_hamilton_spectrum_20261020/proofs.md)
specialized to two actual feedback vertices.

A Hamilton path orders w,p in one of two ways. Its remaining chain
vertices belong to three increasing gaps. For each of the first two
gaps retain three states: empty, nonempty with a valid last-to-next-port
arc, or nonempty with an invalid such arc. Retain empty/nonempty for the
last gap. There are 2*3*3*2=36 states. On reading a chain vertex, append it
to exactly one gap. A first insertion into a gap must respect the preceding
port; the latest insertion updates the right-boundary test. At termination
both right-boundary tests must pass; an empty middle gap requires the
actual feedback arc. Each accepted assignment is exactly one path.

This proves the nonnegative integer four-letter transfer in
`two_feedback.cpp`. It retains the actual source word. Independent subset
DP checks all2730 words of lengths0..5 and both feedback edges.

For a first-chain clone beaten only by w, the prefix is1^m. The requested
coefficients demand H(1 tail)=4277 and H(11 tail)=11735. A strong base has
order<=16 because f(17)=5625>4277, so its tail length is at most13. The
exhaustive search visits2589833 prefixes and prunes1745718 once either
required value is exceeded for both edge choices. Pruning is valid because
canonical insertion injects old Hamilton paths into the enlarged tournament.
There are zero hits; the100-hit safety cap is never reached.

Thus even without requesting readiness, this particular two-feedback
uniform-prefix architecture cannot cancel the1182 constant debt at any
strong base size. This does not cover a clone in an arbitrary chain
position, general split incidences, or three feedback vertices.

## What this establishes about the main spectrum question

The all-order impossibility of7 and21 is inherited, most recently audited
in [the exclusion packet](../exclusion_7_21_20261103/README.md).
The remaining statement is that every other positive **odd** integer is
realizable. The prime/7p reduction is
[THM-4094, matching-deficit and two-prime-lane completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md).
Neither faster class counting nor a ready-return obstruction proves that
existence assertion. The constructive successes here are the five infinite
families, an exact all-N census recurrence implementation, and a sharper
description of which uniform return architectures must be enlarged.

# Positive occupancy responses and finite families of bounded-defect programs

**PROVED; independent root-lane proof review PASS.** The objects are actual
rooted tournaments, optionally carrying a specified finite sequence of birth,
deletion and arc-reversal operations. All class counts below are ordinary
rooted isomorphism classes. They do not forget the root or the operation
history. The fixed branching condition is d=|N^-(r)|>=2.

The inherited [fixed-defect proof](../defect/BOUNDED_TEMPLATES.md) gives finite
transitive-module templates. Here their path responses have an exact positive
signature, their finite search threshold becomes explicit, weak intermediate
graphs are included by a canonical SCC composition, and fixed finite programs
receive their own family theorem. Chan--Pak source-return closure remains open.

## 1. The positive occupancy signature

Fix a finite rooted skeleton. Replace specified vertices by nonempty transitive
modules TT_(k_1),...,TT_(k_s), retaining all external incidence. Suppose these
modules lie in singleton intervals of the residual SCC chain
R=T-({r} union I), where I=N^-(r). Every k_i is a positive integer.

For each Hamilton path, or ordered spanning cover by a fixed number of
nonempty paths, erase the vertices in those modules. Retain the order of all
remaining vertices, the ordered path-bin boundaries, and which of the modules
occupy each intervening run. Within one run the module order is forced by the
residual SCC order, and within a module its vertices appear in their unique
transitive order.

If module i occupies j_i distinct named runs, assigning its k_i vertices to
those runs is exactly a surjection onto j_i named bins. The resulting count is

    L(k_i,j_i)=j_i! S(k_i,j_i).

The assignments for different modules are independent: every external
incidence has already been checked in the retained finite run pattern. Each
path object determines its pattern and assignments uniquely. Consequently
every such response has a finite expansion

    response(k_1,...,k_s)
       = sum_j A_j product_i L(k_i,j_i),       A_j nonnegative integers. (1)

All j_i>=1 because every module is nonempty and every vertex is covered. Empty
modules must be split into separate skeleton cases. This is an exact finite
combinatorial enumeration, not fitting the response at sampled lengths.

The positive defect identity is

    b = sum_(a in I) start_(T-r)(a)
          + sum_(a<c in I) [W_ac(T-r)+W_ca(T-r)].                 (2)

Thus b itself has the nonnegative signature (1). Here W_ac counts two ordered
nonempty spanning paths, the first ending at a and the second starting at c.
The ordinary Hamilton response uses one unconstrained path. Arbitrary fixed
endpoint constraints among retained vertices can be included in the same
finite run-pattern enumeration.

## 2. Only d runs are possible in a positive defect object

Between two visits to a particular residual TT module in the SAME path, there
must be an inport. A path entirely inside R cannot leave its position in the
ordered SCC chain and return. Inserting other retained residual vertices
therefore cannot circumvent the required inport.

For a Hamilton path starting at a specified a in I, the initial inport is not
between module visits. At most d-1 other inports separate visits, so there are
at most d visits. For W_ac, if both paths visit the module, j total visits
require j-2 separating inports. The specified end a of the first path and
start c of the second are two additional, distinct inports. Hence j<=d.
If only one path visits the module, the bound is even stronger.

Therefore every exponent index in the defect signature satisfies

    1 <= j_i <= d.                                             (3)

For the ordinary Hamilton count the corresponding bound is j_i<=d+1. Indeed
T-I is r followed by the residual SCC chain, so it also forbids a return to
the module without an inport. More generally an unconstrained ordered cover
by m paths has j_i<=d+m.

Using the elementary identity

    L(k,j)=sum_(h=0)^j (-1)^(j-h) binom(j,h) h^k,                (4)

all these responses are exact multivariate exponential polynomials. With all
but one module length fixed, H satisfies a constant-coefficient recurrence
whose characteristic polynomial divides product_(j=1)^(d+1)(X-j). This is
for positive module lengths; zero-length boundary cases were separated.
The positive basis (1), rather than the signed form (4), certifies monotonicity
and the threshold below.

## 3. A universal logarithmic finite grid

For 2<=j<=k,

    L(k,j) >= 2^k-2.                                          (5)

One elementary proof uses the recurrence
L(k,j)=j[L(k-1,j)+L(k-1,j-1)]. The j2 boundary is equality. The jk boundary is
k!>=2^k-2. For the remaining indices, induction applies to both summands;
the first nontrivial interior case is k4,j3. This proves all k>=2.

Fix a positive defect b and put

    K_b = floor(log2(b+2))+1.                                 (6)

The sharp inherited bound b>=2^d-2 gives d<K_b. Fix all other module lengths
arbitrarily. After those lengths are substituted into (1), the coefficients
of the remaining L(k,j) are still nonnegative integers. If k>=K_b and any
j>=2 coefficient is positive, then (3)--(5) give

    defect >=2^k-2>b.

It follows that an accepted length k>=K_b has only a j1 contribution. The
actual defect is then unchanged for EVERY positive length in that coordinate,
not merely all sufficiently large lengths. Strongness, when imposed, is also
unchanged by nonempty module inflation.

Thus the complete exact-b test for any finite template uses only the grid

    k_i in {1,...,K_b}, with K_b representing every length >=K_b. (7)

This threshold is uniform in all other lengths and independent of the size
of the residual kernel. Applying it successively to the coordinates proves
the finite disjoint Cartesian-box classification without a DFA. Moreover
each accepted box can be extended down to length1 in all its free coordinates
simultaneously. Those enlarged boxes can overlap; the disjoint grid boxes in
(7) remain the convenient canonical counting convention.

The two-inport exceptional-singleton estimate also simplifies. From

    m0+binom(m10,2)+binom(m01,2)<=floor(b/2)

and m-binom(m,2)<=1 for every nonnegative integer m, one gets

    m0+m10+m01<=floor(b/2)+2.                                 (8)

For d>=3 the existing bound is floor(b/2). The zero vertices in (8) are all
zero-contact vertices in a strong source, or the charged early ones in the
weak extension below.

## 4. Weak graphs have a canonical composition law

Every rooted tournament has a unique SCC decomposition

    T = A -> C -> B,

where C is the strong component containing r, and A,B may be empty. A
singleton counts as a strong component. Write H(empty)=1 and Q_empty(2)=1.
The root response factors exactly as

    H(T)=H(A)H(C)H(B),
    F_r(T)=Q_A(2) F_r(C) H(B).                                (9)

For the second identity, every B vertex must be in the path-bin AFTER r.
The A vertices split arbitrarily between the two bins and contribute Q_A(2).
The C restriction is precisely the original root response. This is a
bijection of actual path objects, respecting their ordered bins.

Consequently

    b(T,r)=H(B)[Q_A(2)b(C,r)+(Q_A(2)-H(A))H(C)].               (10)

If A is empty, d(T,r)=d(C,r), so a branching root comes from a strong
branching source C, followed by an arbitrary tournament B of fixed H.
If A is nonempty, Q_A(2)>=2H(A), and b(C,r)>=0, whence

    H(T)<=b(T,r).                                             (11)

Also |A|<=d<=floor(log2(b+2)). Thus only finitely many A and rooted strong C
occur in this second case: H(C)<=b, and the fixed-H strong-kernel bound makes
their order bounded. B again has a specified fixed Hamilton count.

Let S_b(x) count strong rooted branching sources of exact defect b, and let
A_h(x) count all ordinary unrooted tournaments of exact Hamilton count h.
The latter fixed-H series is the earlier finite SCC-kernel/singleton-gap
grammar. Define D(A,C)=Q_A(2)F_r(C)-H(A)H(C). The exact weak rooted series is

    W_b(x)=sum_(t|b) S_(b/t)(x) A_t(x)
       + sum_(A nonempty,C rooted strong,d>=2,D(A,C)|b)
             x^(|A|+|C|) A_(b/D(A,C))(x).                    (12)

The second sum is finite by (11). Only d>=2 pairs are included. Distinct
terms describe distinct rooted SCC decompositions, so no hidden root-forgetting
or repeated-component quotient occurs. Both factors in the first sum already
have denominators which are powers of 1-x. This proves finite inflation
families and eventual polynomial rooted class counts for ALL rooted tournaments
of fixed branching defect, whether strong or weak.

Combining (10) with the independently reviewed
[dyadic prefix theorem](../arithmetic/DYADIC_DEFECT.md) also extends its scalar
sparsity result to weak graphs. With A empty, H(T)=H(B)H(C) merely scales a
strong dyadic ray by a divisor of b. With A nonempty, H(T)<=b gives finitely
many values. Therefore ALL rooted tournaments with d>=2 and defect at most
a fixed positive B have only O_B(log(X+2)) distinct Hamilton values <=X.
Weak sink padding does not evade the bounded-budget obstruction.

There is a local form useful for operations. In the residual SCC chain, a
zero-contact singleton before some later vertex y beating an inport is charged
by the same tagged five-cycle used in the strong proof. Only the suffix AFTER
the last nonzero-contact residual SCC can contain unboundedly many uncharged
zero singletons. Its all-zero runs, together with the earlier all-one runs,
are homogeneous TT slots separated by finitely many retained objects. One
finite divider records the last-nonzero boundary, including an adjacent
all-one/all-zero pair of runs. The threshold (6) applies to both kinds.

## 5. Complete weak branching enumerators through defect six

Combining (12) with the complete strong classification gives the following
new ordinary ROOTED series. They include strong and weak sources together:

|b|W_b(x)|Eventual count at order n|
|---:|---|---|
|2|(x^4+2x^5+2x^6)/(1-x)^2|5n-21, n>=5|
|3|(x^3+x^4)/(1-x)|2, n>=4|
|4|(x^4+5x^6+10x^7)/(1-x)^2|16n-88, n>=6|
|5|3x^5/(1-x)|3, n>=5|
|6|(4x^5-2x^6+12x^7+12x^8-16x^9)/(1-x)^3|5n^2-29n-12, n>=7|

There are no branching sources of defect1. For b2,4 the only new freedom is
transitive sink padding of a strong source. At b6 there is also a strong
defect2 source followed by an H3 tournament, whose series is x^3/(1-x)^2.

To check every nonempty-A contribution for b<=6, (11) restricts C to a
singleton, C3 or the unique strong H5 tournament S4. If |A|=3 the minimum D
already exceeds6. If |A|=2, the only surviving case is A=TT2 and C={r},
with D3. If |A|=1, the only surviving cases are C3 with its defect0 root,
giving D3, or the indegree1 defect0 root of S4, giving D5. The other indegree1
root of S4 has defect1 and gives D7. Thus the two extra b3 families are

    TT2 -> r -> TT_l,
    singleton -> (C3,r) -> TT_l,

and the one extra b5 family is singleton -> (S4,r_zero) -> TT_l.
All have l>=0. No H2 tournament exists, so D3 contributes nothing at b6.
These account for all terms of the displayed table.

### Sharp growth degrees on infinite families of budgets

Let Omega(u) count prime factors with multiplicity and let u be the odd part
of b. The dyadic prefix grammar gives a substantially sharper degree bound
than the general exceptional-vertex estimate. A nontrivial prefix core has
odd H>=3, and the product of these H values divides b. Hence at most Omega(u)
such cores occur. A prefix with c cores has c+1 free gaps and eventual class
degree c. Therefore the strong fixed-b count has eventual degree at most
Omega(u); an odd-b strong stratum is finite.

For even b, (12) adds a suffix with fixed odd Hamilton count t|b. Its core
count is at most Omega(t). The total pole order is at most Omega(u)+2, so
the weak count has eventual degree at most Omega(u)+1. For odd b>1 the
strong factor is finite, and every nonzero defect factor b/t is an odd
integer at least3. Thus the weak degree is at most Omega(b)-1. Its Hamilton
spectrum is finite even if the weak class count grows: the strong odd-defect
part is finite, and the nonempty-A part has H<=b.
For the second sum of (12), the same observation applies to D(A,C): it is
the defect of A->C with total indegree at least2, and is odd when b is odd,
so D(A,C)>=3.

These bounds are sharp on an explicit infinite set of budgets. For

    b=2^v 3^a 5^c,       v>=1, a,c>=0,

take a prefix P with a triangle cores and c copies of the strong H5 core S4,
in any order, separated by transitive gaps of total length m. Let the tail
W_v be TT_(v+2), with its FIRST and LAST vertices as the two inports. The
v interior vertices lie between them. One ordered pair-cover direction is
impossible; the other assigns each interior vertex independently to either
bin. Thus its capacity is C=2^v. The last tail SCC is an inport, so the
prefix construction is strong. For the corresponding finite rooted kernel,

    q0=1, E=2^(v+1)+1, C=2^v.

The prefix cores have (H,Q(2))=(3,12) and (5,32), respectively. Therefore

    H=(2^(v+1)+1) 12^a 32^c 2^m + (1-2^v)3^a5^c,
    n=m+3a+4c+v+3,             b=2^v3^a5^c.                 (13)

Use m>=0 if a+c>0 and m>=1 otherwise. The one scalar ray (13) has exactly

    binom(a+c,a) binom(m+a+c,a+c)

rooted classes in this specified construction, with rooted automorphism
group C3^a. The first binomial counts the distinct ordered core lists and
the second counts their gap distributions. It gives positive growth of
degree a+c; the general upper bound proves that the COMPLETE strong
stratum has exactly this eventual degree. Adding an arbitrary transitive
sink suffix gives degree a+c+1 for the complete weak stratum.

For odd b=3^a5^c>1, take one of the finite D3 or D5 leading-root families
already identified above, then a suffix of Hamilton count b/3 or b/5 made
of triangle and S4 cores. This achieves degree a+c-1, again matching the
upper bound. Consequently the exact degree table on these infinite budget
families is

|Budget|Strong eventual degree|Weak eventual degree|
|---|---:|---:|
|2^v3^a5^c, v>=1|a+c|a+c+1|
|3^a5^c>1|finite stratum|a+c-1|

This gives genuine all-order enumeration transitions: even budgets6 and10
have linear strong growth; budgets18,30 and50 have quadratic strong growth;
each gains one degree when weak padding is included. These statements concern
degree, not a claimed full leading coefficient or a count of unrooted
tournaments after forgetting r.

There is also an exact bridge back to the prime Chan--Pak frontier. For an odd
prime p, let P_p(x) be the finite polynomial counting ordinary strong cores
with H=p by order. At b=2p, any nontrivial neutral prefix has total H=p and
contains exactly one such core. Its remaining capacity is2, so its bounded
tail is one of the complete five pentagon tails. Hence the entire part of
S_(2p)(x) with a pole of order2 is

    P_p(x) x^3(1+2x+2x^2)/(1-x)^2.                          (14)

The remaining terms have at most a simple pole. Consequently the eventual
linear slope is exactly 5P_p(1). It is positive if and only if p has a
tournament realization: an ordinary tournament of prime H has exactly one
nontrivial strong component, itself of H=p. The excluded value p7 gives no
linear growth at b14, while p3 gives the slope5 at b6. This is an equivalent
enumeration frontier, not a proof that every other prime is realizable.

For cumulative budgets there is a fully explicit degree law at EVERY B>=2.
The number of strong rooted branching sources of defect at most B has exact
eventual polynomial degree

    floor(log_3(B/2));

the weak count has degree one higher. Indeed c prefix cores force
b>=2*3^c. The triangle-prefix construction at b=2*3^c attains the upper
bound. In weak sources, the combined prefix and suffix core count obeys
the same product bound, and the suffix contributes one extra free gap;
transitive sink padding attains it. Thus the degree transitions occur
exactly at budgets2,6,18,54,... . This determines the growth degree, not all
polynomial coefficients or the unrestricted A000568 sequence.

## 6. A finite-family theorem for bounded operation programs

Fix a finite word of operation types: one-vertex births, one-vertex deletions,
and individual arc reversals. A program includes the chosen vertices/arcs and
newborn attachments, not just its sequence of scalar H values. A specified
root persists through the program. Assume every intermediate rooted graph has
d>=2 and defect at most a fixed B. Strongness may be required at selected
stages, but weak stages are allowed.

**Program theorem.** All such programs, up to isomorphisms preserving the root,
step order and operation incidences, are a finite union of computable
transitive-inflation boxes. The free clones persist untouched throughout the
program. Each intermediate H and each fixed path-bin response has the exact
multivariate signature (1). Counts by initial order are eventually polynomial.

Proof. At each stage retain its root, inports, nontrivial residual cores, charged
exceptional singletons, and the finite last-nonzero divider. Retain also every
vertex that is ever born, deleted, or incident to a reversed arc. The union of
these sets over finitely many stages is bounded by B and the program length.

Every remaining vertex persists and is a singleton of the residual chain at
every stage. It has ordinary all-one or permitted terminal all-zero contact.
The orientation between any two such vertices is never changed, so they share
one intrinsic transitive order throughout. Against a retained residual vertex
present at a given stage, their incidence changes at most once along this
order: they lie before or after its residual SCC. Against an inport, their
ordinary contacts can change only at the finite last-nonzero divider. The
root dominates them all. Taking maximal consecutive intervals on which all
these incidence patterns are constant gives the canonical common refinement:
each resulting TT slot is homogeneous across EVERY stage.

There are therefore finitely many operation skeletons. Replacing a positive
slot by any positive TT length respects every operation incidence and leaves
every selected strongness condition unchanged. At each stage its defect has
the positive signature with j<=d_stage. Applying (6) with B simultaneously to
all stages proves membership stabilizes on the finite grid K_B. Exact chosen
defect values can be imposed as well as the upper bounds. The finite skeleton
isomorphism group fixes every ordered slot, so it introduces no periodic
length conditions. This proves the theorem.

This theorem includes homogeneous sink-padding stages such as those in the
uniform reset, provided the stated branching defect bound covers their actual
profiles. It does not require them to remain strong during the padding.

The intermediate bound is essential. Start with r->TT_k->TT2->r, of defect2.
Birth a temporary vertex v with v->r and one fixed inport a->v, give v an
arbitrary binary contact word along TT_k, then delete v. Every intermediate
graph is strong, and the endpoint sources are identical, but these are 2^k
distinct rooted programs. Rooted automorphisms fix the ordered TT_k chain.
An endpoint-only bound therefore does not yield the claimed finite-box
program classification; the intermediate defects are unbounded.

Native path-complement steps are outside this theorem: a path can change
unboundedly many arcs and brings a separate path-position carrier. Nor does
the theorem automatically solve conditions such as H_out=3H_in+2. Those
become explicit exponential-polynomial equations in the slot lengths, and
their solution sets require a separate arithmetic argument.

## 7. Reproduction, controls and relation to the stronger prefix theorem

Run `python -B response_families.py` from this directory. The independent
endpoint-subset DP checks every labeled tournament through order5, every
branching root, and every maximal all-one or terminal-zero slot. It verifies
nonnegative integral occupancy coefficients, the d/d+1 run bounds, the
explicit threshold, and the root-SCC factorization. A separate exact
root-fixed permutation quotient checks the weak rooted series through order5.

The receipt [response_families.json](response_families.json) records 3,654
rooted factorization checks, 1,688 slots and 7,072 literal inflation queries,
including 784 weak-source slots. The 360 nonconstant strong-slot cases are
hostile controls: their defects grow under inflation, so a homogeneous TT
module alone does not imply a free fixed-defect direction.
Another84 literal controls verify the explicit sharp-degree family (13),
including both leading and trailing gap placements and the nontransitive
prefix cores. Normal and optimized Python runs both pass.

For example rows(24,1,3,6,14), root2, slot{0} has Hamilton occupancy
(11,17,4) and defect occupancy(9,2). Its H response contains a genuine
three-run contribution, and its defect is 9+2(2^k-2). This demonstrates why
the unrestricted positive response theorem must allow bases beyond2.

The concurrent [prefix-tail classification](../arithmetic/DYADIC_DEFECT.md)
strengthens the ordinary fixed-defect
description further: late all-one singletons after any bad contact are also
bounded, so every unbounded part is an initial all-inport-beating prefix.
Its dyadic H rays are compatible with, and more specific than, (1). The
occupancy signature here remains useful for bounded finite lengths, arbitrary
path-bin observables, common refinements across programs, and hostile
nonconstant templates. No claim of classifying the unrestricted metagraph
or proving recursive source-return closure is made.

Independent review: the root lane checked the full occupancy, threshold,
weak-composition and finite-program proofs and reported PASS. Its separate
`audit/weak_factor.py` implementation checked 5,405 root-SCC factorizations
over all rooted labeled tournaments through order5, including the same 3,654
branching roots and 1,406 nonempty-leading-region bounds. These are independent
implementations supplementing the structural proofs.
The arithmetic lane independently reviewed the same six sections and also
reported PASS, including the uniform threshold and the finite-program quotient.
It separately reviewed the sharp growth-degree and prime-slope corollaries
and reported PASS.

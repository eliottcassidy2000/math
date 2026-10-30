# Suffix-entry domination and endpoint controllers

**PROVED scoped, same-agent proof review:** for every fixed positive valuation
controller, the largest ordinary five-adic content mean over all finite
positive entries occurs at a suffix of that controller. For a controller
with one exceptional letter, an endpoint marking maximizes every truncated
layer mean simultaneously. **FINITE-EXACT:** all controllers of integer clock
at most14 satisfy the asymptotic 5/8 bound for every entry. **OPEN:** the
all-controller bound, a uniform changing-controller error, finite-word debt,
and convergence of every positive Collatz source.

## 1. Inherited mechanism and preserved data

Use the exact definitions and ordinary-limit theorem of
[Session114, periodic content rates](../session114_continuation_symmetry_20261030/collatz/PERIODIC_CONTENT_RATES.md).
For a positive valuation word x,

    F_x(X)=(3^len(x) X+B_x)/2^A_x,  A_x=sum(x),
    c(x)=v_5(gcd(|3^len(x)-2^A_x|,B_x)).

Concatenation is chronological: F_(xy)=F_y F_x. For fixed entry u and
nonempty controller w of length d, let v_j=w[:j], u_j=w[j:], 1<=j<=d.
The charge after entry, for N controller periods, is

    S(N;u,w)=sum_(k=0)^(N-1) sum_(j=1)^d c(u w^k v_j).

No prefix of u is charged twice or reset; A_w remains the original integer
clock. The inherited theorem proves the ordinary mean R(u,w), with error
O_(u,w)(log(N+2)); it is not merely a Haar-average computation.

Write M=F_w and let H_e=<M mod5^e> inside G_e=AGL_1(Z/5^e).
Its order is T_e=t 5^max(e-r,0), with t,r given in that theorem.
For a phase j define

    D_j(u)=sup{e>=1 : F_u F_(v_j) belongs to H_e},

using zero for the empty set and infinity for membership at all levels.
The inherited exact formula says

    R(u,w)=sum_j sum_(1<=e<=D_j(u)) 1/T_e.             (1)

The geometric tail converges. There is exactly one successful residue of
period index k modulo T_e at every admitted level; its offset is not recorded
by D_j. This distinction will remain essential.

## 2. The finite suffix theorem — PROVED

**Theorem.** For every positive controller w and every finite positive entry
u, there is a j in {1,...,d} such that

    D_i(u)<=D_i(u_j) for every phase i.                (2)

Thus the same suffix entry dominates every truncated mean, and

    sup_u R(u,w)=max_(1<=j<=d) R(u_j,w).              (3)

The empty entry is allowed. The supremum is attained by a legal finite
positive word, rather than an artificial inverse affine map.

**Proof.** The actual prefix is F_(v_i) M^k F_u. It can be identity modulo
5^e for some k exactly when F_u lies in the right coset

    C_(i,e)=H_e F_(v_i)^(-1).

The compatible preimages of these cosets at greater precision are nested.
Since w=v_j u_j, F_(u_j)=M F_(v_j)^(-1), so this suffix entry lies in
C_(j,e) at every level e.

Choose j maximizing the finitely many depths D_i(u). If e<=D_i(u), then
also e<=D_j(u). The same entry F_u lies in both right cosets C_(i,e) and
C_(j,e). Right cosets of a subgroup are equal or disjoint, so they are
equal, and F_(u_j) lies in C_(i,e). This proves (2). Positive layer weights
in (1) prove the upper bound in (3); the suffixes themselves are allowed
entries, giving equality. The proof includes infinite maximal depth. QED.

More generally, the argument works for any nested subgroup chain and a
finite family of compatible cosets, provided each labelled coset has an
actual representative at every depth. It does not require normal subgroups
or a commutative ambient group.

### A finite ultrametric carrier

Set D_ij=D_j(u_i). Then D_ii=infinity, D_ij=D_ji, and

    D_ik >= min(D_ij,D_jk).

Indeed, D_ij>=e exactly means that the two phase cosets agree modulo5^e.
For each e this is an equivalence relation on the d phase labels. The
relations refine as e grows, so the whole mean optimization is a finite
nested partition tree with edge lengths determined by arithmetic valuations.
For a suffix entry, its layer-e mean is its class size divided by T_e.
Thus one can replace an unbounded entry search by d weighted tree rows.
Duplicate leaves, including infinite-depth equivalences, are retained or
merged with their multiplicities; phase weights must not be discarded.

### Cyclic markings

For w=v_j u_j, put z_j=u_j v_j. The identity

    u_j w^N=z_j^N u_j

is an identity of actual words. The inherited fixed-exit estimate makes the
terminal u_j charge O(log(N+2)), while the initial u_j has fixed charge.
Dividing by N shows R(u_j,w)=R(empty,z_j). Thus the right side of (3) is
also the maximum of the empty-entry means over cyclic rotations of w.
A rotation is used to bound a mean; it does not change the original source
word or its finite-time return address.

## 3. Single-outlier endpoint principle — PROVED

Let a,b be distinct positive integers, L>=1, w=(a)^L(b), and G=F_(a).
For each e let ell_e be the least positive integer k such that
G^k mod5^e belongs to H_e. Such an integer exists because G_e is finite.
The preimage of H_e under the homomorphism k -> G^k is precisely ell_e Z.

The phase cosets are represented by G^i for i=1,...,L and by M for the last
phase; the latter has the same coset as G^0. Hence two of the L+1 phase
labels i,j agree at depth e exactly when ell_e divides i-j. The class size
at label i is

    1+floor(i/ell_e)+floor((L-i)/ell_e)
       <=1+floor(L/ell_e).

Both endpoints i=0,L attain this upper bound at every precision. The label
0 is the empty-entry coset, so (2) gives the all-entry theorem

    R(u,(a^L,b)) <= R(empty,(a^L,b))
                 = sum_(e>=1) (1+floor(L/ell_e))/T_e.          (4)

The same domination holds before summing beyond any finite precision.
This is an all-L, all-a,b reduction, not yet an all-L 5/8 clock inequality.
It eliminates entry optimization without claiming that ell_e or T_e obey
uniform estimates when the controller changes.

### Complete maximizing-entry grammar for the same family

The incoming Session115
[all-reset classification](../session115_odd_eighth_classifications_20261030/collatz/GENERAL_RESET.md)
already proves the a=1 case, with a stronger explicit formula for ell_e and
its mixed radices. The endpoint argument above independently recovers that
case and extends the group-theoretic statement to every repeated letter a.
Its equality argument extends too: all maximal entries are exactly those
satisfying

    (a)^j u=w^B,  0<=j<=L, B>=0,
    j mod ell_e <= L mod ell_e for every e>=1.         (4a)

For B=0 only j=0,u empty is possible. For B>=1 this is an actual positive
suffix followed by controller powers, so no inverse-map realization is needed.

**Proof.** Every layer deficit from the endpoint maximum is nonnegative.
If all depths of u are finite, it misses high layers which the endpoint
continues to occupy, so it cannot attain the same full mean. An infinite-
depth phase commutes with w as a rational affine map. The inherited faithful
word centralizer then forces its completed word to be a power of w, since
w=(a)^L(b), a!=b, is primitive. This gives the first equation in(4a).
For such an entry the layer-e class size at j is the floor expression above;
its equality condition is exactly the residue inequality in(4a).

The intersection moduli ell_e form a divisibility chain and tend to infinity.
For otherwise a fixed nonzero power of F_a would commute with F_w at every
precision. Their rational fixed points differ when a!=b (the shift created
by the last letter is nonzero), so this is impossible. A direct check uses
F_a's fixed point 1/(2^a-3), which F_b fixes only if b=a.

Delete repetitions from1,ell_1,ell_2,...; write j,L in this mixed radix.
The residue inequalities are equivalent to coordinatewise digit inequalities.
Thus if L_i are the mixed-radix digits, the exact number of maximizing phase
labels is product_i(L_i+1). The incoming a=1 theorem supplies explicit
2/5/10 radix classifications; this general version does not assume those
same radices for an arbitrary repeated letter. It preserves the completion
word as well as the digit count.

## 4. Hostiles: neither entry invariance nor finite-time domination

For w=(1,3,1,1,3), A_w=9, the empty-entry mean is5/8. The five suffix-entry
means, in phase order, are

    5/8, 369/200, 73/40, 369/200, 5/8.

The depth matrix is

    inf  0    0    0    0
    0    inf  2    3    0
    0    2    inf  2    0
    0    3    2    inf  0
    0    0    0    0    inf.

The maximizing suffixes see a three-leaf cluster at level2 and a two-leaf
cluster at level3. Empty entry misses it. Even a cyclic re-marking changes
this directional mean; a trace average would conceal the largest row.

Pointwise finite-time domination does not follow from (2): for controller
(1), empty entry has first charge c(1)=0, while entry(1) has first future
charge c(1,1)=1. Both phase depths are infinite. The missing coordinate is
the first-return offset, not an extra depth or an averaged source weight.

## 5. Exact finite certificate and limits

`python3 entry_audit.py` (also `python3 -O entry_audit.py`) computes all
16,383 positive controllers with A_w<=14, using the inherited exact rational
rate formula for every suffix. All satisfy

    sup_u lim_N S(N;u,w)/(N A_w)<=5/8.

This finite controller certificate has unrestricted actual entry and
unrestricted repetition length in the limiting statement. Exactly14 tested
controllers attain equality, all ones. The largest tested non-ones rate is
45/104, at a cyclic marking of (1^8,5). Among the controllers,3,305 have
primitive-root order stalls outside the inherited sufficient growth test.
These counts are finite classifications, not universal equality claims.

Controls additionally include8,160 componentwise entry comparisons,22,528
ultrametric triples, and1,120 single-outlier examples. An independent
letter-by-letter finite-group implementation enumerates every affine entry
modulo5 and25 for controllers of clock<=5:16,120 entry tests. Its maximum
truncated mean agrees with the suffix formula. A further3,240 actual-word tests check the complete maximizing-entry
grammar using literal cyclic-group intersections. The two hostiles above are
required controls, and failures remain active under Python -O.

## 6. An explicit entry-clock discrepancy bound — PROVED

The inherited O_(u,w)(log N) can be given a uniform, if coarse, constant
which keeps the actual entry clock. Let U=A_u, A=A_w, d=len(w), N>=1, and
let K be the least nonnegative integer satisfying

    5^K >= 2 * 3^(U+4A) * (U+N A+1).

Then, with R_max(w)=max_j R(u_j,w),

    S(N;u,w) <= N R_max(w)+d K.                     (5)

This holds for **every** entry, controller and N with exactly these
quantifiers. It does not erase U or claim a uniform constant as words vary.

**Proof of the height bound.** For any nonempty word x of clock B,
its integer coefficients obey

    0<B_x<=3^B-2^B,   |3^len(x)-2^B|<3^B.

The translation bound follows by induction from B_(xa)=3 B_x+2^B;
for a>=1, 3^(B+a)-2^(B+a) >= 3^(B+1)-2^(B+1).
For a phase j let H=F_(v_j u). If its commutator numerator

    C=(P_w-Q_w)B_(v_j u)-B_w(P_(v_j u)-Q_(v_j u))

is nonzero, any successful identity return modulo5^e forces5^e|C.
Its absolute value is less than2*3^(U+2A), so the actual content height
at that phase is below K, independently of k.

If C=0, the inherited faithful-word commuting theorem gives w=z^a,
v_j u=z^b, with a,b positive integers. After conjugating by F_(v_j),
the actual prefix has the same content height as F_z^(ak+b). For its root
of clock A_z, the geometric-sum/LTE formulas bound this height by

    4 A_z log_5(3)+log_5(ak+b).

Indeed the residue multiplier order is at most4, its initial valuation is
less than4 A_z log_5(3), and the remaining valuation is at most log_5(ak+b).
The multiplier-one case only improves the constant. Since A_z<=A,
a<=A, b<=U+A and k<=N-1, this is bounded by4A log_5(3)+log_5(U+NA),
again below K. Conjugation is by an invertible integral five-adic affine
map, so it preserves being identity modulo5^e at every e.

**Counting with the original period index.** No observed height exceeds K.
For each phase and each layer e<=K, the successful k form either no
progression or one residue class modulo T_e. Their number for0<=k<N is
at most N/T_e+1. Therefore

    S(N;u,w)<=N R_K(u,w)+d K<=N R_max(w)+d K.

This proves(5) without an exchange of an infinite time and precision limit.
In fact R(u,w) can replace R_max in the upper bound; the suffix theorem
makes the latter finite to optimize. QED.

For controllers of clock<=14, the exact certificate above yields the
concrete finite-time consequence

    S(N;u,w)<=(5/8)N A+d K,

for arbitrary entries and all N. For changing controllers, these overheads
must still be paid, and the displayed K can grow linearly with the already
accumulated entry clock. A useful next refinement would replace that cost
by a paid change of the actual commutator data, not silently discard it.

## 7. The next sufficient target

The all-controller asymptotic conjecture is now equivalent to a finite
weighted-row inequality for each controller:

    max_i sum_(e>=1) |[i]_e|/T_e <= (5/8) A_w.       OPEN

The objects are actual prefix phase classes, not arbitrary ultrametric
trees: their weights, cyclic order and integer clock arise from positive
valuation letters. A generic tree inequality cannot enforce those costs.
Even a proof of this displayed target would still need uniform discrepancy
and entry-debt control for changing controllers before it could imply the
finite-word content bound or any universal Collatz conclusion.

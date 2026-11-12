# A finite budget for exceptional unselected twin-block lengths

**PROVED (scoped).** Fixing the exterior and selected set produces either a
freely inflatable port or only finitely many successful positive lengths.
This note supplies a structural bound on that finite number, an exact finite
certificate, and a physical example where arbitrarily accurate local
agreement conceals literal failure. It does not classify all marked cores.

## 0. Inheritance and novelty boundary

Let `(T,O)` be a tournament with a nonempty selected vertex set O, and choose
v outside O. Replace v by the transitive unselected module TT_m, m>=1, with
uniform exterior incidences. Let H(m) count Hamilton paths. Let W(m) count
ordered spanning covers by two **nonempty** directed paths P,Q such that
last(P) and first(Q) belong to O. Set δ(m)=H(m)−W(m). A unit seed has δ=1;
its positive multiplier U and the tail construction are inherited from
[the five-stress seed note](../../five_stress_connection_atlas_20261112/SEEDS.md).

The run/Stirling substitution is inherited from THM-3121
[path-cover-walk-content-substitution-kernel](../../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md)
and [Session260](../../session260_quasiperfect_chanpak_20261111/constructors/RESEARCH.md).
The exact Hamilton run degree ν+1 is already
[Session196, Theorem2](../../../../session196/spectral-moves/PROOF.md)
and [Session197](../../../../session197/ribbon-clone-matching/PROOF.md),
not new here. The contribution is the selected-endpoint tax bound and its
combination with an elementary exponential zero theorem to classify the
**number and certifiability of allowed lengths**. We include the old matching
proof so the parameter used by the new bound is transparent.

## 1. Return capacity and the marked tax

Put A=N+(v), B=N−(v), and let ν be the maximum matching number of the
bipartite graph with edge a b precisely when a in A, b in B and a→b.
These are disjoint directed-triangle return gates through v. Write

    R_H=ν+1,
    R=min(ν+2, 1+min(|A|,|B|)).                        (1)

For r>=1 put F_r(m)=r! S(m,r), the number of surjections from an m-element
set onto r ordered nonempty boxes. There are nonnegative integers h_r,w_r
such that, for every m>=1,

    H(m)=sum_(r=1)^R_H h_r F_r(m),  h_(R_H)>0,
    W(m)=sum_(r=1)^R w_r F_r(m).                       (2)

**Hamilton part, inherited.** Compress maximal module runs. Each gap between
two runs contains an A→B crossing. Distinct gaps give disjoint crossings,
hence r−1<=ν. Conversely alternate ν matching edges with ν+1 copies of v.
Every unused exterior vertex inserts into this directed word while retaining
the relative order of its existing letters; the distinct module copies stay
separated. This proves exact rank ν+1.

**Tax refinement.** A two-path cover with r module runs has at least r−2
internal gaps, since at most two paths contain runs. Selecting an A→B
crossing in each gap gives r−2<=ν. There is also an endpoint capacity bound.
If both paths contain module runs, P supplies an A vertex after every run
and a B vertex before all but its first run; Q supplies a B vertex before
every run and an A vertex after all but its last. These are distinct
exterior vertices, so at least r−1 of each type are needed. If only P
contains runs it supplies at least r A vertices and r−1 B vertices; if only
Q contains runs the reverse holds. In all cases r<=1+min(|A|,|B|).
Endpoint membership in O is essential to these boundary counts. The bound
is not asserted for arbitrary unrestricted two-path-cover observables.

Compressing each run leaves a finite skeleton, and each skeleton has exactly
F_r(m) lifts, as in the inherited substitution theorem. Thus (2) follows.
The Hamilton coefficients also fit the common bound R, since ν<=min(|A|,|B|).

If w_(ν+2)>0, the defect has leading exponential coefficient −w_(ν+2),
and is eventually negative. Such a port cannot be freely unit-preserving.
This statement applies only when that coefficient actually exists; (1) is
an upper bound, not a formula for the tax's exact rank.

## 2. The length-rigidity theorem

Set a_r=h_r−w_r, padding missing terms with zeros. Inclusion–exclusion gives

    f(m):=δ(m)−1=sum_(j=1)^R c_j j^m,
    c_j=sum_(r=j)^R (−1)^(r−j) binom(r,j)a_r − [j=1].   (3)

Let q be the number of nonzero c_j. Exactly one alternative holds:

1. **Free:** q=0, and δ(m)=1 for every positive integer m.
2. **Rigid:** q>0, and at most q−1<=R−1 positive integers m satisfy δ(m)=1.

Indeed the stronger zero bound holds on the entire real line after interpreting
j^x=exp(x log j), with multiplicities.

**Proof of the zero bound.** A nonzero sum of q exponentials with strictly
increasing real exponents has at most q−1 real zeros counted with
multiplicity. Induct on q. Divide by its smallest positive exponential.
The derivative is a nonzero sum of q−1 exponentials, hence has at most q−2
zeros. If the original quotient had q zeros, Rolle's theorem (including
multiplicity) would give at least q−1 derivative zeros on their finite
spanning interval, a contradiction. The q=1 case is immediate. This is an
elementary Chebyshev-system argument, not a new general zero theorem. ∎

**Finite certificate.** Success at **any R distinct positive lengths** proves
free stability. In particular, an unstable unit seed first fails at some
m<=R. There is no need to test n+1 or n+2 lengths once the endpoint and
matching bounds are retained.

For consecutive lengths one can instead use the exact recurrence with
annihilator product_(j=1)^R(E−j), E the forward shift. Or triangularly
recover the a_r from values at m=1,...,R, since F_r(m)=0 for r>m and
F_m(m)=m!. No zero-length datum is used.

**Delayed failure has a structural price.** If δ(1)=...=δ(L)=1 but the port
is unstable, then R>=L+1. Therefore

    |T|>=2L+1,   ν>=L−1.                             (4)

Long deceptive stability requires a larger exterior and many disjoint
triangle-return gates; it cannot be produced at arbitrary delay inside a
fixed small chart. This is a necessary condition, not a sufficiency theorem
or a claim that equality is realizable for every L.

The first failed length r0 is exactly the first r>=2 with a_r≠0 when
δ(1)=1: all lower coefficients are a_1=1, a_2=...=a_(r0−1)=0. Its failure
is δ(r0)−1=a_(r0) r0!. Thus a delayed failure exposes a previously forbidden
number of module visits, rather than an unexplained numeric anomaly.

## 3. Exact decision, and what deletion loses

Equation (3) is an exact complete algorithm for the allowed positive lengths.
If q=0 return ALL. Otherwise let M=max{j:c_j≠0}, and, if lower modes exist,
put Q=max{j<M:c_j≠0}, S=sum_(j<M)|c_j|. Find an integer K>=1 with

    |c_M| M^K > S Q^K.                               (5)

Integer repeated multiplication finds K and must terminate. For all m>=K
the leading term dominates and f(m) has its sign. Test only 1<=m<K.
If there are no lower modes, no length succeeds. The number of successful
lengths is also bounded by q−1, regardless of the magnitude of K. This
combines an effective stopping rule with a structural root-count bound;
the largest-base cutoff itself was already used in the seed atlas.

These statements concern m>=1 only. A zero-size module changes the graph
rather than merely its positive-length response. For the inherited
Session260 row A,O={1,3},port2, δ(m)=1 for every m>=1 but δ(0)=−1.
Thus deleting an unselected block cannot be inferred from freely inflating it.
The marked endpoint set must remain outside the changed block throughout.

## 4. A sharp physical two-success example and p-adic false comfort

Use the inherited order7 witness

    rows=(90,28,9,16,36,79,30), O={3,6}, v=2.

The port has (ν,d+,d−)=(2,2,4), hence R=3. Its exact defect is

    δ(m)=1−4 F_3(m)=1−4(3^m−3·2^m+3), m>=1.          (6)

Consequently the successful lengths are **exactly {1,2}**. The first six
values are 1,1,−23,−143,−599,−2159. This attains the R−1 length budget
physically and shows why one extra clone is not a certificate. The graph
and formula are inherited; their sharpness for the present zero theorem
is the new interpretation. Algebraically F_R(m) attains R−1 positive integer
zeros for every R; no physical realization of those general sharp examples
is asserted here.

For every prime p>=5 and precision k>=1 set

    m=1+(p−1)p^(k−1).

Fermat/Euler unit congruences imply 2^m≡2 and 3^m≡3 modulo p^k. Thus
δ(m)≡1 mod p^k while m>=5 and δ(m)≠1. This is an actual tournament family
with arbitrarily deep local unit-seed agreement and no literal success
beyond m=2. It complements Session287's p-adic closure results: local
accuracy does not remove the real exponential zero/height obstruction.

## 5. Compensating clocks: every fixed affine ray is decidable

The same mechanism gives a scoped multiblock tool. Suppose s unselected
ports induce a transitive subtournament and k other vertices remain fixed.
Replace the ports by independent transitive modules of positive lengths
m_1,...,m_s. Between consecutive fixed vertices, each port can be visited
at most once: the port directions form an acyclic total order. Across at
most two paths there are at most k+2 such exterior-free intervals. Thus
every port has at most k+2 runs, and the inherited multivariate skeleton
formula is a **finite** tensor of products of F_r(m_i).

Expand δ(m_1,...,m_s)−1 in positive exponentials. On a prescribed affine
ray m_i=a_i+b_i t, with a_i>=1, b_i>=0 integers and t>=0, every term becomes

    C (product_i j_i^a_i) (product_i j_i^b_i)^t.

Collect equal bases before counting them. If q distinct nonzero modes
survive, the same zero theorem permits at most q−1 successful t; if all
cancel, every point of the ray is a unit seed. The dominance cutoff gives
an exact terminating decision for that ray. Success at q distinct t values
certifies the identity provided q is a known bound on the number of modes.
Use the uncollected finite tensor as that bound if cancellation has not yet
been calculated. This is the same elementary zero mechanism, not a new
multivariate zero-set theorem.

The test detects genuinely compensating joint inflation even when individual
axes fail. Such compensation requires **exact collisions and cancellation of
product bases**, a multiplicative relation in the mode registers. It cannot
be inferred from approximate matching, coefficient gcd, or separate local
p-adic inverses. For example the abstract tensor 2^a−2^b vanishes identically
on a=b while neither coordinate is freely variable; this is an algebraic
control, not an asserted physical tournament defect. The theorem does not
cover arbitrary cyclic port sets or nonlinear paths through length space.

## 6. Reproduction and limits

Run `python3 verify.py` and `python3 -O verify.py` from this directory.
The verifier uses exact subset DP for paths and two-path covers, checks all
labelled tournaments of orders2..4 and every nonempty exterior selected set,
plus eighteen deterministic higher-order core samples. It verifies the exact
Hamilton matching degree, the tax/defect degree and an additional literal
length. For orders2..3 a second method replaces clones by mutually
nonadjacent labelled twins in an auxiliary digraph, recovering each skeleton
coefficient after division by r!. Those auxiliary digraphs are not claimed
tournament outputs. The literal order7 sharp witness and general algebraic
zero controls are separate.

This classifies each **fixed one-block length problem** and its rigidity
budget, not the interaction of several moving lengths, all stable-port
births, arbitrary selected marks within a growing block, or Chan–Pak.

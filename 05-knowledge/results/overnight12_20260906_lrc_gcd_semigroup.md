# Exact signed coefficient boxes and a gcd closure criterion for actual LRC entry

**Status: PROVED RELATIVE TO THE STATED INHERITED INPUTS + INDEPENDENTLY
AUDITED; FINITE-EXACT controls.** The all-parameter conclusions are relative to the explicitly inherited
THM-3818 actual-entry hypotheses and cited lower-dimensional LRC. The elementary
signed-box and three-label relation lemmas have no LRC dependency. The incoming
hereditary gcd theorem is used only in Section 7, not in the main closure proof.
The eleventh unit-component theorem remains frozen with its separate lineage.
No theorem identifier or external priority claim is requested. LRC(14) remains
open.

The [independent referee](overnight12_20260906_lrc_gcd_semigroup_audit.md)
accepts the complete analytic proof and both actual-entry controls with no
remaining repair. Its separate engine passes 149,274 normal/optimized gates.
The root independently read and accepted the full analytic proof as well.

## 1. Inheritance and the map between the two threads

The closest proved mechanism is the
[audited eleventh unit-component closure](overnight11_20260906_lrc_unit_component.md),
with its [independent referee](overnight11_20260906_lrc_unit_component_audit.md).
It combines the actual bounded crossing prohibition in **THM-3818**, Section
6.4, (15o)--(15q), with symmetric cyclic gluing in Section 6.5, (15aa):
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`.
The finite-box supplier is **THM-2052**,
`01-canon/theorems/THM-2052-finite-height-forces-high-rank-bounded-relation-code.md`.
Those source arguments were read for the eleventh audit and are retained here.

The recovered, less-used sidecar is the two-generator semigroup coordinate in
[the earlier trinomial-return classification](synthesis_20260905_moments_trinomial.md),
Section 2. That result distinguishes support return from coefficient-weighted
cancellation. Here the source object is an ordinary two-generator nonnegative
semigroup, the target is an integer coefficient box, and the exact map is

    x = Q(a+b) - (a u+b v).

This preserves the integer equality and converts bounds on u,v into signed
coefficient heights. It discards multinomial weights, phases, and return
multiplicities: no moment noncancellation theorem is transferred. The retained
coordinate is the position relative to the first semigroup gap and to the
finite box boundary. The cheapest decisive test is a complete small signed
box, including its first missing integer, followed by a literal three-label
relation test. Both are performed below.

The canonical normalization hostile is the primorial eleven-component from
**THM-4052**, `01-canon/theorems/THM-4052-lrc14-affine-component-width-escape-cones.md`.
Its gcd is one while its maximum and every relevant endpoint gcd can be large.
Section 7 explicitly reclassifies it after the incoming hereditary sieve: it
is not a surviving unsafe candidate. The corrected near miss is the virtual
ratio-four family whose claimed split was excluded by a bounded crossing row;
large pairwise heights never establish W=V_dec. Targeted MISTAKE-486/490 checks
retain the all-scale atlas and physical-box scope. Exact searches for the
signed-box formula, cutoff, and primitive-unit extension found no duplicate
statement on the searched current surfaces; this is not an external priority
claim.

The live board has six objects: semigroup gaps, signed coefficient boxes,
minimal outside coefficients, endpoint gcds, global protected arcs, and
hereditary subset gcds. The principal result is a faithful consumer of the
first five. The sixth sharpens the remaining hypothetical-failure boundary.

## 2. The exact central interval of a signed two-coordinate box

Let `1<=a<b<=Q` be integers with gcd(a,b)=1 and put

    B_Q(a,b)={a r+b s: r,s integers, |r|,|s|<=Q},
    L=Q(a+b),
    R=Q(a+b)-(a-1)(b-1).                              (1)

**Lemma 1.** Every integer in [-R,R] belongs to B_Q(a,b), while R+1 and
-(R+1) do not. Thus R is its exact largest centered interval radius. Moreover

    Qb <= R <= L,                 2R>L.                (2)

The assertion concerns a centered interval. The box can contain integers
beyond its first gap: for Q=3, (a,b)=(2,3), its first positive gap is 14 but
15 belongs to the box.

**Proof.** For 0<=x<=R, set M=L-x. Then

    (a-1)(b-1) <= M <= L.

Choose the unique u in [0,b-1] with au=M modulo b and set v=(M-au)/b.
If v<=-1 then M<=a(b-1)-b=ab-a-b, a contradiction. Hence v>=0. Also

    v<=M/b<=Q(a+b)/b<2Q,

and `0<=u<=b-1<=Q`. Therefore

    x=(Q-u)a+(Q-v)b

has both coefficients of magnitude at most Q. Negative x follow by symmetry.

If a>=2, at x=R+1 the complement is F=ab-a-b. A nonnegative representation
F=au+bv must have u=b-1 modulo b. The smallest such u is b-1, which already
forces v=-1; all larger choices make v smaller. Thus none exists. If x had
a signed-box representation, its complement would be a nonnegative one with
coefficients at most 2Q, impossible. For a=1, R=L and R+1 is outside the
entire support. This proves the first gap and its negative counterpart.

Finally `R-Qb=Qa-(a-1)(b-1)>0` for Q>=b, and

    2R-L >= b(a+b)-2(a-1)(b-1)
          = b(b-a)+2a+2b-2 >0.

This proves (2). The proof is elementary and does not invoke an external
semigroup theorem.

There is also an exact membership test, including all points beyond the
central interval. For any integer x, put

    r0 = a^(-1)x mod b,             0<=r0<b,
    s0 = (x-a r0)/b.

Every solution is r=r0+kb, s=s0-ka. Hence x belongs to B_Q(a,b) exactly when

    max(ceil((-Q-r0)/b), ceil((s0-Q)/a))
      <= min(floor((Q-r0)/b), floor((s0+Q)/a)).         (3)

Any integer in this intersection constructs a witness. Unlike a numerical
search, (3) is a constant-size exact arithmetic test at arbitrary height.

## 3. Minimal outside coefficient gives an exact three-label test

Let A=tDa and B=tDb be positive physical core coordinates, and let Y be one
positive physical outside coordinate. Keep `1<=a<b<=Q`, gcd(a,b)=1, and put

    delta=gcd(tD,Y),       c=tD/delta,       x=Y/delta.

**Lemma 2.** There is a relation

    C Y = r A+s B,
    1<=C<=Q,       |r|,|s|<=Q                         (4)

if and only if

    c<=Q and x belongs to B_Q(a,b).                    (5)

This is an exact existence test on the selected three labels, including
relations whose support drops to two. It does not decide all possible
crossing supports at once.

**Proof.** Divisibility in (4) forces C=jc with j a positive integer; after
division by tD it becomes jx=ar+bs. If c>Q this is impossible. If x is in
the box, j=1 supplies (4). If x is not in the box, Lemma 1 gives x>R>L/2.
Then every j>=2 has jx>L, outside the entire box support, while j=1 is absent
by assumption. Thus no higher outside coefficient can repair the missing
point. This proves the converse and explains why only the minimal cleared
coefficient is needed.

The inherited bound b<=Q is load-bearing. With Q=2, a=1, b=6, x=3, the
minimal point x is missing but 2x=6 is present. The relation 2*3=6 defeats
the minimal-coefficient conclusion outside its stated scope. The c<=Q gate
is independently necessary: with core coordinates 3 and 6, outside coordinate
1, and Q=2, the normalized x=1 belongs to the box but its minimal outside
coefficient is 3, and no bounded crossing exists.

## 4. Actual-entry closure and its strongest stated scalar threshold

Use exactly the actual-entry hypotheses of the eleventh report. In particular
the primitive thirteen-speed row lies in the THM-3818 box `sum speeds<=Q^2`,
has positive distinct physical coordinates, and is written

    tV disjoint-union g(p,q),
    |V|=11, gcd(V)=1, t,g>0, gcd(t,g)=1,
    1<=p<q, gcd(p,q)=1,

with the two actual decoder components of sizes eleven and two, and with
the full support-at-most-three, height-at-most-Q relation span W equal to
the rank-eleven decoder span V_dec. The pair is in the actual 5,855-entry
atlas; in particular `p+q<=356`, `p<=177`, `q<=355`. A displayed split or
atlas ratio alone does not establish these hypotheses. Here V need not
contain one. Put K=max V.

For any u<v in V, define D=gcd(u,v), a=u/D, b=v/D. The inherited internal
pair-height theorem gives b<=Q even when K itself exceeds Q. For either
outside shape coordinate `s in {p,q}`, put

    delta=gcd(tD,gs),    c=tD/delta,    x=gs/delta,
    R=Q(a+b)-(a-1)(b-1).                              (6)

If c<=Q, equality entry forces x to be absent from B_Q(a,b). Otherwise
Lemma 2 gives a literal bounded crossing on labels tu,tv,gs. Its outside
partial sum is nonzero, so it lies in W but outside V_dec. This is precisely
the prohibited actual-entry relation. Therefore

    c<=Q  implies  x>R.                               (7)

This argument retains signs, zero coefficients, the coefficient height,
and the component restrictions. It never infers W=V_dec from a selected
rank calculation. Equation (3), together with c<=Q, is an exact selected
three-label crossing test and can reject points beyond the scalar radius.

The strongest integer lower bound following from (7) and its retained
divisibility is explicit. Put

    z=gcd(delta,s),       d=delta/z,       e=s/z,
    B(u,v;s)=d*(floor(R/e)+1).                         (8)

Because gcd(d,e)=1 and delta divides gs, d divides g. Writing g=dh gives
x=eh, so (7) is exactly h>=floor(R/e)+1. Consequently any eligible pair
with c<=Q forces

    g>=B(u,v;s).                                      (9)

The bound (8) is optimal for the information `x>R` and `delta divides gs`:
its endpoint attains those two conditions. It is not claimed to optimize
the defining gcd equality for delta, the full decoder constraints, the exact
membership gaps, or all relations.
Those retain more information. A simpler sufficient comparison is
`delta R/s>=42K`, which forces the strict inequality g>42K.

**Theorem 3 (native closure criterion).** Such an actual entry is safe at
clearance at least 1/14 if either

    11t>=21q,

or there is a core pair u<v and s in {p,q} for which

    c<=Q and B(u,v;s)>=42K.                           (10)

Thus one may maximize (8) over all eligible core pairs and both outside
coordinates. There are 110 selected supports, each requiring only exact
integer arithmetic. These are two-core/one-outside supports. They do not
include the eleven supports using both outside labels and one core label;
this collection alone is not a complete crossing test across orientations.

**Proof of the phase conclusion.** The first branch is inherited symmetric
cyclic gluing: a pair one-third-clearance phase gives a closed 1/14-safe arc
of length at least 11/(21q), hit by the complete t-grid. In the second branch,
cited eleven-speed LRC supplies a 1/12-clearance phase for V. The full core
clearance is K-Lipschitz, so its 1/14-safe set contains a closed arc of length
1/(42K). By (9)--(10), g>=42K. The g lifts of any pair-safe phase become a
complete shifted g-grid in the core clock, because gcd(t,g)=1. A closed arc
of length at least one grid spacing contains a point of the grid, including
equality. This gives a safe physical phase. The theorem claims the weak
target, not strictness at every boundary.

The external supplier is the same one read for the eleventh referee:
Sungkawichai--Trakulthongchai, *Eleven, twelve, and thirteen lonely runners*,
arXiv:2604.23906v2, Theorem 1.3, printed page 2. It gives LRC(k) for k<=12
nonzero speeds and is used here at k=11. Its proof is **CITED**, not
independently reproduced here. [Primary PDF](https://arxiv.org/pdf/2604.23906v2).

## 5. A uniform endpoint-gcd cutoff, without a unit assumption

Set

    H=floor(Q/(42*177))=76,388,115,
    Q=567,869,252,041.                                (11)

**Corollary 4.** If V contains a coordinate v<K with

    gcd(v,K)<=H,

then the actual entry is safe at clearance at least 1/14. In particular
this holds if a coordinate is coprime to the maximum. The primitive core
need not contain a literal unit. A minimum-coordinate bound min V<=H is a
weaker sufficient condition, not the intrinsic invariant.

**Proof.** If 11t>=21q, Theorem 3 applies. Otherwise the atlas gives
t<=677. Choose the endpoint pair (v,K), put D=gcd(v,K), and use s=p.
Then

    c=tD/delta<=677D<=677H=51,714,753,855<Q.

By Lemma 1, `R>=Qb=QK/D`. Equality entry therefore forces

    g>delta R/p>=QK/(Dp)>=QK/(Hp)>=42K,

which triggers the second phase argument. The exact R and the exact
integer threshold (8) can close additional rows beyond the uniform cap;
no sharpness claim for H as an LRC boundary is made.

Consequently any still-unsafe actual entry must have gcd(v,K)>H for every
v<K. This is a necessary coordinate obstruction with the global maximum
retained. It is stronger than merely excluding normalized unit components.

## 6. Arbitrary pairs improve a score, but their maximum does not replace K

The following unitless shape provides a fully typed nonvacuity control:

    V=(2,3,4,5,6,10,12,14,15,20,30),
    t=1, (p,q)=(1,3), g=60Q+1=34,072,155,122,461.      (12)

Its core is safe at 1/17 with clearance at least 2/17>1/12. The complete
physical row is primitive, has thirteen distinct coordinates, and lies
inside the actual box. A literal reconstruction of the atlas graph gives
exactly the connected eleven-component and the pair component. Since
g>2QK, no bounded crossing of support at most three exists: one pair term
dominates two core terms, and the nonzero sum of two pair terms is a
multiple of g and dominates one core term. A zero pair sum cannot balance
a nonzero core term. Internal connected decoder spans already have dimension
eleven, so W=V_dec. This establishes actual entry independently of a small
gcd slogan. The source also constructs and checks a literal full safe phase.

At s=p=1, g is coprime to every core coordinate, so delta=1 throughout.
Among pairs containing K=30, the largest radius is at (14,30):

    R_maximum=22Q-84.

The best unrestricted pair is the interior pair (14,15):

    R_interior=29Q-182 > R_maximum.                   (13)

All 55 core-pair radii are checked exactly. Thus arbitrary pairs can
strictly improve the relation-based lower bound on g in an actual-entry
example. This example already closes by the uniform endpoint-gcd criterion;
it does not establish a separate class that all maximum-pair tests miss.
Finding such a class under all current necessary conditions remains an
opportunity, not a claimed theorem.

The protected core arc still depends on K, the maximum of the entire core.
Replacing it by the selected pair maximum is false. Let

    V=(1,2,...,10,85),       y=1/12.

Every coordinate has clearance at least 1/12 at y. Nevertheless 7/85 lies
at distance 1/1020<1/840 from y and has zero clearance for the last speed.
Thus the radius `1/(84*10)` suggested by a pair with maximum ten does not
protect the full core. Arbitrary pairs help the integer relation score;
the present proof does not transfer their Lipschitz constant to the whole
eleven-speed target.

## 7. Incoming hereditary gcd restrictions and honest hostile status

During this work, the fully proved and independently audited incoming result
[recursive hereditary gcd sieve](lrc14_recursive_gcd_empty_core_next_sep06.md)
and its [independent audit](gcd_pair_hunter_audit_empty_core_next_sep06.md)
were read in full. They gave subset-gcd caps 1,2,4,9,32,96 at subset sizes
12,11,10,9,8,7 for a primitive thirteen-speed strict failure. During the
present independent review, the subsequent proved and audited joint-shadow
result was also read in full:
`05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.md`.
It sharpens the current caps to **1,2,4,9,30,90**, by excluding the previous
maximum clocks 32 and 96. These finite classifiers are inherited; this source
does not rerun or replace them. The corresponding namespace reservation is
not used as a proved dependency: the source is the completed, audited result
note itself.

In a hypothetical unsafe actual entry, the eleven physical core has gcd t,
so t<=2. This is a restriction on a strict failure, not on every decoder
entry: safe entries can have larger t. It may be used to simplify the
remaining native coefficient gate to c<=2D/delta. The main proof of
Corollary 4 did not require this new input; it used the older, weaker t<=677
in the unresolved branch. The full hereditary table further restricts every
large subcollection of the physical core.

For the canonical primorial normalization hostile, let

    primes=(37,43,61,67,73,79,97,103,127),
    P=15 product(primes),
    V={2P/r:r in primes} union {P/3,P/5}.

It has gcd one, no unit, maximum `237907127334685115>Q`, and every primitive
internal pair height at most 127. Those facts still disprove the inference
from primitive gcd one to K<=Q. However it is already excluded as a strict
failure by the new hierarchy: the first seven, eight, and nine displayed
coordinates have gcd 392430, 3810, and 30, respectively; adjoining P/3 to
the first nine gives ten coordinates of gcd five. Each exceeds the
corresponding inherited cap. Scaling this core by any t only increases
these subset gcds in a primitive physical completion. The example is
therefore retained solely as an algebraic normalization hostile, not as
a live residual unsafe core.

The remaining structural problem is precise: an unsafe actual entry would
need t<=2, all large-subset gcd caps, no eligible exact crossing in the
two-core/one-outside orientation tested here, and no native pair scale bound
reaching 42K. The pair and subset conditions
retain different coordinates. No implication that their conjunction is
empty is claimed here.

## 8. Exact controls and reproduction

[The standalone source](../../04-computation/overnight12_20260906_lrc_gcd_semigroup.py)
uses standard-library integer and rational arithmetic, no producer imports,
and always-active gates. Its bounded universes are:

* all 576 primitive pairs 1<=a<b<=Qtoy for Qtoy=2,...,17, constructing full
  signed coefficient sets and comparing all 223,530 test points with (3);
* all 7,830 specified three-label coefficient boxes for Qtoy=2,...,6,
  primitive a<b<=Qtoy, D in {1,2,5}, t in {1,2,3}, and outside coordinate
  1,...,30, comparing every bounded outside coefficient with Lemma 2;
* 4,320 complete divisibility controls R=1,...,30, delta,s=1,...,12,
  verifying that (8) is the first possible integer scale under its two inputs;
* exact inherited-height cutoff arithmetic, both typed actual-entry controls,
  all 110 selected-support tests for each, a literal safe phase for each,
  and all 55 pair scores in the unitless control;
* the coefficient-gate, missing-height, first-gap, global-arc, and primorial
  hierarchy hostiles described above.

The actual entry controls use a freshly generated atlas and graph traversal,
plus a proved dominance exclusion of every crossing support; the generic
toy coefficient bank is not mislabeled as actual entry. The finite bank
tests formulas and failure boundaries. The all-height claims follow from
the displayed proofs.

    python3 -B 04-computation/overnight12_20260906_lrc_gcd_semigroup.py
    python3 -B -O 04-computation/overnight12_20260906_lrc_gcd_semigroup.py

Both runs pass **673,482 always-active gates**. Normal and optimized outputs
are frozen byte-identical LF. Source SHA256:
`15a3fa8f511640c2fee3404cb5cc5b377af4164f8db3dcec7aa1702c2f15e4c3`.
Output and optimized-output SHA256:
`b3f47f0f979c9b126dd6a67b2a83c0dc8dc30fe46b9834d389d3425f82f66ab8`.
No Git, shared navigation, or prior frozen file was modified.

**Filing:** root integrated these independently audited artifacts in the twelfth
checkpoint. Reproduction commands are relative to the repository root.

# An arbitrary-source wrapper and a new family of neutral preparations

**PROVED IN NOTE:** the all-source marked profile, the new all-parameter controller-entry family, and the stated traps. **FINITE-EXACT:** literal graph controls and the named H847 replacement. Chan-Pak and the next ready 405*2^k+71 family remain OPEN. This note does not use the word readiness for two different defects.

## Inheritance and the actual target

The inherited [H17 preparation](../../eight_nine_twentyseven_20261112/constructors/SIX_RETURN.md) changes a named internal root's response while preserving H17. Incoming [Session 298](../../session298_exact_target_geometry_20261113/README.md) extends that to every 2^m+1, and its [two-width chart](../../session298_exact_target_geometry_20261113/TWO_WIDTH.md) treats C3[TT_m,1,TT_b]. The inherited [native controller](../../source_selection_closure_20261103/NATIVE_SELECTION.md#7-closed-primitive-cycles-and-an-eight-state-branching-controller) accepts an actual balanced rooted source (b,t,t,b) with gcd(b,t)=1, 3|t, 3 not dividing b and then renews it physically.

The new step is to replace TT_b by an arbitrary tournament B, retaining its first three path-cover evaluations. This gives a source family outside the two-width transitive chart. A strong four-vertex B then supplies infinitely many exact initial targets where a fixed neutral root preparation unlocks the native guard using the same three-birth word RNR.

Closest mechanism: Session 298's marked two-width formula. Canonical hostile: a source module with a surviving common profile factor. Corrected near miss: full H and even rooted coordinate gcd 1 do not ensure that a chosen continuation reaches its guard. Least-used sidecar: the arbitrary source triple(H(B),P_B(2),P_B(3)), together with the position of the named root outside B.

The local board is: exact unrooted target; internal-root phase; three source evaluations; integral return defect; prime-phase exclusion; physical balanced reentry. The Anchor is reusable marked preparation. The Niche is an arbitrary-source extension of the inherited chart. The Wildcard is the unexpected 199/847 family and its arithmetic phase selector.

## P1. The arbitrary-source marked profile

Let B be any nonempty finite tournament. Let P_B(j) count ordered spanning covers by j directed paths, with empty bins allowed. Put

    h=H(B), p=P_B(2), t=P_B(3).

Build a directed triangle of parts C=TT_m, X={x}, B, with C→x→B→C. Assume m>=2. Mark a named resource vertex with l predecessors and r successors in C, l+r=m-1. The graph is strong for every nonempty B, whether or not B itself is strong. Changing only the resource order keeps the same unrooted tournament type and its full Hamilton count.

For the marked state (h0,s,q,F), let s count full paths starting at the root, q count full paths after its deletion, and F count ordered spanning pairs(P,Q) with P nonempty ending at the root and Q possibly empty. Then

    h0=(2^m-1)(p-h)+2h,
    s=2^r h,
    q=(2^(m-1)-1)(p-h)+2h,
    F=2(2^r3^l-2^l)(t-p)+3*2^l p-2*2^l h-[r=0](p-h).       (1)

All statements retain m>=2; after deleting the only resource at m=1, the cyclic quotient disappears and the displayed q formula would be wrong.

**Proof of the extension.** A Hamilton path can enter B from outside only at its own start or from x. Therefore it has at most two B-runs. A spanning ordered pair of paths can have two initial B-runs, one per bin, and at most one extra B-run entered from x: at most three total. The same uniform-module skeleton decomposition used in the [previous compiler](../../cycle_opening_six_return_20261112/counts/PROOF.md#3-exact-run-coefficients-and-arbitrary-modules) expresses these observations as finite nonnegative linear combinations of ordered nonempty path-cover counts K_j(B), with j<=2 for h0,q and j<=3 for F.

On B=TT_b, the inherited two-width theorem supplies (1) for every b>=1. Since K_j(TT_b)=j!S(b,j) and the first three such functions are linearly independent, those all-b identities identify every possible coefficient and hence prove (1) for every B. This is a bounded-degree skeleton proof using an inherited all-b formula, not a three-point numerical extrapolation.

Explicitly, with A=2^r3^l, D=2^l and I=[r=0], the F coefficients in the ordered nonempty basis are

    F=(2A+2D-I)K_1(B)+(4A-D-I)K_2(B)+2(A-D)K_3(B),          (2)

because p=2K_1+K_2 and t=3K_1+3K_2+K_3. All coefficients are nonnegative. For the unrooted count,

    h0=(2^m+1)K_1(B)+(2^m-1)K_2(B).

A root-starting path traverses an increasing subset of the r successors, then x, then one full path of B, then the remaining resources. This proves s directly. Deleting the root leaves TT_(m-1), proving q. QED.

Moving the physical root across one adjacent resource reverses exactly their mutual arc. Each intermediate resource order is transitive, so h0 and q stay fixed at every step. The inherited [adjacent-port implementation](../../source_selection_closure_20261103/NATIVE_SELECTION.md#4-select-one-port-with-a-linear-native-clock) also compiles each such arc swap into actual native moves with one temporary vertex. Ordinary adjacent arc swaps preserve H at their endpoints; temporary native stages need not have that same H.

## P2. A fixed preparation and fixed exit at every 27*2^m-17, m odd

Choose the actual strong four-vertex source

    B=(6,12,8,1)               (outgoing bit masks),
    (h,p,t)=(5,32,117).

Let m>=3 be odd, let r=m-2 be odd, and put z=2^r. The initial unrooted target is

    h0=27*2^m-17=108z-17.                                  (3)

Begin with the named root last in C. Move it to the second position by m-2 adjacent swaps. This preserves every vertex name, the whole unrooted type, q and h0; only arcs incident to that root within C change. In its prepared phase l=1, the exact state is

    (h0,s,q,F)=(108z-17,5z,54z-17,510z-168).                (4)

The first targets are 199, 847, 3439, 13807, 55279.

Use the inherited physical operations R and N: add a vertex y with root→y→all other old vertices, and either retain the root (R) or mark y (N). Their exact state maps are

    R(h0,s,q,F)=(h0+F-s,q,q,2F-q),
    N(h0,s,q,F)=(h0+F-s,h0-s,h0,h0+F).

Define the return defect

    Delta=F-h0-q+s,
    d=F-q.

This is distinct from the indegree-one readiness defect q-s. For any k>=1, R^kNR returns the actual balanced state (b_k,t_k,t_k,b_k), with

    t_k=2^k d+q-Delta,  b_k=3t_k+2Delta.                   (5)

For (4),

    Delta=353z-134, d=456z-151,
    t_k=(456z-151)2^k-299z+117.

In particular the **single fixed word RNR** gives

    t_1=613z-185,       b_1=2545z-823.                    (6)

**Theorem.** For every odd m>=3, (6) satisfies

    gcd(b_1,t_1)=1,       3|t_1,       3 not dividing b_1.

Thus the same exact-target preparation and same three-birth exit feed the inherited repeatable native controller at every target (3).

**Proof of primitivity.** z is even, so t_1 is odd. Since b_1-3t_1=2Delta,

    gcd(b_1,t_1)=gcd(t_1,2Delta).

Modulo Delta, t_1 is 2F-q=966z-319. The linear resultant is

    353(966z-319)-966(353z-134)=16837=113*149.

Any common odd prime therefore belongs to {113,149}. Divisibility of Delta would require z=58 modulo 113 or z=81 modulo 149. The following finite arithmetic identities rule both out:

    2^28=1 mod113,       58^28=98 mod113;
    2^74=-1 mod149,      81^74=1 mod149.

The first forbids 58 among every power of 2. The second forbids 81 among odd powers of 2, which is exactly the r restriction. Both moduli are prime (trial division suffices). Therefore the gcd is 1. Finally z=2 modulo3, so t_1=0 and b_1=1 modulo3. QED.

The source is strong. Every R/N birth preserves strongness. The final R seals an actual root-child pair, not only a formal state. The coordinate gcd before this final integer state map divides the output gcd 1, so that pair meets the inherited primitive-content guard as well. The native controller therefore applies to the concrete output graph and retained interface.

## P3. Why the neutral preparation matters, and infinitely many exits

At the original last-root phase l=m-1, r=0, equation (1) gives

    F=170*3^l-84*2^l-27,
    Delta=170*3^l-165*2^l+12.

Here l>=2, so Delta is divisible by 3. Consequently every returned b_k in (5) is divisible by 3. No choice of k reaches the declared native guard from this phase. The preparation repairs a specific continuation capability at the **same** full target; it does not merely improve an unrooted scalar count.

At m=3 the unprepared and prepared states are

    (199,5,91,1167)  ->  (199,10,91,852),
    Delta=882        ->  Delta=572.

Only one old arc changes. RNR then gives (4267,1041,1041,4267). The useful preparation **decreases** Delta. Thus maximizing a defect or a convex response is not equivalent to optimizing the arithmetic controller-entry predicate; the residue and gcd conditions remain necessary sidecars.

There are also infinitely many successful exit lengths at every fixed source target (3). Let Q be the odd part of Delta and put

    L=lcm(2, ord_Q(2)),

with ord_1(2)=1. Every k>=1 with k=1 modulo L has 2^k=2 modulo Q. The preceding resultant proof gives gcd(2F-q,Q)=1, hence gcd(t_k,2Delta)=1. Such k are odd, and substituting z=2 mod3 in (5) gives 3|t_k and 3 not dividing b_k. This supplies a whole exact progression of physical primitive exits, while k=1 already gives a short uniform choice.

The repeatability is the separately proved C3/C5 controller after entry. Reapplying R^kNR to an arbitrary balanced output is not asserted to be the same renewal. Nor does (6) equal the desired 3h0+2 target. The indegree-one readiness defect of these balanced outputs is zero, not one; this lane does not construct the missing ready 405*2^k+71 family.

## P4. The 141 hostile identifies what must be selected

For B=C3, the source triple is (3,12,33). At m=4, every root position in C has full H=141. Its four states are

    (141,24,69,324), (141,12,69,480),
    (141,6,69,708),  (141,3,69,1029).

Every coordinate is divisible by 3. More generally, gamma=gcd(h,p,t) divides every coordinate in(1), for every resource phase. Integer R/N state maps preserve this common divisor. If gamma>1, no neutral rearrangement of this particular resource block followed by R^kNR can give a primitive balanced state.

This is a complete obstruction for the specified arbitrary-source wrapper and operation menu. It does not exclude selecting a port **inside B**, changing B, or a count-neutral operation that destroys B's uniform module boundary. The example explains why blindly replacing the transitive source by a cyclic source can make an apparently richer preparation permanently trapped.

The source triple is only a necessary arithmetic interface. gamma=1 does not by itself ensure a useful phase or branch. The positive family above proves its extra conditions using an actual source and two explicit finite-field phase tests.

## P4a. Complete unit-ready port boundary of this wrapper

For a root inside C, equation (1) gives

    q-s=(2^(m-1)-1)K_2(B)+(2^(m-1)+1-2^r)H(B).             (7)

Both terms are nonnegative and the second coefficient is at least one.
Thus q-s>=1, with equality exactly when B is a singleton and the root is
first in C. Indeed K_2(B)=0 holds exactly for a singleton: any nontrivial
bipartition contributes Hamilton paths when |B|>=2. This statement concerns
the resource-root chart; other physical ports need separate treatment.

Here is that treatment, derived with the coordinator. A C vertex with l
predecessors has indegree |B|+l, while x has indegree m>=2. A vertex b in B
has indegree 1+indegree_B(b), so the only possible indegree-one port in B
is a global source. Write B=b->D. If D is nonempty, deleting b leaves
U=C3[TT_m,{x},D], and putting b back is exactly the N birth at x. Hence

    q_b-s_b = s_x(U) = H(D).                              (8)

The last equality follows because a path starting at x must traverse one
full path of D and then the increasing C chain. H(D)=1 exactly when D is
transitive: a nontransitive tournament contains a directed triangle, whose
three Hamilton paths extend distinctly by vertex insertion. If D is empty,
the deleted graph is C->x and s_x(U)=0, so b is not unit-ready.

Consequently every C3[TT_m,{x},B], m>=2, has the following exact classification
of ports satisfying both indegree one and q-s=1:

- B singleton: exactly the first vertex of C.
- B transitive of order at least two: exactly the first vertex of B.
- B nontransitive: no unit-ready port anywhere.

For transitive B of order b, the target is
(2^m-1)(2^b-1)+2. This scalar family was already the rectangular one-apex
maximum in [Session 26](../../session26_random_bridges_20261022/session26-tournament/CUT_ADDENDUM.md#maximum-and-equality),
and is the inherited Session 298 two-width chart. The port classification
makes explicit why the nontransitive source in P2 supplies balanced native
entries while leaving unit readiness. A selective extension seeking a ready
return must leave this wrapper or change its source/port architecture.

## P5. A literal changing-core bridge at H847

The earlier H17 preparation followed by R^3NR gives an eleven-vertex source at

    (h0,s,q,F)=(847,269,269,847).

Its t=269 is 2 modulo3, so it misses the direct native-entry guard. Our new m=5 prepared carrier has ten vertices and state

    (847,40,415,3912).

Identify the old physical root with the new carrier's marked root, retain nine other names in their order, delete the one unmatched vertex, and set the remaining arcs to the new target carrier. The receipt gives both literal adjacency matrices, the retained-name map, the deletion, and every changed arc. Both endpoints have full H=847 and the root's physical name is retained. Its incident arcs and indegree are allowed to change. This is an endpoint-neutral replacement macro; no claim of H-neutrality at each intermediate individual edit is made.

RNR on the replacement gives

    (19537,4719,4719,19537),

a primitive native-controller entry. Thus an explicit count-preserving change of core repairs a known awkward branch of the H17 family. The two initial-count families intersect only at 847: solving 105*2^k+7=27*2^m-17 for positive k,m reduces to 35*2^k+8=9*2^m; its 2-adic valuation gives k=3, m=5. No uniform intersection or arbitrary-target replacement is inferred.

## Verification and connection contract

Run `python -B 05-knowledge/results/changing_core_returns_20261113/preparation/check.py` and repeat with `-O -B`. Both receipts agree. The exact universe comprises all labelled B of orders 1..4, resource widths m=2..4, and every named-root position: 675 physical contexts. A separate literal permutation/cut enumeration checks 35 of the small contexts. Selected actual RNR and RRNR outputs are counted independently by subset paths; the all-parameter theorem is the proof above. Prime-power phase identities, the 141 hostile and the 847 replacement are explicit controls.

Source -> target: arbitrary B and a transitive resource block -> an exact same-H marked phase -> a physically constructed balanced controller entry. The map preserves the initial unrooted type and H during preparation. It changes rooted endpoint data and then deliberately changes H during the exit. The required sidecars are (h,p,t), the root phase, Delta, its odd-part clock, and the literal root-child interface. The cheapest hostile is B=C3 at H=141; the cheapest positive case is H=199. The family is sparse and supplies no general exact-target selector. The separate `port_readiness.py` checks P4a at all5,415 physical ports in the same675 contexts, in normal and optimized modes. The [881 changing-core study](NEUTRAL_881.md) treats different endpoint operations and proves their exact minimum arc widths. Chan-Pak and the ready ternary return remain OPEN.

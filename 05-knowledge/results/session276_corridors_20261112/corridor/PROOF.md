# A central switch is not yet a corridor: exact fixed-core islands

**Session276. PROVED in this note for the isolation theorem and charge-fiber
classification. Finite experiments are separately labelled. Chan–Pak OPEN.**

## 1. The carrier and the inherited opportunity

Session275's `landscape/PROOF.md` defines a count corridor using actual
framed moves with changes0 or2, neutral padding, and reframing. Its initialized
component contains23. Its elementary intermediate-value criterion says that
unbounded **count**, rather than unbounded size, in that component suffices.

The Session276 digits lane identifies two rooted six-vertex tournaments
(root0, outgoing rows)

    T41=(56,13,25,48,34,6),
    T43=(52,13,24,49,34,6).

Their common root-deletion tournament is the cyclic five-vertex tournament
C5, represented on0,...,4 by outgoing rows(6,12,24,17,3). They have the same
full fusion interface except for H, with difference2E14. This central-switch
identity is established in the digits lane, not claimed again here.

That observation suggests a cheap corridor: surround C5 with arbitrarily
many transitive vertices, retain the apex, and use small changes of its
incidences. The following all-order theorem determines exactly why that
particular suggestion fails and which extra operation is needed.

## 2. Exact isolated island, at every padding length

Fix l,r>=0. Let U=TT_l -> C5 -> TT_r and add apex a. Write its outgoing
incidences as a binary word. A **neutral-tail cut** has zero on the left
transitive block, one on the right block, and an arbitrary five-bit core cut.
The outside vertices are actual global sources/sinks and can be deleted
without changing H.

Let I_(l,r) be the20 neutral-tail cuts whose core cut has size2 or3.
On the whole apex-cut cube retain edges that change at most two incidences
and satisfy |Delta H|<=2. This is a RELAXATION of the same-root framed
star-move graph: a common Hamilton frame is not required in its definition.

**Theorem (fixed-C5 island).** For every l,r>=0, I_(l,r) is exactly one
connected component of this relaxed graph. Its count set is {41,43}.
Consequently it is also closed under all legal framed one/two-star moves.
Neutral padding, its inverse, and apex-preserving core automorphisms do not
let this island escape. The conclusion does not concern arbitrary reframing
that changes the chosen apex, internal-core edge changes, or general cut gauges.

The minimum count jump to leave I_(l,r) by one/two star changes is at least8.
If any tail incidence is changed, the jump is at least14. The8 bound is
attained by a core-only move41->33. The14 bound is attained by43->57 when
one tail vertex is available on the appropriate side.

### Proof

The32 core cuts of C5 have the following complete table by cut size k:

    k                 0    1     2       3     4    5
    possible H       15   33   41,43   41,43  33   15.

This is a finite lemma, reproducible by independent permutation counting of
six vertices in `check.py`. It can equivalently be obtained from the explicit
quadratic cut form. All20 middle cuts are connected: the Johnson graph on
2-subsets uses two-bit exchanges, the Johnson graph on3-subsets does likewise,
and a2-subset joins any3-subset containing it by a one-bit move. Every such
edge has count difference0 or2.

For the quantitative tail bound use the inherited seam coordinates of
Session275 `seams/PROOF.md`. C5 is vertex transitive, has h=15 and q=120,
and its single-port masses are s_v=e_v=3 and R_v=L_v=21. Thus a core cut of
size k has s=e=3k and R=L=21k.

Reverse the apex incidence of one left tail vertex, leaving all other tail
incidences neutral. Delete the earlier neutral sources and every neutral
sink. If d neutral source vertices remain between this vertex and C5,
the resulting count is

    H(core cut)+120-24k,                   d=0,
    H(core cut)+2^d(120-21k),              d>=1.                 (L)

For a reversed right-tail incidence, with d neutral sink vertices between it
and C5, the count is

    H(core cut)+24k,                       d=0,
    H(core cut)+2^d(21k+15),                d>=1.                 (R)

For (L), neutral source padding transforms (q,L,s) to(2q,2L,0),
while a source seam adds q-L-s. For(R), neutral sink padding transforms
R to2R+h and e to h; the final sink seam adds R+e. These are positive
physical operations, including their d=0 boundaries.

Starting in the island, a move changing one core bit and one tail bit has
new core size1,2,3,or4. Equations(L),(R) and the finite table give final
H>=57: the smallest case is k=4 on the left, or k=1 on the right,
with33+24=57. The starting count is at most43, so the difference is at least14.
A move changing one or two tail bits but no core bits has final H>=H_old+48.
For two changed tail bits, delete one of them, then every other neutral tail
vertex: Hamiltonian-path count of a tournament is monotone under induced
subtournaments, so the single-anomaly lower bound still applies. For clarity,
monotonicity follows by inserting a new vertex at the first lawful location
in each old Hamilton path; deleting it recovers that path, giving an injection.

Finally, a core-only move outside sizes2/3 lands at count15 or33, so its
absolute jump from41/43 is at least8. This excludes every possible exit under
at most two incidence flips and proves the theorem. The literal examples
attaining both minima are in the replay. QED.

## 3. A genuine escape changes the deletion core

The obstruction identifies a missing operation rather than a global failure.
The single INTERNAL core edge reversal{3,4} in the displayed T43 produces

    T45=(52,13,24,33,42,6), H=45.

Both tournaments have common Hamilton frame(0,2,4,1,3,5), which does not
use{3,4}. This is therefore an actual one-free-tile corridor move43->45,
not merely a scalar adjacency. Together with the shared frame for41->43
it supplies a literal41->43->45 corridor of three values.

But its transported interface is different. The old/new root packets
(q,beta,gamma,delta) are

    T43: (120,-78,-57,73),
    T45: (104,-60,-45,57).

With positive source/sink fan lengths m,n and x=2^m,y=2^n, its count
change becomes

    -16xy+18x+12y-16.

For x,y>=2 this is at most-20: the partial differences are negative there,
and the maximum at x=y=2 is-20. Thus the core-changing escape is NOT
preserved as a0/2 move by positive bifans. Single source continuation preserves
it once (91->93), but its difference after m>=1 source seams is2^(m+1)-2,
so it too stops being a low-jump escape after that first step.

This is a precise compatibility obligation for a corridor proof: it needs
both central motion along one interface AND an interface-changing escape
whose admissibility can be restored recursively. Keeping only either one
of those is insufficient. No universal impossibility is claimed.

## 4. The exact charge fibers over a transitive deletion core

Let U=TT_n in increasing order and let a_i indicate apex->i. Then

    R_i=2^(n-1-i), L_i=2^i,
    c(A)=R(A)-L(A)=sum_i[2^(n-1-i)-2^i]a_i.

For i<n/2 define d_i=a_i-a_(n-1-i) in{-1,0,1}, and
w_i=2^(n-1-i)-2^i. Then c=sum_i w_i d_i.

**Theorem (charge-fiber cube).** The integer c uniquely determines the whole
ternary vector(d_i). Its fiber consists EXACTLY of choosing00 or11 in every
reflected pair with d_i=0; pairs with d_i=±1 are fixed. When n is odd the
middle bit is also free. Consequently there are exactly3^floor(n/2) charge
fibers; the fiber dimension is the number of zero d_i plus(n mod2).
For n>=2, z=e(A)-s(A)=a_(n-1)-a_0=-d_0, so c already determines z.

**Proof.** The positive weights satisfy v2(w_i)=i and
w_i>sum_(j>i)w_j. Suppose two ternary vectors have the same c and i is
their first differing coordinate. Reduce the difference modulo2^(i+1):
the coefficient difference is even, hence±2. Its term has absolute value
2w_i, whereas all later terms together have absolute value at most
2sum_(j>i)w_j<2w_i. Contradiction. The fiber description follows. QED.

These fibers are precisely the cut families to which Session275's balanced
translation theorem applies when the deletion core is transitive. This is
more information than knowing its unique individually balanced middle port.
It retains the paired incidence choices needed to assess moves.

For tests use the inherited odd-cycle expansion, or derive it directly by
choosing the spine vertices of each cycle through the apex:

    H(a)=1+2 sum_(i<j) a_i(1-a_j) f(j-i),
    f(1)=1, f(d)=2^(d-2) for d>=2.

All cycles contain the apex, so no two nontrivial odd cycles can pack.
The number of even-size selected spine sets with endpoints i,j is f(j-i).
The one-apex cycle-packing formula is inherited, not asserted novel.

`charge_probe.py` enumerates every cut and every natural fiber-cube edge
through n=14. The0/2 components all have count width at most2, despite
increasing cube dimensions. This is **FINITE-EXACT evidence only**, not an
all-order width theorem. Arbitrary simultaneous flips, changes of deletion
core, or other framings are outside that experiment.

## 5. Connections and a repaired next target

The less direct source was the occupied-helper stabilization work in
`05-knowledge/results/compatible_charged_helpers_20261023/fibres/stabilization.py`.
Its useful lesson is that a context can preserve an old charge yet fail to
connect its old endpoints; a new occupied/unoccupied helper changes the live
fiber. Here the source is a same-interface central switch, the map is neutral
transitive padding, the preserved predicate is the exact41/43 count pair,
and the lost opportunity is movement of the deletion core. The explicit
43->45 edge is a native core-changing helper, but its continuation packet
shows exactly what it costs. No theorem from the helper problem is imported
as a tournament theorem.

The ternary charge digits are also an example of combining an Archimedean
bound with a2-adic pivot: dominance alone is too weak to decode differences
in{-2,...,2}, while valuations force the leading difference to have magnitude2.
This is a precise mechanism, not a numerical analogy between unrelated gaps.

A repaired constructive target is therefore a **renewable core-change gate**:
a finite or growing controlled menu of interface-changing moves that can be
rearmed after central translation, together with an exact initial attachment
to23 and increasing count height. The island theorem specifies a family where
that extra gate is indispensable. The displayed gate supplies its first
physical example. Section6 then finds its correct unbounded context language,
while keeping the missing inter-scale return map explicit.

## 6. A positive repair: classify every one-sided lossless escape context

The bifan failure does not end the investigation. In fact the interface-changing
43->45 gate has an unbounded family of exact transports, and the admissible
one-sided contexts can be completely classified.

Use the inherited fusion matrix M of Session264/Session275. For the two gate
endpoints its difference is exactly

    D=M(T45)-M(T43)
      =[-2,-4,2,4; 0,-16,0,2; 0,0,0,2; 0,0,0,-2].

Let C be ANY pointed tournament with h_C=H(C-root)=1. Its deletion is
transitive (a directed triangle would already give at least3 paths), so C is
a one-apex transitive-spine tournament. The singleton identity is included.
Direct matrix multiplication gives the following exact scalar formulas:

    H(C*T45)-H(C*T43)=2(U_C+d_C-H_C),                         (8)
    H(T45*C)-H(T43*C)=6-2H_C-4V_C+2e_C.                     (9)

Here * is physical rooted fusion, not formal subtraction. The matrices are
only an evaluator for the actual endpoints.

**Theorem (lossless escape contexts).** The left fusion in(8) retains difference
+2 iff the root of C beats at most one deletion-spine vertex. The right fusion
in(9) retains difference+2 iff C's root is a global source (including the
singleton). The right-side contexts are count-neutral padding. The left-side
family contains genuine unbounded, non-neutral transports.

**Proof of the left classification.** Label the transitive spine0,...,n-1,
and put a_i=[root->i]. The exact one-apex formulas give

    U_C=1+sum_i a_i2^(n-1-i),   d_C=1-a_(n-1),
    H_C=1+2 sum_(i<j) a_i(1-a_j)f(j-i),

where f is as in§4. For each i<n-1,
2 sum_(j>i)f(j-i)=2^(n-1-i). Therefore

    U_C+d_C-H_C
      =1+2 sum_(i<j) a_i a_j f(j-i).                         (10)

Every f is positive, so(10) equals1 exactly when at most one a_i is1.
The n=0 singleton has H=U=d=1 and satisfies the conclusion separately.

For the right classification, all endpoint path masses are nonnegative,
|e_C|<=V_C, and H_C,V_C>=1. Equality(9)=2 requires
H_C+2V_C-e_C=2. It forces H_C=V_C=e_C=1. The weighted root-ending mass V_C
contains its singleton contribution h_C=1. Any incoming root arc would supply
an additional positive length-two path with its remaining complement's
positive Hamilton count, using H(empty)=1 at the two-vertex boundary. Thus V_C=1 iff the root has no incoming arc.
Conversely a source root gives H=V=e=1 when its deletion is transitive. QED.

### Exact translations of a whole three-value physical corridor

The41->43 central switch also survives every h1 context. Thus the classification
above transports the whole41->43->45 chain, including its core-changing edge.
The three inputs all have receiving coordinates V=57,e=9. If C's root beats
no spine vertex, left fusion is just neutral source padding. If it beats a
single vertex with t later spine vertices, its exact translation is

    shift=48                  when t=0,
    shift=57*2^t              when t>=1.

Preceding spine vertices are neutral global sources. For t>=1,
(H_C,U_C,d_C)=(1+2^t,1+2^t,1), so the first row of M(C) is
(1,2^t,0,2^t); multiplying M(T) gives H(C*T)=H(T)+2^t V_T.
For t=0 it is(1,1,-1,0), giving H(C*T)=H(T)+V_T-e_T.

Consequently the complete set of retained three-value translated corridors
under arbitrary h1 left contexts is

    {41,43,45}, {89,91,93},
    {41+57*2^t,43+57*2^t,45+57*2^t},       t>=1.             (11)

These are not merely scalar pairs. Both old moves have common frames starting
at the root, recorded in CHECKS.json. If t>=1, prepend the entire context
spine to that common frame; its last vertex beats the root, so it is a shared
Hamilton frame of both fused endpoints. If t=0, put the sole outgoing context
vertex immediately after the root and its preceding global sources first.
For no outgoing vertex prepend the whole spine. Every changed old edge stays
free. Thus(11) is an all-length family of actual three-value low-jump paths
with a deletion-core-changing middle gate.

This is **renewal as a family of physical constructions**, not a claim that
successive corridors are connected by a low-jump return path. Their shifts
are exactly{0,48} union{57*2^t:t>=1}; their widths are4 and all intervals are
disjoint. Hence even allowing EVERY h1 left context that preserves the entire
gate cannot make these translated ranges overlap. This is an exact
classification/no-overlap result, not a finite failure to find a context.

The escape therefore can be preserved indefinitely once the correct context
language is chosen, but its absolute landing remains lacunary. To go further
one must widen the seed corridor, permit controlled non-preservation of this
gate, use two-sided correlated contexts, or vary the receiving interface.
The successful native repair and the exact remaining obstruction are both
retained; neither renewable local gates nor their unbounded heights prove
Chan–Pak.

### Receiving-level transport, and its exact continuation leakage

There is a useful generalization beyond the three literal inputs. For t>=1,
let C_t have a transitive deletion spine of t+1 vertices, with its apex beating
only the first. Its full relevant data are

    h_C=1, Q_C=2^(t+1), H_C=U_C=1+2^t,
    V_C=2^(t+1)-1, d_C=e_C=1.

For every pointed T, physical left fusion gives

    h'=h,    H'=H+2^t V,    V'=2^(t+1)V-h.                  (12)

Thus an arbitrary family with common receiving mass V is translated by
2^tV, even when its deletion counts and other interface coordinates differ.
Every common old Hamilton frame survives by prefixing the context spine;
the final spine vertex beats the root and every other old vertex. So every
old free-tile move remains that same free-tile move. This is a receiving-level
transport theorem, strictly less restrictive than requiring full equal
interfaces. Its proof is the first and second rows of the fusion product,
and the literal frame construction.

But the last coordinate in(12) records a genuine leakage. The41,43,45 inputs
have h=15,15,13. At t=1 their new counts are155,157,159 and their receiving
masses are213,213,215. The interface-changing gate remains low-jump, but
moves off the common-V fiber needed for another indiscriminate application.
A second C_s context makes its count difference2+2^(s+1), not2.

Nor can further h1 LEFT contexts alone equalize these receiving masses.
For such a context D, Q_D>=V_D and Q_D>=1 (its transitive deletion has
Q_D=2^n and V_D=1 plus a sub-sum of1,2,...,2^(n-1)). WithDelta h=-2
andDelta V>0 the next difference is

    Delta V'=Q_D Delta V +2(Q_D-V_D) >0.

This last obstruction is one-sided and preserves the old mark. Rerooting,
correlated right contexts, or a new deletion-changing gadget are not excluded.
Keeping (12), not merely its scalar translation, is essential to an actual
return-map proof.

### Sharp two-sided nonclosure under all unit-deletion fusion programs

The one-sided leakage obstruction can be completed for arbitrary same-mark
h1 sandwiches, not just repeated C_t contexts.

**Theorem (unit-context nonclosure of the core-changing gate).** Let L,R be
any pointed tournaments with transitive root deletion, including singleton
identities. If the two physical outputs L*T43*R and L*T45*R have the same
receiving mass V, then L and R perform only neutral global-source/global-sink
padding. In particular no such sandwich can both increase their counts and
restore the common receiving coordinate. The theorem concerns a common
context program applied to both endpoints with the same retained root.
It does not exclude graph-dependent contexts, rerooting, or h>1 contexts.

**Proof.** The original differences are

    Delta h=-2, Delta Q=-16, Delta V=0.

Write q=Q_L. Since L has transitive deletion,

    V_L=q-A,      0<=A<=q-1,
    V_R=1+B,     B>=0.

Indeed V_L=1+sum_(i: i->root)2^i=q-sum_(i: root->i)2^i.
The singleton has q=1,A=0. Fusion first on the left gives
Delta V=2A andDelta Q=-16q. Fusion on the right then gives

    Delta V_final=2A-16qB.                                  (13)

If B>=1, the expression is strictly negative because2A<=2q-2<16q.
If B=0, it vanishes iff A=0. A=0 says the left root beats no deletion
vertex, so the whole left block is global-source padding; B=0 says the
right root receives from no deletion vertex, so the whole right block is
global-sink padding. Both are count neutral. Conversely those neutral
contexts preserve the original equal receiving masses. QED.

An arbitrary finite sequence of left/right h1 fusions without changing the
root collapses to such a sandwich: fusion is associative and the root-deleted
parts of h1 contexts concatenate transitively, so their products again have
h1. Thus(13) is a program-level obstruction, not merely a one-step failure.
It pinpoints why the transported three-value gate cannot close a reusable
same-receiver loop inside the unit-context language.

The search then allowed actual root changes. For the six transported triples
with t=1,...,6, the three complete root menus have no common V. This is only
a **finite negative control**, recorded in `reroot_receipt.jsonl`; no all-order
root-menu obstruction is inferred. The remaining positive direction is a
larger corridor or a context-dependent core mutation that changes the
receiver defect while preserving exact small count steps.

## 7. One finite core repair beyond the unit-context obstruction

The sharp obstruction suggested changing internal edges rather than adding
more unit contexts. At t=1 the transported T45 has H159,V215. Flip pairs
{2,3},{3,7}, in the labels used by `escape_check.fuse`, to obtain

    R159=(210,252,241,181,96,129,168,16),
    (H,h,Q,U,V,d,e)=(159,17,448,312,213,8,17).

This is a real zero-count two-edge repair, with common frame
(1,2,5,7,4,6,3,0). It restores the receiving mass213 of the transported
T41/T43 endpoints, but changes other interface data. The same one/two-edge
repair search finds no such V repair for the next two transported T45 inputs
(H273 and501); that statement is limited to the exact local neighborhood.

There is an important distinction between this finite reset and a repeatable
corridor. The path157->159->R159 has receiving values213,215,213. Applying
receiving-level translation to the whole path does NOT preserve its first
edge. The direct157->R159 transition flips THREE old arcs,
{2,3},{3,7},{5,6}. It has common frame(0,1,2,4,6,3,5,7) and count difference2,
so it is a legal free-tile move only in an explicitly ENLARGED graph permitting
three simultaneous free-arc flips. It is not silently inserted into Session275's
one/two-tile move set.

This produces an exact same-V repaired pair in that enlarged move language
and a finite way around the fusion-only obstruction. It does not prove an
unbounded rearming law or an inter-scale overlap. Its surviving question is
whether controlled larger moves can repair the whole interface in a recurring
pattern. `repair_probe.py`/`repair_receipt.jsonl` preserve the search universe;
`escape_check.py` independently checks the displayed positive repair.

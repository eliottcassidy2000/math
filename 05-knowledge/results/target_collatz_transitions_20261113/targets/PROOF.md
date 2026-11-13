# A ternary fork enters an actual reusable two-clock chart

**PROVED:** the arbitrary-module core replacement, retained-name map,
complete clock-root response, and guarded clock continuation below.
**FINITE-EXACT:** the small literal controls in `receipt.json`.
**OPEN:** exact all-target coverage and a return between all the source
charts. The count recurrence is not a proof of either ordinary Collatz
or Chan–Pak.

## 1. Inherited interfaces and the new handoff

Let H count Hamilton paths, and let P_K(d) count spanning directed-path
covers in d ordered slots, allowing empty slots. For nonempty K write

    h=H(K), p=P_K(2), w=P_K(3).

The incoming [CORE_REPLACEMENTS §2](../../prime_ternary_lift_20261113/CORE_REPLACEMENTS.md)
has a common source frame S and two ternary target frames:

    S=(2,20,9,35,45,7), module at0,
    H(S[M])=15P_M(2)−H(M),
    H(B2[M])=45P_M(2)−H(M),
    H(B4[M])=45P_M(2)+H(M).

For transitive M these are the exact +2 and +4 children. For arbitrary
M the digits are2H(M),4H(M); that qualification is inherited and retained.

The older two-clock Π from
[CLOCKS §§3–4](../../clone_clocks_negative_collatz_20261113/CLOCKS.md)
has vertices `(a,v,c,u,d0,d1,x0,f0,f1)` in that order. All arcs point
forward except

    d0→a, d1→a, f0→c, f1→c.

Its rows are `(462,508,120,496,481,449,384,260,4)`, and its universal
two-module identity is

    H(Π[u←M,v←N])=60P_M(3)P_N(2)+H(M)H(N).         (1)

Both the two-clock and its one-arc transfer are inherited. The new
connection is a physical core rebase from the incoming +4 fork into it:

    S[TT2→K]  ->  Π[u←TT1,v←K],
    60p−h     -> 180p+h =3(60p−h)+4h.              (2)

For every nonempty K the move below retains its entire induced module
and the old ready-root **name**, while changing the interface type. It
uses two fixed deletions, three new vertices, three fixed old-arc
reversals and three whole module-incidence reversals. Its order grows
by exactly one. This is **not** an induced one-vertex extension.

The intermediate B4[TT2→K] also has count180p+h. Thus the last rebase is
count-neutral for every module, not merely at a transitive sample. This
is the exact source-family handoff absent from a comparison of the two
scalar fork values alone.

## 2. All physical labels and arc changes

Keep the source frame labels1,...,5. Its module0 is the ordinal sum
`A→B→K`, where A,B are the two singleton prefix vertices. Delete B and
source vertex5. Retain A, every vertex of K, and source vertices1,2,3,4.
Add three new singleton vertices x0,f0,f1. Identify the retained names
with the target Π as follows:

| target name | retained source name |
|---|---|
| a | A |
| v module | K, with its internal labels unchanged |
| c |1|
| u |4|
| d0 |3|
| d1 |2|

Reverse exactly these pairs of retained fixed vertices:

    {2,3}, {1,3}, {4,A}.

Also reverse the three uniform incidences

    2→K, 3→K, 4→K

into K→2,3,4. Give the three new vertices exactly their Π incidences.
Every remaining retained-old arc already agrees with Π. Reading the six
rows of S proves this directly; the script checks the complete named arc
set for each literal control. The number of reversed arcs between old
vertices is3+3|K|. The number of changed incidence **classes** is constant;
these are different notions of cost.

This map can also be seen through B4. In its singleton reference picture
the apex word is `1001 00 K 110`, with K in the zero block. For |K|=1
label its spine0,...,9 and apex10, so K=6. Delete2,5 and use

    (a,v,c,u,d0,d1,x0,f0,f1)=(4,K,10,0,1,3,7,8,9).

Reverse fixed pairs{0,4},{1,10},{8,10} and the0,1,3→K stars. This is
the neutral B4-to-Π step. It is a quotient-core replacement whose only
variable object is the retained K; it does not rebuild K from its count.

Strongness holds throughout the two specified endpoint constructions.
The source and target frames are strong; inflating a vertex of a strong
tournament by a nonempty tournament leaves the result strong, since
exterior directed paths connect every pair of vertices in the module.
No claim that arbitrary intermediate single flips remain strong or
count-neutral is needed.

## 3. What happens to the ready root and its missing sidecar

Source vertex4 has unique predecessor1. Its full source state after
substituting TT2→K is

    (H,s,q,F)=(60p−h,12p+2h,12p+3h,60p).           (3)

Here q=H after deleting the root, s counts starts there, and F counts
ordered spanning covers(P,Q) with P nonempty ending at the root and Q
possibly empty. Thus it is h-ready and is unit-ready for transitive K.
This profile is inherited from
[H9_INFLATION §5](../../h17_exact_repairs_20261112/constructors/H9_INFLATION.md).
To compare its frame G with S explicitly, the isomorphism is

    (G0,G1,G2,G3,G4,G5) -> (S1,S0,S5,S2,S4,S3),

preserving its module and ready-root marks; then P_(TT2→K)(2)=4p.

The handoff retains this physical vertex as u and its old predecessor
as c. It reverses the incidences from u to A and K, so u now has
indegree |K|+2. The old unit-ready root star is deliberately surrendered.
Its complete new state is

    (180p+h, 57h, 60p+h, 186w+48p−53h).             (4)

In particular, for transitive K the target is primitive at this retained
root: H−3q=−2 and both H,q are odd, so gcd(H,s,q,F)=1. The graph therefore
has both an actual two-clock module interface and a primitive named-root
interface. It is not asserted to remain unit-ready. At K=TT2 the concrete
transition is

    (239,50,51,240) -> (721,57,241,1813).            (5)

The new F depends on the third-bin response w. Thus a claimed general
closure on the two numbers(h,p) would discard a coordinate that future
root operations really observe.

### Complete finite coefficient proof of(4)

The K module in Π has exactly one exterior predecessor, a. In a Hamilton
path it has at most two runs: a can supply at most one reentry. In a
two-path cover it has at most three: both paths may start in K, with at
most one additional entry through a. This bounds the entire universal
run expansion, even when K is not transitive.

With K_i counting ordered covers in i nonempty slots, the exact exterior
skeleton coefficients are

| observation at u | K1 | K2 | K3 |
|---|---:|---:|---:|
| H |361|180|0|
| s |57|0|0|
| q |121|60|0|
| F |601|606|186|

`universal_u_certificate()` counts these fixed finite skeletons directly.
Its states are the used exterior subset, number of K placeholders and
last symbol. Equal adjacent placeholders are forbidden. For F it sums
over two complementary exterior subsets, with the first word ending at
u; the second may be empty. The ordered assignment of K paths to these
placeholders is exactly the definition of K_i, with no extra factorial.
The at-most-three proof makes the table complete; no polynomial fit is
being extrapolated. Substituting `K1=h`, `K2=p−2h`, `K3=w−3p+3h` proves
(4). The literal subset-DP controls are an independent implementation.

## 4. Continued exact transitions after the handoff

For transitive K=TT_k, k≥1, the source in(2) is60·2^k−1. Its +4 child
has been rebased into the actual two-clock

    H(j,k)=1+60·3^j·2^k,      j,k≥1,                (6)

starting at j=1. Ternary and dyadic module growth give the inherited
exact maps3H−2 and2H−1, and the inherited single-arc transfer gives

    (j,k)->(j+1,k−1),  H->(3H−1)/2,  k≥2.         (7)

The transferred vertex is the last vertex of the v block; reverse its
arc to c and reinterpret it as the first u vertex. Every other vertex
and every other arc stays fixed. Thus, for example,

    239 -> 721=7·103 -> 1081

consists of the explicit core handoff followed by one arc reversal.
This is a proved change of construction language, not simply a sequence
of attained scalar values. The k=1 boundary is treated separately in
the sibling clock packet: deleting the entire v module cannot be
justified by substituting the empty tournament into(1).

There is a useful retained-core generalization. Let K=R→TT_k, where R
is any fixed **nonempty** tournament, with h=H(R),p0=P_R(2). Then

    H(j,k)=h+60·3^j·2^k p0,      j≥1, k≥0.         (8)

The same transfer is valid for every k≥1, including the last tail
vertex, because the v module R remains nonempty. Its actual law is
`H→(3H−h)/2`. The moved vertex is dominated by every other v vertex,
including all of R, so its only required change is still its arc to c.
All R vertices and internal arcs are retained. This supplies a graph
interpretation of the anchor h; it is not a free numerical shift.

## 5. An exact finite target selector, and the separate return gate

For a prescribed positive odd target Y, membership in the **destination**
chart(2) is equivalent to the following finite test. Enumerate

    1≤h≤floor(Y/361), h odd, h≡Y (mod360),
    p=(Y−h)/180.

For each candidate enumerate actual nonempty tournaments K with

    |K|≤min(floor(p/h)−1, floor(log2 p)),
    (H(K),P_K(2))=(h,p).                            (9a)

Both order bounds are necessary: the inherited path-cut injection gives
p≥(|K|+1)h, and every one of the2^|K| induced-subset summands in p is
positive, giving p≥2^|K|. The first also gives Y≥361h. Slot-swapping is
fixed-point-free for nonempty K, so p is even and h≡Y modulo360.
Conversely, any accepted actual K supplies both the source60p−h and
the physical target transition(2), with all labels computable by§2.
This is an exact inverse for this constructor, not an assumption that
every arithmetic candidate is a physical response.

Its endpoint primitivity test is especially simple:

    gcd(H,q)=gcd(180p+h,60p+h)=gcd(h,15p).           (9b)

Subtract three times q and use that q is odd, or reduce gcd(q,120p)
by its power-of-two factor. Thus gcd(h,15p)=1 guarantees the full rooted
state is primitive. Every prime target in this chart is automatically
primitive because0<q<H. For composite targets, (4) supplies the full
test; failure of(9b) alone is not an obstruction. For instance K=Q4
has(h,p,w)=(5,32,117), target state(5765,285,1925,23033), and full gcd1
even though gcd(H,q)=5.

The same finite selector covers the full reusable chartΠ[TT_j,K]. For
each j≥1 with120·3^j+1≤Y, enumerate

    h≤floor(Y/(120·3^j+1)), h≡Y (mod120·3^j),
    p=(Y−h)/(60·3^j),

and apply exactly the physical-module test(9a). Every accepted K supplies
the source60p−h, the handoff(2), and then j−1 literal u-module growths
to the prescribed target. The bound on j follows from h≥1,p≥2, so the
whole search is finite. Its output is an actual graph and construction
address; no efficiency claim is made for unrestricted target sizes.
The supplied smaller source has the exact value
`S=(Y−(3^j+1)h)/3^j=60p−h>0`. Thus an accepted target has a genuine
decreasing parent with a physical transition certificate. Iterating this
descent still requires the next source to pass an appropriate chart test.

`select_target()` implements this finite inverse. For H(K)=1 it uses
the inherited equivalence K transitive, so p must be a power of two.
For other candidates it enumerates the bounded labelled graph universe.
The controls include Y=1081: its j=1 arithmetic candidate(h,p)=(1,6)
has no physical module, whereas j=2 gives K=TT1 and the actual path
119→361→1081. Thus the selector distinguishes the response-set test from
an integer solution. At every j the first u vertex has deletion count
`q=60·3^(j−1)p+h`; for j≥2 this deletes a transitive clone and for j=1
it is the separately certified formula(4). Hence H−3q=−2h, preserving
the sufficient primitive guard gcd(h,15p)=1, and prime targets again
have0<q<H.

The handoff is universal over its stated input family and the clock
really can continue. It does not automatically return to the handoff's
source family. For a prescribed Y, such a return requires an actual
nonempty module L with

    Y=60P_L(2)−H(L).                               (9)

Writing h_L=H(L), the path-cut injection gives
P_L(2)≥(|L|+1)h_L. Thus(9) is a finite physical-source test:
enumerate positive odd h_L≤Y/119, set p_L=(Y+h_L)/60, and test actual
modules of order at most floor(p_L/h_L)−1 with response(h_L,p_L).
Integrality and evenness of p_L are necessary, not substitutes for the
physical module test.

There is an immediate hostile to an automatic return. Every count in
(6) is1 modulo120. Because P_L(2) is even, (9) would force
h_L≡119 modulo120. Hence h_L≥119, so L has at least three vertices
and P_L(2)≥4h_L. Consequently

    every return from(6) through(9) has Y≥239·119=28441. (10)

In particular neither721 nor1081 can be silently reclassified as a
member of the handoff source family, even using a completely different
module. This is a restriction of this chart, not of Hamilton-count
attainment. Larger returns remain an explicit response-set problem.

The inherited finite ready deficit and fixed homogeneous-cloning
obstructions motivated changing the core. The underused positive object
was Π's two-module interface. Their common scalar response180p+h reveals
the exact rebase, while(4) identifies the information needed beyond that
scalar equality. These are the preserved and lost objects:

| operation | preserved | changed / guard |
|---|---|---|
| core handoff(2) | all K vertices/arcs; ready-root and predecessor names | two fixed deletions, three insertions; old root star changes |
| neutral B4 rebase | exact H; K itself | fixed core, interface becomes(u,v) |
| clock transfer(7),(8) | every vertex and all but one arc | one tail vertex moves between module roles; v must remain nonempty |
| attempted return(9) | prescribed scalar target | requires an actual bounded module response, not an arbitrary integer pair |

## 6. Reproduction

    python3 -B 05-knowledge/results/target_collatz_transitions_20261113/targets/handoff.py

The receipt contains five literal module handoffs (TT1,TT2,TT3,C3,Q4),
the exact old-to-new names and changed arcs, complete source/target root
states, the universal run table, and three one-arc transfers. The C3
control has digit12 rather than4, explicitly guarding the arbitrary-
module qualification. The R=C3,k=1 transfer tests that the final tail
vertex can move when a nonempty core remains. No large census or all-
target verification is claimed. Inherited sources are pinned to the
snapshot recorded in the script.

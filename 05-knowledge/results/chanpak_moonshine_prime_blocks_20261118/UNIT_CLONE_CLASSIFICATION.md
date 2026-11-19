# Exact classification of unit transitive clones, and the marked descent boundary

**All-order:** a transitive vertex-clone adds exactly two Hamiltonian paths only in the padded-triangle family. After that move, rerooting cannot make another unit clone. **Inherited:** general reader substitution and odd-cycle deletion laws. **Not asserted:** impossibility of other physical +2 moves, changing old arcs, selective contacts or universal Chan–Pak coverage.

## Fresh random inspiration

A new independent Random(191137) draw from sorted `07-reflections/*.md` files, excluding filenames containing tournament, chanpak, hamilton or metagraph, selected:

- `qr-resonance-principle.md`;
- `lrc14-p5-owner-permutation-apex-obstruction-20260902.md`;
- `the-mod-6-acceleration-law-and-a-complete-small-K-tight-census-kps-S137.md`.

The LRC owner note distinguishes total capacity from where owners actually meet; the acceleration note distinguishes a locally promising operation from reusable closure. They motivated testing whether a successful unit clone remains ready after rerooting. Their numerical or broad conjectural assertions are not premises of the theorem below. No claims from the QR-resonance reflection are imported.

## 1. Physical definition and exact reader

Let T be a tournament with marked vertex v. A transitive clone replaces v by adjacent twins v1->v2 with identical contacts to every old outside vertex. More generally let h(m)=H(T[v<-TT_m]). The physical Stirling expansion from `random_reader_interpolation/EXACT_READER_INTERPOLATION.md`, a specialization of THM-3121/3134, is

    h(m)=sum_r c_r r! S(m,r),

where c_r counts outside-vertex orders with r v-slots forming a directed word. In particular c_1=H(T), and h(2)-h(1)=2c_2.

When an outside vertex is inserted into a directed path in a tournament there is always a legal slot. Choose one deterministically. Distinct original paths remain distinct after insertion, since deleting the added vertex recovers the originals. This remains true for paths containing both distinguished twins and preserves which old vertices lie between them.

## 2. Complete unit-clone theorem

The following are equivalent:

1. The clone adds exactly two paths.
2. c_2=1.
3. T is the ordinal composition TT_a -> C3 -> TT_b, with v in the triangle, for some a,b>=0.

In this case H(T)=3 and the child has H=5.

**Proof.** If c_2>0, its two-visit word contains a segment v->u_1->...->u_k->v. Its first outsider is beaten by v and its last beats v. The first transition from the former side to the latter supplies a directed triangle v->u->w->v. Conversely such a triangle supplies the path v1,u,w,v2, which extends to a full Hamiltonian path after inserting all remaining vertices, and hence supplies a two-visit skeleton.

Assume c_2=1 and fix that triangle. If an outside vertex t has nonuniform contacts with its three vertices, then t has at least two legal insertion slots in the path v1,u,w,v2. To see this, write the four bits recording whether t beats each path vertex. The first and last bits coincide. If they are1, the initial slot is legal and nonuniformity gives a0-to-1 transition, another legal slot. If they are0, the terminal slot is legal and nonuniformity again gives a0-to-1 transition. Both resulting paths retain separated twins, and inserting every other outsider preserves their distinct restrictions. This gives two different two-visit skeletons, a contradiction. Thus the triangle is a module.

Let I be the outsiders beating the entire triangle and O the outsiders losing to it. If some o in O beats some i in I, two distinct directed paths in the two-twin inflation of {v,u,w,i,o} are

    i,v1,u,w,v2,o
    w,v1,o,i,v2,u.

The original u,w lie between the twins in the first and outside the twins in the second. Both paths are valid from the declared directions. Extend with remaining vertices to contradict c_2=1. Hence I->O and T=I->C3->O.

Ordinal multiplication gives h(m)=H(I)(2^m+1)H(O), so c_2=H(I)H(O). It equals1 only if both factors have count1. A tournament has count1 exactly when it is transitive (for example by the inherited odd-cycle formula), so I and O are transitive. The converse follows from the same formula. Empty I/O have count1. QED.

**Rerooting consequence.** The child has count5, so no choice of its marked vertex can admit a +2 transitive clone: the theorem requires its source count to be3. This is an all-order failure of renewal, not merely failure at the inherited mark. It does not prevent a different one-vertex attachment, a coordinated arc edit, or a change of source.

## 3. Exact marked transport, not just loss of readiness

Keep either clone as the distinguished vertex after one cloning step. Inflating it by k vertices gives the same transitive module as inflating the original mark by k+1. The Stirling recurrence therefore gives the full profile update

    c'_j = j c_j + (j+1)c_(j+1),       j>=1,

with absent coefficients zero. This is a physical all-height marked return rule. In particular c'_2=2c_2+3c_3. A ready c_2=1 becomes at least2, so the retained port's next increment is at least4. Repeated cloning cannot hide the growing increment in an unchanged endpoint scalar. The unit-clone theorem strengthens this to show that no alternative root repairs the first unit step either.

## 4. What exact deletion descent really guarantees

THM-026 already gives, from the odd-cycle formula,

    H(T)-H(T-v)=2 sum_(directed odd cycles C containing v) H(T-V(C)).

The cycles are counted up to cyclic rotation, with their directed orientation retained; empty complements have H=1. Each packing containing v has a unique cycle containing it, yielding this identity. The descent is strict exactly when v belongs to a nontrivial strong component, equivalently a triangle. It is not a target-first source-selection theorem.

A unit deletion therefore has an exact physical certificate: there is precisely one directed odd cycle through v, it is necessarily a triangle, and its complement is transitive. The uniqueness is of all odd directed cycles, not merely triangles. This unit-deletion criterion is less restrictive than transitive cloning: arbitrary attachments may change contacts rather than duplicate them. Thus the clone obstruction must not be promoted to a general +2 obstruction.

## Replay and scope

Standalone `unit_clone_replay.py` checks every labelled rooted tournament through order5 by independent subset DP. The numbers of unit-clone roots at orders3,4,5 are6,48,360; all have source count3, child count5, the claimed padded-triangle structure, and no ready child root. It also checks the full marked profile transport at both retained twins for every labelled rooted source through order4 using direct slot-word enumeration. Normal and optimized outputs agree. The finite census supports implementation; the classification proof is all-order.

A targeted search of available local canon found the inherited substitution and deletion mechanisms but no exact unit-clone classification under that description. This is not a comprehensive historical novelty claim.

# Receiver rearming by core surgery, then a six-count carry

**Session280. PROVED in this note from the inherited physical fusion law;
literal finite controls are separate. Chan–Pak remains OPEN.**

## 1. Inheritance and the construction that must be repaired

[Session276 corridor §6–7](../../session276_corridors_20261112/corridor/PROOF.md)
transports the physical 41/43/45 corridor but exposes a receiver leak: its
upper endpoint has deletion count13 rather than15. The old proof excludes
repair by common unit-deletion fusion sandwiches. We first strengthen that
obstruction to arbitrary physical fusion contexts, then deliberately leave
that operation language by changing internal arcs. The resulting repaired
bank supports an infinite six-count family, not merely a finite numerical run.

The inherited rooted-fusion carrier from
[Session264](../../session264_constructive_20261111/fusion/PROOF.md) is

    M(T)=[h,U-h,d-h,H-h; 0,Q,0,V-h; 0,0,k,e-h; 0,0,0,h].

Here H is Hamilton count, h the root-deletion count, Q the ordered two-path
cover count of the deletion (empty slots permitted), V the weighted root-ending
path mass, and e its odd-minus-even counterpart. H(empty)=1. Physical fusion
identifies roots and directs all nonroot left vertices to all nonroot right
vertices. Matrices multiply in that order. The singleton is the identity.

The source endpoints are

    T41=(56,13,25,48,34,6), T43=(52,13,24,49,34,6),
    T45=(52,13,24,33,42,6),

outgoing row masks, root0. Their receiving data are

    (H,h,Q,V,e): (41,15,120,57,9), (43,15,120,57,9),
                  (45,13,104,57,9).

No novelty is claimed for these seeds, the fusion law, or the elementary
endpoint-cover injection used below.

## 2. A universal receiver-gap lemma

For any root-deleted tournament D, let E_i be the sum of e_i(S)H(D\\S)
over all nonempty S containing i, where e_i counts paths ending at i.
Every Hamilton path of D splits uniquely immediately after i, yielding one
such two-slot cover. Hence E_i>=h=H(D). Summing over all i gives Q-h.
If I is the set of vertices beating the root,

    V=h+sum_(i in I)E_i,      Q-V=sum_(i not in I)E_i.       (1)

Consequently:

* V=h iff the root is a source; otherwise V>=2h;
* Q=V iff the root is a sink; otherwise Q-V>=h;
* always 0<V<=Q.

The empty deletion is included: h=Q=V=1, and the singleton is both source and
sink. These are exact positive combinatorial bounds, not signed estimates.

## 3. Complete same-root sandwich obstruction, without the unit restriction

Let L,R be ANY actual pointed tournaments. Fusion gives

    V(L*T)=Q_L V_T+(V_L-Q_L)h_T.

For the difference T45−T43, Δh=−2, ΔQ=−16 and ΔV=0. Therefore

    ΔV(L*T*R)=2(Q_L−V_L)h_R−16Q_L(V_R−h_R).              (2)

If R is not a source, (1) gives V_R−h_R>=h_R, while Q_L−V_L<Q_L.
Thus (2)<−14Q_L h_R<0. If R is a source, it vanishes iff Q_L=V_L,
which means L is a sink. Thus receiving equality is restored **iff L has a
sink root and R has a source root**.

Those exceptional contexts do not perform a hidden repair: the physical
output is the ordinal sum (L−root) -> T -> (R−root), so

    H(L*T*R)=h_L H(T)h_R.

The output count difference is2h_Lh_R. To retain a +2 gate, both positive
integers h_L,h_R must be1; their deletions are transitive and the contexts
are neutral global-source/global-sink padding. Thus allowing arbitrary
nonunit contexts does not evade the old obstruction. A sequence of common
same-mark left/right fusions collapses to a sandwich by associativity, so
this is a program-level statement. Root changes, endpoint-dependent contexts
and internal-arc changes are not excluded.

More generally, for a same-V gate with Δh=−a, ΔQ=−b, 0<a<=b,

    ΔV=a(Q_L−V_L)h_R−bQ_L(V_R−h_R),

and the same extremal-root classification holds. The conclusion for
b<a is not asserted. This threshold retains the magnitude of the structural
change rather than treating all nonzero defects as equivalent.

## 4. Positive repair: restore the continuing coordinates, not just H

Let C_t (t>=1) be a pointed one-apex tournament with transitive deletion
spine of length t+1, root beating only its first vertex. It has inherited data

    h=1, Q=2^(t+1), H=U=1+2^t,
    V=2^(t+1)−1, d=e=1.

For every pointed T,

    h'=h, H'=H+2^t V, V'=2^(t+1)V−h, e'=h.             (3)

Take t=2. The three outputs C_2*T41,C_2*T43,C_2*T45 have counts
269,271,273. The first two have(h,V,e)=(15,441,15), but the upper endpoint
has(13,443,13). This is precisely the leakage, not an erased coordinate.

In C_2*T45 reverse the four arcs on unordered pairs

    {0,5}, {2,4}, {4,7}, {4,8}.

The result is the actual nine-vertex tournament

    R273=(386,508,489,497,485,193,257,320,32),
    (H,h,Q,U,V,d,e)=(273,15,768,441,441,9,15).          (4)

Thus the bank

    B0=C_2*T41, B1=C_2*T43, B2=R273

has counts269,271,273 and SHARED(h,V,e)=(15,441,15). The prior repair at159
restored only V and leaked again; this repair restores the entire continuing
receiving triple needed below.

This is a finite certified seed, not a claim of uniqueness or minimality.
The deterministic search considered every four-edge reversal of C_2*T45;
there are58,905 candidates and three matching outputs. Only the displayed
witness is needed by the theorem.

There is also a physical count-corridor interpretation in an explicitly
expanded move language. B0→B1 reverses two old root-incidence arcs {0,5},{0,6}, and B1→B2
changes five arcs {0,5},{2,4},{4,7},{4,8},{6,7}. A common Hamilton frame for
the latter is(0,1,2,3,4,6,8,5,7), avoiding every changed arc. Hence this is
an actual +2 free-arc move if five simultaneous free-arc flips are allowed.
It is NOT a one/two-tile corridor edge. The four-edge273→R273 repair also
has a common frame(0,1,2,3,4,5,7,6,8) and zero count change.

## 5. Rearmed transport under arbitrary left unit contexts

For any common pointed L with h_L=1, matrix multiplication yields

    H(L*T)=H_T+(H_L−1)h+(U_L−1)(V−h)+(d_L−1)(e−h).

It depends on T only through H,h,V,e. The receiving updates are

    h'=h,
    V'=Q_L V+(V_L−Q_L)h,
    e'=k_L e+(e_L−k_L)h.                               (5)

Therefore every common left unit-deletion context translates the three bank
counts by EXACTLY the same amount and retains common h,V,e. This can be
iterated without limit. It is no longer the old common-V-but-different-h
family whose very next transport fails. No equality of Q,U,d is required.

For repeated C_1 write the resulting bank B_i^(k), k>=0. Equations(3) give

    h_k=e_k=15,
    V_k=436·4^k+5,
    H(B_i^(k))=269+2i+(872/3)(4^k−1)+10k.              (6)

The displayed quotient is integral. Every B_i^(k) is a literal tournament
with9+2k vertices. Common old frames and free arcs survive these C_1 transports:
prefix the two new spine vertices to the old frame; the second beats every
old vertex including the root. Thus the enlarged five-free-arc corridor
also persists. Counts grow without bound but width remains4; no connection
between successive scales is claimed.

## 6. Combine rearming with an exact carry: six consecutive odds infinitely often

For P>=1 let L_P be the one-apex transitive-deletion context whose root-incidence
word is the usual binary expansion of P, read along the increasing spine.
Set L_0 to the singleton. The inherited one-apex recursion is

    g(0)=1, g(2v+1)=g(v),
    g(2v)=g(v)+v+(v mod2), v>=1.

Here H(L_P)=g(P), U(L_P)=P+1, h_L=1. For the bank in(6), e=h=15, so

    H(L_P*B_i^(k))=H(B_i^(k))+T_k(P),
    T_k(P)=15(g(P)−1)+(V_k−15)P.                       (7)

Take any EVEN k>=0. Since436·4^k+5≡21 modulo30 at every such k (including k=0),

    z_k=(V_k−21)/15=(436·4^k−16)/15

is a positive EVEN integer. Let P_k=2z_k. The exact even-to-odd carry gives

    T_k(P_k+1)−T_k(P_k)
      =(V_k−15)−15z_k=6.                              (8)

The two three-count banks with contexts P_k and P_k+1 therefore realize

    N_k, N_k+2, N_k+4, N_k+6, N_k+8, N_k+10,
    N_k=H(B0^(k))+T_k(P_k).                            (9)

This is an explicit infinite family of SIX consecutive odd Hamilton counts,
constructed by a reusable receiver repair and a paid binary carry. At k=0,
z=28,P=56, and the six counts are

    24875,24877,24879,24881,24883,24885.

They are independently counted on the six literal fifteen-vertex tournaments
in the verifier. Larger even k use the proved recurrence, not inference from
a numerical range. N_k tends to infinity because V_k and P_k do and all terms
in(7) are nonnegative.

The two translated banks do NOT share their receiving coordinate with each
other. If n is the common bitlength of P_k and P_k+1, their context root words
differ only at the last bit. Their two-cover counts agree, but V of the
context decreases by2^(n−1). Therefore the receiver of the second triple is
exactly15·2^(n−1) smaller. For k=0 the two receivers are28119 and27639,
difference−480. Both retain h=e=15.

No common further physical fusion sandwich, even with nonunit contexts, can
equalize this defect. Compare corresponding seed index i in the two banks:
their h and Q agree, while ΔV=−15·2^(n−1). Left fusion multiplies ΔV by
Q_L and right fusion by h_R, so the final difference is Q_L h_R ΔV≠0.
This is a smaller obstruction than requiring the whole six-bank to match,
and already excludes that match. Nor can independent neutral
global-source/sink padding do so. Source padding multiplies V by a power of2;
sink padding leaves V fixed. The ratio of the two positive receivers lies
strictly between1 and2: writing q=2^n, the larger is at least q(V_k−15)+15,
and its difference from the smaller is15q/2, less than half the larger
since V_k>30. A nontrivial power of2 cannot lie strictly between1 and2.

Thus(9) is an output interval, not a renewable six-way shared-interface
bank. Reapplying the three-bank theorem to this union would be unsupported. The missing constructive obligation is a receiver repair for that
widened bank, or another operation that increases width while retaining the
continuing state. This six-wide family does not prove an eventual tail.

### A general paid width-doubling criterion

More generally, suppose r actual pointed tournaments have counts
H0,H0+2,...,H0+2r−2 and common(h,V,e) with e=h>0. If

    V>h+2r,     V≡h+2r (mod 2h),

then z=(V−h−2r)/h is a positive even integer. The same contexts P=2z,P+1
shift the two banks by exactly2r, realizing2r consecutive odd counts.
The carry incurs receiving defect−h·2^(bitlength(P)−1), with Δh=ΔQ=0
between corresponding outputs. This explicitly exposes both the positive
width-doubling gate and the quantity a new surgery must repair. The lemma is
conditional on actual shared-interface seeds; it cannot supply those seeds
by itself. The construction above meets its hypothesis at every even k with
h=15,r=3. No generic repeated doubling is inferred.

## 7. Connections and failed alternatives

The change of representation is the same discipline seen in the distant
[recognizable-frame/guarded-composition work](../../support_image_composition_20261026/COMPOSITION.md):
retain exactly the coordinates consumed by future operations, not merely
an invariant scalar. Here the map is physical core surgery, the preserved
predicate is H=273, the repaired sidecar is(h,V,e), and the cheap decisive
test is one more C_1 transport. This is a concrete transfer of method, not
an assertion that unrelated models are equivalent.

The second connection is to the exact carry methodology of
[Session276](../../session276_corridors_20261112/carry/PROOF.md): choose a
carry with known correction, then solve its integer equation with a positive
physical parameter. The new source bank changes which equation is solvable.
No residue-only or gcd-only inference is used.

Exploration controls: graph-dependent unit-left contexts of spine length<=10
produced no repaired adjacent T43/T45 pair; three-edge receiver+h repairs at
C_2*T45 failed; the fourth edge supplied the displayed positive seed. These
bounded failures are not all-order impossibility claims. The theorem in§3
has its separate proof and much wider common-context scope.

## Reproduction and scope

Run `python3 verify.py` here (also with `python3 -O verify.py`). The verifier
uses independent literal subset dynamic programming for counts and root
masses, exact integer fusion arithmetic for long recurrences, explicit
exception guards, and hostile one-step leakage controls. `RESULTS.json`
records the resulting counts. Optional `search_repair.cpp` reruns the bounded
surgery search. No floating decision or numerical cutoff proves(9).

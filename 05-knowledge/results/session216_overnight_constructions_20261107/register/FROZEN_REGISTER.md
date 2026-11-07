# Fixed-H5 register completeness and a reset-created adjacent pair

**PROVED elementary construction; separate from the frozen three-core packet.** No Hamilton-count completeness or source-menu width amplification is claimed.

Use the strong core C=(8,1,3,6), with H(C)=5, and a nonempty cut. Let α count ordered spanning covers(P,Q) with Q starting in the cut and P optional, τ require additionally P nonempty ending in the cut, and e count Hamilton paths ending in the cut. A selected/unselected new sink ε gives

    (h,α,τ,e) → (h,2α+εh,τ+ε(α+e),εh).

The finite core-cut values of(α,e), in cut order1..15, are

    (5,1),(6,2),(11,3),(9,1),(14,2),(15,3),(20,4),
    (7,1),(12,2),(13,3),(18,4),(16,2),(21,3),(22,4),(27,5).

These are directly recomputable from the80-state table of the three-core packet.

## 1. An exact unbounded register spectrum

Among all bodies obtained from this core by finitely many sink appends, retaining its nonempty initial cut, the attainable α values are exactly

    {integers α>=5} minus {8}.

The lower exclusion follows because all seed values are at least5 and the maps x→2x and x→2x+5 never decrease x. To reach8 from a non-seed would require predecessor4, impossible.

For positive existence, the seed set already contains all integers5..16 except8 and10, while10=2·5. It also contains21. For an integer α>=17, use α/2 if α is even and(α−5)/2 if it is odd, except that21 is a direct seed. The predecessor is a smaller attainable integer at least5; the only possible descent into8 is α=16 or21, both already seeded. Induction proves the claim. The inverse gives an explicit body on at most4+floor(log2(α/5)) vertices.

The joint(α,e) support is also exact:

- e=0: every even α>=10 except16;
- e=5: every odd α>=15 except21;
- e=1: α in{5,7,9};
- e=2: α in{6,12,14,16};
- e=3: α in{11,13,15,21};
- e=4: α in{18,20,22}.

Indeed any nonempty append word ends with e=0 or5, with the last bit determined by the parity of α. The predecessor α is unique, so the first two lines follow from the register spectrum. The intermediate endpoint counts occur only at the original core.

This realizes an arbitrary register at fixed body Hamilton count5. The closing count2α+5−τ still depends on τ, so register completeness is not Hamilton-count completeness.

## 2. A physical ambiguity creates a persistent adjacent pair

Core cut6 has(α,τ,e)=(15,6,3). Append one selected sink to obtain

    (h,α,τ,e)=(5,35,24,5), on5 body vertices.

Core cut1 has(α,τ,e)=(5,0,1). Append two selected sinks, obtaining successively(5,15,6,5) and

    (h,α,τ,e)=(5,35,26,5), on6 body vertices.

The first construction and the intermediate second construction agree in h,α,τ but differ in e. The final selected sink synchronizes e and turns that lost endpoint difference into a τ difference of2. After every common further append word, both bodies retain equal h,α,e and τ difference2.

For any nonempty tournament U, close either body by the R=empty star with root module U. The two counts are

    α Q_U(2)+(5−τ)H(U),
    α Q_U(2)+(3−τ)H(U),

where τ denotes the first body's value. For transitive U these are adjacent odd counts. At the empty continuation and U a singleton they are51 and49. All such closed outputs are strong: the initial cut meets the strong core, and appended sink vertices remain reachable from that core and reach the closure root.

This is a returning construction on cut-labeled bodies with different orders, not an induced +2 extension and not a scalar operation on an arbitrary realized count. It needs only two finite core states, not a census or a universal-pair assumption.

## 3. What remains

For a fixed α, stripping append bits is deterministic until one chooses to stop at a core-cut value. Beyond α=27, all ambiguity in the suffix-only τ menu therefore comes from finitely many terminal core cuts. Prefix freedom in the three-core grammar changes that picture and remains the useful additional coordinate for inverse counting.

The four-choice H5 menu discovered with prefix/suffix padding has(h,α,e)=(5,3807,5), τ=2982,2984,2986,2988. Its current status is a finite exploratory witness family; the separately frozen star packet's four-choice H9 menu is already fully controlled. Neither menu gives a five-state returning digit bank.

## 4. A same-order pair stable under both-sided padding

There is a stronger exact pair, using the same core cut10 in both bodies:

    L_A = TT_6(cut111011) ⇒ C(cut10) ⇒ TT_1(cut1),
    L_B = TT_3(cut111) ⇒ C(cut10) ⇒ TT_4(cut0011).

Both bodies have11 vertices and share the full coordinates

    (h,q,α,γ,s,e)=(5,4096,3807,3219,5,5),

but τ_A=2984 and τ_B=2986. The two closed counts are4635 and4633. These numbers follow exactly from the ordinal composition law, not interpolation.

For any common cut-labeled transitive prefix and suffix, all six displayed coordinates remain equal and the τ difference stays2. More generally arbitrary common outer bodies scale the τ difference by the product of their Hamilton counts. Thus arbitrary two-sided transitive padding preserves an adjacent closed-count pair. The underlying two bodies have the same order; no statement about the distinct O/Z/J interface is implied.

The local mechanism is particularly small. Before the common prefix111, the pair is

    TT_3(cut011) ⇒ C(cut10) ⇒ TT_1(cut1),
    C(cut10) ⇒ TT_4(cut0011).

Both have(h,q,α,γ,τ,e)=(5,512,223,398,168,5), but their start-cut counts are0 and2. The common selected prefix converts this hidden endpoint difference into a τ difference while synchronizing the total start count. This is an endpoint-sensitive rebase, not a repeatable unit-width enlargement. Its direct local rewrite changes the necessary suffix pattern, so iteration cannot simply be asserted.

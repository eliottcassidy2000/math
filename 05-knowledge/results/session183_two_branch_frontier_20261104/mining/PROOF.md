# Terminal-block return laws with arbitrary preceding marked context

**PROVED:** context-uniform terminal-block +2 construction and exact paired +4 obstruction. **FINITE-EXACT:** independent literal graph controls. This is a scoped extension of the newest common-entry construction, not a closed two-digit grammar or an actual continuation of H1537279.

## Sources and what is new

At pinned commit `3f2fa002801868c280fafed0c05a1d9060df6c86`, the concurrent packet `05-knowledge/results/multiport_defect_grammar_20261103/transport/TWO_DIGIT_ENTRY.md` constructs both digits from B0=TT3 with every vertex marked. It has three actual suffix signatures:

- B0=TT3, W0={0,1,2}: (1,8,7,7,6,1,1).
- B2=TT2 -> C3[2,3,4], W2={0,1,2}, cycle 2->3->4->2: (3,48,21,39,16,1,3).
- B4=TT5, W4={0,2,4}: (1,32,21,21,12,1,1).

The signature order is (h,p,f,g,j,e,s) from Session182 `induction/PROOF.md`: Hamilton count, unrestricted ordered two-slot covers, marked first endpoint, marked second start, their joint, marked Hamilton endpoint, marked Hamilton start. These are actual port sets. The new result is the complete effect of prepending an arbitrary marked context to these three fixed blocks. It identifies exactly how far this pair of local recipes extends, while giving a much larger +2 source family.

## The context laws

Let L be an arbitrary nonempty tournament with arbitrary actual ports U and signature (h,p,f,g,j,e,s). Put D_i=L->B_i, with ports U union W_i. Let A be any nonempty tournament. Form the usual attachment M(A;D_i,U union W_i): A beats D_i, each port beats the root, and the root beats A and the nonports. Then all three attachments are strong, and

    H(M_2) - 3 H(M_0) = 2 H(A) h,                         (1)
    H(M_4) - 3 H(M_0) = H(A)(4h+2j+2e) - 2 P_A(2) f.     (2)

The empty context is permitted in the displayed h/f/j formulas with formal coordinates (1,1,0,0,0,0,0), yielding the original common entry. This tuple is not claimed to be an identity for all seven assembly coordinates; no empty actual tournament is required in the full attachment.

Proof. Every D_i has a marked vertex in its terminal SCC, so the inherited terminal-port criterion gives strongness. For an ordered sum L->B with second signature (H,P,F,G,J,E,S), the inherited exact assembly is

    f_total=f H+p F,
    j_total=j H+h J+e S+g F,
    h_total=h H.

Subtract three times the B0 signature. The B2 difference is

    Delta(H,P,F,G,J,E,S)=(0,24,0,18,-2,-2,0),

and the B4 difference is

    (-2,8,0,0,-6,-2,-2).

Consequently B2 has Delta f_total=0 and Delta(h_total-j_total)=2h; B4 has Delta f_total=-2f and Delta(h_total-j_total)=4h+2j+2e. Substitute these into the inherited attachment identity

    H(M(A;D,W))=H(A)(h_D-j_D)+P_A(2) f_D.

This proves (1)-(2) at every order with no strongness condition on A or L. QED.

## A broader literal +2 family

Take A=TT_t, t>=1, and let L be any transitive tournament, with **arbitrary** port markings U. Since h=1, (1) gives exactly 3H+2. Thus every transitive suffix whose last three vertices are ports has a uniform +2 construction, with no restriction on any earlier port bits.

The same construction actually needs only the final two suffix vertices to be ports: choose the shorter local source TT2/all-marked and replace its terminal singleton by a singly marked C3. The local signature difference from three times the source has Delta h=Delta f=Delta s=0 and Delta j=-2, so the identical proof applies. We retain the three-vertex blocks above to compare exactly with the inherited +4 recipe. The short source and target signatures are respectively (1,4,3,3,2,1,1) and (3,24,9,15,4,1,3), independently included in the controls. A marked terminal singleton instead gives zero offset under this recipe, so the final two marked vertices are a substantive guard, not just a naming convention.

For the three-vertex version, retain A, the root r, L and old block names 0,1,2. Birth 4 with old in-neighbors A union L union {0,1,r}; then birth 3 with in-neighbors A union L union {0,1,2,r}, with 3->4. The old edge 2->3 is part of the second birth. All old names and all old edges remain unchanged. After the first birth the terminal suffix vertex is still marked 2; after the second it lies in the marked terminal triangle. Therefore both completed states are strong. These are direct adjacency births, not asserted Hamilton-complement native programs.

## Exactly when this pair of fixed recipes supplies both digits

Fix L and U and vary the common outer transitive clock t. From (2),

    H(M_4)-3H(M_0)=4h+2j+2e-2^(t+1) f.

This equals the constant 4 at two distinct clocks, hence at all clocks, if and only if f=0 and 4h+2j+2e=4. Because f counts two-slot covers ending at the ports, f=0 iff U is empty: for any port v, the singleton first path v and a Hamilton path on L-v give at least one cover. With U empty, j=e=0. Thus +4 is uniform iff U is empty and h=1. Equation (1) then also gives +2.

For a nonempty tournament, h=1 iff it is transitive. One elementary proof of the only-if direction: a nontransitive tournament has an induced directed triangle with three Hamilton paths. Any directed path can be extended by inserting a new vertex without changing its existing relative order. Iterating a deterministic insertion rule gives distinct full Hamilton paths from the triangle's three distinct orders. Hence h>=3. A transitive tournament has exactly one Hamilton path.

Therefore these **particular paired local recipes**, after an arbitrary fixed preceding marked context, give both clock-uniform digits 2 and 4 exactly when the preceding context is an unmarked transitive block (or is empty). This does not prohibit a different +4 transformation on the larger +2 source family.

The obstruction is a real transported port debt: each earlier marked context contributes the nonzero slope -2f to the proposed +4 offset. Ignoring the joint correction does not remove this debt, and selecting a favorable single clock cannot establish uniformity. The formulas allow isolated-clock successes outside the uniform locus, which are not excluded here.

## Independent finite controls

`verify.py` computes path endpoints by directed subset recurrence, starts by the converse recurrence, and all seven signature coordinates by explicit subset partitions. It does not use the ordered-sum formula to obtain the counted signatures. It checks all labelled preceding contexts through order 3 with every port set, the three local suffixes, literal full Hamilton counts for outer TT1, TT2 and C3, and strongness by reachability. Explicit exception gates remain active under Python -O. The receipt records the exact finite universe; these checks support, but do not replace, the proofs.

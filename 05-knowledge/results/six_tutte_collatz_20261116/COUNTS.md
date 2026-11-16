# The 3→10→5 carry, and the two missing tournament counts

## The exact graph map

Let C_s(n)=n/2 on even n and 3n+s on odd n (s=±1). Set E(n)=2n+1. In the odd-count coordinate h=E(n),
`E C_s E⁻¹(h) = (h+1)/2 on h≡1 mod4; 3h+2s−2 on h≡3 mod4.`
This is a conjugacy of arithmetic graphs, not a tournament construction. The plus sheet sends 3→10→5 to **7→21→11**. It identifies the odd source and the carried even state with the holes 7 and 21. The other endpoint 5 itself encodes 11, which is attained. Thus {3,5} is the endpoint pair of a two-step segment, not a pair of forbidden counts under an isomorphism.

For odd x=2k+1, x under the rising rule 3x+1 and 3k+2 under doubling have common image 6k+4. In count coordinates their two incoming legs have coordinates
`h=4k+3, j=6k+5, y=12k+9=3h=2j−1.`
So the exact repair objective is to get a witness at 2j−1 from a witness at j in this progression, with the necessary exception j=11. Ignoring the direction of actual Collatz time, the diagram has one physically automatic leg: ordinal sum with a directed triangle multiplies the Hamiltonian count by three. The other leg cannot be universally implemented, because its input 11 is attained and its output 21 is forbidden. Applied to a realized count with h≡3 mod4, the two-step arrow can always realize its *intermediate* 3h; realization of the landing state (3h+1)/2 is a different obligation. Conversely, 63 is attained and ordinal sum proves every 7·3^j for j≥2 is attained; division by 3 is not sound (63/3=21). Iteration of pure multiplication by 3 is not iteration of Collatz: 21 takes the halving branch.

On the known component of 1 in the ordinary positive graph, removing 3 and 10 splits it into exactly three pieces, with boundary vertices 6,20,5. These pieces are respectively the ray of multiples 6·2^i, the full backward basin of 20, and the piece containing the 1–4–2 cycle. Indeed the only outside neighbors of {3,10} are those three, and the connected functional component has only its original cycle. The 20 side branches, e.g. 13→40→20, so the two incoming chambers cannot be exchanged by a rooted symmetry. In count coordinates the ports are 13,41,11. This statement only concerns the actual component of 1 and needs no assumption about other components.

## An unconditional three-branch descent criterion

Overlay the plus and minus sheets. On even n halve; on odd n choose the unique sign making 3n+s divisible by 4 and divide by 4. Then every n>1 decreases and 1 is fixed. This proves connectivity for the union without proving it for either single sheet. In the count chart, the decreasing map D is
* h≡1 mod4: (h+1)/2;
* h≡3 mod8: (3h+3)/4;
* h≡7 mod8: (3h−1)/4.

The complete incoming boundary of {7,21} for D is {13,27,41}: D(13)=7, D(27)=21, D(41)=21. Check even inverse 2d−1 and rising inverses (4d−3)/3 and (4d+1)/3, keeping their displayed residue guards. Thus **Chan–Pak would follow from an inverse realization lemma for the three guarded branches, with only the target exceptions {7,21}**: whenever both d and its proposed inverse target are admissible and there is a physical witness at d, construct one at the target. Seed H=1 and H=3 and the three immediate boundary counts; then ordinary strong induction applies. Start masks for 13,27,41 are (12,21,24,18,1), (30,52,56,50,32,1), and (14,28,56,48,33,3), respectively; bit j in row i means i→j. This is a sufficient repair program, not a solution: it asks for just three branch types but does not supply their universal physical lifts.

Compare [four-vertex Collatz sheets](../procgen_tourn_20260924_four_vertex_sheets_redei.md) and the [recent cubic-radix classification](../decimal_1024_extensions_20261116/README.md). The arithmetic transform here does not replace the structural nonexistence proofs of 7 and 21.

## A literal converse pair + fixed class, and a completed symmetry ray

The [fixed-count multiplicity archive](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md) proves there are just two strong H=63 classes globally, exchanged by converse; it proves there is no self-converse H=63 at any order. Paley-7 is a strong self-converse H=189 class (multiplying labels by a nonsquare reverses its arcs). Together these form a genuine 2+1 converse action on three chosen classes. They do not exhaust H=189, and they do not turn forbidden counts into objects. The two 63 classes and the Paley class can even be placed at order 9 by padding the former with a source or a sink respectively, and the latter with one on either side.

New strong self-converse center H=567 on nine vertices: outgoing masks
`(50,252,225,501,420,128,177,257,103)`.
The anti-automorphism is i↦8−i; the starting counts are (33,81,53,224,50,18,51,18,39). The accompanying checker verifies the count, converse and strongness independently. This supplies the missing even-exponent center in the archive and closes the growth-degree question there:

* H=7·3^a is self-converse-realizable exactly for a≥3 (a≥0). For odd a use center 189 and a mirrored pair of triangle chains; for even a use 567.
* For a≥3, the degree of the eventual self-converse multiplicity quasipolynomial as vertex order grows is exactly floor((a−3)/2), in the sense of maximum degree over the two parities. To see this, in a self-converse SCC word the component containing 7 has to be the center, and all other strong components pair with their converses. Each nontrivial component in the first half uses at least one factor 3, while the center exponent is at least 3 (because 7,21 are absent and 63 is not self-converse), and has the parity of a. Therefore the number of nontrivial components in the half is at most floor((a−3)/2). The archived word formula counts singleton insertions into t+1 gaps of a half having t nontrivial components by a binomial polynomial of degree t, separately on permitted order parities. The explicit centers and mirrored triangle words attain the bound. This is an all-order symmetry theorem on the entire ×3 ray, not full Chan–Pak coverage.

The odd-exponent ray actually attains that degree on both parities. Alongside Paley-7 there is a strong self-converse H=189 tournament on eight vertices, outgoing masks (148,125,120,241,96,193,1,86), reversed by i↦7−i. Starting counts are (7,68,29,33,13,13,7,19). Using either center in the mirrored triangle construction supplies the two parities. The even-exponent ray occurs only at odd orders: by the inherited even-order self-converse congruence (repo THM-791), H≡1 mod4 at even order, whereas 7·3^(2m)≡3. For the congruence, a self-converse tournament admits an involutory arc reverser: automorphisms have odd order as no involution can preserve the edge joining a swapped pair, so an odd power of a reverser yields an involution. At even order it has no fixed vertices (it cannot have two). It pairs directed odd cycles by reversing their directions; an odd vertex set cannot be invariant. The [odd-cycle packing formula](https://arxiv.org/abs/2307.05569) then gives H≡1 mod4, since packings with two or more cycles contribute multiples of four.

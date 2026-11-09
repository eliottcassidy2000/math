# Restriction-stable Hamilton-path gaps and the cost of small symmetry

Session 232, 2026-11-09. **PROVED, conditional only on the named inherited canon (including its finite certificates).** No new proof of the inherited full-automorphism gap theorem is claimed. The accompanying independent small census is **FINITE-EXACT** and is a control, not the all-order proof. Chan–Pak completeness remains **OPEN**.

## 1. Inheritance and scope

For a finite tournament T, let P(T) be its set of directed Hamiltonian vertex words, H(T)=|P(T)|, G=Aut(T), and for any subgroup A≤G put

    q_A(T)=|P(T)/A|=H(T)/|A|,       q(T)=q_G(T).

The free action is inherited: fixing an ordered spanning word fixes every vertex. Rédei oddness makes |A| and q_A positive odd integers. We include the one-vertex tournament; the empty case can be added separately with H=1.

The load-bearing inherited inputs, pinned to commit `4b23190e9c39`, are:

1. [THM-002, odd-cycle collection formula](https://github.com/eliottcassidy2000/math/blob/4b23190e9c39/01-canon/theorems/THM-002-ocf.md):
   H(T)=Σ_C 2^|C|, over vertex-disjoint collections of directed odd cycles of length at least three, including the empty collection.
2. [THM-4812, parts (b), (c), (d), with its audit repairs](https://github.com/eliottcassidy2000/math/blob/4b23190e9c39/01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md): q(T) is never 7 or21; every vertex-transitive tournament with q<23 is one of K1,C3,R5,P7, having (H,|Aut|,q) respectively (1,1,1),(3,3,1),(15,5,3),(189,21,9).
3. [LEM-003, the free action](https://github.com/eliottcassidy2000/math/blob/4b23190e9c39/01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md), and [THM-4747, explicit seven-vertex controls](https://github.com/eliottcassidy2000/math/blob/4b23190e9c39/01-canon/theorems/THM-4747-hamiltonian-cycles-under-automorphisms-fermat-wilson-congruences-for-circulant-and-paley-tournaments.md).

We do **not** use THM-4812(g), whose extension to q=23 has a different audit scope. Our small-symmetry theorem has the strict threshold q_A<23.

Inheritance board: the closest mechanism is THM-4812's orbit-product inequality; the hostile is Aut(P7) of order21 although H(P7)=189 and q(P7)=9; the corrected near miss is the false transfer from an excluded H-value to an excluded automorphism order; the retained sidecar is the subgroup index [G:A] and the actual vertex-orbit partition. The anchor is subgroup-spectrum inheritance, the niche is quantitative prime cost and equality, and the wildcard is transport through a faithful undirected graph encoding (owned by the parallel graph packet).

## 2. Orbit products with their missing index

Let O_1,…,O_s be the vertex orbits of the **full** group G, and put T_i=T[O_i], G_i=Aut(T_i), P=∏_i G_i, d=[P:G], using the injective restriction G→P. Every T_i is vertex-transitive. Positivity of OCF gives

    H(T) ≥ ∏ H(T_i),
    q(T) ≥ d ∏ q(T_i).                         (1)

The first inequality retains just the odd-cycle collections whose individual cycles are contained in one orbit. The second includes the index d, not merely the product of quotient counts. Restricting the free action to A gives the exact formula

    q_A(T)=[G:A]q(T).                          (2)

These formulas retain actual group actions. Knowing only an abstract odd integer or only the number of orbits would not suffice.

## 3. Prime cost below23, with complete equality structure

**Theorem A.** If A≤Aut(T) and q_A(T)<23, every prime divisor of |A| is in {3,5,7}. If

    α=v_5(|A|),       β=v_7(|A|),

then

    q_A(T) ≥ 3^(α+2β).                         (3)

Equality holds if and only if A=Aut(T) and T is an ordinal sum of blocks from {K1,C3,R5,P7}, containing exactly α copies of R5 and β copies of P7. Here the total is still assumed to satisfy q_A<23. There may be arbitrarily many K1 and C3 blocks, in any order.

**Proof of the inequality.** By (2), q(T)≤q_A(T)<23. From (1), each q(T_i)<23, so the inherited vertex-transitive lemma gives only the four listed orbit types. Let c,r,p count C3,R5,P7 orbits. Their full product group has order

    |P|=3^c 5^r 21^p,

and the quotient product is 3^(r+2p). Since A≤G≤P,

    v_5(|A|)≤v_5(|G|)≤r,
    v_7(|A|)≤v_7(|G|)≤p.

The prime-support conclusion follows immediately. Combining (1) and (2),

    q_A ≥ [G:A] d 3^(r+2p) ≥ 3^(α+2β).

**Equality.** Equality forces [G:A]=d=1, r=α, p=β, and H(T)=∏H(T_i). In OCF every term is nonnegative; a directed triangle meeting more than one orbit would contribute a further positive term. Thus there is no cross-orbit directed triangle.

Each T_i is strongly connected or a singleton: strong components of a tournament are linearly ordered, and every automorphism preserves that order, so a vertex-transitive tournament cannot have multiple components.

If B is a strongly connected orbit block and x is outside B, x cannot see B in both directions. Otherwise choose a directed path in B from a vertex u with x→u to a vertex v with v→x. At its first change from a vertex beaten by x to a vertex beating x, an arc y→z gives x→y→z→x, a forbidden cross-orbit triangle. Thus every exterior vertex sees B uniformly. This is automatic for a singleton B. Consequently every pair of orbit blocks has uniform incidence: the preceding argument applied to both blocks rules out variation inside either one.

The resulting quotient tournament on the blocks has no directed triangle, since any quotient triangle lifts to a cross-orbit triangle. A tournament without directed triangles is transitive. Therefore T is an ordinal sum of its orbit blocks, of the stipulated types.

Conversely, every automorphism of an ordinal sum fixes each strong component, since their order is intrinsic. Its group is the direct product of the component groups, and H multiplies because a Hamiltonian word must visit the strong components in their order. Hence A=G and the specified ordinal sum has |G|=3^c5^α21^β and q_A=3^(α+2β). QED.

**Sharp consequences.** In the range q_A<23:

- |A| has no prime divisor≥11.
- 7|A| costs at least9 path orbits; 5|A| costs at least3.
- 5·7∤|A|, 7²∤|A|, and 5³∤|A|, because their respective lower costs are27,81,27.
- If q_A=9 and7|A|, then A is the full group and T is a K1/C3/P7 stack with exactly one P7. In particular its group is F21×C3^c, not just C7.
- If q_A=3 and5|A|, then A is full and there is exactly one R5 among K1/C3 blocks.
- If q_A=9 and25|A|, then A is full and there are exactly two R5 among K1/C3 blocks.

The boundary is genuinely action-sensitive: for P7, the subgroup C7 has q_A=27, while the full group F21 has q_G=9. The subgroup C3 has63 orbits. A denominator chosen without its action/index would lose these distinctions.

## 4. The two gaps survive every automorphism subgroup

**Theorem B.** For every finite tournament T and every A≤Aut(T),

    q_A(T) ∉ {7,21}.                           (4)

This includes actions by an arbitrary finite group through automorphisms: replace it by its image if its action has a kernel. It does not include anti-automorphisms or an action on cycles instead of ordered paths.

**Proof.** Let G=Aut(T). Equation (2) makes q(T) a divisor of q_A(T). If q_A=7, the inherited full-group exclusion leaves only q(T)=1 and [G:A]=7. If q_A=21, it leaves q(T)=1 or3 and [G:A]=21 or7. In every case q(T)≤3 and7 divides |G|. Apply Theorem A to A=G: its lower bound forces q(T)≥9, a contradiction. QED.

The mechanism is not the invalid claim that multiplying an allowed odd number cannot produce7 or21. It is that the potentially dangerous subgroup indices7 and21 are unavailable precisely in the low-q group strata where they would be needed.

If S_H is the ordinary H-spectrum, S_full the q(T)-spectrum, and S_sub the spectrum of all q_A(T), then

    S_H ⊆ S_sub,
    S_full ⊆ S_sub ⊆ positive odd integers minus {7,21}.

The first inclusion uses A={1}, the second A=Aut(T). Consequently ordinary Chan–Pak completeness would imply completeness for the all-subgroup spectrum. No converse and no unconditional completeness are proved here. In particular, the subgroup-spectrum corollary is an omission theorem, not a solution to Chan–Pak.

## 5. Automorphism orders7 and21 are present, and why that is compatible

Let Q7=Cay(Z/7,{1,2,3}) and P7=Cay(Z/7,{1,2,4}). Exact independent counting gives:

- Q7: H=175, Aut=C7 of order7, q=25. Its subgroup path-orbit counts are25 and175.
- P7: H=189, Aut=F21 of order21, q=9. Its literal subgroup counts by order are1:1,3:7,7:1,21:1. The corresponding path-orbit counts are189,63,27,9.

Thus the two missing path-count values are both **present as automorphism-group orders**, already in order7. The ordinary spectrum, full-orbit spectrum, subgroup-orbit spectrum, and automorphism-order spectrum are four different typed objects.

Sharp positive controls also include R5, with (H,|Aut|,q)=(15,5,3); R5→R5 with (225,25,9); and P7→C3 with (567,63,9). The latter illustrates the arbitrary C3 padding in the equality theorem.

## 6. Independent finite certification and graph handoff

Run:

    python3 verify.py > normal.out
    python3 -O verify.py > optimized.out
    cmp normal.out optimized.out

The program is standard-library-only and uses explicit exceptions, not Python assertions. It enumerates every tournament with the canonical Hamiltonian word0→1→…→n−1 for n=1,…,7, then canonicalizes under all score-respecting permutations. Every tournament has a Hamiltonian path, so this universe covers every ordinary isomorphism class. Exact multiplicity of each representative is separately checked to be H/|Aut|; the aggregate mass is2^((n−1)(n−2)/2). Full automorphism groups are found by all degree-preserving vertex permutations. Every literal subgroup is enumerated by closure, and every subgroup orbit count is checked by literal Hamilton words, independently of endpoint subset DP. Known class counts1,1,2,4,12,56,456 are controls, not substitutes for the construction.

The program also independently checks the parallel graph worker's smaller arc encoding: u→v becomes u–L–R–v, with a pendant path of length1 at L and2 at R. The 2-core is the subdivision of K_n for n≥3; original/L/R roles have attached pendant lengths0/1/2. For each seven-vertex control the resulting uncoloured graph has112 vertices,126 edges,42 leaves. A generic paired refinement/individualization graph-automorphism enumerator, with no gadget colors supplied, finds exactly7 and21 automorphisms. All have positive Royle edge-reversal sign. Decoding recovers the tournament, and175/189 admissible original-traversal paths are checked literally. Ordinary graph Hamilton paths are0, because of42 leaves. This check is an independent validation of that packet's graph bridge, not a claim that graph Hamilton-path counts inherit the gaps.

The new statements above use the existing canon rather than the user-supplied preprints. No external preprint theorem is being imported without inspection. The natural conceptual link to restriction/transfer is equation(2) and the prime-cost obstruction to the dangerous indices; any wider comparison belongs to the separate source-intake packet.

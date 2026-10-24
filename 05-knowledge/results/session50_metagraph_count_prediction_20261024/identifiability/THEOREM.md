# What the metagraph determines, and what its aggregates do not

## Status and scope

**PROVED:** the complete characterization of the abstract odd-fibre C,R layer in Section 2; the finite-moment obstruction in Section 4; the merged-degree congruences in Section 5. **FINITE-EXACT:** the actual-cube countermodels and native support controls in Section 3. No countermodel here is a second tournament census. Each is explicitly certified not to be the tournament-isomorphism partition. The results do not rule out native all-order identities or an n-to-n+1 prediction theorem.

The exact count formula, reflection fixed-point characterization, odd class fibres, and black twins are inherited, not new claims: THM-643, THM-796, and THM-4731. The pinned historical copies are in `archive/` as byte-preserved `.txt` files. The earlier distance-matrix constraints also already distinguish cube partitions from native tournament fibres. This note contributes a sufficiency theorem for the relaxed structure, constructive hostile controls, and a clearly normalized arithmetic extraction step.

## 1. The exact native object and the information level

For n >= 3, let m=(n-1)(n-2)/2, f=floor((n-1)/2), N=2^m and B=2^((m+f)/2). A word in X={0,1}^m specifies the free arcs outside the fixed Hamiltonian path. C complements every free bit. R permutes coordinates by (i,j) -> (n-1-j,n-1-i). Thus C and R commute; C has no fixed word. R has f fixed coordinates, including (0,n-1), so CR has no fixed word either. Exactly B words are R-fixed.

Ordinary tournament isomorphism fibres have odd size H(T)/|Aut(T)|. R induces converse. A merged node is one self-converse ordinary fibre or a converse pair of ordinary fibres. Write V for the number of merged nodes, S for the number of self-converse nodes, A for the number of ordinary classes, d_u for the literal degree, and b_u for its blue degree. Loops count twice in both degrees.

Then:

- d_u is odd at an SC node and is 2 modulo 4 at an NS node.
- b_u is positive odd at an SC node and zero at an NS node.
- S = #{u : d_u odd} = #{u : b_u > 0}; A = 2V-S.
- The black multigraph is 2K for an integer multigraph K, loops included.
- At an NS node, deg_K(u)=d_u/2 is odd; at an SC node, deg_K(u)=(d_u-b_u)/2.

Consequently the full literal uncolored multigraph already determines A by degree parity. The full colored support graph also determines A, provided the full vertex set is retained: S is the number of vertices incident to blue support. Multiplicity is unnecessary for that second statement. An uncolored simple support graph does not retain these parity degrees. A node-signature histogram with its multiplicities also determines S; merely allowing signatures from an observed palette does not fix their multiplicities.

These are reconstruction statements at order n. If the order-n graph or the relevant histogram has already been completely enumerated, applying A=2V-S is not an independent prediction of A(n), and still less a prediction of A(n+1).

## 2. Exact characterization of the relaxed C,R layer

An **abstract odd-fibre partition** of the actual X is a partition into nonempty odd ordinary fibres that R permutes with orbits of size one or two. Merge each such orbit. No tournament-isomorphism condition is imposed.

**Realization theorem.** A finite colored multigraph G with vertex set U is the complement-line graph of such a partition if and only if:

1. Its blue degree b_u is either zero or positive odd.
2. Its black multiplicities are even: G_black=2K, including loops.
3. If b_u=0, then deg_K(u) is positive odd.
4. The blue and black totals are |E_blue|=B/2 and |E_K|=(N-B)/4.

Empty graph vertices are excluded by conditions 1 and 3. SC vertices are exactly the positive-blue vertices.

**Necessity.** The odd-fibre/reflection argument gives conditions 1 and 3. Since CR has no fixed words, R acts freely on black C-lines. The two lines in an R-orbit have the same merged endpoints, yielding condition 2. The cube orbit totals give condition 4.

**Sufficiency on the actual cube, not just an invented action.** Its Klein action has precisely B/2 orbits of size two, with R fixing both words, and (N-B)/4 free orbits of size four. Match the first kind bijectively to blue edge instances. For an edge u-v, assign x to u and Cx to v. Match the second kind bijectively to K edge instances. For an edge u-v assign x,Rx to u and Cx,CRx to v. These assignments exhaust the actual X and commute with its existing C,R action.

At a positive-blue vertex put all assigned words into one ordinary fibre. Its size is b_u+2deg_K(u), positive odd. At a zero-blue vertex split each R-pair between two sheets, using x,Cx in sheet zero when those words occur at that endpoint and Rx,CRx in sheet one. For a K-loop both endpoint occurrences are included. Each sheet has size deg_K(u), positive odd, and R exchanges the sheets. Thus the required partition exists. A blue loop contributes two R-fixed words; a K-loop contributes two words per ordinary sheet, or four words to an SC fibre. This also verifies the loop conventions.

This theorem is a complete description of the declared relaxation only. It does not characterize tournament isomorphism, Hamiltonian-path normalization, score profiles, strong components, or extension compatibility.

## 3. Explicit hostile controls with no false native claim

The saved `certificate.json` contains actual word-to-merged and word-to-ordinary maps for four partitions of the n=6 cube. `verify.py` checks every word, the genuine coordinate reflection, actual complementation, odd ordinary fibres, all edge multiplicities and the count formula. It does not import the producer.

### 3.1 Bounded, native-signature, native-total hostile

Both bounded models have:

- N=1024, B=64, V=34, and M2=sum_u d_u^2=48724.
- Exactly two singleton pure-blue vertices.
- Every individual (d_u,b_u) drawn from the actual native n=6 signature palette.
- Exactly two blue loops and 24 literal black loops, the native totals.
- Connected blue support on its nonzero vertices, connected black support on its nonzero vertices, and connected full support.

The first model uses the native histogram, hence S=12 and A_abstract=56. The second has S=20 and A_abstract=48, with this histogram:

| (d,b) | multiplicity |
|---|---:|
| (1,1) | 2 |
| (2,0) | 1 |
| (10,0) | 1 |
| (15,3) | 16 |
| (29,7) | 1 |
| (33,7) | 1 |
| (38,0) | 2 |
| (50,0) | 1 |
| (58,0) | 1 |
| (62,0) | 1 |
| (66,0) | 7 |

**Exact aggregate ambiguity set.** Under these declared n=6 constraints, the two values above are the entire possible set, not merely two examples. The native palette has blue degree one only at the two singleton vertices; every other positive blue degree lies from three to nine. Thus

    3S-4 <= B=64 <= 9S-16.

Since S is even, 10<=S<=22. The second-moment congruence with V=34 and M2=48724 gives S=4 modulo8. Therefore S is either12 or20, and A is either56 or48. Both possibilities are realized by the saved actual-cube partitions. The exact identification set in this relaxation is consequently {48,56}. Native tournament realization, a separate stronger condition, selects56.

The maximum merged degrees are 86 and 66. Their SC masses and NS half-masses are all actually observed native fibre sizes. The alternative uses an allowed subset of the native signature palette; it does **not** preserve the palette's full set or multiplicities. Full colored incidence and full literal degree histograms are also not claimed equal.

The native palette was independently re-enumerated by all fixed-path words and score-block permutation canonicalization, for n=3..6. Its degree histogram separately agrees with the byte-pinned historical fibre inverse census. The bounded hostile therefore does not depend on an obviously impossible huge fibre. Nevertheless both lifted partitions are non-native: the checker exhibits two words in one ordinary fibre with different tournament outdegree multisets. Those explicit failures are part of the certificate, not a caveat left unchecked.

The result establishes insufficiency of the stated totals, second moment, local allowed palette, loop totals and connectivity, **within the abstract C,R relaxation**. It does not establish insufficiency of any specified collection of native order-n graphs to predict their next native order.

### 3.2 Equal uncolored support hostile

A separate pair has literal degree multisets

    {2 repeated 8 times, 5, 1003},
    {1 repeated 3 times, 3 repeated 6 times, 1003}.

Both have V=10, total degree 1024, blue mass 64, and M2=1006066. Both have exactly the same uncolored simple support: a nine-leaf star and a loop at its center. They have S=2 and S=10, respectively, hence A_abstract=18 and 10. The giant center is deliberately an abstract construction and violates native mass bounds at n=6. It is used only for the exact support-collision claim, not the bounded claim above.

## 4. Every fixed finite degree-moment list loses information in the relaxation

Fix r>=1 and the odd interpolation nodes a_j=2j+1, j=0,...,r. Let L_j be their Lagrange basis polynomials and choose D to clear the denominators of L_j(2). Then the integer signed measure

    D delta_2 - sum_j D L_j(2) delta_(a_j)

annihilates d^q for every 0<=q<=r. Its positive and negative parts therefore have the same number of vertices and the same first r degree moments. Their numbers of odd vertices differ by D, because sum_j L_j(2)=1. D is even: the first-moment identity modulo two forces the sum of the odd coefficients to be even. For r=2 the relation is

    8 delta_2 + delta_5 - 3 delta_1 - 6 delta_3.

If the common leaf mass N0 is even, add one shared degree-one leaf, making N0 odd. Choose a sufficiently large native order n so B>2 max(S_leaf) and N>B+2N0. Add one common center of degree N-N0, odd. Each odd leaf receives one blue edge to the center; each leaf receives (d-[d odd])/2 K-edges to the center. Allocate

    (B-2S_leaf)/2 blue loops at the center,
    (N-B-2N0+2S_leaf)/4 K-loops at the center.

These are positive integers. The resulting blue/black degrees satisfy the realization theorem. The two graphs have the same star-plus-center-loop uncolored support, same actual cube totals, same V and same moments through r, but different S and A. Matching Klein orbits realizes both on the same actual fixed C,R cube. Thus no order-independent finite truncation of literal-degree moments identifies A throughout this abstract relaxation. The finite certificates check r=1,...,8; the interpolation proof covers all r.

This does not contradict finite-degree exact recovery when n or the support range is bounded: sufficiently many moments recover the bounded degree distribution, and Section 5 gives earlier arithmetic extraction in some intervals.

## 5. Merged-degree arithmetic extraction

Write M_q=sum_u d_u^q, so M_0=V. For odd d, d^2=1 modulo 8; for d=2 modulo 4, d^2=4 modulo 8. Therefore

    M2 = 4V-3S (mod 8),
    S = 4V-3M2 (mod 8),
    A = 3M2-2V (mod 8).

This explains why the second-moment hostile separates S by eight rather than two. For j>=2, odd-unit arithmetic gives d^(2^j)=1 modulo 2^(j+2) at an SC vertex; an NS degree is divisible by two and its power vanishes modulo that modulus. Hence

    M_(2^j) = S (mod 2^(j+2)),
    A = 2V-M_(2^j) (mod 2^(j+2)).

If an independently justified interval for S contains only one integer of the required residue, this determines S and therefore A exactly. For example 0<S<=min(V,B), with S even, is always available. This is an arithmetic extraction principle, not a method for producing the moments at the next order. The companion mass-prediction lane develops a stronger factorial-divisibility hierarchy for a different normalization: positive power sums of **ordinary odd fibre masses**, not these merged degrees. Those formulas must not be substituted without converting normalization.

## 6. The actual predictive gap

The missing native ingredient is a mechanism generating enough order-(n+1) incidence, signature multiplicities, or moments from independently computable smaller-order data. The odd-fibre/Klein action imposes real restrictions, including congruences, but the realization theorem shows exactly how much freedom that layer still permits. Extra native relations must cut down this relaxation. A full next-order graph already contains its own count; an abstract rewiring is not evidence that no native recurrence exists.

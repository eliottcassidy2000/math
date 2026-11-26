# A Valiant coefficient cut on the actual marked contact board

2026-11-26. **FINITE-EXACT** result on the 1,741-option marked board. This is a physical follow-up to [the state-sum and separator calculation](VALIANT-CUT-RANK-AND-SEPARATORS.md), and does not require the preprint's Euclidean pricing theorem or the existence of a full Moore graph. It concerns one specified branch cut, not arbitrary orders or arbitrary circuits.

## The board and its switch submatrix

Use the [anchored action board](FROZEN-SPECTRUM-AND-ANCHORED-ACTION.md) on source (1,1), fixing marker (2,7) and the three centers (56,0,28). Run the existing board generator and filter contacts whose center images use the marker destinations, as in anchored_wilson.cpp. The remaining 52 branches have 1,741 options. Each selected contact puts one unit in each of three center slots. With the three reserved destinations and the marker destination removed from each center, there are 52 usable slots per center. The square count in this board is

    C = sum_s binom(load[s],2).

A pair of contacts agreeing in multiple center slots contributes one distinct square per agreement. The marker and known base destinations cannot contribute squares to the remaining choices by construction. The square count is the marked-source score, not a final global girth certificate.

Introduce one independent variable per branch choice and form the polynomial with coefficient z^C on each complete selection. Flatten by the two sides of a branch cut. Fix the two vacancy contacts to (13,16),(33,19), and keep unswitched contacts at reference. The [switch list](marked_cut_switches.tsv) records 17 old/new choices on each side, in 34 distinct branches; put all other branches on either side with only the indicated choices retained. A single-switch difference vector has +1 at its new slots and −1 at its old slots, canceling any repeated coordinates. On the actual board, the 289 cross-side inner products form I₁₇. If u,v are the two bit vectors, this gives the exact mixed effect

    C(u,v) = 3 + a(u) + b(v) + u·v,
    a(u)=C(u,0)−3,   b(v)=C(0,v)−3.

The mixed finite difference of binom(load,2) is the product of the left and right signed increments in that slot. Consequently the cut submatrix, up to invertible row and column monomials and a common z³, is the 17-fold Kronecker product of [[1,1],[1,z]]. Its rank is **131,072** over any field of rational functions in z and at every nonzero specialization z≠1. A single bilinear message at this cut needs at least this many coordinates. At z=1 the matrix on this box is all ones and has rank one. Keeping an extra invertible indeterminate for each side's edit count does not affect the generic rank.

## Which lower bounds survive conditioning on the score?

The generic rank need not survive taking a single coefficient. Here there is a separate field-independent minor. Among u of weight 9, take those with a(u)=24 and b(complement u)=20. There are exactly **900**. Use those u as rows and their complements as columns. Every entry then has 17 reference edits, excluding the two fixed former vacancies. Its exponent is 47+u·complement(u'), which equals 47 exactly on the diagonal (equal weight forces u=u' if u⊆u'). Thus the coefficient at **47 squares and 17 edits has rank at least 900**. This includes no statement that those high-square configurations are good Moore candidates.

For the low-score contrast, the replay exhausts the same 34-switch box up to 12 squares. There are no selections up to cap 2 and exactly one at cap 3. At cap 12 there are 7,667 selections, and the restricted cumulative 0/1 matrix has rank 52 over F₂. The full table is in [the receipt](marked_cut_rank.out). Branching counts existing occupancies before placing a contact, so the current score cannot decrease and the cap pruning is exact for this box. Neither the generic minor nor this box exhaustion rules out a low-score choice outside the box.

## Exact contraction without storing every exterior contact

For an arbitrary cut retain only counts on slots attainable on both sides. If the boundary count signatures are L and R and the two internal scores are C_L,C_R, the exact total is C_L+C_R+L·R. For the same signature and same edit spend, a higher-cost partial choice is dominated in a minimum search. At square cap two a shared-slot entry can take three values (0,1,2), since three hits already cost three squares. This gives an exact cut-dependent meet-in-the-middle rule, not a guarantee that its space is small.

For the zero-square coefficient, discard any internally colliding partial selection. Each side's signature is now a subset in each center; the full 52 hits must occupy all 52 available slots there. A row signature therefore joins only its componentwise complement, subject also to the already checked private slots. Equivalently keep full slot masks (or reject any partial choice leaving uncovered private slots before projecting to shared slots). The zero-square matrix is a disjoint union of all-ones blocks, so over any field its rank is exactly the number of realized signatures on one side whose complements are realized on the other.

Run **python3 marked_cut_rank.py** in this directory. It reconstructs the marked board, verifies all mixed products and direct physical scores, and replays both the fixed-score minor and the low-cap search.

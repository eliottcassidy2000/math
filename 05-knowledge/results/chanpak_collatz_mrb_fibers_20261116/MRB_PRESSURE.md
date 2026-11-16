# MRB and tournament returns: control the total tail, not its largest term

**PROVED, scoped to these return languages.** The mixed dyadic/decimal controller with 111 distinct actions has O(X^(18/23)) distinct absolute outputs through X. The full bank of 168 eight-letter returns has O(X^(37/43)). This is an absolute-count obstruction despite complete residue steering; it does not constrain arbitrary span-three ribbons or arbitrary tournaments.

## Transfer idea

Session31's MRB audit distinguishes an exact residue from reconstruction of the alternating functional. Its entire perturbations preserve every branch jump, any fixed initial segment and all inverse-power positive-axis asymptotics, yet change the sum by an explicitly nonzero amount. Its positive result uses an Euler transform with a certified *whole-tail* estimate, not agreement to many digits or just the size of one contribution. See [the contour audit, §§2,6](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/session31_mirror_mrb_20261022/mrb/MRB_CONTOUR_AUDIT.md) and [the entire deformation note](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/session31_mirror_mrb_20261022/ENTIRE_MRB_DEFORMATION.md). The present analogue is an exact Dirichlet tail bound controlling an entire physical return bank, despite its local residue flexibility. Euler's alternating-series identity is not being applied to tournaments.

## Grouped pressure theorem

Fix a physical ribbon seed T and a finite return bank W acting on its common continuation space. A stop may also append any partial next block. For a word u of q complete blocks let B(u) be H of its isolated induced tournament with no seed, not the H increment at the live seam. For a positive rational s define

    P_q(s) = sum_{u in W^q} B(u)^(-s).

If P_q(s)<1 then the absolute-value counting function is O(X^s). In fact the sum of h^(-s) over distinct values is finite. To prove it, group the return schedule chronologically into q-block groups and a final stub. By the odd-cycle packing formula, independently selected packs from disjoint induced groups, and one pack in T, inject into packs of the full tournament. Cross-group arcs are unused in this injection. Thus H ≥ H(T)∏B(u). A given length of leftover suffix is one of finitely many choices independent of the number of complete groups. Sum the resulting majorant over the words of k groups; the product of their weights sums to P_q(s)^k. Summing this geometric series and then multiplying by the tail-menu size bounds the whole Dirichlet sum, including partial stops. Counting distinct integers rather than words only decreases it. The same proof works uniformly for arbitrarily many admissible physical seeds if their sum of H(T)^(-s) is finite.

For the mixed bank one can replace each complete block by an arbitrary representative with the same action on the common continuation space; this preserves the state exactly at the boundary and hence the final readout for any chosen partial stop. There need not be equality of the *isolated* B(u) between different representatives. Compute B from the specific physical representatives used in the pack proof. In this replay the last representative in the inherited ordering is used. The full bank uses all 168 words.

Ordinal sums multiply Hamiltonian path counts. After dropping h=1, the closure under these products still has finite Dirichlet mass at the same s: the mass over distinct values is bounded by the multiset product ∏_h (1−h^(-s))^(-1), finite since Σ_h h^(-s)<∞ and all h>1. No unique factorization in the bank is assumed.

A changing-core version is immediate when finitely many seam states with labeled q-block transitions encode a proposed language. Replace P_q(s) by the nonnegative matrix with (i,j) entry equal to the sum B(u)^(-s) over transitions i→j. Over all length-k walks the weight is a nonnegative endpoint-vector pairing against the matrix power. Spectral radius below one certifies an obstruction: it implies finite Dirichlet mass, ruling out coverage at any s<1. Radius at least one removes this obstruction only; it is not a coverage criterion, since word collisions and sparse output intervals remain possible. An infinite changing-core proposal could instead seek positive weights v_i with row bound Σ_j P_ij(s)v_j ≤ κv_i, κ<1, provided the terminal costs admit a uniform upper bound times v_i.

## Exact replay

The literal transfer matrices from the eight-braid packet give each isolated three-block count; the inherited matrix products check that their nonnegative int64 arithmetic cannot overflow. Bucket each H down to L(H)=2^15 floor(H/2^15)>0. With D=10^12 and s=p/q, use the least integer m satisfying m^q L(H)^p ≥ D^q as a per-term upper numerator, then count multiplicities of L(H). All exponent comparison and sum arithmetic here is integer and directed upward.

| Bank | Triple terms | s | Directed upper bound |
|---|---:|---:|---:|
| Full physical bank | 4,741,632 | 37/43 | 985140783540/10^12 |
| Invertible-at-five actions | 1,367,631 | 18/23 | 986083961573/10^12 |

The prior minimum-only bounds used exponents about 0.93644 and 0.86070, respectively. Aggregate growth captures the fact that the overwhelmingly prevalent profiles expand faster than the worst triple. Every large dyadic interval [X,2X] therefore contains runs of Ω(X^(5/23)) consecutive odd sites absent from the mixed bank. Over an even modulus M for which the mixed bank hits all odd classes, the maximum of the least positive representatives it realizes is Ω(M^(23/18)); a uniform selector bounded by cM is impossible. By ordinary pigeonhole this abundance of holes survives any finite repertoire of fixed compatible seeds and all their ordinal products. It does not rule out families with genuinely more mass, different seams, or an exact signed-cancellation construction outside this positive-growth architecture.

Replay in the working tree: `python certify_pressure.py` (uses the sibling eight/braid and decimal-rank packets and their generated return automaton; exact pressure rows in `certify_pressure.out`).

A sharper asymptotic follows without changing the certificate: finite Dirichlet mass implies o(X^s), by splitting the sum at a fixed cutoff and sending that cutoff to infinity after X. Thus the maximal missing-odd run divided by X^(5/23) actually tends to infinity along dyadic bands, and the worst least representative divided by M^(23/18) tends to infinity along allowed moduli.

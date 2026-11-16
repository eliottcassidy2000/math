# The seam with a growing spacer is an exact divisor problem

**Proved, all widths.** In a single transitive tape with one head, insert an unmarked transitive spacer of length t≥1 at a fixed cut. The two now-separated marks may change on both sides, and their high contributions need not be frozen individually. The complete readout reduces to a product of two bounded factors and a separate bit-reversal remainder.

## Setup and factorization

Fix left and right transitive blocks of widths k,r≥1. Index each from zero in forward tape order. Mark any positions, but leave the very first and last tape vertices unmarked. Keep the entire tape TT_k → TT_t → TT_r, the sentinel t₀ and head x from the earlier [head construction](https://github.com/eliottcassidy2000/math/blob/6d38e8ceb0ccc8a061270e52507820bd4ecf53e6/05-knowledge/results/chanpak_computational_tapes_20261115/HEAD.md). The completed tournament is strong. Let b encode the left mark **backward** and a encode the right mark **forward**. Thus

    0 ≤ b < 2^(k−1),        0 ≤ a < 2^(r−1).

Write rev_w for reversal of w bits, with leading zeroes retained, and J(z) for the within-block pair energy

    J(z) = Σ_{i>j, bits i and j of z both set} (2^(i−j−1)+[i=j+1]).

The exact count is

    H = P·2^t + C,
    P = (2^k−b)(2^r+a),
    C = 1 + rev_k(b) − rev_r(a) − J(b) − J(a).             (★)

Proof: let A_R be the right block's forward marked weight before inserting the spacer, and B_L the left block's backward marked weight. Then A_R=2^k a and B_L=2^r b. In the whole-tape one-head formula, the signed linear weights are rev_k(b)−rev_r(a) + 2^t(A_R−B_L). The tape baseline is 1+2^(k+r+t). Within-block J is unaffected by insertion. Across the cut the forward contact is impossible since the spacer is nonempty, so the pair interaction is exactly

    Σ_{u left,v right marked} 2^(t+v−u−1)
    = 2^t · 2^(−k−r) · B_L A_R.

Here u,v are pre-insertion tape positions, and reversing the left encoding preserves its within-block J. Gather growing coefficients: 2^(k+r)+2^k a−2^r b−ab factors as displayed. All expressions count physical paths; the equality is not an approximation or a congruence.

## A genuine three-state bank at the same growing seam

Take k=r=10. Mark the masks encoded below, without changing the cut, outer guards, internal orientations or sentinel connections:

| forward right a | backward left b | left factor | right factor | C |
|---:|---:|---:|---:|---:|
|153|334|690|1177|−391|
|494|489|535|1518|−389|
|241|382|642|1265|−387|

Every product is **812130**, although both factors change on every transition. For **every t≥1**, these give three consecutive odds

    812130·2^t −391,  812130·2^t −389,  812130·2^t −387.

At t=1 the counts are **1,623,869 → 1,623,871 → 1,623,873**. Only head contacts change, and the two transitions change different pairs of factor coordinates; this is not a repeated constant drift. The same reversal convention at every height defines one fixed pair of outside mark-address spaces. Unlike simply increasing exterior guards, the unbounded growth is inserted directly into the seam between the two active mark sets.

A small control is k=5,r=4: (a,b)=(5,13),(3,11) gives P=399 and C=3,5. Its counts are 801,803 at t=1 and fit independent full subset-DP on the literal tournaments. The larger fixture is checked separately by directly summing the pre-head and post-head path options by prescribed marked endpoints; the proof of (★) gives the all-height assertion.

## What a wider bank now requires

At a fixed cut, two configurations have the same growing coefficient iff the corresponding bounded factors have the same exact product. The list of available residues at that product is not decided by the divisor count:

- enumerate divisors d within (2^(k−1), 2^k]; retain e=P/d in [2^r, 3·2^(r−1));
- for each, set b=2^k−d and a=e−2^r and calculate the bit-reversal remainder in (★);
- only exact adjacency by 2 between those remainder values makes a repair path.

This is a complete decoder within this family. In particular a prime leading coefficient is impossible when both k,r≥2. The width is at most the number of admitted factor pairs; many divisors alone cannot certify a window. Also J(z)≤z−[z odd]: each occupied upper bit i can contribute at most 2^i to the pair sum, and i=0 contributes nothing. Thus |C|<2^(k+1)+2^(r+1), so at large spacers a prescribed target anchors P to at most one integer as soon as 2^t exceeds twice this bound. The remaining exact obstruction is then product factorization *and* digital remainder, rather than only a local residue.

The discovery program exhausts all masks and cuts through the chosen size bound; compile with `c++ -O3 -std=c++17 cut_product.cpp -o cut_product` and call `./cut_product 22`. Through that bound the maximum is three, without any global optimality claim. This was motivated by the older [arithmetic seam ledger](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/arithmetic_seams_20260921_dynamics.md), which keeps a second sheet beyond the projected return. Here the two sheets are explicit as P and C: a divisor collision on the first sheet alone is not an actual repair.

## Same factorization for arbitrary nontransitive cores

Put fixed outer guards TT_a, TT_b (a,b≥1) around D₁ → TT_t → D₂ with t≥1. Use profiles (h_i,q_i,E_i,L_i,J_i) for the two actual marked cores and put X=2^a,Y=2^b. Then the exact polarized factorization is

    H = (q₁X−L₁)(q₂Y+E₂)·2^t
         + h₁h₂ + h₂(XE₁−J₁) − h₁(YL₂+J₂).             (★★)

Indeed the nonempty spacer kills the forward contact, while the reverse contact contributes −2^t E₂L₁. Combine it with the three other growing terms q₁q₂XY·2^t, q₁XE₂·2^t and −q₂YL₁·2^t. The remaining two one-sided profiles combine exactly as stated. The proof only uses uniform ordinal transitions and endpoint-weighted paths.

This gives a new exact separator obstruction: on a fixed-product fibre, all changes in H are multiples of gcd(h₁,h₂), because the remainder is a sum of two terms with those multipliers. Hence no unit repair with a common odd factor can survive on the same high product after the spacer is made nonempty. With a fixed finite bank of masks and fixed cores/outer guards, any persistent bounded difference at increasing spacer length eventually has identical high product since its bounded remainder cannot absorb a nonzero multiple of 2^t. Both assertions fail to describe the contact term at t=0: the deleted term −T₁S₂ need not be a multiple of the shared factor. This is exactly how the earlier adjacent h₁=h₂=13 move obtained +2 despite its common multiplier. The route to larger growing-spacer windows on two fixed nontransitive cores therefore needs coprime h's, or genuine product carries with an expanding bank, not merely more masks on the 13/13 pair.

One more restriction if the *outer* guards grow too: for fixed two cores, suppose two fixed mask states keep identical spacer coefficient at infinitely many synchronized outer heights X=x·2^s,Y=y·2^s. The degree-one and constant terms in 2^s must both agree. Thus q₁x ΔE₂=q₂y ΔL₁ and Δ(L₁E₂)=0. Nonzero changes in the first equation have the same sign, which strictly changes the product of the two nonnegative coordinates. Consequently E₂ and L₁ must both be unchanged. A genuine factor-pair swap of this type needs at least one bounded/changing outer address; synchronously growing every guard is a different problem.

### Dyadic scaling cannot by itself grow the bank

For each odd divisor d₀ of P there is at most one admissible left factor d=d₀·2^j: the interval (2^(k−1),2^k] meets a dyadic orbit in at most one point. In the transitive two-block case each admitted factor pair determines both masks, so the full width is at most the number of divisors of the **odd part** of P. The power of two in the growing product cannot increase this budget.

More precisely, start from k,r,P and replace them by k+u,r+v,P·2^(u+v). For a fixed odd divisor d₀, the unique normalized factor (allowing rational dyadic values temporarily) changes from d to d·2^u and its paired factor from e to e·2^v. The interval tests are therefore unchanged; each candidate becomes integral once both increments are sufficiently large. At that point the associated b and a also scale by 2^u and 2^v, respectively. Appending zero low bits while increasing a bit width keeps the reversal value fixed, and translating every set bit by the same amount keeps J fixed. Thus the whole exact remainder table stabilizes as both widths grow. Any unbounded bank in this architecture must change the odd product (or the architecture), not just apply larger dyadic dilations to a fixed product.

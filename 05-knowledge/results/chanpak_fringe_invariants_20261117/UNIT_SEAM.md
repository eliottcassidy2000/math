# An exact two-count seam with a growing, target-invertible source

Here is a physical absolute repair, independent of the modular fringe compiler. A one-apex spine word has a marked central bit, written P b S, with nonempty left and right words P,S. Let F be the usual binary value of P, Z the value of the reverse-complement of S, and e,z their least significant bits. From the odd-pack inversion formula,

    H(P1S)−H(P0S)=Z−F+z−e.

To prove it, concatenate two words: the crossing inversion sum is the left binary value times the right reverse-complement value, with one extra for an adjacent 1→0 at the concatenation. For P0, the new last zero adds F+e, and the new left binary value for the final concatenation is 2F. For P1 there are no newly finished inversions and the binary value is 2F+1; the new adjacent contribution against S is z. Subtracting gives the identity.

The invariant is the **difference of two physical addresses**, not a residue. For any F≥1 and d≥1, set Z=F+2d, use the ordinary binary expansion of F as P and the reverse-complement of the ordinary binary expansion of Z as S. These are literal words; the two apex tournaments differ only on the arc at the middle bit. Both are strong because their first bit is 1 and their last bit 0. Their count difference is exactly 2d. In particular the d=1 repair survives arbitrary growth of the address.

For d=1 the lower read is

    H₀(F)=2F(F+2)+F+H(P)+H(S)−1+e.

For every 1-bit in P, its total possible later-zero contributions are at most its binary weight, so H(P)≤F+1. The analogous bound summing by zero bits in S gives H(S)≤Z+1. Hence both reads lie in

    [2F²+5F+1, 2F²+7F+6].

The intervals are disjoint for consecutive F≥1. This yields Ω(√X) distinct adjacent-odd windows up to X and an exact inverse for this source family: from the prescribed target, locate the only candidate interval by integer binary search, then evaluate the explicit inversion sum on its two candidate words. No candidate interval or a mismatch in the literal read rules out this family at that target. This is not interval coverage; the separated envelopes and their one address each make clear why the repair alone does not fill the spectrum.

Replay: `python replay_unit_seam.py` checks the exact windows, including literal subset Hamilton counts on small direct controls.

## Carry-free address: preserve all longer odd cycles

There is extra structure at the unit step. Restrict F to have binary bit 1 equal to zero (positions start at zero), so Z=F+2 changes just that one address bit. In the central arc flip, losses of inversions ending at the central bit are indexed by the binary bits of F, with distance d corresponding to place d−1. Gains of inversions starting at the central bit are indexed by Z in exactly the same distance coordinates. Thus all the inversion counts by distance agree except for **one extra distance-two inversion** on the upper side.

In a one-apex word, a distance-d inversion gives β_{d,k}=binom(d−1,k) odd cycles of length 3+k for even k. Hence the +2 switch adds exactly one directed triangle to the cycle-length census and leaves every longer odd-cycle count unchanged. All directed cycles go through the apex, so no two can coexist in a spanning cycle pack. Writing G(x,y) for the spanning odd-cycle and singleton packing polynomial at a common spine length L, the **whole** difference is the single monomial

    G_upper(x,y)−G_lower(x,y) = x y^(L−2).

In particular the entire ordered cover read changes by ΔP_q=2q^(L−1); the Hamiltonian difference is 2 and ΔP₂=2^L. Thus the unit repair preserves all L low bits of P₂ at growing size, not merely a fixed residue. It cannot preserve P₂ exactly. The same no-carry argument with Z=F+2^j (j≥1 and F's bit j zero) changes just the distance-(j+1) inversion census, leaving all odd cycles of length greater than j+2 unchanged. The restricted unit family still has a positive fraction of the addresses and still gives Ω(√X) separated windows.

## A different exact switch preserves the full two-cover count

Use the same central flip with addresses F=20+64u and Z=32+64u for arbitrary u≥0. The two physical strong tournaments have ΔH=12. The low six address bits differ by removing places 2,4 and adding place 5 (zero-based), and all higher bits are identical. For slot variable q, the inversion-distance formula assigns address place s the weight

    q^(L−s−1)[(q+1)^s+(q−1)^s].

At q=2, the normalized weights at places 2,4,5 are 5/2,41/8,61/8, and the last is the sum of the first two. Thus **ΔP₂=0 exactly**, at every u. At q=3 the same difference gives ΔP₃=−100·3^(L−5), so the invariant has a concrete scope: a third-slot reader can see the replacement. The nonlinear wing reads cancel in the response; they are free to grow in the baseline. The exact H windows have quadratic baseline growth and occur at unbounded orders, not just after neutral padding.

A reader using only P₂ therefore cannot distinguish these two actual growing sources. Any affine reader of (P₂,H) sees only the H change with its own coefficient. For example the [physical sevenfold reader](../chanpak_changed_contacts_20261117/MERSENNE.md) 7(64P₂−H) transfers this to an exact absolute 84-count gap. This remains a gap construction rather than a universal choice of baseline.

## Sharpness and a limit on serial neutral contacts

The 12 step is the smallest positive H change of a single apex-arc switch preserving P₂ in this word model. At fixed order, P₂≡(−1)^order H mod 3 and H is odd, so the only smaller possible positive step is 6. Orient the comparison from central 0 to central 1 (or exchange the two wings to reverse the sign). Then H changes by Z−F+z−e. For step 6 the integer address difference Z−F is 5, 6, or 7. To get ΔP₂=0, the two addresses must have equal digit weight sum with weights λ_s=(3^s+1)/2^s.

If adding that small difference propagates a carry into positions ≥3, the altered high bits are a block of ones from position 3 through t−1 changed to zeros and a zero at position t changed to one. Its signed weight is

    λ_t−Σ_(s=3)^(t−1) λ_s = 13/2−(3^t−3)/2^t.

The low three positions can offset at most 13/2, so t≥7 is impossible. With no high carry or t≤6 all differing bits lie below position 7. Enumerating their 128 possible low addresses gives no exact equality with H gain 6; the replay verifies this finite last step. Thus a one-arc +6 with exact P₂ is impossible here.

There is also at most one usable P₂-neutral apex contact in a fixed strong spine of length at least six. The normalized ΔP₂ for flipping position i from 0 to 1 is the sum of λ_s over following zeros minus the sum over preceding ones, where s is distance minus one. As i moves right the latter sum increases and the former decreases. For adjacent internal indices the decrease is strict: if the first bit 1 is already at distance at least two, its weight increases strictly, and otherwise the last bit 0 is far enough to strictly decrease its weight. Here λ₀=λ₁=2 and all subsequent weights strictly increase. Flipping the one neutral index leaves its own normalized ΔP₂ unchanged; the same monotonicity prohibits a new neutral index in that word. So serial single-contact neutral changes can only flip back. More than one usable digit in the same exact P₂ fiber requires changing the wing or using a different kind of move.

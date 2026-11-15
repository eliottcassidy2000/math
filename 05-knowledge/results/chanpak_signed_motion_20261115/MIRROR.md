# Mirror the live block, not each large coordinate

The two-copy repair with a common divisor has a general norm form. It is a signed cancellation at the level of path geometry, not a combination of independently amplified +2 shifts.

Let a physical core D have a contact mask V and profile (h,q,E,L,J,T,S) as in [the two-block formula](signed_move.md#two-active-nontransitive-blocks-cancel-the-amplification-or-bypass-it). Choose the adjacent right block to be the converse D* with the reflected copy of the very same mask (same labels before any isomorphism). Reversal interchanges E with L and T with S and preserves h,q,J. The forward cross-contact term is therefore T²; the reverse cross-contact term is L². Substitution in the two-block identity gives, for independent positive transitive guard lengths a,b,

    H = q²·2^(a+b) + (h E + q L)(2^a−2^b) + h² − (2h J + L² + T²).

Hence equal guards cancel the entire large contact response, even if E and L both change. The exact degree-two energy for the selected mask is `2hJ+L²+T²`; separate coordinate equality is unnecessary. With asymmetric guards it suffices to preserve the single combination hE+qL.

This can break a common divisibility barrier: when a prime p divides both h and q the residue is exactly `−L²−T² mod p`. In particular when p≡3 mod 4 this residue is zero exactly when L and T both vanish modulo p. That is the surviving two-contact information, not the divisibility of either component's unmarked path count.

## An exact shared-factor kernel

For lexicographic forward code 955 on five vertices the reversal map i↦4−i is an anti-automorphism. Thus the same physical internal arc set can be used in both halves. Moving the left mask 21→7 and the right mask 21→28 implements the mirror rule. For the two masks V=21 and V=7 the profiles give h=13, q=104, hE+qL=6786 both; (L,T,J)=(58,8,26) or (59,7,22). The energies are 4104 and 4102. Algebraically, if L increases by one and T decreases by one, the exact energy difference is

    Δenergy = 2 (h ΔJ + L−T+1).

Here -52+50+1=-1. This identifies the arithmetic relation a new mask bank would need to preserve or traverse. Each isolated edit contributes a multiple of 13 if it changes only J; the joint move gives +2 anyway.

An unmarked transitive block of length c≥1 between the live pair changes the formula: the forward contact is impossible, the reverse term becomes 2^c L², and the leading coefficient is scaled by 2^c. At equal guards the exact difference for this **same** mask move becomes 104−117·2^c, a multiple of 13. So this unit does not come from “presenting two copies” alone: the head must be able to traverse the seam in the forward order. No small shared-factor correction can be inferred after inserting even one intervening vertex. For the no-separator pair, all 32 mirror masks have longest step-two run two; this is not a repeatable decoder.

Finally, appending more copies of the same common-factor bank to a one-control ordinal tape cannot by itself repeat this construction. If r blocks each have p dividing h and q, then each term of the exact two-slot expansion still contains one divisible unmarked factor from all but at most two of them, so p^(r−2) divides every readout for r≥2. In particular all such outputs are multiples of 13 once a third copy of code 955 is added, regardless of its own mask. To extend the repair one must share the same two-contact geometry among many states, change the bank, or add a second interacting control. Copying a successful repair as an opaque block is not a valid repetition.

## Necessary arithmetic for a repeatable drift

Suppose two successive +2 repairs with the same fixed h use the same signed drift (u,v) in (L,T). Each energy drop is −2, so starting at (L,T),

    2(Lu+Tv) + u²+v² ≡ −2 (mod h).

At the next step the left side gains 2(u²+v²). Since h is odd, repetition requires `u²+v² ≡ 0 mod h`. Thus the successful (1,−1) move cannot simply be repeated on the h=13 mirror bank, even if another suitable mask were available. For h=13 a nonzero constant drift must have u≡±5v (mod 13), besides satisfying the first affine condition. More generally, if any prime p≡3 mod 4 divides h, a repeated constant drift is impossible: isotropy forces u≡v≡0 mod p, making the first condition read 0≡−2. This does not rule out changing drift, changing h, or adding a control. It says exactly what a constant-drifts-only extension would need to avoid, and suggests testing isotropic signed directions before searching more masks.

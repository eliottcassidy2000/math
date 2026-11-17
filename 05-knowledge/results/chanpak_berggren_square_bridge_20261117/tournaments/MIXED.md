# Break the diagonal norm: a literal taxicab sevenfold family

The [diagonal construction](PROOF.md#an-exact-physical-sum-of-two-squares) uses a weak migration followed by its transpose and reads a sum of squares. To break its even 7-adic valuation restriction, choose **two different** words. For any nonempty binary w,z pad them as a=0w1 and b=0z1. From the rooted transitive pair v=(1,2), apply a and then the reverse-complement of b, followed by the same strong singleton reader. The identity M(b*)=M(b)^T now gives the **mixed** count

    v^T M(b)^T M(a) v = h_w h_z+R_w R_z.

All moves are literal, and the proof of strength for the final reader is unchanged. There is no subtraction of tournaments.

Fix w=10, so (h_w,R_w)=(15,38). Choose z=1^j00, j≥1, and let A=2^j. The padded word 0 1^j 00 1 has d=16A, u=8A−7, rev=10A−2 and F(u)=3A−1, from the same digit recurrence. Hence

    (h_z,R_z)=(19A−13,38A+2).

Reading the dot product gives a strong tournament on j+14 vertices with exactly

    1729·2^j−119 = 7(247·2^j−17)

Hamiltonian paths. The taxicab number 1729=7·13·19 is **the actual coefficient** of the two-port read, not an assumption about cube representations. Moreover 247A−17 ≡ 2A+4 modulo 7, which never vanishes for A in the orbit {1,2,4}. Thus the mixed family has **exactly one factor of seven at every height**. For j=2 the quotient is the verified prime 971, giving H=6797 on 16 vertices, also checked by ordinary literal subset DP. This does not imply infinitely many prime quotients or prove that they attain first order at these sizes.

In this construction the distinction between a norm and a mixed pairing is operational: keeping the two words identical loses the entire v₇=1 corridor, while a fixed left word and a systematically different right word recover an infinite physical lane. The useful remaining challenge is to choose the two physical words from a prescribed prime quotient, not to show that such a factor pattern is possible at one scale.

Run `python mixed.py` for the digit formulas, partial inverse arithmetic, and literal 6797 control.

# A fixed rejected triple has sharp cubic height, not just an average obstruction

The [joint residue theorem](JOINT_WORD.md) gives simultaneous congruences
in the strong one-apex transitive-spine family but, following the
local/height distinction in [Sun §§6–7](https://arxiv.org/html/2611.14605v1#S6),
does not imply an exact solution. For the fixed triple (63,65,65)
the least ordinary count of a physical match at modulus M=2^k is
**Θ(M³)**. Here H is the ordinary count and the other two coordinates
are the same apex's positive response in the tournament and its
converse. Both the lower and upper bounds below hold for every k≥15.

## Freeze two endpoint words, then read the middle

More generally fix any odd response r≥3 in both readings. Put
x=r−1, l=bitlen(x), ρ=rev_l(x), and L=2^(k−l). At width n≥k,
the response congruences R≡R_c≡r mod 2^k fix the last k bits B of
the contact word to x padded by zeros, and the first k bits A to
the reverse-complement of B. Their addresses are

    P_A=Q_B=2^k−1−ρ L,      g(A)=g(B)=g(x).

For sufficiently large k these end pieces cannot overlap. With
k−l>l, an overlap of at most k−l positions puts ones against
zeros; a longer overlap would require the suffix of A to have
at least k−l zeros, although A has at most l. Write the full
word A Z B, with t bits in the free middle Z, and let Y be the
bitwise complement of Z at its declared width.

Split the weighted inversions from the [binary reader](JOINT_WORD.md)
into the three internal terms, the three cross products, and the
adjacency correction. For t>0 the ordinary count is

    H = 2g(x)+g(Z)−2 + P_A Q_Z + P_Z Q_B
           + 2^t P_A Q_B + 1−Z_first+Z_last.

The [complement identity](../chanpak_absolute_reflection_20261114/COMPLEMENTARY_READER.md)
is g(Z)−P_Z+rev_t(P_Z)=g(Y)+Z_first−Z_last. With
Q_Z=2^t−1−rev_t(P_Z), the middle readout simplifies to

    H ≡ 2g(x)+g(Y) + ρL (2^(t+1)−P_Z−Q_Z) + ρ²L²2^t   (mod 2^k).

In particular H≡2g(x)+g(Y) mod L, also for t=0:
the direct cross term is P_A Q_B and the adjacent correction is
one, giving 2g(x)+1 modulo L as required.

If a fixed desired count N differs from g(x), widths n<k are
excluded for sufficiently large k because they would force
P=x exactly and hence H=g(x). Now put D=N−2g(x).
If D is negative or is not a count of any one-apex binary word,
then the congruence forces a middle count of order L or larger.
Since g(Y)≤2^t, the free middle needs at least k−l bits for all
sufficiently large k. This is a third logarithmic block in the
size of one common witness, in addition to the two prescribed ends.
It applies in particular when D=7 or 21 (the scalar gaps), without
asserting anything about whether N itself is attained.

## At (63,65,65), both bounds are explicit

Here x=64, l=7, ρ=1, g(x)=65. For k≥15 the required end words are
A=1111110 followed by k−7 ones and B=1000000 preceded by k−7 zeros;
they cannot overlap. Smaller width than k forces H=65. The middle
must satisfy g(Y)≡−67 mod L, hence g(Y)≥L−67>L/2. Therefore
t≥k−7, n≥3k−7, and the mandatory first-to-last inversion gives

    H ≥ 2^(3k−9)+1.

For a matching upper bound, use fixed short words
U=1101101 (address 109), V=101011 (address 43), and take

    Y = U 0^(k−7) V,      Z = complement(Y).

For arbitrary fixed nonempty U,V and s≥1, the same inversion
split gives
g(U 0^s V)=P_U(1+Q_V)·2^s+g(U)−P_U+(P_U mod 2)+g(V)−1.
In this case g(U)=33, g(V)=9, Q_V=10, and therefore
g(Y)=1199 L−67. Also P_Z+Q_Z≡rev(U)−P_V−1=91−43−1=47
mod 128, since the long zero segment isolates the two ends.
As 1199≡47 mod 128, the unreduced middle identity gives
H≡130+1199 L−67−47 L≡63 mod 2^k. The fixed ends supply
both rooted congruences. This is a literal strong word of
width n=3k+6 with H≤2^(3k+6), so it achieves the same cubic
order as the lower bound. It does not give H=63 at these
rooted responses.

The two-short-word calculation also gives a target-first search
rule for other fixed equal-response triples: seek U,V with
constant term D=N−2g(x) and with
P_U(1+Q_V)≡ρ(rev(U)−P_V−1) mod 2^l. If found, the same middle
padding produces literal simultaneous matches of cubic height.
This is a sufficient construction for the local lift, not for the
original scalar target. The one-apex chart's height loss has now been
measured on a prescribed example, rather than just inferred by
counting all residue classes.

The fixed-63 lower obstruction originated in a background
pass. The equal-response reduction, explicit saturated-order
upper construction and replays are derived here. Run
**node replay_fixed_triple_height.js** for 35,770 general middle
checks and 82 exact modular witnesses at the fixed triple.

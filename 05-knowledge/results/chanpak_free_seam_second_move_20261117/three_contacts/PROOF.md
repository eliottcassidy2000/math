# A near-balanced three-contact repair at the one-third scale

This extends the [free reflected seam](../../chanpak_wheel_square_selector_20261117/FREE_SEAM.md)
and the [one-arc obstruction](../one_arc/PROOF.md). A proposed three-contact
construction from the background research was independently checked, then
extended here to every complementary suffix, with an endpoint count for the
entire enlarged bank.

Use the same literal seven-vertex deletion core and apex. For opposite
parity positive addresses, the lower count is

    H(u,w) = 81 + 142(u+w) + 256uw + 5(F(u)+F(w)),

with the odd address placed first. Reflecting the flanks gives H+2.
Choose a dyadic G>=4 and any odd r in [1,G-1], and let s=G-1-r.
For positive integers x,y with x even, take

    u=G(16x+5)+r, w=G(16y+4)+s,
    u'=u+4G,     w'=w-4G.

The transitions in the four bits above the suffixes are 0101→1001
and 0100→0000, so these are exactly three literal apex arc reversals
on each of the two reflected hosts. The sum is unchanged, equals
-1 modulo G, and both parities and leading bits remain intact.
All four hosts are strong by the original seam construction.

## Exact displacement and inverse

Let C_r and C_s be the reverse complements of r and s in log2(G)
bits. Their leading-zero terms in the outer F formula are unaffected,
because both prefix changes are even. The inner four-bit formula
gives F(16x+9)-F(16x+5)=x+2-[x odd] and
F(16y)-F(16y+4)=2y-4. Writing z=x-y, the common displacement on
the lower and upper tracks is exactly

    Delta = 15y - A z - B_r,
    A = 16384G^2-5,
    B_r = 5120G^2 - 1024G(s-r) + 10 - 20(C_r-C_s).

Indeed the product difference is 4G(w-u)-16G^2 and the outer F
linear contributions are 4(C_r-C_s). At r=1 this specializes
to B_r=4096G^2+3072G+30.

A is a unit modulo 15: G^2 is 1 or 4 mod15, so A is -1 or -4.
Given any even Delta with |Delta|<=G^2, take the unique z in
1,...,15 making Delta+A z+B_r divisible by 15, then set y to
that quotient and x=y+z. Since A is odd and B_r and Delta even,
x is automatically even. These are positive, uniformly
Theta(G^2), and z is bounded, independently of the chosen suffix.

Consequently the lower source count is Theta(G^6) and |u-w|=O(G).
At that scale the repaired shift can be any prescribed even integer
of order up to G^2, i.e. the scale of the one-third error, without
sacrificing the exact +2 or the complementary suffix. This is a
shift selector, not a selector for a prescribed absolute count.
For G=4, r=1, Delta=4, the lower source is 22047013898432061
and the lower destination 22047013898432065, with both +2 partners.

Within the fixed grammar of one apex and these binary words, two
parity- and sum-preserving edits can only transfer a bit between
the two words in the same column: two signed powers of 2 sum to
zero only if their magnitudes coincide and their signs differ.
The [paired-digit estimate](../paired/PROOF.md) applies to the
whole words (pad with leading zeros if needed), giving minimum
nonzero displacement 512. Thus the three-contact scheme is
genuinely beyond same-column trades.

## Even the enlarged bank does not prescribe the absolute center

Fix X>=2 and 0<=E<=sqrt(X). Count endpoints at most X from this
template with |Delta|<=E, allowing **every** dyadic G and **every**
complementary suffix. There are O((E+1)(sqrt(X)+E)) such values.

To prove this, note each address on either side is at least a
constant times G times its corresponding positive prefix x or y.
The count is at least 256 times their product, so G^2 xy=O(X),
and G=O(sqrt(X)). Uniformly B_r=O(G^2) and B_r>0. When z>=0
we have x>=y, hence y=O(sqrt(X)/G). The displacement equation then
forces z=O(sqrt(X)/G^3+E/G^2+1). When z<0, positivity of y and
the same equation give |z|=O(E/G^2+1). For fixed (G,r,z) the
number of y satisfying the displacement bound is O(E+1). There
are O(G) choices of r. Sum over dyadic G: the three terms are
O(sqrt(X) sum G^-2), O(E sum G^-1), and O(sum_(G<=const sqrt X)
G), respectively. Multiplying by four also counts the reflected
partners. In particular E=O(X^(1/3)) still gives density zero.

This is a constructive gain over a one-bit repair: the *same
geometric regime* now supports fine moves. It isolates the real
remaining condition: to cover every target, the small displacement
must be chosen together with the starting count, not in place of it.

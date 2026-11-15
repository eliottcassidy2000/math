# Exact cancellation with a hidden cut graph

The [contiguous-control obstruction](../chanpak_graceful_ternary_20261115/RECTANGLE.md)
does not prohibit interleaving. Write a one-apex tournament on an
ordered transitive spine by bits x_i; 1 means the apex beats site i.
Its count is

    H=1+2 Σ_(i<j) c(j−i)x_i(1−x_j),
    c(1)=c(2)=1, c(d)=2^(d−2) for d≥2.

This is the [prior odd-cycle formula](../chanpak_graceful_ternary_20261115/GRACEFUL.md).
For disjoint signed changes δ,η the mixed count is exactly
minus twice the pairing with kernel c(|i−j|). The prior
[Boolean continuation dividends](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-2259-boolean-continuation-hasse-field-and-signed-interaction-dividends.md)
give the general language: independence is the vanishing of a
mixed difference, not a visual separation of the supports.

## A literal four-count unit rectangle

Using the 11-bit center integer 1158, flip either of the disjoint
xor masks 32 and 220, or both. The exact ledger is:

| Choices | Word integer | Count |
|---|---:|---:|
| Neither | 1158 | 721 |
| 32 | 1190 | 723 |
| 220 | 1114 | 725 |
| Both | 1146 | 727 |

The first control flips bit 5 (from least significant numbering).
The other changes bits 2,3,4,6,7 with signed changes −1,+1,+1,+1,−1.
Its interaction with bit 5 is −2+1+1+1−1=0 in the kernel.
All four tournaments are strong because the word starts at 1 and
ends at 0, closing a directed spanning cycle through the apex.
The replay checks actual Hamilton paths by subset DP. The additive
minimum length is 11 within this model (the finite smaller-span
check is separate and is not needed for construction).

## The same switch, an exact boundary ledger

Let a be the ordinary binary left address and b the reverse
complement right address. For a fixed-width center with ordinary
address u and reverse complement v, the two-cut seam formula is

    H(LWR)=H(a)+H(u)+H(b)−2+av+bu+2^m ab
              +[a odd][v odd]+[u odd][b odd].

This follows by applying the [two-block seam law](../chanpak_graceful_ternary_20261115/GRACEFUL.md)
twice. Both original controls fix the end bits of the center, so
the change from the exterior is exactly a Δv+b Δu. Their steps are

    A: 2+32(b−a),        B: 4+56a−44b.

The mixed term stays zero for every exterior. The original +2,+4
steps survive simultaneously only at a=b=0. For a=b=t they give
two +2 pairs, with starts separated by 4+12t. The lower start is
2048t²+2428t+721+2p(t), where p(t)=H(t)−1 on a binary word.
This is an infinite family of strong exact +2 pairs, not arbitrary
placement of their absolute baseline.

## A mirror-neutral fiber is an exact weighted cut graph

Two center words have the same correction relative to their
interior count in every mirror exterior a=b=t iff their ordinary
minus reversed binary addresses agree. Equal signatures imply
that their signed bit change is palindromic. To see this, let
P(z)=Σ_i(δ_i−δ_(m−1−i)) z^i. It has roots 2 and 1/2 and
coefficients in [−2,2]. If its nonzero degree is d, multiplying
P(1/2)=0 by 2^d shows the leading coefficient is even, hence has
magnitude two. The leading term at 2 is then too large for all
lower terms combined (2^(d+1) vs at most 2^(d+1)−2).
So P vanishes. This also verifies that the bounded parity guards
of the seam agree. See the [older two-sided radix proof](../chanpak_collatz_norm_clock_20261115/WALSH.md).

Each mirror pair is therefore either fixed to 01 or 10, or free
to choose 00 or 11. Any middle bit is free. Put s_i=2x_i−1:
the count summand is c(d)(1+s_i−s_j−s_i s_j)/2.
The linear term is antisymmetric under reflection and the quadratic
form commutes with reflection. The free symmetric choices have no
interaction with the fixed antisymmetric coordinates; other fixed
symmetric coordinates on mixed pairs are zero. Thus the excess
over the all-aligned state is exactly the weight of a cut of the
complete graph on free pairs and the possible middle bit. Pair
representatives i,j give edge weight 2[c(|i−j|)+c(|i−mirror(j)|)].
A pair and the middle give weight 2c(distance). Every edge
weighs at least two, so with f free vertices the cheapest positive
rise is at least 2(f−1). Interior unit differences are not ruled out.

## Why the singleton cannot be prolonged on any positive wall ray

Suppose two independent switches are neutral along a positive
rational ray of exteriors a:b, and one flips just bit i.
The coefficient of its external response forces b/a=2^(m−1−2i)
(with i numbered from the least significant bit). Hence the
second switch must satisfy equality of its two signed dyadic
moments about i. After zero padding around i, the two-root
argument shows its signed changes symmetric about i. If its mixed
difference with the singleton vanishes, its radial signed sum
against weights 1,1,2,4,… vanishes. At the largest active
distance D this requires D≥2, and every smaller radius active
with the sign opposite to that at D, because the sum of all
smaller weights equals the last one. At these radii the initial
bit pairs are equal within each pair and opposite between the
inner group and outermost pair. Thus their contribution
Σ_(d≤D)c(d)(1−x_(i−d)−x_(i+d)) to the singleton half-step
cancels. At greater distances every weight is divisible by
2^(D−1), also on one-sided padding. Its actual step is divisible
by 2^D and cannot be ±2.

This pinpoints the extra requirement beyond integer orthogonality:
to transport a unit step with the baseline one needs two genuinely
multisite switches, a new reader or a target-dependent interface.

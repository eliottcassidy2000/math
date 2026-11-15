# Mirrored endpoints: a sharp obstruction and an eight-site exact triple

**Proved:** complete classification of unit moves for the simultaneous
endpoint toggle, and a different three-choice move producing three
consecutive odd counts in every transitive mirror frame. **Finite-exact:**
all binary-center buckets through the orders in shelf.out. No exact
selection of an arbitrary target is claimed.

Use the [inherited mirror count](https://github.com/eliottcassidy2000/math/blob/9afc28f177c5cd7dc7872c089565298fb76382d2/05-knowledge/results/chanpak_tangent_survey_20261115/MIRROR_NORMS.md).
An apex beats the 1-sites of a transitive center C. Write n_C=1+2I(C), where

    I(C) = sum over i<j of 2^max(0,j-i-2) C_i (1-C_j).

This follows by taking the first and last spine vertices of each odd
cycle through the apex, then choosing an even subset of the intermediate
vertices. The positive odd-cycle expansion counts both the empty
packing and those single cycles. A left painted exterior L and its
reverse-complement counterpart R are held fixed throughout. In this
mirror, the readout depends on the center via its length m, its integer
address skew w−rev(w), and endpoint difference C_last−C_first, together
with n_C. Indeed s_C+t_C=2^m−1+w−rev(w), and e_C+f_C=1+C_last−C_first,
in the inherited notation. Any changes in n_C within one such bucket
are multiplied by g², where g=H(L), regardless of the other exterior
coordinates.

## 1. Lengthening both-end toggles never gives a new unit

Flip both endpoints from 0 to 1, leaving the L=m−2 interior bits held.
Let b be the interior address (most significant bit first), rev the
reversal in exactly L bits, and e the sum of its endpoint bits. Write
D=I(new)−I(old). For L≥2,

    2D = 2^L − b − rev_L(b) − e.

Only pairs incident to the flipped endpoints contribute; the pair
between the two flipped sites contributes zero both times. The
one-sided weight sum is 1+1+2+⋯+2^(L−2)=2^(L−1).
For the marked interior, doubling its weight reproduces its binary
address with one extra copy of its nearest digit. Adding the reversed
side gives the formula.

If e=0, then b and rev are each at most 2^(L−1)−2, so D≥2.
Equality holds when all internal sites of the interior are 1.
If e=2, both addresses are at least 2^(L−1)+1, so D≤−2.
If e=1, removing the two opposite bits leaves a word c of length
d=L−2 and gives D=2^d−1−c−rev_d(c). For d=0 this is zero, and
for d=1 it is +1 or −1. For d≥2, D=±1 would require an even sum
c+rev(c), so c must have equal endpoint bits. With both endpoints
0, c+rev(c)≤2^d−4; with both 1, c+rev(c)≥2^d+2. Neither meets
2^d−2 or 2^d. Length three gives D=1−2b directly; length two gives
zero. Thus all **upward** unit states are

    000    → 101
    01000  → 11001
    00010  → 10011

for the 00→11 toggle. The corresponding downward unit states are
010, 01100 and 00110. At every other m≥4 the minimum strictly
positive difference is +4 in the unscaled apex count (choose e=0
as above). The statement is about this particular endpoint toggle,
not every switch, and is independent of the frame.

## 2. A different move makes three consecutive odds

The three centers

    10100101   10011001   10000001

are palindromic with endpoints 1, hence have the same skew 0 and
endpoint difference 0. Their I values are 30, 31 and 32 respectively.
They induce the three exact outputs

    2g n_L + (59+2j)g² + 255gS + gd_L + 256S²,  j=0,1,2.

In particular the transitive mirror with a marked singleton on the
left gives the literal strong order-11 counts **573, 575, 577**.
Every transitive mirror with nonzero address X gives three consecutive
odd counts, with leading term 256 X². For centers of *any* even
length at least 8, absorb extra complementary singleton pairs into
the mirror; the three centers still share their readout bucket and
keep the +2 gaps since the extra ordinary counts are 1.

This is an exact family of triples, not a coverage theorem. For
transitive left exterior the other terms are 255X+61+(X mod 2)+2p(X)
in the j=0 branch, with 0≤p(X)≤X as in the inherited formula.
The triples are disjoint and give only Θ(sqrt(T)) values up to T.
The missing operation must move their common baseline, rather than
just making the same pair closer.

### Sensitivity to a mismatched exterior

The same transfer applies with independent left and right painted
exteriors, whose ordinary counts are g_L,g_R, start and end weights
S_L,T_R respectively. Set A=g_Lg_R, D=g_R S_L−g_L T_R. The three
addresses are 165,153,129, with s_C=w and t_C=255−w. Both end bits
are held. The two varying crossing terms are g_R S_L t_C for the
left/center seam and g_L T_R s_C for the center/right seam; crossings
from left to right use the fixed optional two-slot count 2^8 of the
center. The two signed differences in the displayed order are
2A+12D and 2A+24D. In particular the three outputs can be consecutive
odds in some order only if D=0 and A=1: both signed differences
are congruent modulo 12; in a permutation of three consecutive odds
they can share that residue only when they are both +2 or both −2.
This forces D=0 and, by positivity, A=1. Thus the triple is exact
for matched trivial-H envelopes but does not also let their two
responses drift independently.

## Replays

    python verify.py
    c++ -O3 -std=c++17 shelf.cpp -o shelf
    ./shelf 23

The second command groups every word by its length, skew and
endpoint difference, then computes the longest unit run of I. The
finite run lengths are not asserted globally, except for the
explicit existence and endpoint classification above.

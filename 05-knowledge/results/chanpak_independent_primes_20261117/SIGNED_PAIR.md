# Let two prescribed targets choose the normalization defect

Use [the two-variable cyclic scaffold](../chanpak_prime_dichotomy_20261117/RECEIVER_GATE.md)
F=C₃[TT₃→B,C,a] with one apex a and marks at the first and third vertices of TT₃. Define g=H(F), W as the weighted two-cover start count, and D=U+E−g as the previous normalization defect. The two displayed bilinear forms in that result give T=3W+2g and

    D=H(B)(8H(C)−P₂(C)).

Do **not** impose D=0. Adding the temporary root y gives h=W−D and w=g+W, so the ordinary menu count is p=h+2w=T−D. Taking the labels {y,s} in the return construction gives exactly

    H_return=7p+4D=7T−3D.

The same source therefore realizes p and 7q precisely when it has those two physical reads and

    D=7(q−p)/4.

The coefficient 4 is not a rounding tolerance; it is the absolute-target defect. The constructions remain strong after discarding the neutral ordinal sink for the ordinary count, or after closing the return at y and s for the sevenfold count. No normalization is needed for this strength argument.

## Complete upward prime branch of this scaffold

Suppose p,q>3 are primes and q>p. Then r=8H(C)−P₂(C)>0. Because P₂(C)≥(|C|+1)H(C), this implies |C|≤6. The exact Hamilton-frame replay gives the following exhaustive positive receivers and ordinary count forms (same B):

| (H(C),P₂(C),P₃(C)) | r | coefficients of (H(B),P₂(B),P₃(B)) in p |
|---|---:|---|
| (1,2,3) | 6 | (−7,46,78) |
| (1,4,9) | 4 | (−15,48,390) |
| (3,12,33) | 12 | (−45,−36,1638) |
| (5,32,117) | 8 | (−123,−828,6630) |
| (15,108,477) | 12 | (−417,−5172,28782) |

One may fix the frame since every tournament has a Hamilton path; all its other arcs vary independently. The last four rows make p divisible by 3. Thus **the receiver has to be a singleton** on this branch and the target difference determines the source:

    m=(q−p)/24 must be a positive odd integer,
    H(B)=7m,
    39P₃(B)+23P₂(B)=(p+49m)/2.

These equations are sufficient as soon as an actual B exists. They also give a finite inverse from the target pair: enumerate positive P₃ up to (p+49m)/78 in its unique residue modulo 23 and set P₂ by the linear equation. The possible order is bounded by 3^|B|≤P₃, from the positive all-singleton term in the odd-pack formula. For prime p>3, the source must have odd order and 3∤m: modulo 3, P₂(B)≡(−1)^|B|H(B), so p≡P₂(B)−H(B).

The pair (p,q)=(731761,732073) gives m=13. An actual compatible B has outgoing masks (48,89,35,85,4,26,53) and profile (91,1256,8649). Thus an actual ordinary prime and the distinct count 7·732073=5124511 land together. Both quotients are prime by exact trial division. The smaller required source count is H(B)=91=7·13; generally this is a descent constraint on the difference, not a proof that every independently chosen pair has that form or a physical source.

## An exact small-gap obstruction in both directions

If |q−p|=16, the required |D|=28. Since H(B) is odd and is never 7, it must be 1. Thus the receiver must have signed error ±28. Positive 28 is absent in the catalog above. For negative 28, P₂(C)=8H(C)+28. If |C|≥8, the split bound gives (|C|−7)H(C)≤28 and the all-singleton term gives 2^|C|≤P₂(C). This contradicts 2^|C|≥256>8·28+28. A Hamilton-frame scan through order seven finds no receiver with negative error −28. So even allowing the normalization to change cannot pair primes of difference 16 by this particular sevenfold branch at any order.

Replay: compile [positive_receivers.cpp](positive_receivers.cpp) in C++20; it prints the positive profiles and would also print any profile at error −28. [replay_pairs.py](replay_pairs.py) checks the positive fixture by physical subset DP.

# Ternary branch witnesses and a closure obstruction

Status: finite certificates and an all-order necessary bound, not a proof of the spectrum or of a closed state bank. Counts are numbers of Hamiltonian paths. This is a test of Session 151's ternary compiler, not evidence that its offset gadgets exist generally.

For order-five codes, use lexicographic pairs i<j with bit 1 for i→j. An ear mask m adds x→i iff bit i is set; each new vertex receives the next index. Starting from code 759, mask 2 gives a strong parent of H=31. Its masks 50, 24, 28 produce H=93, 95, 97, i.e. all 3H, 3H+2, 3H+4. Independently, code 8 followed by mask 10 gives H=29, whose ear masks 19 and 33 produce H=89 and 91. Each is strong and all path counts agree between subset DP and literal permutation enumeration (`verify.py`).

Nevertheless, enumerating every nonconstant mask of the H=93 and H=95 children finds none of their respective 3H, 3H+2, 3H+4 as a *single-ear* child. For the H=97 child, masks 2 and 78 give H=293 and 295. This rules out closure for this particular one-ear state selection, not the possibility of different representatives, multiple ears, or the Session 151 bank. The 3H branch always exists via the ordinal product T→C3, but its state membership is a separate question. `chains.py` prints complete target-mask lists for these examples and further comparison states; `verify.py` independently confirms the specified witnesses.

## A strengthened complementary-cut bound

Let S be a nonempty proper vertex set of a strong n-vertex tournament T, H its path count, and F(S) its ear response as in the cut-induction note. If T has a directed Hamiltonian cycle with r edges crossing the unordered cut (S,V−S), then

    F(S)+F(V−S) ≥ 2H + 2 ceil((H−n+(n−1)r)/2).

Indeed each of the n rotations of the cycle is a distinct old Hamiltonian path and together these rotations cross the cut (n−1)r times. Each other old Hamiltonian path crosses it at least once. Thus C(S)≥H−n+(n−1)r. The marked-break identity from `cut_realizability` says F(S)+F(V−S)=2H+K(δS), where K(δS)=C(S)+B_cross(S) is even and B_cross≥0. The stated bound follows. Since r≥2, it improves or matches the previous strong-parent bound, and it is strictly stronger when a cycle with more crossings is known. It is a necessary condition and does not prove that complementary target responses are realizable.

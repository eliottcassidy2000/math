# Why the certified interval starts at 69: the seam misses 67

The [finite range certificate](FERMAT_RANGE.md) starts at 69 for a real structural reason. **No rooted tournament of any order** can reach H=67 by one one-way overlap with a near-transitive R_j, even when j need not be a power of two. This is not a statement about all tournaments realizing 67.

Write the source as the ordinal join P→C→D, where C is the strong component containing the shared root v (possibly a singleton), and P,D can be empty. A post-v path in the relaxed L count cannot use a vertex of P, which v loses to, since any such vertex would precede the entire later part. All of P must be in the pre-piece. The two pieces on D are independently arbitrary except for their joint coverage, because v and all of C dominate D. Therefore

    H_A=H(P)H(C)H(D),   L_A=H(P)L_C Q_D(2).
    H_out=H(P)[H(C)H(D)+2^j L_C Q_D(2)].

The bracket is at least three, so if H_out=67 (prime), H(P)=1; the prefix can be discarded. Since L_C≥H(C)≥1, the budget forces Q_D(2)≤33 and H(C)≤22. By the exact small-clock profile catalogue the only nontransitive D possibilities at that Q budget have (H,Q)=(3,12),(3,24),(5,32). The last imposes an output at least 69. In the first two, if C is nontrivial then H(C)≥3 makes the output at least 81. If C is a singleton, H(C)=L_C=1 and the output is 3+12·2^s, never 67.

It remains that D is transitive, of any size compatible with the bound. Its Q is a power of two and H(D)=1. Absorb it into the exponent: 67 would have to equal H(C)+2^k L_C for some k≥1. The [Busch floor](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/) excludes strong C of order seven or larger at count ≤22. The literal complete order-at-most-six enumeration has no matching rooted profile; k≥7 is already too large. The finite replay checks every allowed (C,root,k), together with the small D catalogue. This proves the claimed failure for 67 even before imposing the nine-vertex source bound or the true Fermat exponents.

For a strong source and the three true Fermat exponents 1,2,4, the entire lower window below69 misses exactly 1,7,11,13,21,47,67. Thus the new finite contiguous interval has the sharp immediate lower endpoint for its seam grammar, while the older selected-ear packets and the full spectrum are different families.

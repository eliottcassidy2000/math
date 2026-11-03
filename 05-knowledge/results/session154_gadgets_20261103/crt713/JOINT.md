# Exact q2 joint image modulo 7·13²

Use the actual TRANSPORT seed and pure-q2 prefix: z_k=(25,34)A^k, A=[[11,18],[4,7]], h_k=Y_k+6. This is an exact statement about one source family, not about all strong tournaments or mixed prefixes.

**Theorem.** Write u=h_k mod7, v=h_k mod169. The exact simultaneous image consists of the following 934 pairs (u,v):

* If u≠6, then v mod13 is one of 0,1,...,11. When v≡6 mod13, retain only v∈{32,45,136,149}; all 13 lifts of each of the other eleven allowed residues occur.
* If u=6, then v mod13 belongs to {4,9,11,12}, with all 13 lifts allowed.

In particular the image is strictly smaller than the product of the individual images, which has 7·160=1120 pairs. The first six rows have 147 pairs each and the last has 52. Each pair occurs for infinitely many nonnegative k, yielding infinitely many distinct h_k.

**Finite exact certificate.** The integer checks prove ord(A mod7)=48 and ord(A mod169)=728, whence 4368 is a common period. For each k in 0,...,4367, compute z_k modulo 1183=7·169 and retain ((h_k mod7),(h_k mod169)). The included code verifies membership equals the displayed condition for every one of the 1183 possible pairs, verifies exactly 934 attained pairs, and produces the fibre-size histogram

    size:    2   3   4   5   6  12  13
    count: 234  78 156 156 247  39  24.

This is an all-index proof, not a sample extrapolation: A^4368=I modulo both 7 and 169, so every k has an equivalent index in that finite range. The bound on 13² uses 728=56·13; its order is independently checked. The generalized CRT shows why compatibility can fail: gcd(48,728)=8, not 1, and a mod7 phase and mod169 phase can coexist iff their indices agree mod8. The code's complete cycle tests all compatible phases directly. Positive growth of h_k supplies infinitely many distinct counts in each attained residue class.

The two actual companion families have exact images obtained by applying v↦3v+2 or v↦3v+4 modulo 1183 to this set, since 3 is a unit modulo 1183. Equivalently invert 3 in each local ring before applying the displayed criterion. They differ by 2 for the same k, so their residues are not independent. This result does not assert that the same forbidden pattern persists at higher 7- or 13-powers: their index-period gcd can grow with the 7-power.

Source: TRANSPORT.md at f6d5304849e50b29df6e368239e92fffa92e57e3; the joint enumeration and phase compatibility are deductions here. Run `python check.py` for the certificate.

# Online binary and decimal steering on the same rank-four image

Fix the strong span-3 seed ribbon `(0,0,3,5,2,4)`, with 25 Hamilton paths. The prior [eight-mask and braid packet](https://github.com/eliottcassidy2000/math/blob/3a2423f0bb11b8a63dc81b7b3616a881c101b3f6/05-knowledge/results/chanpak_eight_braid_20261116/THREE_LEASHES.md) supplies 168 eight-letter return blocks, with 136 distinct integer matrices in a common four-coordinate lattice. A block is named by its first three letters; its final five are always `(0,3,5,2,4)`. Write (B_w) for that matrix. Every block maps the exact four-dimensional image onto itself over the rationals, so every sequence of these blocks preserves rank four at every letter (ranks cannot increase between endpoints). The same lattice reads out (H=fv), where (f=(32,14,10,27)) and the seed vector is ((-43,-5,4,53)^t).

**Theorem.** For every fixed (k,dge1), the very same seed admits online steering to any sequence of odd residues modulo (2^k5^d), with finitely many literal extension letters per request and exact image return at every eight-letter boundary. In particular this works for every fixed decimal precision (10^d). The number of letters in this all-power proof is not small. For modulus 800, a separate short certificate gives **exactly 40 letters per request after a one-time 24-letter initialization**. The initial ribbon and every eight-letter return boundary are strong, although individual letters between those boundaries need not be. To check strongness, contract the old strong ribbon to one vertex, which has an outgoing arc to each new vertex because some old site is more than three places back. The new vertex at one-based block position \(j\) can reach the old vertex directly precisely when its mask has a set bit for some offset \(\ge j\). The replay checks every listed block in this nine-vertex graph; all new vertices reach the old one. Thus the property persists at all return boundaries. These are congruence results, not absolute count selection.

## All precisions: the cheap coordinates split

Take (A=B_{(4,3,5)}). The exact matrix checks give
```
A^620 = I mod 25,
A^4   = l r mod 32,
l=(28,24,24,19)^t,  r=(24,10,14,27),
f l = r l = 1 mod 32.
```
In particular (A^4mod32) is idempotent and its reduction mod 2 has its unique nonzero entry in position (4,4).

We may work at (kge5); for a smaller requested precision use 32 internally and reduce the result. Choose an exponent (E) divisible by 1240 with (E/4ge lceil k/5
ceil), and let (S=A^E). The exterior-square matrix of (A^4) is zero modulo 32. Consequently every 2-minor of (S=(A^4)^{E/4}) vanishes modulo (2^k). Its fourth column has unit readout because the reduction mod 2 of (S) remains the projector onto the fourth coordinate. Dividing that column by its readout gives (L), with (fL=1) modulo (2^k). All columns of (S) are multiples of (L): use the vanishing minors with the unit pivot in position (4,4). Thus every post-reset state has dyadic part (tL), (t=Hmod 2^k), and (Lequiv lmod32). Also (Sequiv Imod25). This collapses hidden dyadic coordinates without reducing rational rank. The other coordinates are retained modulo (5^d); they stay primitive (not all divisible by 5) whenever we use only (B)'s invertible modulo 5. The seed is primitive and its dyadic state after initialization by (S) is odd.

Use macros (X_w=SB_w) (append (w)'s return block, then the reset). Each maps the dyadic line to itself by a scalar (alpha_w=fX_wL). The scalar is odd for each of the 111 distinct (B_w)'s invertible mod 5, as checked mod 2. On the (5^d) side every such macro is invertible. It remains to program two channels without interfering with one another.

### First program the factor (5^d)

For every nonzero vector (umod5), the outputs (fB_w u), over those 111 blocks, run through all five residues. The finite check needs only the 156 projective representatives with first nonzero coordinate 1; unit rescaling permutes the outputs. Since (S=Imod25), one macro can therefore set the first 5-adic digit.

The following four (C_i=X_{w_i}^{o_i}) are identity modulo 5. The displayed rows are (f(C_i-I)/5mod5); because (Sequiv Imod25), their value is independent of the chosen (E).
```
w_i         o_i     row
(2,1,6)     24      (4,2,4,2)
(2,1,7)     24      (3,1,1,0)
(3,1,7)     20      (0,0,4,4)
(2,3,3)     24      (0,0,4,2)
```
The rows form a basis over (mathbb F_5). To lift from a selected residue modulo (5^{j-1}) to one modulo (5^j), (jge2), choose a row nonzero on the current primitive vector mod 5. If (C_i=I+5D_i), then
```
C_i^(2^k 5^(j−2)) = I + 2^k 5^(j−1) D_i    mod 5^j.
```
The matrix binomial expansion proves the identity by induction. Its first five powers, exponents 0 through 4, choose every lift of the readout. The factor (2^k) kills any odd dyadic multiplier, so this action also preserves the whole dyadic line state. Primitivity at 5 is unchanged.

### Then program the dyadic factor

For (w=(3,3,5)) and (w=(4,1,6)), respectively, the macro scalar is 5 and 7 modulo 8, and (B_wmod5) has orders 15 and 31. Therefore raising these macros to (15cdot5^{d-1}) and (31cdot5^{d-1}), respectively, makes the matrix identity mod (5^d). This follows by repeated fifth powers of a matrix that is identity mod 5. These exponents are odd, so their dyadic scalars are still 5 and 7 modulo 8. Such units generate the entire odd unit group modulo (2^k): the first generates all units that are 1 modulo 4, since (v_2((1+4u)^{2^m}-1)=m+2) for (u) odd; the second supplies the other coset. Use positive powers (or zero) to multiply the present odd dyadic readout into the target. This leaves the programmed (5^d) residue intact. CRT finishes the request. The same primitive post-reset conditions hold for the next request. This proves online steering, not merely a one-time menu.

## A shorter mod-800 controller

Use the different reset (P=B_{(4,5,5)} B_{(4,3,7)} B_{(2,1,7)}), in chronological order (217,437,455). It is invertible modulo 25, and mod 32 it equals (l_0 r_0) with
```
l_0=(7,2,21,9)^t, r_0=(16,30,24,27), f l_0=1.
```
Initialize by (P). From the normalized state (t l_0mod32) and primitive (umod25), choose **two** blocks (B_a,B_b) invertible at 5, then (P). Each choice has a signature
```
g = f P B_b B_a l_0 mod 32,
lambda = f P B_b B_a mod 25,
```
and the readouts on the two channels are (g t) and (lambda u). There are 12,308 distinct signatures from the 111-by-111 choices. For each of the 16 odd (g), its rows (lambda) evaluate onto all 25 residues at every primitive (u). The replay normalizes the first unit coordinate of (u) to 1, leaving 19,500 vectors, and checks the full 25-bit coverage mask for every (g,u). Unit rescaling completes the proof for all primitive vectors. Each request is five eight-letter blocks, i.e. 40 letters.

### Replay and limitations

Run `python check_decimal_controller.py`; it reads `return_automaton.json` from the sibling eight-mask packet (or an explicitly supplied path). The higher-power proof uses just the small matrix identities and finite mod-5 check; the stronger fixed latency uses a larger NumPy census. Integer matrix arithmetic is within int64 for these finite checks; the general proofs use binomial divisibility, not a search at growing moduli. The powers may increase rapidly with precision, and image return is not vector return. Nothing here places (H) at an absolute target or proves density one.

## Why this control language alone cannot settle Chan–Pak

The [earlier three-block growth certificate](https://github.com/eliottcassidy2000/math/blob/3a2423f0bb11b8a63dc81b7b3616a881c101b3f6/05-knowledge/results/chanpak_eight_braid_20261116/RETURN_GROWTH.md) gives a minimum of (gamma=13,460,883) Hamilton paths on any isolated three-block suffix in this bank. Thus after (m) whole blocks, (Hge25gamma^{lfloor m/3floor}) by packing on disjoint three-block supports. Every extension used above belongs to the 111 different matrices invertible at 5. For (Hle X) there are at most
[
 O(111^{3log_{gamma} X}) = O(X^{3log(111)/log(gamma)}) = O(X^{0.860697})
]
different outputs in this full repertoire, including the unrestricted choices and therefore any selected policy. Allowing an incomplete final block adds at most a constant factor, since after exact image return its observable state determines each continuation count, and count is monotone while adding vertices. The larger full 136-matrix repertoire also has an exponent less than 0.898. Thus even a repertoire that can follow arbitrary increasing decimal precisions remains sparse in the absolute integers. Modular universality and an absolute spectrum covering theorem are different requirements.

A precision need not be frozen. To increase or change it before a request, apply the reset for the newly requested (k) to the current physical ribbon. Its mod-5 vector stays primitive, and the new dyadic line condition is restored. With (E=1240lceil k/1550ceil), the all-power algorithm needs (O((k+1)2^k5^d)) letters per request (replace (k) by 5 internally when smaller); this is a constructive upper bound, not a claim of optimal latency.

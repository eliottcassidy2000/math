# Arbitrary independent pair residues away from 2 and 3

For every modulus m coprime to 6 and every independently prescribed (p0,q0) mod m, there is a **strong one-apex** B such that the unit-determinant mixed lift reads p≡p0 and q≡q0 (its upper graph count is 7q). The claim is only a local compiler. In particular it does not place either prescribed absolute prime.

## A finite physical word action

On a binary spine word of length n, let b_i=1 mean the apex beats spine i, indexed from zero. Write F for its ordinary binary value, D=2^n, and normalized registers
    G = (Σ_i b_i 2^i 3^(n−1−i))/D,
    J = (Σ_i b_i 2^i)/D.
Let H be its full Hamiltonian path count and P=P₂ its two-cover count, and use barred registers
    X=H−F,   Y=P/D−2G+2J.
Every denominator is a unit mod m. For the empty word H=1,P=2 (one apex).

Use two available physical suffix commands at a word ending in 0:
    A = append 0,
    B = append 10.
For the bit0 addition, the one-apex inversion formula gives H+=F+last_bit and P_new=2P+2D(G+J). Appending 1 instead fixes H and doubles P. These formulas follow by separating all old inversion contributions and adding those ending at the new zero. Normalizing and using that both A and B finish in 0 gives
    A: (F,G,J) -> (2F,3G/2,J/2); D->2D; X,Y fixed,
    B: (F,G,J) -> (4F+2,(9G+3)/4,(J+1)/4);
       D->4D; X->X−F, Y->Y−G−J.
In both cases the modular actions are invertible: gcd(m,6)=1. A negative power is physically available by repeating the corresponding forward permutation until taking an inverse step on its cycle. Thus any finite group expression in A,B has a positive suffix-word implementation. Initialize the actual spine word 10, so any continuation by these symbols has first bit 1 and last bit 0 and its tournament remains strong.

## Independently program the two read coordinates

Use standard composition right-to-left. Put L=diag(2,3/2,1/2) and W=B A^(−2). It fixes D, translates (F,G,J) by v=(2,3/4,1/4), and has vertical shears
    X -> X−F/4,
    Y -> Y−4G/9−4J.
For W_j=A^j W A^(−j), its translation is L^j v and its shears are the old rows multiplied by L^(−j). The operation V_j=(W_j W)^(−1) (W W_j) fixes F,G,J,D. Its change in X, or in Y, is l L^j v−l L^(−j) v, for the corresponding W shear row l. At j=1,2 this gives
    V₁=(-3/4,11/9), V₂=(-15/8,85/27)
in the (X,Y) plane. Their determinant is −5/72. Hence they independently set X,Y modulo every m coprime to 30, without shifting the registers determining H−X, P/D−Y or D.

The same independence holds at every 5-power with one extra word. Modulo 5 let Θ=Y+2X and Q(F,G,J)=3 F J+G². Then A preserves Θ−Q because F and J scale inversely and (3/2)²=1 mod5. Under W the Θ shear is 2F+4G+J, whereas Q's shift under v=(2,2,4) mod5 is 2F+4G+J+3. Thus W increases Θ−Q by 2, as do all W_j. The concatenation with powers
    W₀^(−6), W₁^19, W₂^(−16), W₃^4
fixes F,G,J,D exactly at every permitted modulus. Indeed its translation polynomial 4t³−16t²+19t−6 annihilates L. The sum of its exponents is 1, so it increases Θ by 2 mod5. Since V₁ has Θ increment zero and X increment nonzero mod5, the new word together with V₁ is invertible on the vertical plane modulo each 5^a. Their changes may be evaluated at the fixed starting trans registers and then are independent of X,Y. CRT combines the two generating pairs at coprime factors, so we can set arbitrary H,P modulo any m coprime to 6.

Finally the associated exact mixed lift is p=44P−13H, q=71P−21H with determinant −1. Set H≡71p0−44q0, P≡21p0−13q0, and use the compiler above. Actual magnitudes can be very large: inverse implementation by permutation powers is not optimized and constrains no absolute count. At modulus 3 there is a genuine physical restriction P₂≡(−1)^|B|H, so both residues cannot be varied arbitrarily; this result makes no claim at moduli divisible by 2 or 3.

## Preserve the prime-compatible modulo-three phase and even source order

We can require the actual tournament B to have even order and prescribe either of the two nonzero phases p≡−q (mod 3), while still placing arbitrary independent residues modulo every m coprime to 6. In the physical implementation enlarge the finite state by the parity of spine length and by (F,X) mod3, even if m does not include 3. On this component A flips the parity and sends F→2F, X→X; B preserves parity and sends F→F+2, X→X−F. Both are permutations. Each formal vertical generator has net length zero and fixes F mod3; for the extra 5-word this follows from the same translation polynomial at 2. Its X shift is independent of X, so cubing every vertical generator preserves both mod3 registers. Cubing is a unit on the m coordinates, leaving their span unchanged. Start from spine 100 or 10100, of odd length and with count 5 or 19. Those give the two H mod3 choices, and even tournament order gives P₂≡H mod3. The two prime-compatible output orientations follow from the unimodular reads. This does not prescribe residues at powers of 3 or 2.

Replay `python replay_joint_phase.py` compiles literal positive suffix words at moduli 5,7,11,25,35, independently reads each word with the unnormalized inversion recurrence, and checks that even order and the selected phase at 3 survive. Some words in this explicit implementation already have hundreds of thousands of bits. This emphasizes the remaining gap between phase and absolute magnitude.

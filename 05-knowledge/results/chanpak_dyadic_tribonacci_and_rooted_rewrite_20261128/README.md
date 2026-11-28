# The moving-root grammar at the prime 2

2026-11-28. This extends the [all-odd-modulus result](../chanpak_joint_crt_and_dyadic_boundary_20261128) and the [moving-root construction](../chanpak_moving_root_and_module_interfaces_20261128). The parity-compatible chart union is now completely reachable modulo each power of 2, and also modulo every modulus with exactly one factor of 2. A separate rooted isomorphism removes one entire positive-integer predecessor type. These results concern register states of this grammar, not exact scalar surjectivity for Chan–Pak.

The chronological letters, starting at (1,1,1,1), are
\[
P:(h,e,u,g)\mapsto(h,0,2u,g),\quad
K:(h,e,u,g)\mapsto(h-e+u,g,2u-g,g),\quad
M:(h,e,u,g)\mapsto(h-e+u,h-e,h+u,h).
\]

## 1. Exact dyadic reachability

For \(N=2^a\), \(a\ge1\), we claim exactly the following union is reachable, always with h and g odd:
\[
\mathcal P:\ e=0,\ u\text{ even};\qquad
\mathcal K:\ e=g,\ u\text{ odd};\qquad
\mathcal M:\ u=h+g-e.
\]
Necessity follows from the displayed maps and the preservation of e,u having the same parity. P and K are disjoint and each has \(N^3/8\) states; M has \(N^3/4\); it meets each of the other charts in \(N^2/4\) states. Thus there are exactly \(N^2(N-1)/2\) states. The proof below constructs positive words of length \(O(N)\).

**A bit-clearing Tribonacci action on M.** In coordinates \((h,g,d)\) with \(d=u-h=g-e\), the append M has matrix
\[
A=\begin{pmatrix}2&-1&2\\1&0&0\\0&1&-1\end{pmatrix},
\]
with characteristic polynomial \(X^3-X^2-X-1\) and determinant 1. It preserves odd h,g and switches the parity of d. The last rows of \(A^4,A^8\) are \((2,0,1)\) and \((24,-4,17)\); \(A^4-I\) is divisible by 2 and \(A^8-I\) by 4. For \(j\ge2\), squaring a matrix \(I+2^jB\) preserves B modulo 2 at the next level, so the last row of \((A^{2^{j+1}}-I)/2^j\) is \((0,1,0)\) modulo 2. First make d even by at most one M; then use \(A^{2^{j+1}}\) when needed to clear bit j, starting at j=1. The parity of h,g supplies the bit toggle at each stage. Lower bits stay fixed, and the total exponent used is less than \(2N\). Since \(A^{2N}=I\bmod N\), we can reverse the bit-clearing orbit by a positive power. Thus it is enough to construct every wall state \(W(H,G)=(H,G,H,G)\).

**Wall lifts by a one-call induction.** Modulo 2 the only wall state is the seed. Modulo 4 the four \((H,G)=(1,1),(3,1),(1,3),(3,3)\) are given by the empty word, PK, PKPM, PPKMM. For \(N\ge8\), assume all K states modulo \(N/2\) have been constructed. For a prescribed odd G modulo N and any odd w modulo \(N/2\), lift a K state with
\[
g=w,\quad u=w+N/4,\quad h=G-N/4\pmod{N/2}.
\]
On K write \(x=u-g\); K sends \(x\) to \(2x\) and h to \(h+x\). One extra K therefore gives \(x=N/2\bmod N\), while \(h=G,\ u=w\bmod N/2\). An optional second K toggles the top bit of h without changing u modulo \(N/2\). This realizes prescribed h modulo N and any prescribed odd u modulo \(N/2\). To obtain an even u modulo \(N/2\), use its odd part followed by enough P's to double it to that residue (for zero use \(a-1\) P's). In particular arrange \(h=G,\ u=(H-G)/2\bmod N/2\), and then PM produces \(W(H,G)\).

**P and K from M.** Appending P to some M state with the desired h and g and half the desired u fills P. On K use \(x=u-g,\ y=h-u\); both are even. Appending K to P can fill any x congruent to 2 modulo 4 and any even y at fixed g: solve \(x=2u_{\rm old}-2g\), \(y=h_{\rm old}-u_{\rm old}+g\), with \(u_{\rm old}\) even. Further K's double x and fix y, reaching every even x, including zero. Modulo 2 the only K state is the seed. The induction now closes with one smaller-modulus call and \(O(N)\) extra moves. No interpolation or extrapolation from a finite census is used.

## 2. A mixed-modulus case follows without assuming independent words glue

**Theorem.** The same parity-compatible chart union is the reachable set modulo \(2m\) for every odd \(m\ge1\).

For K, use the earlier theorem to access a prescribed odd-modulus K preimage and then apply K, which is invertible on that chart modulo m. The final K forces the unique mod-2 K point regardless of the parity of the chosen word. P follows in the same fashion from K by P (invertible as a map of chart coordinates modulo m). On M, the map from K to M is unimodular in chart coordinates. It produces the mod-2 M point with d odd. To produce the other mod-2 M point, select the odd-modulus target one inverse M back, and append one extra M. This flips d parity while reaching the required odd part. Thus the *same physical word* serves both moduli. There is no claim here about simultaneous reachability when the modulus has both an odd factor and a factor of 4 or higher.

## 3. A literal isomorphism simplifies the exact backward search

For every starting pointed tournament the two-letter words MK and PM give isomorphic pointed tournaments, fixing every old vertex. Write a for the old root and b,c for the newcomers. For MK the triangle is \(a\to c\to b\to a\), with root b. For PM it is \(a\to b\to c\to a\), with root c. All old nonroots point to both newcomers, so swapping b,c is the required isomorphism. Subsequent moves depend only on the pointed tournament and preserve the equivalence. Each rewriting MK to PM removes one K, so it terminates in an equivalent word with no MK.

**Exact inverse consequence.** In this presentation, if the last letter is K then its unique possible predecessor *state* has
\[
u_{\rm old}=(u+g)/2,\qquad
e_{\rm old}=\begin{cases}
g,&u_{\rm old}\text{ odd},\\
0,&u_{\rm old}\text{ even},
\end{cases}\qquad
h_{\rm old}=h+e_{\rm old}-u_{\rm old},\qquad g_{\rm old}=g.
\]
Its preceding letter is K or P, or it is the seed, since MK is absent. The parity distinguishes these chart choices. This extends the earlier off-wall uniqueness across the wall. It is a necessary predecessor test, not automatic feasibility: the recovered state can violate positivity or fail to have an earlier word, and a terminal state can also admit a different last letter. It targets the exact descent rather than just increasing the congruence bank.

## Replay

Run node replay.js for independent constructive target checks through modulus 32, independent dyadic breadth-first searches through 64, literal rooted swaps, and finite checks of the mixed \(2m\) theorem. The proof is uniform. The results do not show that congruence compatibility implies positive integer compatibility, and they do not give any new obstruction on the scalar height alone.

# Simultaneous local inverses for a physical Farey flip

**PROVED.** Fix any of the eight canonical spine lengths \(4,8,\ldots,32\) with left tags 12 and middle tags 35. For every modulus \(m\) and every odd integer \(N\), a binary left address of \(O(\log^2(2m))\) bits and an empty right flank give a physical flip with heights \(N,N+2\bmod m\). Either prescribed left parity works. The same word carries all prime-power conditions for the given modulus. On this very axis at lengths 28 and 32, the exponent-61 absolute crossing is impossible by the exact one-sided inverse. [Physical Farey family](BERGGREN-RECEIVERS.md) · [Axis exclusions](FERMAT-RECEIVER-SEARCH.md)

The local reading also distinguishes the eight lengths. At the product of the 15 supersingular primes, the length-12 chart has relative target density **49/34** for two long independent words. The other seven have relative density below one. Only 17 and 29 affect these factors. The same densities hold if the primes are raised to any positive powers.

## A complete local test for arbitrary two-flank charts

The inherited binary read is \(f(0)=0, f(2u+1)=f(u), f(2u)=f(u)+u+[u \text{ odd}]\). Put \(F=f/2\). It is integer-valued. Fix a physical unit chart
\[
 L=K+t x+c y+qxy+r[x \text{ odd}]+e[y \text{ odd}]+h(f(x)+f(y)).
\]
A prescribed parity pair \(a,b\) means \(x=2u+a\), \(y=2v+b\). Put
\[
 K_{ab}=K+(t+r)a+(c+e)b+qab,\quad
 A=t+qb,\quad B=c+qa,\quad g_{ab}=\gcd(h,A,B,2q).
\]
For every positive \(M\), the actual lower heights in that branch modulo \(2M\) are **exactly** the classes \(K_{ab}\bmod 2\gcd(M,g_{ab})\). The odd-modulus version follows by ignoring the factor two. In particular, locality at every modulus is equivalent to \(N=K_{ab}\bmod 2g_{ab}\) for that branch.

Indeed, the automaton below realizes \(u,F(x),v,F(y)\) freely modulo \(M\). After subtracting \(K_{ab}\) and dividing by two the image is that of
\[
                 A u+B v+2q uv+hZ .
\]
Work one prime power at a time. Divide out the greatest common power of the prime in the coefficients and the modulus. If the resulting modulus is one, there is nothing to solve. Otherwise a remaining coefficient is a unit. If \(h,A\) or \(B\) is a unit, use the corresponding linear variable with the others zero. In the remaining case \(2q\) is a unit and \(A\) a multiple of the prime, so set \(v=1\) and solve in \(u\). CRT combines the residues of all three variables. One may take \(Z=F(x)\) and \(F(y)=0\). This proves the converse as well as the gcd obstruction. The chart coordinates and the parity pair remain fixed throughout the lift.

For the original exponent-15 tags \((1,35,451)\), all four gcds are 5. The even-left branches have lower height 7 modulo 10 and the odd-left branches 9 modulo 10. These two progressions are its **complete local image**. Consequently 5 supplies its only odd local obstruction to the exponent-61 lower count. At tags \((12,35,451)\) all four gcds become one. For the canonical tower with left tags 12 the stronger one-sided calculation gives \(\gcd(h,t)=1\) at each of the eight lengths. For \(y=0\) solve the single integer equation
\[
 t u+h z=(N-K-(t+r)a)/2 .
\]
The same fixed pair of integer coordinates \((2u+a,z)\) can be realized by an actual address modulo \((2M,M)\) for every \(M\). Only the address changes with the precision. [Exact profiles and physical seed](E15-FAREY.md) · [Controller and certificate generator](local_word_inverse.py)

## A short word between any two states

**PROVED.** The binary automaton with state \((X,Z)=(x\bmod 2M,F(x)\bmod M)\) has a single length \(L_M=O(\log^2(2M))\) such that every ordered pair of states is joined by some word of exactly that length. Its update on bit \(b\) is
\[
       X'=2X+b,\qquad Z'=Z+(1-b)\lceil X/2\rceil .
\]
The ceiling is taken on any representative of \(X\bmod 2M\). The previous reset proof already gave reachability, equidistribution and a short dyadic controller. The bound here is polynomial in bit precision for the odd part too. [Prior reset theorem](../chanpak_mirimanoff_unit_lifts_20261122/RESET-PORT.md)

The following switch writes a chosen difference into the second register. Take two nonnegative word values \(U,V\) of equal padded width \(b\) with \(U=V\bmod 2M\). Append either word after loading an even address \(X=2d\bmod 2M\). Both options end at the same address residue. The difference in the second register, in the direction \(V-U\), is
\[
 s d+C,\qquad s=\operatorname{rev}_b(U)-\operatorname{rev}_b(V),
                 \qquad C=F(V)-F(U).
\]
In fact concatenating \(W\) with an even prefix \(X\) gives \(f(2^bX+W)=f(X)+f(W)+(2^b-1-\operatorname{rev}_b(W))X\). This follows by the defining bit recurrence. Loading any desired address residue after any current state requires just \(w=\lceil\log_2(2M)\rceil\) bits by appending the \(w\)-digit representative of \(X_{\rm wanted}-2^wX_{\rm current}\bmod 2M\). A loader changes the height register by a quantity that depends only on the address register.

**Three switches suffice for every modulus.** Write \(2M=2^r m\), \(m\) odd. Choose \(u\in[0,7]\) congruent to \(5-m\bmod8\), and put \(v=u+m\). For the first two switches take
\[
  (U_0,V_0)=2^r(u,v),\qquad (U_1,V_1)=2^r(u+1,v+1)
\]
at the same width \(b\), the bit length of \(2^r(v+1)\). Take \((U_2,V_2)=(0,2M)\) at the bit length of \(2M\). Each switch has width \(w+O(1)\). In passing from \(u\) to \(u+1\), its least bit turns on. Since \(v=5\bmod8\), passing from \(v\) to \(v+1\) turns off its least bit and turns on the next. Reading those changes in reverse gives \(s_1-s_0=3\cdot2^{b-r-2}\). The third slope is odd. Therefore \(g=\gcd(M,s_0,s_1,s_2)\) divides three.

When \(g=3\), at least one of \(C_0,C_1\) is a unit modulo three. To see this, write \(u=2i, v=4j+1\), with \(j\) odd. Using
\[
 f(2^r z)=f(z)+(2^r-1)z+[z \text{ odd}]\qquad(r\ge1)
\]
gives \(2(C_1-C_0)=i+j+[i \text{ odd}]-[j \text{ odd}]\). Here the linear shift in \(z\) cancels in the double difference. The recurrence gives the unshifted increments \(f(v+1)-f(v)=j+2-[j \text{ odd}]\) and \(f(u+1)-f(u)=-i-[i \text{ odd}]\). Since \(3\mid m=v-u\), \(i+j=-1\bmod3\), which proves the claim. Exchange the two options of one unit-constant switch if necessary to ensure \(C=\sum C_i\) is nonzero modulo three.

**Synthesis.** For each \(j<k=\lceil\log_2 M\rceil\) make a group of these switches. Select \(\epsilon_j\in\{1,-1\}\) with \(\epsilon_j 2^j=C\bmod g\). Bézout's identity chooses even loading addresses \(2d_i\) so that
\[
             \sum_i s_i d_i=\epsilon_j 2^j-C\pmod M .
\]
Use a loader before every switch. In this group prescribe either all its first options or all its second options. The end addresses and all intermediate loader inputs agree in the two runs. The register difference is exactly \(\epsilon_j 2^j\). All groups can therefore be installed independently, followed by one final loader to the target address. Compare the desired \(Z\) with the baseline using every first option. Subset sums of the signed binary powers form an interval of \(2^k\) consecutive integers, covering every residue. Ordinary binary digits select the switches. This uses at most \(3k(2w+O(1))+w\) bits regardless of the incoming state, target, or prime factors of \(M\). A single switch suffices when its slope already has gcd one with \(M\), or gcd three with unit constant modulo three. The optional shortcut in the replay uses this observation.

The construction works after an arbitrary prescribed word. To prescribe terminal bits as well, append the fixed terminal word of width \(j\) and value \(v\) after steering. A final address \(X'\) is possible precisely when \(2^j X+v=X'\bmod 2M\) has a solution; choose any such \(X\) and subtract the known terminal-word increment from the desired \(Z\) before steering. In particular prescribed boundary words can coexist with every modular height-register value whenever the address congruence permits them.

At \(M=2^s\), each step has a cellular implementation. Shift the \(s+1\) address bits, and write the input at bit zero. For the height update use the old bits \(x_i,z_i\), carry \(d_0=(1-b)x_0\), and addend \(a_i=(1-b)x_{i+1}\). Bit \(i\) writes \(z_i\mathbin{\rm xor}a_i\mathbin{\rm xor}d_i\) and sends \({\rm majority}(z_i,a_i,d_i)\) to cell \(i+1\). With the fixed input distributed to the register, one ripple pass suffices. Finite binary tapes also implement the general-modulus loader, reversal, addition and extended Euclidean algorithm. The proof bounds the number of bits a deterministic steering algorithm writes. The finite-state observer at fixed precision has \(2M^2\) states; reaching every state from one origin at the same length already needs at least \(\lceil\log_2(2M^2)\rceil\) bits.

## What the supersingular primes detect

Duncan and Ono explain the agreement between the prime support of the Monster and the 15 characteristics in which all supersingular elliptic \(j\)-invariants lie in the prime field. Here they supply a specified list of local tests. The tournament density is determined by the coefficients of the physical flip. [Duncan and Ono, *The Jack Daniels Problem*](https://arxiv.org/abs/1411.5354)

Let \(P=2\cdot3\cdot5\cdot7\cdot11\cdot13\cdot17\cdot19\cdot23\cdot29\cdot31\cdot41\cdot47\cdot59\cdot71=1{,}618{,}964{,}990{,}108{,}856{,}390\). Let two independent binary addresses with leading one grow to arbitrary lengths, with all other bits independent and fair. The modular state has uniform stationary distribution because every state has two bit-labelled predecessors. The fixed-length steering theorem makes the finite chain primitive. Consequently its distribution converges to uniform, also conditionally on the terminal bit. For the two addresses the CRT components are independent after the two parities are fixed.

For the E61 lower target and **any** modulus \(D\) whose prime divisors are exactly these 15 primes, normalize the limiting probability of a match modulo \(D\) by the uniform odd-residue probability \(2/D\). The resulting ratios are

| Spine length | Right word even | Right word odd | Unconditional |
|---|---:|---:|---:|
| 4, 20, 28 | 448/493 | 448/493 | 448/493 |
| 8, 16, 24, 32 | 28/29 | 28/29 | 28/29 |
| 12 | 16/17 | 33/17 | 49/34 |

Both left parities give the same entries. For an odd prime \(\ell\nmid h\) the height register makes the observer uniform. If \(\ell\mid h,q\), the unit left slope makes the reduced bilinear law uniform at the prime. At higher powers, if \(q\) is divisible by \(\ell\), then \(t+qy\) stays a unit for every \(y\); the uniform left-address residue gives uniform output with the independent register readings held fixed. At 2, \(h\) is odd so the even nonlinear increment is uniform among even residues, conditional on the parities. The only remaining primes in the list are 17 and 29, each to exponent one in \(h\) when it contributes. At \(\ell\mid h, \ell\nmid q\), complete the product:
\[
        (qx+c)(qy+t)=q(N-K-ra-eb)+ct\pmod\ell .
\]
A nonzero right side gives \(\ell-1\) solutions in the two residues. Zero gives \(2\ell-1\). Dividing by \(\ell\) gives the local relative density. Since \(v_\ell(h)=1\) here, the independent height register supplies a uniform lift above the first power. At \(n=12\) and odd \(y\) the right side is zero modulo 17. The two lines are exactly \(x=4\) or \(y=16\bmod17\). All other contributing product tests have a nonzero right side. This proves the table, including arbitrary positive exponents at all 15 primes.

The radical \(P\) has 61 binary digits. Its square is about 98.59 percent of the prescribed positive lower count \(N=2^{60}(2^{61}-1)-1\). The replay constructs left addresses on every spine of precisely 7,810 and 30,131 bits, giving the two exact endpoint residues modulo \(2P\) and \(2P^2\) respectively. Although \(2P^2>N\), these words produce much larger tournaments and counts. A height bound \(0<L<2P^2\) on that same construction would turn its modular success into equality. This is a concrete missing size condition.

From this directory run `python3 -B local_word_inverse.py`. It checks the unconditional controller without the shortcut on moduli through 128, prints the density factors, and checks all 16 actual addresses against the original integer height formula. For a full binary address use `python3 -B local_word_inverse.py --n 12 --modulus 1618964990108856390`. The modulus option specifies \(M\), giving a certificate modulo \(2M\). The separately constructed E15 pair itself supplies one common ordinary address pair across every modulus for its own target \(E_{15}-1\). [Exact E15 certificate](E15-FAREY.md)

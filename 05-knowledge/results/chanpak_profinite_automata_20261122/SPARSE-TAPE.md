# Sparse tape programming and an absolute-reader obstruction

Use the one fixed strong physical seam in [RESET-PORT](../chanpak_mirimanoff_unit_lifts_20261122/RESET-PORT.md), whose two actual Hamilton counts are
\[
 L(y)=149+162y+25f(y),\qquad L(y)+2.
\]
The common marked arc and core are unchanged. The binary reader is \(f(0)=0\), \(f(2u+1)=f(u)\), \(f(2u)=f(u)+u+[u\text{ odd}]\). In particular \(f\) is always even. The all-modulus output theorem for this seam is inherited; the new results constrain the addressing words, prescribe the whole two-register state, and separate modular control from exact finite automata.

## 1. Sparse two-register programming

Let \(M=2^r m\) with \(m\) odd, \(r\ge1\); for an odd modulus one may adjoin a factor of 2. Fix \(L_0\ge r\) positive and divisible by \(\operatorname{ord}_m(2)\) (take order 1 at \(m=1\)). Start at a positive even address \(y\equiv0\pmod{2^r}\); from any existing word one can ensure this by appending \(L_0\) zeros. The following chronologically appended macros all end in at least \(r\) zeros and contain no consecutive ones:
\[
 P=0^{L_0+1},\quad Q=1\,0^{L_0},\quad R_s=1\,0^s\,1\,0^{L_0}\quad(1\le s<r).
\]
Modulo \(2^r\), \(P,Q\) preserve \(f\) and \(R_s\) adds \(-2^s\) to \(f\); all keep the address zero there. To prove the latter shift, for an even address \(u\) appending \(1\,0^j\) changes address to \(2^j(2u+1)\) and changes \(f\) by \(2(2^j-1)u+2^j\). The first half of \(R_s\) therefore changes \(f\) by \(2^s\) modulo \(2^r\), and the second by \(-2\cdot2^s\). Choose \(R_s\) at the set bits of the negated desired even difference to set the dyadic read in at most \(r-1\) macros.

Modulo \(m\), \(P\) acts on \((y,f)\) as \((2y,f+y)\), and \(Q\) as \((2y+1,f+1)\). Their state actions are finite permutations. In fact \(P^{L_0}\) is the identity there, so the positive word for \(T=Q\circ P^{-1}\) is \(P^{L_0-1}\) followed by \(Q\). It acts by
\[
 T:(y,f)\longmapsto(y+1,f+1-y/2)\pmod m .
\]
Since \(m\) is odd, summing the increments over \(m\) successive values of \(y\) shows \(T^m\) is the identity. Apply at most \(m-1\) copies of \(T\) to set any prescribed odd-modulus address. At this fixed address, translate it to \(-1\) by a power of \(T\), apply any number \(j\) of copies of \(Q\), then translate back using the inverse power. At \(-1\), \(Q\) fixes the address and adds 1 to \(f\); because every constituent map has coefficient 1 on the old \(f\), the conjugate adds \(j\) at the original address. All these operations preserve the dyadic read. This sets arbitrary residues of \(y\bmod m,\ f\bmod m\), and arbitrary even \(f\bmod2^r\) simultaneously, with \(y=0\bmod2^r\) and with no consecutive ones in the appended tape.

The program uses \(O(r(r+L_0)+mL_0^2)\) appended bits: at most \(r-1\) dyadic macros and \(O(m)\) uses of \(T,Q\), each \(T\) being a positive word of length \(O(L_0^2)\). Thus for fixed odd primes and fixed odd exponents, the dyadic precision cost remains quadratic in \(r\), despite the sparse alphabet constraint.

To set an arbitrary \(a=y\bmod M\), finish with its padded \(r\)-bit suffix \(b=a\bmod2^r\). Program the temporary prefix \(z\equiv 2^{-r}(a-b)\pmod m\), \(z=0\bmod2^r\), then append \(b\). The inherited suffix identity gives
\[
 f(2^rz+b)=f(z)+(2^r-1-\operatorname{rev}_r b)z+f(b)
\]
because \(z\) is even. The fixed correction has fixed residue modulo \(M\), so precompensate it in the read. Only this optional suffix may contain consecutive ones. This proves that the whole graph of \(f\) is dense in \(\widehat{\mathbb Z}\times 2\widehat{\mathbb Z}\), not merely that a linear output is dense. In particular \(f\) has no continuous profinite extension on ordinary addresses: even at any one positive integer address there are convergent sequences of addresses with arbitrary prescribed, hence different, read phases. These are finite congruence statements; a word in a residue class need not attain its least representative as an absolute count.

### Sharp cost of wider symbol separation

If every pair of successive ones is separated by at least \(K\ge1\) zeros and a checkpoint ends in at least \(K\) zeros, then both its address and read are zero modulo \(2^K\). This follows from the isolated-one block formula \(f\mapsto f+2(2^j-1)y+2^j\) with \(j\ge K\), inducting on the ones (the first block has read \(2^j\)). This is the only loss in the controllable boundary residues: choose \(L_0\ge\max(K,r)\), restrict \(R_s\) to \(s\ge K\), and leave all odd-address operations unchanged. At dyadic precision \(2^r\) this reaches exactly the read subgroup \(2^K\mathbb Z/2^r\mathbb Z\) whenever \(r\ge K\), with the temporary address zero. Consequently the sparse strong seam still attains every lower residue compatible with \(L(y)\equiv149\pmod{2^K}\), simultaneously with arbitrary odd-prime phases. The exclusion \(K=1\) is the largest uniform minimum gap admitting all odd residues.

## 2. One append-only tape for simultaneous coherent local inverses

For the fixed physical seam, the closure of the even-address states is \((2\widehat{\mathbb Z})^2\). The exact linear map from state to lower height has a continuous section there: since \(-2\cdot162+13\cdot25=1\), for any odd profinite target \(\tau\) choose
\[
 a=-2(\tau-149),\qquad b=13(\tau-149).
\]
Thus local inverse phases can glue as two-register state, although they cannot define a continuous absolute forward read from the ordinary address register alone.

For a section that is writable as **one sparse infinite tape**, instead choose at each prime
\[
 (a_p,b_p)=
 \begin{cases}
 ((\tau_p-149)/162,\,0),&p=5,\\
 (0,\,(\tau_p-149)/25),&p\ne5 .
 \end{cases}
\]
Here the divisions use units, and \(b_2\) is even. In particular \(a_2=0\), so no nonsparse final suffix is needed. Pick any increasing sequence of precision moduli with each dividing the next. At each stage extend the current sparse word by zeros (resetting the stricter dyadic address), then run the finite program from §1 for the new residues. No earlier bits are changed. The state and the literal pair of counts converge along the chosen checkpoint prefixes at every prime whose requested precision tends to infinity. This is equally valid for any chosen finite set of primes, including all fifteen Ogg primes, or for all primes using a cofinal modulus sequence. There is no claim that the read converges along *all* prefixes, nor that one finite checkpoint realizes the exact integer target. A Turing machine with access to the requested local digits constructs the tape by the explicit loops above.

## 3. The exact reader is not a synchronized automaton

Say that a binary integer function is synchronized if the graph of its value is a regular language of padded same-length input/output columns. For \(y=(1^a0^b)_2\), \(a,b\ge2\), the recurrence gives
\[
 f(y)=2^{a+b}-2^a-2^b+2.
\]
If \(f\) were synchronized on these inputs, so would be the graph of \(g(y)=2^{\operatorname{bitlength}(y)}-f(y)+2\). Indeed padded binary addition/subtraction and multiplication by fixed constants have finite digit-carry relations, and the power one place to the left of the leading input bit is also synchronous. Regularity survives composing these relations by existentially projecting auxiliary tracks. Now \(g(y)=2^a+2^b\) is a power of 2 exactly when \(a=b\). Intersect its output track with the regular language of powers of 2 and project to obtain the nonregular language \(\{1^a0^a:a\ge2\}\), a contradiction. For any fixed \(K,c,h\) with \(h\ne0\), a synchronized graph of \(K+cy+h f(y)\) would give one for \(f\) by the same fixed-coefficient arithmetic relation. Hence neither \(f\) nor the literal seam height \(149+162y+25 f(y)\) is synchronized even restricted to those even addresses. This is stronger than failure of a particular deterministic streaming algorithm.

The synchronized obstruction persists under the very same sparse language of §1. On \((10^a10^b)_2\), \(a,b\ge1\), the isolated-one formula gives \(f=2^{a+b+1}-2^a+2^b\). The regular relation “output is the bit of the input's leading position” means output \(2^{a+b+1}\), and it holds exactly when \(a=b\). Projection would force the nonregular language \(\{10^a10^a:a\ge1\}\). Thus excluding consecutive ones helps controllability without turning the exact reader into a synchronized one.

A fixed-radius cellular automaton with input and output aligned in ordinary binary also needs \(\Omega(n)\) time at bitlength \(n\), in either tape orientation. For the most-significant end, \(2^n\) and \(2^n+2\) share all but the final two positions, but their reads differ by \(2^{n-2}-2\), as their reads are \(2^n\) and \(3\cdot2^{n-2}+2\), respectively. The two output heights therefore differ by magnitude \(\Theta(2^n)\); with a common \(n+O(1)\)-bit padded output, at least one of the first \(O(1)\) output bits differs, at linear distance from the changed input. For the least-significant end, \(4\cdot2^n\) and \(5\cdot2^n\) have the same length and share at least their lowest \(n\) bits, but the zero-run formula gives reads differing by \(2\bmod4\) for \(n\ge2\). The literal seam heights do so as well. A radius-\(R\) CA cannot transmit influence farther than \(Rt\) in \(t\) steps. This lower bound is for aligned binary exact output without a precomputed nonlocal auxiliary track; modular finite-state reading and the append-only checkpoint construction are compatible with it.

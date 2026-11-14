# The physical 3+5 bits, the missing size, and a reciprocal obstruction

A terminal-ready rooted tournament has \(L=H\). Attach the near-transitive \(R_j\) by the inherited ordered seam and reroot at its new terminal. This multiplies the Hamiltonian-path count by \(2^j+1\), adds \(j+1\) vertices, and returns to \(L=H\). Starting at a singleton gives literal strong tournaments as soon as the word is nonempty. With the three genuine Fermat blocks \(j=1,2,4\) the exact counts at added order \(M\) are
\[
H=3^a5^b17^c,\qquad 2a+3b+5c=M,\qquad a,b,c\ge0.
\]
The order of those blocks does not change this scalar. Below the size and the residues are both physical, not two independent freedoms.

## The two elementary equal-size substitutions

Trade one factor 17 for \(3\cdot5\) at cost five. Put \(A=a+c\), \(B=b+c\), so \(2A+3B=M\), \(0\le c\le\min(A,B)\), and
\[
H=3^A5^B\lambda^c,\qquad
\alpha=27/25,\quad\lambda=17/15.
\]
The other equal-size trade \(3^3\) for \(5^2\) costs six and changes the result by \(\alpha\) or its inverse. All ratios here are taken only in rings where their denominators are units.

Modulo \(2^n\), \(n\ge3\), \(\alpha\equiv3\pmod8\) and \(v_2(\alpha^2-1)=3\). Repeated squaring gives order \(2^{n-2}\) (with the same direct check at \(n=3\)); its subgroup consists exactly of the units congruent to 1 or 3 mod 8. Since \(\lambda\equiv7\pmod8\), these two equal-cost substitutions generate the entire unit group. With only factors 3 and 5, at a fixed order all outputs lie in one coset of \(\langle\alpha\rangle\) and cannot cover it and its complement.

In the zero-based odd address \(x=(H-1)/2\), the low three bits are \(H\bmod16\); a full eight-bit address is \(H\bmod512\). Thus the physical fixed-size problem really is
\[
1\longrightarrow \langle\alpha^4\rangle \simeq C_{32}
\longrightarrow (\mathbb Z/512)^\times
\longrightarrow(\mathbb Z/16)^\times\simeq C_2\times C_4
\longrightarrow1 .
\]
The group in the middle is \(C_2\times C_{128}\), not a direct product \(C_{32}\times C_2\times C_4\): their exponents are different. The eight choices in the quotient are the three visible bits, the cyclic kernel supplies 32 choices per low address, and the extension retains a carry. This is also the familiar restriction map for the cyclotomic Galois actions on 512th roots and on 16th roots, where multiplication by a unit acts by raising the root to that power. Neither finite multiplication nor this restriction is an additive ring structure on the set of odd integers.

## Exact size frontier, including the one missing size

This argument works for every modulus below for which \(\alpha\) generates a cyclic index-two subgroup and \(\lambda\) is outside it. Write the subgroup as \(C_D\) with \(\lambda^2=\alpha^w\). At odd \(M=2u+1\) put \(B=1+2t,A=u-1-3t\); at even \(M=2u\) put \(B=2t,A=u-3t\). Write \(c=2s+\epsilon\) with \(\epsilon=0,1\). In its coset, after a fixed common unit factor, the output address is \(t-ws\) in \(\mathbb Z/D\) (reversing the coordinate exponent if needed). For odd size the respective allowed \(t\)-intervals are
\[
[s,\lfloor(u-2s-1)/3\rfloor],\quad
[s,\lfloor(u-2s-2)/3\rfloor],
\]
and for even size they are
\[
[s,\lfloor(u-2s)/3\rfloor],\quad
[s+1,\lfloor(u-2s-1)/3\rfloor].
\]
Ignore empty intervals and use \(s\ge0\). The \(\epsilon=0\) interval always contains the \(\epsilon=1\) interval at the same \(s\), so saturation of the allowed whole group or coset is equivalent to saturation of the latter interval family modulo \(D\).

For even \(2u\), subtract one from \(t\) in the second interval. The result is **exactly** the odd \(2u-3\) family, up to the harmless common translation of cyclic addresses. This proves an exact parity law:
> Size \(M\) even is full if and only if size \(M-3\) is full.

A full size stays full at sizes two and three higher (append the corresponding multiplier), relative to its expected coset when coupled. Consequently the first full added size \(f\) is odd, \(f+1\) necessarily fails, and every \(M\ge f+2\) is full. This is an all-order explanation of the isolated subsequent hole in the finite experiments.

For \(3\le n\le6\), \(\lambda^2\equiv1\) because \(17^2-15^2=64\). The odd second interval with \(s=0\) contains all others, and the criterion gives \(f=6\cdot2^{n-2}-1\). For the actual eight-bit modulus 512, \(\alpha\equiv83,\lambda\equiv479,\lambda^2\equiv\alpha^{16}\). Thus the limiting intervals at odd \(M=2u+1\) are
\[
[-15s,\lfloor(u-2s-2)/3\rfloor-16s]\pmod{128}.
\]
At \(M=157\) their union misses exactly address 37; at \(M=159\) their union is all 128. Adding two vertices propagates fullness upward through odd sizes. The parity law proves **159 is the first full added size** (literal order 160), **160 fails** with the unique missing Hamilton residue 87 mod 512, and **every added size at least 161 is full**. This gives a fixed-size realization of every hidden five-bit payload at every three-bit observation, not merely a count of all possible bit strings.

## Beyond the finite size checks

There is a geometric upper bound which improves the naive order proportional to the modulus. Work first mod \(2^n\), and take the lattice
\[
\Lambda=\{(t,s)\in\mathbb Z^2:\alpha^t(\lambda^2)^s\equiv1\pmod{2^n}\}.
\]
Its determinant is \(D=2^{n-2}\). The two positive rationals are multiplicatively independent: only \(\lambda^2\) has prime factor 17, then \(\alpha\) is nontrivial. For a nonzero lattice vector, clearing the odd denominators gives a nonzero multiple of \(2^n\) whose absolute value is less than \(289^{|t|+|s|}\). Thus every nonzero vector has length \(\Omega(n)\). In dimension two a shortest vector \(v\) is primitive. Extend it to a basis and subtract its integer multiples from the second basis vector so that the parallel component has length at most \(\|v\|/2\), with perpendicular length \(D/\|v\|\). Minimality of \(v\) bounds \(\|v\|=O(D/\|v\|)\) as well. Rounding coefficients in this basis gives covering radius \(O(D/\|v\|)=O(2^n/n)\).

In particular all cosets of \(\Lambda\) have representatives in \([0,T]^2\) for \(T=O(2^n/n)\): approximate the center of the square after choosing \(T\) twice the radius and rounding it up. To realize these simultaneously at any added size \(M\ge16T+7\), choose \(A_0\ge2T+1,B_0\ge4T+1\) with \(2A_0+3B_0=M\) (the remaining integer is at least two and is a nonnegative combination of 2 and 3). Put \(A=A_0+3t,B=B_0-2t,c=2s+\epsilon\). They satisfy every nonnegativity restriction. First choose \(\epsilon\) to reach the right coset and then choose \((t,s)\) in the square to reach its exact residue. This proves an \(O(2^n/n)\) conductor. Counting parameter pairs gives the distinct obstruction \(M=\Omega(2^{n/2})\) for any full size. The proof leaves a real gap between these bounds; it does not assume random \(2\)-adic digits.

### Infinitely many precisions have optimal square-root cost

In fact the lower square-root scale is attained along infinitely many precisions, with the *whole* stable tail above that cost. Let \(v_n\) be a shortest vector in \(\Lambda_n\) and \(D_n=2^{n-2}\). The next lattice has index two in this one and contains \(2v_n\), hence \(\|v_{n+1}\|\le2\|v_n\|\). The direction of the shortest vector must change infinitely often. Otherwise eventually all \(v_n=m_n v\) on a fixed primitive integer direction \(v\). The associated positive rational \(z=\alpha^{v_1}(\lambda^2)^{v_2}\) is neither \(1\) nor \(-1\). For odd rational numerator and denominator the factorization of \(z^m-1\), repeatedly separating \(z^j+1\) after a square, gives
\[
v_2(z^m-1)\le C_z+v_2(m).
\]
Thus \(|m_n|\gg 2^n\). But the lattice has a vector of length \(O(\sqrt{D_n})\): pigeonhole the integer points in a square with more than \(D_n\) points by its \(D_n\) lattice cosets and subtract two in the same coset. Contradiction.

At a direction change, \(v_n,v_{n+1}\) are independent vectors of \(\Lambda_n\), so
\[
D_n\le|\det(v_n,v_{n+1})|\le2\|v_n\|^2 .
\]
Hence \(\|v_n\|\) is of order \(\sqrt{D_n}\) at infinitely many \(n\). The covering-radius argument then supplies all phases for every \(M\ge C\sqrt{2^n}\) for a uniform constant \(C\) on those levels. Combined with the counting lower bound, the first/full stable order cost is \(\Theta(2^{n/2})\) along infinitely many precisions. This does not settle the worst-case precision cost; it uses forced changes of shortest direction, not a conjecture on ordinary digits of a \(2\)-adic logarithm.

## A quadratic coupling invisible in either single modulus

Now impose \(2^n\) and \(7^k\) together, \(n\ge3,k\ge1\). **Every** genuine Fermat block satisfies
\[
\left(\frac{7}{2^j+1}\right)=(-1)^{j+1}
\]
as a Jacobi symbol. For \(j=1\) both sides are \(+1\). For genuine even \(j\), \(2^j+1\equiv1\pmod4\), and modulo 7 it is either 3 or 5 since \(j=2^a\) is not divisible by 3. Both are nonresidues. The Jacobi symbol therefore gives, at fixed added order,
\[
\left(\frac{H}{7}\right)=(-1)^M\left(\frac{-1}{H}\right).
\]
Here \((H/7)\) is the Legendre sign; \((-1/H)\) is just the mod-4 sign. This is an exact order-dependent coupling, not a condition discovered only in the finite table. No larger genuine Fermat block breaks it. In fact \(7\nmid 2^j+1\) for every \(j\), genuine or not, because 2 has order three modulo 7.

For the three smallest genuine blocks the coupled half is **exactly** the eventual joint image at each size. Modulo 49, \(\alpha\equiv5\) and \(5^6\equiv43\ne1\), while 5 generates the units modulo 7. Thus \(\alpha\) generates modulo all \(7^k\). The power \(\alpha^{2^{n-2}}\) fixes the \(2^n\) component and generates every square unit modulo \(7^k\). Since \(\alpha,\lambda\) already cover the \(2^n\) component, the joint group has exactly the prescribed half. More specifically \(\langle\alpha\rangle\) in the joint ring has order \(D=2^{n-2}\cdot3\cdot7^{k-1}\) and index two inside the allowed subgroup of equal-size ratios; the other coset is supplied by \(\lambda\). The method of intervals applies with this \(D\). The lattice argument above also applies using the joint modulus and gives a conductor \(O(2^n7^k/(n+k))\) for the allowed half.

### The cheapest repair by a single different near-transitive exponent

Replace three count-3 blocks (six vertices, scalar multiplier 27) by **one \(R_5\)** (six vertices, multiplier 33). This is not a genuine Fermat exponent. Since \((7/33)=-1\) while \((7/27)=+1\), it flips the quadratic gate without costing any additional vertices. The four shorter exponents \(j=1,2,3,4\) all obey the gate (the extra control \(j=3\) is a square multiplier 9, at even cost). At sufficiently large common \(M\), the original alphabet attains the proper half, and a word with a single \(R_5\) followed by an original word of size \(M-6\) attains the other half. Thus the new alphabet obtains every joint unit at a common order, not through a change to the target's size.

For the concrete modulus \(512\cdot7=3584\), independent finite replay first attains all 768 permitted units at added size 319 with the three Fermat blocks; adding \(R_5\) first attains all 1536 units at added size 173. These are literal phase coverage theorems, not absolute Chan–Pak coverage. All counts in each named alphabet still factor over finitely many primes. Up to \(X\) their number grows only as a power of \(\log X\), even though entire quotient spaces are filled at common orders.

## A rank law for several odd phases

For any odd prime \(p\nmid3\cdot5\cdot17\), let \(\chi_p\) denote the quadratic character of nonzero residues. There is a unique character \(\psi_p\) on \((\mathbb Z/8)^\times\) with
\[
\psi_p(3)=\chi_p(\alpha)=\chi_p(3),\qquad
\psi_p(7)=\chi_p(\lambda)=\chi_p(17)\chi_p(3)\chi_p(5).
\]
Indeed 3 and 7 are the two independent nonidentity generators mod 8. In particular \(\psi_p(5)=\chi_p(17)\chi_p(5)\). Check on the three literal multipliers that
\[
\chi_p(H)=\chi_p(17)^M\psi_p(H\bmod8)
\]
for every original word at size \(M\). Thus for any list of \(r\) such odd primes there are \(r\) independent sign constraints once the binary address and order are fixed.

Extra types of blocks with multipliers \(F\) and costs \(b\) (all coprime to the moduli) contribute one vector each of sign violations, whose coordinate at \(p\) is
\[
\chi_p(F)\,\psi_p(F\bmod8)^{-1}\,\chi_p(17)^{-b}.
\]
The reachable sign syndromes lie in the span over \(\mathbb F_2\) of these vectors. Conversely, for every sufficiently large common \(M\) and every specified binary phase, every syndrome in that span actually occurs: use at most one of each selected extra type, then adjust the remaining base word at its lower order using full binary saturation. This is an exact rank test for independent **quadratic** phases, not a claim about higher odd residue coordinates in general.

There is a second full-product repair. For \(2^n\), \(7^k\), and \(13\) together, the original alphabet's eventual image is exactly the quarter determined by the two quadratic signs. To see there are no other constraints here, let \(D=2^{n-2}\). The element \(\alpha^D\) fixes the binary and the 13 coordinates (since \(\alpha\equiv-1\bmod13\)) and spans the square units modulo \(7^k\). Choose the even \(w\) such that \(\lambda^2\equiv\alpha^w\bmod2^n\) (even since the square is 1 mod 8). The element \(\lambda^2\alpha^{-w}\) fixes the binary coordinate and is 4 mod 13, generating the squares there because \(\lambda\equiv2\bmod13\). Its seven coordinate is a square, removable using \(\alpha^D\). Thus the subgroup over binary identity contains every square pair, exactly the required kernel. Positive words at common large order reach the full coset: as above allow all integer trade exponents over the finite group, choose nonnegative representatives by adding their finite orders, then choose the baseline \(A_0,B_0\) large enough for all of them. No conjecture about primes is used.

Now \(R_5\) of multiplier 33 contributes sign vector \((-, -)\) at \(7,13\); \(R_7\) of multiplier 129 contributes \((-,+)\). Their costs 6 and 8 are even and their multipliers are 1 mod 8, so these signs follow by reducing the multipliers themselves modulo 7 and 13. They are independent. At sufficiently large fixed order the alphabet enlarged by **both** types gives every joint unit modulo \(2^n7^k13\), whereas either one alone still leaves at least half missing. This is a true two-axis phase repair without extra vertices relative to the assigned size; the price is changing the composition.

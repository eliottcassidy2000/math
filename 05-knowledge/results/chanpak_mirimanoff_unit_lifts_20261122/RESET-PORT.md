# One fixed arc has every compatible residue, with a short dyadic controller

**PROVED.** Use the fixed [exponent-31 internal arc and binary receiver](../chanpak_e31_physical_20261122/README.md) with receiver address \(v=1\). This is a tagged two-vertex transitive receiver. With a common right suffix address \(y\ge 0\), the two strong tournaments have counts
\[
            L(y)=149+162y+25 f(y),\qquad L(y)+2,
\]
where \(f(0)=0,\ f(2u+1)=f(u),\ f(2u)=f(u)+u+[u\ {\rm odd}]\).
The rows of the fixed core are \((4,5,24,19,35,15)\), the initial tag is vertex 0 and the flip is \((1,3)\). The joined receiver adds a tag at vertex 6. Vertex 0 reaches the entire core in either orientation, and the untagged receiver vertex is reachable from the entire deletion and beats the outer apex. The apex reaches 0. Thus each host is strong for any suffix (in the ordinal deletion each suffix vertex sees the receiver in the forward direction, and the last suffix vertex when present beats the apex). The physical proof of the constant gap is the previously established transfer \(h_V\Delta k-R_V\Delta z=2\#\{\text{full paths starting tagged in }V\}\).

## Exact reset law and a constructive all-modulus selector

Let \(q=2^r m\), \(m\) odd, with \(r\ge1\); for a purely odd modulus we may impose an extra factor 2. Choose \(t\ge r\) a positive multiple of \(d=\operatorname{ord}_m(2)\), taking \(d=1\) for \(m=1\). A word \(B_a=1^a0^t\) of binary digits, \(a\ge1\), appends to the address. The intermediate state after the ones is \(A=2^a(y+1)-1\). It is odd and it preserves \(f\). By applying the \(t\) zero steps,
\[
  y'=2^t A,\qquad f(y')-f(y)=(2^t-1)A+1 .
\]
From \(y=0\bmod 2^r\) this keeps \(y'=0\bmod 2^r\) and moves \(f\) by \(2-2^a\bmod 2^r\). Modulo \(m\), it multiplies \(y+1\) by \(2^a\) and adds 1 to \(f\), independently of the incoming address. In particular \(a=t\) leaves the whole address fixed modulo \(q\) and adds the CRT lift of \((2\bmod 2^r,1\bmod m)\). This generates the even residues modulo \(q\); the addend \(25\) in the read is not a unit at 5, but that causes no loss: fix \(y\) first modulo \(g=\gcd(25,m)\) using \(\gcd(162,g)=1\).

Here is a selector with polynomial address cost in \(r\) when \(m\) is fixed, avoiding up to \(q\) rounds of naive modular lifting. Let \(N\) be any odd target residue and select a positive multiple \(y_0\) of \(2^r\) with \(149+162y_0=N\bmod g\). Choose target residues for \(f\) modulo \(2^r,m\) so that \(149+162y_0+25 f=N\bmod q\); the first is even, the second exists by the choice of \(y_0\). Call the required differences from \(f(y_0)\) \(b_2,b_m\). They are computed by two linear congruences, with \(b_2\) even. Reserve a final block with \(a\ge r\) to add 2 to the dyadic component. Write \(-(b_2-2)/2\) modulo \(2^{r-1}\) in binary. For each set bit \(j\), append \(B_{j+1}\) then \(B_2\). The pair of blocks contributes \(-2^{j+1}\) modulo \(2^r\), including \(j=0\). Append \(0\) to \(m-1\) copies of \(B_1\) so that the total number of blocks, including the reserved one, is \(b_m\bmod m\). These copies do not change the dyadic component. At the final block choose \(a\) in \([r,r+d-1]\) making the sum of **all** \(a\)'s divisible by \(d\). This restores the address modulo \(m\), and still supplies the reserved \(+2\) modulo \(2^r\). No choice depends on intermediate growth.

The number of appended bits is \(O((r+m)(r+d))\), and the initial address can have \(O(r+\log(2m))\) bits. In particular all odd residues modulo \(2^r\) have explicit literal unit seams of order \(O(r^2)\) in this single-core family. Arbitrarily many more addresses with the same phase follow by appending the repeat block \(B_t\) a multiple of its translation order. This gives an infinite strong +2 family in every compatible phase at **every** modulus, not every absolute count. For instance for any perfect center one can specify as many low bits as desired. Increasing the precision changes the resulting tournament and can increase the absolute count; it is not an inverse at fixed height.

## A physical base-three probe

In the same six-vertex seam use receiver \(v=2^k\), \(k\ge1\), and no added right suffix. In the receiver's binary word \(10^k\) the new receiver apex beats only the first spine vertex. Its ordinary path count is \(1+2^k\); the endpoint mass from the published formula is \(2^{k+1}-1\). Its two-slot count is
\[
           Q=2(3^k+2^{k+1}-1).
\]
Indeed choose the apex slot. If it has a nonempty after-apex part that part starts at the first spine vertex, leaving three independent choices per other vertex. If its after-part is empty, the before-apex part can be any subset except the lone forbidden first vertex; the other path is forced. Plugging into the exact word rule gives the actual centers
\[
     M_k=56\cdot3^k+106\cdot2^k-12,\qquad
     I_k=M_{k+1}-2M_k-12=56\cdot3^k.
\]
Both measurements for \(I_k\) are on literal single-arc pairs. Thus for every prime \(\ell\ge5\), \(\ell\ne7\), the probe \(I_{k+\ell-1}\equiv I_k\bmod\ell^2\) is equivalent to Mirimanoff's condition. The subtraction is a measurement across two hosts, not a new tournament construction. Indeed \(M_k=4\bmod8\) for \(k\ge2\), so these uncorrected centers cannot equal any \(E_p\) for \(p\ge4\). This is one explicit example where the arithmetic modes are visible and the absolute perfect target is nevertheless missed. The right-suffix controller restores all local phases, not that missing absolute equality.

## The modular labels are equidistributed at any fixed address scale

For fixed modulus \(M\), set \((X,Z)=(y\bmod 2M,f(y)/2\bmod M)\). Appending bit \(e\) acts by
\[
     (X,Z)\mapsto(2X+e,\ Z+(1-e)\lceil X/2\rceil).
\]
The factor of two is harmless since \(f\) is always even. Under fair independent bits each state has exactly two predecessors, each of probability \(1/2\): the output address specifies \(e\) and leaves two old addresses; the old \(Z\) is then forced. Uniform measure on the \(2M^2\) states is stationary. It is also the limiting measure, exponentially fast for each fixed \(M\). Indeed a long repeat block as above, with \(t\ge v_2(2M)\), first puts the address in a class fixed under more copies and then translates \(Z\) by a unit modulo \(M\). From this class a prescribed final address is reached by a suffix of length at least \(\lceil\log_2(2M)\rceil\); the increment in \(Z\) from that suffix depends only on the starting \(X\), so it can be precompensated by block repetitions. This proves irreducibility; the loop at \((0,Z)\) proves aperiodicity.

For uniformly sampled positive even \(n\)-bit addresses, the initial and final bits are fixed 1 and 0, and the interiors are independent. After applying the last zero to uniform measure, \(X=2U\) and \(Z\) are independent uniforms on \(\mathbb Z/M\). Therefore
\[
   L(y)=149+324U+50Z\bmod M
\]
is uniform on the classes compatible with oddness, up to a fixed-\(M\) exponential error, since \(\gcd(324,50)=2\). The statement is over addresses; distinct addresses can give the same scalar count.

One can still deduplicate at a nontrivial scale. For \(B=2^{\lfloor n/2\rfloor}\), split \(y=Bz+u\), with \(0\le u<B\). With \(\mathrm{rev}(u)\) taken in the padded width, unrolling the recurrence gives
\[
 f(y)=f(z)+(B-1-\mathrm{rev}(u))z+f(u)+[u<B/2][z\ {\rm odd}].
\]
For a given suffix \(u\), the coefficient of \(z\) in the height is at least \(162B\), while the variable remainder has total width at most \(25(2^n/B+1)<162B\). Consequently any height occurs for at most one prefix per suffix, and at most \(B/2\) even suffixes. There are \(\Omega_M(2^{n/2})\) distinct lower counts in each compatible phase among the even \(n\)-bit addresses. These counts all lie at order of magnitude \(2^n\). This is still far short of absolute surjectivity. Separately, breadth-first search in the same finite graph finds shortest addresses for every phase in \(O(M^2)\) states. Neither argument proves a uniform \(O(\log M)\)-bit controller.

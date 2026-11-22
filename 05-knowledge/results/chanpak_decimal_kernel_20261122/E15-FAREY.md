# A Farey neighbor selects a physical arc across the fifteenth center

**PROVED.** There are two strong tournaments on 24 vertices with Hamiltonian path counts 536,854,527 and 536,854,529 differing in one arc. Their midpoint is \(E_{15}=2^{14}(2^{15}-1)\), the formal Euclid–Euler center at the composite exponent 15. The same physical change has absolute size two for every choice of the two exterior flank addresses. A Farey determinant specifies exactly which exterior blocks preserve this property for the given middle edge.

## Which two boundaries preserve an arc?

Let a nonempty block \(S\) have tags at the vertices beaten by the apex. Write \(h_S\) for its Hamiltonian path count, \(q_S\) for its two-slot path-cover count with empty slots permitted, and \(c_S,t_S\) for the unsigned two-slot readings starting at a tag and ending at a nontag in the first slot, respectively. Write \(e_S\) for full paths ending at a tag, \(r_S\) for full paths starting at a nontag, and \(k_S\) for the Hamiltonian count of the block with its apex. The unsigned endpoint counts and these two full-path counts are the total and signed-parity quantities in the [shared-apex ordinal law](../chanpak_e31_physical_20261122/README.md).

Take two versions of a tagged nonempty middle block \(D\) differing in one physical arc. Suppose their changes in \((h,q,c,t)\) vanish. For any nonempty tagged blocks \(A,V\), give \(A\to D\to V\) one common apex. Every reading of the joined block except \(k\) is preserved, and
\[
\Delta k_{ADV}=h_Ah_V\,\Delta k_D+h_V e_A\,\Delta r_D+h_A r_V\,\Delta e_D .
\]
Indeed, in an ordinal join \(X\to Y\) the readings \(h,q\) multiply, \(c'=q_Yc_X+h_Xc_Y,\ t'=q_Xt_Y+h_Yt_X,\ e'=h_Xe_Y,\ r'=h_Yr_X\), and \(k'=h_Yk_X+h_X(k_Y-h_Y)+c_Xt_Y+e_X r_Y\). Apply this formula twice. In the exterior-flank formula only the constant \(k\) changes. This proves all-address persistence, including either empty flank.

When \((\Delta k_D,\Delta e_D,\Delta r_D)=(-2,2,2)\), the absolute change is two precisely when
\[
                     \big|h_A r_V-(h_A-e_A)h_V\big|=1.
\]
The endpoint proportions \((h_A-e_A)/h_A\) and \(r_V/h_V\) are then reduced Farey neighbors. Indeed every fraction strictly between two reduced fractions of determinant one has denominator at least the sum of theirs, by expressing its numerator-denominator vector in their unimodular basis. Thus the right endpoint count is selected by a Bézout equation, even when the right receiver has hundreds of full paths. For fixed \(A\), this gives a constructive filter on any ordinary candidate \(V\). Retain \(V\) only when \((h_A-e_A)h_V\equiv\pm1\bmod h_A\), then solve a subset-sum problem on its vector of full paths starting at each vertex. The nontags must select precisely the forced numerator \(r_V=((h_A-e_A)h_V\pm1)/h_A\). The remaining tag choices can change the height at a given address while keeping the unit gap.

## The order-24 certificate

A row integer has bit \(j\) set when its vertex beats vertex \(j\), starting at zero. Use these blocks and apex tag masks, in the indicated order.

| Block | Outgoing rows | Tags |
|---|---|---:|
| \(A\) | \(2,8,3,5\) | 1 |
| \(D\) | \(16,5,9,19,38,15\) | 35 |
| \(V\) | \(310,492,504,497,226,448,385,257,16\) | 451 |

In \(D\), reverse the arc on vertices 1 and 4. Its full profile changes from \((29,320,162,129,16,11,113)\) to \((29,320,162,129,18,13,111)\). The actual fractions are \(3/5\) and \(376/627\), giving \(5\cdot376-3\cdot627=-1\). The lower profile of the whole join is
\[
(h,q,c,t,e,r,k)=(90915,147456000,40180840,50982110,62060,72732,14077887).
\]
The upper profile adds two to \(k\). Choose an empty left flank and the right address \(y=13\), so \(f(y)=4\). The lower count is
\[
14{,}077{,}887+62{,}060+13\cdot40{,}180{,}840+4\cdot90{,}915=536{,}854{,}527.
\]
In the 24-vertex hosts the changed arc joins vertices 5 and 8. The independent count splits the actual graph at its apex, enumerates positions of the last pre-apex and first post-apex occupied blocks, and evaluates subset Hamiltonian-path counts within each block. Every path has one such split and traverses its occupied blocks in order on each side. Run `python3 -B e15_pair.py` for both heights, both strong-connectivity checks and the physical arc check. [Certificate and independent subset recurrence](e15_pair.py)

## What reaches the exponent-61 scale?

**PROVED.** For the displayed tags and unrestricted \(x,y\), the lower count modulo five is \(2+2[x\text{ odd}]\), by the complete profile above. Thus this Farey neighbor excludes \(E_{61}-1\) at every address. More generally every odd exponent \(p\equiv1\bmod4\) has \(E_p-1\equiv0\bmod5\). An odd exponent \(p\equiv3\bmod4\) would require \(x\) even in this chart.

**FINITE-EXACT.** Fix the same three ordinary blocks and the same marked middle arc. Among all \(2^{19}\) ways to tag their vertices, exactly 48 have absolute gap two at every pair of exterior addresses. Only four middle masks preserve both unsigned readings. This preservation is necessary because variation of either address separately detects a change in its linear coefficient, while \(h,q\) are fixed by the underlying blocks. Enumerating left and right masks and applying the joined constant formula gives 48 distinct lower profiles.

For prescribed address parities \(a,b\), set \(x=a+2u,\ y=b+2v\). Since \(f\) is even, every height is congruent to \(k+(t+r)a+(c+e)b+qab\) modulo
\[
               \gcd\bigl(2h,2(t+qb),2(c+qa),4q\bigr).
\]
For the right-only face take \(a=0\) and retain just \(2h,2c\) in the gcd. The left-only face is dual. At \(E_{61}-1\) this necessary test excludes 18 of the 48 profiles at every address. For a right flank alone it excludes 46. The two surviving right-only charts have \(h=90{,}915\) and the following data.

| \((A,D,V)\) masks | \(c\) | Even intercept \(k\) | Odd increment \(e\) |
|---|---:|---:|---:|
| \((12,35,451)\) | 86,260,840 | 30,016,615 | 62,060 |
| \((12,35,211)\) | 86,299,990 | 42,552,735 | 46,400 |

Both need a 95-bit address beginning with \(1100011\), since \(0\le f(y)\le y\) bounds \(y\) between \((E_{61}-1-k-be)/(c+h)\) and \((E_{61}-1-k-be)/c\). Exact solutions on either chart remain **OPEN**.

**PROVED with finite coefficients.** Every exponent-61 solution in these 48 permanent-unit charts has construction order at most 116. With two nonempty flanks of bit lengths \(i,j\), the height is at least \(147456000\cdot 2^{i+j-2}\), giving \(20+i+j\le115\). For the right-only charts the width gives order 115. The 30 left-only survivors have slope at least 49,045,698, giving order at most 116 from \(t x\le E_{61}-1\). These bounds apply to the three fixed ordinary blocks and their normalized binary flanks. Run `python3 -B e15_retags.py` in this directory. [All retaggings, parity congruences and order bounds](e15_retags.py)

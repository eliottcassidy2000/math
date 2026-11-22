# A short complementary controller for a growing modulus

**PROVED.** The fixed-sum controller of COMPLEMENTARY-CONTROLLER can be run with a suffix of common length \(O(\log^2(2D))\) at every \(D=2^aM\), \(M\) odd, instead of length \(O(D^2)\). Its other hypotheses and its physical outputs are unchanged. In particular the phase-matched cube-root brackets work when \(\log D=o(\sqrt{\log T})\), including growing powers of any fixed prime bank. This is not an exact crossing; in particular the separate exhaustive E61 cube-root-window search rules out exact equality in that prescribed window.

## Two observations for complementary words

Write \(F=f/2\). Load an even left address \(x=2d\) when the current sum \(s=x+y\) is odd. For words \(U,V\) of the same padded width, compare appending \(U,\bar U\) to \(x,y\) with appending \(V,\bar V\), where a bar denotes bitwise complement. Write \(R=\operatorname{rev}(V)-\operatorname{rev}(U)\) and \(K(W)=f(W)+f(\bar W)+\operatorname{msb}(W)\). The exact two increments, taken as second minus first, are
\[
 \Delta F(x)=-Rd+\tfrac12(f(V)-f(U)),\qquad
 \Delta(f(x)+f(y))=R(s-4d)+K(V)-K(U).
\]
These identities follow from the binary split formula: the two reversal coefficients have opposite signs, and the cross-boundary parity correction arises only on the odd incoming right word, where its difference is the difference of the leading bits of \(V,U\). If \(U\equiv V\bmod 2D\), both outgoing address residues agree in the two runs. Thus any later fixed word adds the same future increments to the tracked registers modulo their respective factors.

A loader can set \(x=2d\bmod2D\) in a fixed width \(w=\lceil\log_2(2D)\rceil\): append the padded representative of \(2d-2^w x\bmod2D\). Always append its complement on the right. All sums then follow the deterministic rule \(s'=2s+1\), independent of the loaded digits; in particular each switch knows its current \(s\bmod M\).

## Three short switch pairs

For \(M>1\), choose \(u\in[0,7]\) congruent to \(5-M\bmod8\), set \(v=u+M\), and choose the common base width \(b=\operatorname{bitlength}(v+1)\). The first two switches have unshifted options \((u,v),(u+1,v+1)\) in width \(b\). The third has options \((0,M)\) in width \(\operatorname{bitlength}(M)\). In each case shift both options left by \(a+1\) bits, adding this many low zeros and increasing their padded width accordingly. The resulting options differ by \(2D\) and hence have equal outgoing address residues. Every reversal difference \(R\) is odd: \(v,v+1\) have leading bit 1 and \(u,u+1\) have leading bit 0 (when \(v=5\) check \(M=3,5\), otherwise \(v\ge13\)), and this is also clear in the third switch. Passing from \(u\) to \(u+1\) turns on just the low bit; passing from \(v=5\bmod8\) to \(v+1\) turns it off and turns on the next. Thus the difference of the first two reversal slopes is up to sign \(3\cdot2^{b-2}\). In particular \(g=\gcd(M,R_0,R_1,R_2)\mid3\).

When \(g=3\), at least one of the first two constants \(K(V)-K(U)\) is a unit modulo 3. To see this, use \(P=2^{a+1}\), \(W=2^b\), and write the constant read on a shifted option \(Pz\) as
\[
 K(Pz)=f(z)+f(W-1-z)+(P-1)z+[z\text{ odd}]+\operatorname{msb}_b(z).
\]
Put \(u=2i,\ v=4j+1\) with \(j\) odd, \(J=W/2-1-i\), and \(J'=W/4-1-j\) even. For the increment under \(u\to u+1\), the two \(f\)-terms increase by \(J-i+[J\text{ odd}]-[i\text{ odd}]\); for \(v\to v+1\) they increase by \((j+1)+(-J'-2)\). The linear \((P-1)\) increments cancel in the double difference; the parity increments are \(+1,-1\). Neither move crosses the leading-bit threshold \(W/2\). Consequently the second switch constant minus the first is
\[
 2(i+j)-3W/4-1-[J\text{ odd}]+[i\text{ odd}].
\]
Since \(3\mid v-u\), \(i+j\equiv-1\bmod3\); the remaining two indicators differ by 1 because \(J\) and \(i\) have opposite parity. This proves the unit-constant claim.

## Program both precision factors independently

For each switch choose \(d\bmod2^a\) to make \(\Delta F(x)=0\bmod2^a\), possible because \(R\) is odd. Independently, the odd factor can use \(d\bmod M\) arbitrarily. In a group of three switches the total odd change is
\[
 \sum_i s_i R_i+\sum_i (K(V_i)-K(U_i))-4\sum_i R_i d_i \pmod M.
\]
Here the sums \(s_i\) are known before any choice, since every group has fixed widths. If \(g=3\), reverse the two options of one switch with unit constant if necessary to make the total constant nonzero modulo 3; the terms \(s_iR_i\) do not affect it. Then for each \(j<\lceil\log_2M\rceil\) choose a sign so the group can contribute \(\pm2^j\): choose the sign in the required class modulo \(g\), then solve the remaining congruence by Bézout for \(d_i\bmod M\). CRT combines each loading address with its already-selected dyadic congruence. This simultaneously makes every group's change of \(F(x)\bmod2^a\) zero. One either toggles all switches in a group or none. Because the options have equal address residues, the groups can be freely combined; subset sums of signed binary powers cover all odd-modulus register targets. When \(M=1\) omit the bank.

After building these groups, run the inherited short single-address controller on the baseline, appending the complement of its word on the right. Use it to choose the required final \(x\bmod2D\) and \(F(x)=0\bmod2^a\). Every toggle leaves the dyadic conditions unchanged and has its already-computed odd effect after this common tail. Apply the fixed-sum proof unchanged: if \(S+1\) is divisible by \(2^{\ell+a+1}\), its invariant \(C\) fixes the dyadic phase by a permutation quadratic in \(x\bmod2^a\). On the odd factor choose \(x=S\bmod M\), \(y=0\bmod M\), then prescribe the required sum-reading residue (after dividing by \(\gcd(h,M)\), as in the earlier proof). The required progression of sums and the actual marked arc are unchanged.

There are \(O(\log(2D))\) groups and each loader and switch has width \(O(\log(2D))\); the final controller also has quadratic length. Therefore \(\ell=O(\log^2(2D))\). In the geometric bracket proof \(2^{\ell+a+1}=T^{o(1)}\) whenever \(\log D=o(\sqrt{\log T})\). This remains smaller than the chosen \(B\asymp T^{1/6}\), and the smooth step bound and the complementary approximation yield the same \(T^{1/3+o(1)}\) bracket error, at the stipulated phase. There is no inequality asserting that this error is smaller than \(D\).

Replay with `python farey-short-complement-check.py` from this packet. It verifies 208 full literal-chart phase controls and the fixed pairwise sum.

# An obstruction to extending the exponent-15 marked arc

**PROVED.** Preserve the middle blocks \(A\to D\to V\) and the actual marked arc of the [exponent-15 certificate](e15_pair.py). Allow arbitrary tagged tournament blocks in **all four ordinal gaps**, with as many new vertices as desired, and replace the exterior suffix by an arbitrary new tagged tournament. Such extensions cannot straddle \(E_p\) at any \(p\equiv1\pmod4\), including 61. The theorem includes inflating individual vertices of the original transitive suffix, inserting a constant binary run at any of its cuts, and arbitrarily many insertions there.

Write \(h,e,r,k\) for a tagged block's path count, signed start reading, signed end reading and apex path count in the certificate's seven-coordinate join law. For nonempty blocks \(L,R\) surrounding \(D\), the two completions have identical \(h,Q,c,t\) in the middle. The upper-minus-lower differences in its remaining coordinates are \((\Delta e_D,\Delta r_D,\Delta k_D)=(-2,-2,2)\). Substitute twice in the join law. The signed half-gap at the apex is
\[
                    h_Lh_R-e_Lh_R-h_Lr_R.
\]
The same law gives \(h_{B\to C}=h_Bh_C\), \(e_{B\to C}=h_Be_C\), \(r_{B\to C}=h_Cr_B\) for nonempty factors. For the original two flanks one has \((h_A,e_A)=(5,2)\), \((h_V,r_V)=(627,376)\), so this expression equals \(3135-1254-1880=1\).

Take \(L=P\to A\to X\) and \(R=Y\to V\to Q\), permitting empty insertions. If \(X\) is nonempty, both \(h_L\) and \(e_L\) are divisible by 5. If \(Y\) is nonempty, both \(h_R\) and \(r_R\) are divisible by 627. Either condition rules out a signed half-gap of \(+1\) or \(-1\). If both interior insertions are empty, the original half-gap is multiplied by \(H(P)H(Q)\), with the empty block assigned count one. The unit is therefore preserved in its original direction exactly when both exterior blocks are transitive. Their contact words may still be arbitrary.

For those two exterior words use the binary addresses \(x,y\) of the [star construction](../chanpak_gap_star_ports_20261122/README.md). The lower completion has profile
\[
(h,Q,c,t,e,r,k)=(90915,147456000,40180840,50982110,62060,72732,14077887).
\]
Its height is \(k+h(f(x)+f(y))+Qxy+tx+cy+r(x\bmod2)+e(y\bmod2)\). Thus every possible unit lower height in this entire extension family equals \(2+2(x\bmod2)\pmod5\). At \(p\equiv1\pmod4\), the target \(2^{p-1}(2^p-1)-1\) is zero modulo 5. For odd exponents \(p\equiv3\pmod4\), it is two modulo 5, so the left address would have to be even.

For a direct single-run check the ordered tape of the certificate is \(A,D,V,0,1,0,0\). Count cuts starting at zero. The first internal cut (after \(A\)) gives signed half-gaps 1,255 and −1,880 for a nonempty run of zeros or ones, respectively. The second (after \(D\)) gives −1,254 and 1,881. The other six cuts retain the unit. The half-gap is independent of the positive run length. Run `python3 -B e15_insertions.py` for the exact chart and all sixteen cases. [Replay](e15_insertions.py)

# Integer-only Farey corridors at spine four

**Finite-exact.** At \(N=2^{60}(2^{61}-1)-1\), the canonical spine-four chart has no solution in the following specified collection of rational-tangent corridors. This is not a full spine-four exclusion, a contiguous interval of excluded addresses, or a claim of coverage disjoint from every earlier screen.

The chart parameters \((h,q,c,t,e,r,k)\) are
\((2465,2293760,1344270,712944,1740,986,423385)\).
The left congruence permits only physical addresses \(x=10X+b_0\), \(b_0=0,1\); write \(y=Y\), and set
\[
 Q=10q,\quad C=c+qb_0,\quad T=10t,\quad
 K=k+tb_0+r[b_0=1],\quad M=N-K,\quad \Delta=QM+TC.
\]
Then \(\Phi(X)=(M-TX)/(QX+C)\) is the upper envelope for \(Y\), with \(\Phi'(X)=-\Delta/(QX+C)^2\) and \(\Phi''(X)>0\). On a corridor \(X\in[A,B]\), \(0\le f(u)\le u\) gives \(\Phi(X)-G\le Y\le\Phi(X)\) for every solution, where
\[
 G=\frac{h(10B+1+\lfloor\Phi(A)\rfloor)+e}{QA+C}.
\]
Reject immediately if \(\Phi(A)<0\).

For a reduced positive slope \(a/b\) and anchor \(X_0\), put \(F(X)=b\Phi(X)+aX\), with derivative \(D_0=a-b\Delta/(QX_0+C)^2\). Every solution lies on one of the integer lines \(z=bY+aX\) satisfying
\[
 F(X_0)+\min\{D_0(A-X_0),D_0(B-X_0)\}-bG
 \ \le z \le \max\{F(A),F(B)\}.
\]
The lower bound is the supporting tangent and the upper bound is the chord bound from convexity. They remain valid when clipping leaves a corridor not containing its original anchor. There is no floating-point rounding.

The replay takes \(L=115000000000000\) and \(R=200000000000000\). Enumerate all reduced slopes of denominator \(1\le b\le128\) between \(\Delta/(QR+C)^2\) and \(\Delta/(QL+C)^2\). Set
\[
 X_0=\left\lfloor\frac{\left\lfloor\sqrt{b\Delta/a}\right\rfloor-C}{Q}\right\rfloor,
 \qquad
 w=\min\left(\left\lfloor X_0/8\right\rfloor,
 \left\lfloor\sqrt{\frac{(QX_0+C)^3}{128b\Delta Q}}\right\rfloor\right).
\]
Keep anchors in \([L,R]\) with \(w>0\). Sort windows \([X_0-w,X_0+w]\cap[L,R]\) and clip overlaps before testing.

On any surviving line \(z\), solve \(X\equiv a^{-1}z\pmod b\); this is vacuous when \(b=1\). Parameterize its integral points by \(X=X_1+bj,\ Y=Y_1-aj\) and require \(Y\ge0\). Split \(j\) by parity so the right parity is constant. For \(j=2u+d\), the remaining height congruence modulo 17 and modulo 29 is quadratic in \(u\), because both primes divide \(h\). Enumerate roots modulo each prime, combine by CRT, and check the full heights on those progressions. This handles linear and identically-zero congruences without a division hypothesis.

| Left residue \(b_0\) | Disjoint corridors | Normalized \(X\) covered | Corridors without an integer line | Total lines | Literal pair readings | Hits |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 29,457 | 3,204,528,847 | 24,990 | 4,467 | 8,234 | 0 |
| 1 | 29,457 | 3,204,528,847 | 24,952 | 4,505 | 7,817 | 0 |

Run standard-library Python file farey-corridors.py with default settings. Before the target run, it compares several complete exact slanted strips against Cartesian brute force on a planted pair.


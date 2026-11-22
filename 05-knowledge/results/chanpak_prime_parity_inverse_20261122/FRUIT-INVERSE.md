# Dyadic lifting of the fruit cubic and exact parity edges of tournament words

**Status: PROVED** for the local classifier, exact reconstruction algorithm, and physical word bounds below. The equivalence of the fruit's local and rational isogeny images uses the repository's cited rank-and-generator result. No elliptic conjecture, RH, or identification of a fruit coordinate with a Hamilton count is used.

## 1. The missing explicit two-adic guard

Use the notation of the [fruit quotient and full image](../session97_overnight_arithmetic_carriers_20261029/fractions/FULL_IMAGE.md). For a normalized triple of sum one, set $q=ab+ac+bc$. The inverse cubic is
$
 f_q(t)=t^3-t^2+qt-(6q-1)/7,
 \qquad 49\operatorname{disc}(f_q)=(q+1)(1+29q-196q^2).
$
A rational target on the orientation cover supplies a rational square root of this discriminant. The repeated-root value $q=-1$ is the excluded rational boundary. The previous repository theorem proves that having a lift over $\mathbb Q_2$ is equivalent to having one over $\mathbb Q$ for these rational quotient points, but its proof explicitly does not turn arbitrary $q$ into a simple dyadic test.

**Dyadic classifier.** For arbitrary $q\in\mathbb Q_2$, the cubic has a $\mathbb Q_2$ root exactly when
$
 q\notin\mathbb Z_2 \quad\text{or}\quad q+1\in16\mathbb Z_2.
\tag{1}
$
In the square-discriminant nonboundary case a root implies full splitting over $\mathbb Q_2$. On a rational point of the orientation cover, by the prior image theorem this is an exact test for rational reconstruction.

- If $v_2(q)<0$, divide $f_q$ by $q$. Its coefficients are integral and its reduction is $t$, so there is a simple, integral Hensel root.
- If $q$ is integral and even, any local root must be integral by monicity, and $f_q\bmod2=t^3+t^2+1$ has no root.
- If $q$ is integral and odd, its reduction is $(t+1)^3$, so any root must have $t=1+2u$ with $u\in\mathbb Z_2$. With $e=q+1$ the exact translated equation is
  $
    f_q(1+2u)=8(u^3+u^2)+2eu+e/7.
  $
  If $v_2(e)=1,2$, the last term has uniquely lowest valuation. If $v_2(e)=3$, division by 8 again gives the rootless polynomial $u^3+u^2+1\bmod2$. If $v_2(e)\ge4$, including $e=0$, division by 8 gives $u^2(u+1)\bmod2$, with a simple root at 1.
- For a simple cubic root $r$, the remaining quadratic's discriminant differs from the cubic discriminant by the square of the quadratic evaluated at $r$; hence a nonzero square cubic discriminant makes any root a full lift.

When $q$ is integral and odd and the discriminant is a square, $v_2(q+1)=1,2,4$ in fact cannot occur. To see this, write $e=2^h u$ with $u$ odd, use
$
 49\operatorname{disc}(f_q)=e(-224+421e-196e^2),
$
and divide by $2^{2h}$ for $h<5$. For $h=1,2$ the remaining odd unit is 5 modulo 8, and for $h=4$ it is $5-6u\bmod8$, never 1. For $h=3$ the remaining unit is always 1 modulo 8, so the repeated root in the first reduction is a genuine hostile phase unless the full guard is checked.

For reduced rational $q=a/d$ with $d>0$, (1) is simply $2\mid d$ or $16\mid a+d$. Thus on nonboundary rational quotient points it agrees with the previous ternary guard $3\mid d$ or $3\mid a+d$. These are not conclusions for arbitrary rational numbers without the orientation-cover equation. In particular, a valid primitive integer fruit triple cannot have exactly one odd coordinate: its sum would then be odd and its $q$ even. Two odd coordinates give nonintegral $q$; if all three coordinates are odd, the normalized $q+1$ has valuation at least 4.

## 2. Recovering the exact rational triple without a Mordell–Weil index

The orientation-cover point can also be inverted in deterministic bit-polynomial time, not just accepted by a local residue. This does not compute its discrete index on the elliptic curve.

For $q=a/d$ reduced, the integral polynomial $7d f_q$ gives denominator bound $B=7d$ by the rational root theorem. The Cauchy root bound gives numerator bound $M=14(d+|a|)$ for any reduced rational root. Choose $P=2^\ell>2MB$. In the nonintegral passing phase, lift the simple root of $f_q/q\bmod2$ by Hensel to mod $P$. In the integral passing phase, lift the simple root of $f_q(1+2u)/8\bmod2$, and return $t=1+2u\bmod P$. Denote either returned residue by $s$. Each lifted polynomial has 2-integral coefficients with odd denominators after simplification, so each lift step is an exact integer parity calculation.

Any rational root in the chosen simple branch is 2-integral and has odd reduced denominator $b$. For that root $p/b$, write $p=sb-Pj$. Since $|p|\le M$ and $b\le B$,
$
 \left|\frac{s}{P}-\frac{j}{b}\right|
 =\frac{|p|}{Pb}<\frac{1}{2b^2}.
$
The fraction $j/b$ is reduced by $\gcd(p,b)=1$, so it is a convergent of $s/P$ by the continued-fraction approximation criterion. Enumerate the convergents and check the recovered rational candidates in the exact cubic. At most one candidate with the height bounds gives the residue, since two would have a nonzero cross product of size at most $2MB<P$ divisible by $P$.

If any rational split exists, all three roots are rational and the chosen simple local branch is one of them; the algorithm finds it. This reasoning gives an elementary exact inverse on the accepted cases independently of the global isogeny theorem, if one allows the exact check to reject. That theorem guarantees success on the accepted nonboundary rational quotient points. Divide out the first root, take the exact rational square root of the residual quadratic discriminant, clear denominators of the three roots, divide their gcd, and fix the positive-sum sign. The supplied orientation chooses three cyclic labels; omitting it permits all six distinct labels. Testing positivity is then exact and independent of local acceptance. The bit length of $P$ is linear in the input height and all modular lifts, Euclidean convergents and exact square-root checks are polynomial in that bit length.

This gives an explicit local-to-global reconstruction with a height budget: finitely many dyadic digits determine one rational witness uniquely, but only **after** an Archimedean root bound and an exact replay. It does not say that an unrestricted tournament residue selects an absolute count.

## 3. A literal parity-sensitive extremal estimate for physical tournament rows

Use the [physical one-tape row law](../chanpak_local_rows_span3_20261115/ROW_LIMIT.md), which concerns
$
 G(n)=c n+k p(n)+\phi[n\text{ odd}]+b_0,
 \qquad p(0)=0,\quad p(2v+1)=p(v),\quad
 p(2v)=p(v)+v+[v\text{ odd}]
$
with positive $c,k$. Our estimate is at the upper endpoint on the $\sqrt N$ scale, not a normalized presentation density. It applies in particular to the inherited mask-zero physical row $c=189,k=11,\phi=11,b_0=11$.

**Exact sorting lemma.** Fix the last bit and the number of ones in the other $L-1$ bits. Among $L$-bit words with leading zeros allowed, the unique maximizer has its upper bits sorted as ones then zeros.

To prove this, replace an adjacent $01$ by $10$ strictly above the fixed last bit. Before suffixing, the prefixes are $4v+1$ and $4v+2$, so the increase in $p$ is $v+2-[v\text{ odd}]$. If the common suffix has length $j\ge1$, appending it adds its nonnegative reverse-complement to this difference, and removes at most one parity correction from the first suffix position. Thus $p$ does not decrease. The integer word $n$ increases strictly, and the fixed terminal bit keeps the $\phi$ contribution fixed. Repeating removes every inversion; any unsorted word was strictly worse.

Let $N=2^L$. The even sorted candidates are $n=N-2^j$, $1\le j\le L-1$, along with $n=0$ if permitted. For these positive candidates,
$
  G(N-2^j)=(c+k)N-(c+k)2^j-k2^{L-j}+2k+b_0.
\tag{2}
$
Indeed the prefix $2^{L-j}-1$ is all ones and has $p=0$, and appending $j$ zeros adds $(2^j-1)(2^{L-j}-1)+1$.

The odd sorted candidates are $n=N-2^j+1$, $1\le j\le L$; for $2\le j\le L-1$,
$
 G(N-2^j+1)=(c+k/2)N-(c+k/2)2^j-k2^{L-j}+
                  2k+c+\phi+b_0.
\tag{3}
$
This follows by stripping the last $1$ and using the even computation in the shorter word. The two boundary candidates are $n=N-1,1$, both with $p=0$. Thus both parity maxima are **exact finite minimizations of two positive monomials**, with only the listed boundary checks.

For either parity let $A=c+k$ or $A=c+k/2$. For large $L$ the maximum is interior: the optimal cost $A2^j+k2^{L-j}$ is $O(\sqrt N)$, whereas the odd boundary $n=N-1$ loses $kN/2+O(1)$ against its leading edge. Write $j_*=L/2+\tfrac12\log_2(k/A)$, and choose a nearest integer $j$ to $j_*$. The cost is exactly
$
 2\sqrt{AkN}\cosh\!\bigl((j-j_*)\log 2\bigr).
$
Hence the top deficit is of order $\sqrt N$ with an explicit oscillation depending on **word-length parity**, in addition to the different leading edges for **input parity**. Only two adjacent $j$ need be tested. These are sharp missing-height bounds for this row and bit width, not missing counts for the full tournament spectrum.

For the mask-zero row and $L\ge8$, put $D=2^{\lfloor L/2\rfloor}$. Both maxima occur at the respective adjacent addresses $N-2^{\lfloor L/2\rfloor-2}$ and $N-2^{\lfloor L/2\rfloor-2}+1$. The exact maxima are:
 
| word length | even input | odd input |
|---|---:|---:|
| $L$ even | $200N-94D+33$ | $(389/2)N-(741/8)D+233$ |
| $L$ odd | $200N-138D+33$ | $(389/2)N-(1093/8)D+233$ |

These choices can be certified without logarithms: the cost difference between $j+1$ and $j$ is $A2^j-k2^{L-j-1}$, whose sign changes once. At $L\ge8$ the resulting odd interior count also exceeds both boundary counts (at $L=8$ the comparison with $n=N-1$ is already positive; for larger widths the linear advantage dominates). In particular any value in the upper interval above the odd maximum must have an even address in this physical row. The condition is nontrivial because the two maxima have different leading slopes, despite every resulting Hamilton count being odd.

The transferable point from the fruit inverse is the separation of a local parity label, a height bound, and a final exact witness. Here a terminal bit really selects a different maximal height, but no modular fact can manufacture the word beyond its sharply bounded height. The fruit equation is used as arithmetic structure and not as a numerical analogy about exceptional tournament counts.

Exact finite checks: `python3 check_fruit_inverse.py`.

### The actual pair-sum valuations in the integral fruit phase

Let $e=q+1\ne0$, $h=v_2(e)\ge4$ in a splitting integral branch. In the normalized $u$ equation, exactly one root is odd, because its reduction is the unique simple residue 1. Each of the two even roots satisfies $2v_2(u)=h-3$: the constant term has valuation $h-3$, the term $(e/4)u$ has strictly greater valuation, and $u^2(1+u)$ has valuation $2v_2(u)$. Thus $h$ is odd and at least 5. In terms of the three normalized roots $t_i=1+2u_i$,
$
 \{v_2(1-t_i):i=1,2,3\}=\{1,(h-1)/2,(h-1)/2\}.
$
For a primitive integer lift these are exactly the pair-sum valuations: all individual coordinates and the sum are odd, since the three normalized roots are 2-adic units and are reduced with odd denominators. This is a more precise local compatibility condition than calling such a triple merely “all odd.”

### Exact scalar decoder and presentation count at the tournament edge

Fix input parity $\epsilon\in\{0,1\}$. Write $A_\epsilon=c+k(1-\epsilon/2)$, $R(n)=\operatorname{rev}_L(n)$, and $R_0=R(n)-\epsilon N/2$. Directly summing the zero-bit recurrence gives
$
 \big|p(n)-n(1-R(n)/N)\big|\le L+1.
$
Indeed replacing each floor by its argument in the zero-bit sum changes it by at most one per bit, including its adjacent-bit correction; the unfloored coefficient is $1-1/N-R(n)/N$. Thus, with $E=k(L+1)+|\phi|+|b_0|$,
$
 \left| A_\epsilon N-G(n)-
   \left(A_\epsilon(N-n)+k(n/N)R_0\right)\right|\le E.
\tag{4}
$

Given an exact target $Y$ with $|Y-A_\epsilon N|\le K\sqrt N$, put $T=A_\epsilon N-Y+E$. Reject if $T<0$. Any preimage has $N-n\le T/A_\epsilon$. When $T\le A_\epsilon N/2$, we also have $n\ge N/2$ and hence $R_0\le2T/k$. These integer bounds determine both ends of the word:

1. Choose $j$ minimally with $2^j\ge\lfloor T/A_\epsilon\rfloor$; a zero deficit bound rejects because $n<N$. Then the top $L-j$ bits of $n$ are ones.
2. Choose $i$ minimally with $2^i>\lfloor2T/k\rfloor$. Then bits 1 through $L-i-1$ of $n$, counted from the least significant end, are zero; bit 0 is fixed to $\epsilon$.

Inconsistent forced bits reject. Otherwise the number of free bits is at most $\max(0,i+j-L)$, so the exact target can be tested by at most $O_{c,k,\phi,b_0}((1+K)^2)$ words for large $L$ and $K\le \eta\sqrt N$, where $\eta>0$ depends only on the row. Enumerating those middle bits and checking the original recurrence returns **all** actual preimages at the specified width and parity. In particular it is a bit-polynomial inverse for a fixed $K$; there is no factoring or modular-residue promise.

This count is of the right order as a presentation bound. For a sufficiently large row-dependent constant $K_0$, and $K_0\le K\le\eta\sqrt N$, choose powers of two $P,Q$ within a factor two below $K\sqrt N/(8A_\epsilon)$ and $K\sqrt N/(8k)$, respectively. For large $L$, $P,Q\le N/2$ and $PQ\ge N$. The two simultaneous conditions $N-P\le n<N$ and $0\le R_0<Q$ fix respectively a high all-one cylinder and a low zero cylinder with terminal bit $\epsilon$. Their bit sets do not overlap and they contain exactly $PQ/N=\Omega(K^2)$ words. By (4) these lie within $K\sqrt N$ of the leading edge; the exact max theorem places all counts below it for large $L$. The number of presentations in this edge window is therefore $\Theta(K^2)$, even on the $\sqrt N$ scale on which a macroscopic local limit has no such accuracy. Distinct values need not equal the number of presentations; the decoder reports any collisions in the word list rather than suppressing them.

### The edge collision problem stabilizes to finitely many affine lines

For fixed $K$, the forced-end proof allows a fixed width $w\ge1$, depending on $K$ and the row only, such that all edge candidate words have the form
$
 L=2m+\lambda,\quad \lambda\in\{0,1\},\quad
 Q=2^{m-w},\quad B=2^{2w},\quad P=2^\lambda Q-1,\qquad
 n=(PB+b)Q+\epsilon,\quad 0\le b<B.
$
Here $P$ is the high all-one prefix and $b$ is the entire padded middle; all bits below it except the terminal one are zero. Take $m$ large enough that the high prefix is nonempty and that $m-w\ge2$. For $V=\operatorname{rev}_{2w}(b)$, the exact two suffix recurrences give
$
\begin{aligned}
G(n)={}&A_\epsilon N-
 Q\bigl(A_\epsilon(B-b)+k\,2^\lambda(V+1)\bigr)\\
 &+k\bigl(V+1-b+p(b)+[b<B/2]+[b\text{ odd}]\bigr)
 +(c+\phi)\epsilon+b_0. \tag{5}
\end{aligned}
$
Indeed $p(P)=0$, so $p(PB+b)=p(b)+(B-1-V)P+[b<B/2]$. Appending the low zero run, and stripping the final one first if $\epsilon=1$, adds $(Q2^{-\epsilon}-1)(PB+b)+[b\text{ odd}]$. Inserting $P=2^\lambda Q-1$ and $Q(PB+b)=N-Q(B-b)$ yields (5).

For fixed parities and $K$, the entire edge is thus contained in an explicit finite family of affine deficit lines of slope and intercept independent of $m$. Two different slopes eventually separate; two equal slopes are compared by their exact constant offsets. This proves eventual stabilization of all collision classes and edge membership, including ties at the edge-window boundary: when the slope times $Q$ equals $K\sqrt N$ identically for the fixed $\lambda$, the constant offset decides membership. The result is a finite exact description at the parity-specific $\sqrt N$ edge, not an estimate in the full integer range.

### Uniform local law for every even right-hand side

The dyadic analysis does not depend on the constant $4$. For any even integer $n$, replace the right-hand side of the three-fraction equation by $n$. If $s=a+b+c$, $q_0=ab+ac+bc$, $r_0=abc$, the numerator and denominator after summing the fractions are $s^3-2s q_0+3r_0$ and $s q_0-r_0$. In sum-one normalization this gives the inverse cubic
$
f_{n,q}(t)=t^3-t^2+qt-\frac{(n+2)q-1}{n+3}.
$
Since $n+3$ is odd and $n+2$ even, $f_{n,q}/q$ again reduces to $t$ when $q$ is nonintegral, and the integral reduction and all valuation cases in (1) are identical: the translated constant is now $(q+1)/(n+3)$ instead of $(q+1)/7$. Thus the dyadic root criterion, and the pair-sum valuation profile in the integral split case, hold for all these even right-hand sides. This does not transport the global $\mathbb Q_2$-to-$\mathbb Q$ isogeny-image theorem beyond the stated $n=4$ curve.

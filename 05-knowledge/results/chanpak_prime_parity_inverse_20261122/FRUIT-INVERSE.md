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

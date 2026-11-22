# Charging forbidden pairs: an exact sibling-height obstruction

**PROVED** for the stated positive alphabet. The linked prime gap of 180 is **CONDITIONAL**, and no prime-gap estimate is used below.

The [Allahyar–Yianni preprint](https://github.com/sahibzada-allahyar/math/blob/main/preprint/prime-gaps-180/prime_gaps_at_most_180.pdf), §2, keeps a 39-shift tuple of diameter 182 and charges the only pair farther than 180, its two endpoints. If $N$ is the number of prime hits and $B$ is the indicator of this marked endpoint pair, a positive $N-1-B$ forces a permissible pair. This is an exact detection interface. The analytic estimates through Theorem 1.2 remain conditional, in particular on their unestablished Type III Hypothesis A4; the Lean theorem checks the implication from its premises, not A4. The analogue below measures the **actual pair** instead of inferring short gaps from marginal distributions.

## 1. Sibling-height gap and fruit cap

Use the preceding affine sieve's positive alphabet
$
A=\begin{pmatrix}3&3\\6&8\end{pmatrix},\quad
B=\begin{pmatrix}9&9\\30&32\end{pmatrix},\quad
C=\begin{pmatrix}3&3\\12&16\end{pmatrix}
$
on marked states $(h,R)$, seed $(1,2)$, and finish $h+2R$. The three same-source next-letter plus finish outputs are literal ordinary counts:
$
N_A=15h+19R,\qquad N_B=69h+73R,\qquad N_C=27h+35R.
$
Put $t=R/h$. The three projective maps preserve $I=[2,16/3]$, so $t\ge2$. We have $N_A<N_C<N_B$, and the adjacent differences $12h+16R$ and $42h+38R$ are at least $44h$ and $118h$, respectively. This ruled-out pair is very different from the same-source constant-difference bilinear switches: applying a finish to each possible single next letter does not implement a bounded gap along increasing history lengths.

There is also an exact obstruction to these three literal counts solving the fruit equation. The normalized shares $x_i=N_i/(N_A+N_B+N_C)$ have denominator $111+127t$. Their derivatives have numerators $204,-660,456$, respectively, over $(111+127t)^2$. Since $x_B$ is the largest share and $f(x)=x/(1-x)$ has strictly increasing derivative,
$
\frac{d}{dt}\sum_i f(x_i)
=\frac{204(f'(x_A)-f'(x_B))+456(f'(x_C)-f'(x_B))}
 {(111+127t)^2}<0.
$
Thus the maximum is attained at the seed $t=2$, and equals
$
\frac{53}{312}+\frac{215}{150}+\frac{97}{268}
=\frac{205397}{104520}<4.
$
This is a physical one-source impossibility, not a limitation on unrelated triples of tournaments.

More generally, for three positive reals of spread $r=\max/\min$, scale them to $1,s,r$. As a function of $s\in[1,r]$, the fruit sum is convex, and its values at $s=1$ and $s=r$ differ by $(r-1)^3/(2r(r+1))$. Its maximum is therefore $r/2+2/(r+1)$. Solving for when that upper bound reaches 4 gives the necessary cutoff $r\ge(7+\sqrt{65})/2$. This positive-real cap does not remove the extra rational-root and chamber guards of the [prior fruit symmetry quotient](../session97_overnight_arithmetic_carriers_20261029/fractions/README.md).

## 2. Short differences of any two fixed linear finishes are exponentially rare

Let two distinct fixed linear finishes have difference $\alpha h+\beta R$. Fix a bound $K>0$. For $\beta=0$, no history of length $n$ with $9^n>K/|\alpha|$ gives a difference within $K$, since every letter multiplies $h$ by at least $9$. For $\beta\ne0$, each eligible length-$n$ history has ratio $t$ in the fixed interval centered at $-\alpha/\beta$ of radius $K/(|\beta|9^n)$.

The three ratio images of $I$ are separated by at least $\delta=4/9$: $A(I)\subset[2,8/3]$, $B(I)\subset[10/3,32/9]$, and $C(I)\subset[4,16/3]$. Each projective map has derivative at least $d=2/361$ on $I$, by $(a+bt)^{-2}\det$ for its matrix; the minimum lower bound comes from $B$. If two histories have distinct suffixes of length $k$, read those suffixes newest-to-oldest up to their first distinct letter. The outputs just after that letter are separated by $\delta$ and their subsequent $j<k$ common letters contract this distance by no more than a factor $d^j$. Thus their final distance is at least $\delta d^{k-1}$. For
$
k=\left\lfloor \frac{n\log 9}{\log(361/2)}-L\right\rfloor
$
and a sufficiently large fixed $L=L(K,\beta)$, the eligible window is strictly narrower than this separation. At most one length-$k$ suffix can occur for eligible histories. Consequently there are at most $3^{n-k}$ histories at depth $n$ with short difference, an exponentially small proportion. This is a history count, not a theorem ruling out infinitely many exceptional histories. It is also limited to two fixed linear finishes on the **same full state**: it does not address separate growth histories or the bilinear exact-$+2$ switches.

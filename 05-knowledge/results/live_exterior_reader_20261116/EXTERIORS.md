# A physical parabolic crossing with all exterior phase choices live

**Theorem.** Fix the literal synchronized reader with parameters
$(I,Z,g)=(-11,82,9)$. For every sufficiently large integer $N$, let
$n=\lfloor\log_2(N/88)\rfloor$, $p=\lceil n/2\rceil$,
$q=\lfloor n/2\rfloor$ and $S=2^{\lfloor n/10\rfloor}$.
For **each** of the $S^2$ pairs $0\le x_0,y_0<S$, there are
physical addresses strictly below and above $N$, each within $O(N^{2/5})$.
They satisfy $b=x2^q+\operatorname{rev}_q(y)$,
$x\equiv x_0\bmod S$, $y\equiv y_0\bmod S$.
The estimates are uniform in these phases. In particular this leaves
$\Omega(N^{1/5})$ distinct addresses on each side. The synchronized
seed swaps still give exactly +2 and +4 at every address. Neither exact
placement nor a comparable lower bound on distinct **counts** follows.

## Split and phase compression

Write $P=2^p$, $Q=2^q$, $D=PQ$, and abbreviate bit reversal of the
specified width as $\bar x,\bar y$. The quadratic kernel is
$J(w)=\sum_{i<j} (2^{j-i-1}+[j=i+1])w_iw_j$.
It is invariant under reversal at any width and $0\le J(w)\le w$
(sum at each upper bit $j$). The cross term of the split address is
$xy+(x\bmod2)(y\bmod2)$. In fact the exponential of the distance
factors as $2^t2^{q-1-i}$ for lower position $i$ and upper position
$q+t$, and only the two seam bits contribute the adjacent surcharge.
Thus the exact reader count is

$$
H=-11+(82P+9x)(2Q-y)-82\bar x-9J(x)
       +18\bar y-9J(y)-9(x\bmod2)(y\bmod2).
$$

Put $k=\log_2 S$, $R=P/S$, $T=Q/S$, with $k\ge1$.
Let $x=x_0+Su,\ y=y_0+Sv$, and $X_0=\operatorname{rev}_k(x_0)$,
$Y_0=\operatorname{rev}_k(y_0)$. The same split gives
$J(x)=J(x_0)+J(u)+uX_0+\eta_x$, with
$\eta_x=(u\bmod2)[x_0\ge S/2]$, and similarly $J(y)$.
Also $\bar x=RX_0+\operatorname{rev}_{p-k}(u)$ and
$\bar y=TY_0+\operatorname{rev}_{q-k}(v)$.
Set $a=82P+9x_0,\ b_0=2Q-y_0$, and define

$$
A=(a+9Y_0/S)/S+9u,\qquad B=(b_0-X_0/S)/S-v.
$$

Then $H=S^2AB+K+E(u,v)$ where

$$
\begin{aligned}
K={}&-11-82RX_0-9J(x_0)+18TY_0-9J(y_0)
         -9(x_0\bmod2)(y_0\bmod2)\\
   &+aX_0/S-9Y_0b_0/S+9Y_0X_0/S^2,\\
E={}&-82\operatorname{rev}_{p-k}(u)-9J(u)
       +18\operatorname{rev}_{q-k}(v)-9J(v)-9(\eta_x+\eta_y).
\end{aligned}
$$

In fact the large phase terms cancel. Define
$\mathcal E_k(w)=w\operatorname{rev}_k(w)/S-J(w)$. Then

$$
K=-11+9\left(\mathcal E_k(x_0)+\mathcal E_k(y_0)
        +Y_0X_0/S^2-(x_0\bmod2)(y_0\bmod2)\right).
$$

The elementary bound $|\mathcal E_k(w)|=O(k)$ follows by expanding
the reversed product: the terms above the diagonal cancel the
exponential part of $J$, while the diagonal and below-diagonal
geometric sums and the adjacency surcharge each have size $O(k)$.
Thus $K=O(1+\log S)$ while $|E|\le M=91R+27T+18$, uniformly. The same
phase compression reproduces the product plus a bounded-bit correction
after a further split.

## Physical walk across the prescribed target

Let $w=(N-K+M+1)/S^2$, $z=w/(RT)$. Then
$z=N/D+O(D^{-1/2})\in[87.9,176.1]$ for large enough $N$.
Choose $\alpha_0=85$ for $z\le165$, and $90$ otherwise.
Let $\ell$ be the nearest integer to $\alpha_0^2/z$ (either
choice at a tie). Put $\alpha=\sqrt{\ell z},\ \beta=z/\alpha$.
These are uniformly interior to $[82,91]\times[1,2]$.
Indeed the square-root displacement from $\alpha_0$ is at most
$88.05/(84+85)<0.53$; hence the two charts respectively keep
$\beta$ in $[87.9/85.53,165/84.47]$ and
$[165/90.53,176.1/89.47]$. The allowed grids are these scaled
intervals up to $O(1)$ shifts/endpoints.
Take $(A_*,B_*)=(R\alpha,T\beta)$. Round both **up** to their allowed
grids, yielding $A_1,B_1$. The rounding is less than 9 and 1 respectively,
the excess in their product is $O(R)$, and the actual count is
at least $N+1$. With $c=R/T\in\{1,2\}$, the ratio
$A_*/B_*=c\ell$. Take positive integer steps $d=9c\ell$, $f=9$;
for integral $t$ set $A_t=A_1+dt,\ B_t=B_1-ft$.
The grids stay legal for all $t=O(\sqrt R)$, since the starting point
has a margin $\Omega(R)$ and the steps are bounded constants.
The product is exactly

$$
A_t B_t=A_1B_1+(dB_1-fA_1)t-dft^2,
\qquad |dB_1-fA_1|=O(1).
$$

For sufficiently large constant $C$, at $t=\lceil C\sqrt R\rceil$
the concave term dominates the $O(R)$ initial excess and bounded
residual, so actual $H<N$. On this entire walk the maximum one-step
jump is $O(S^2\sqrt R+M)$. To get strict witnesses on both sides
even when equality occurs, perform two bisections, maintaining
boundary points at which the predicate is false and true.
For predicate $H\le N$ return the last point in the bisected pair
with $H>N$. For predicate $H<N$ return the point in its pair with
$H<N$. Each pair is adjacent, so the step bound proves the desired
error. This bisection does not need a monotone residual or select a
globally first crossing.

Since $S=\Theta(D^{1/10})$ and $R=\Theta(D^{2/5})$, the bound is
$O(D^{2/5})=O(N^{2/5})$. Concavity offsets the large individual linear
changes between the two actual exterior coordinates; the uncancelled
variation of bit reversal is limited by fixing the low bits in each
phase rather than selecting one global exterior pair. The square-root
walk is inside the lattice of legal addresses, not a signed linear
combination of counts.

## A simultaneous mod-9 count constraint

For each phase pair and each $r\bmod9$ one can require both strict
witnesses to have $H\equiv r\bmod9$, with the same asymptotic error.
The cutoff and implicit constants here are effective but large.

Use $\alpha_0=85,89.5$ in the two charts and choose the nearest positive
$\ell$ not divisible by 3. Its displacement from $\alpha_0^2/z$ is at
most 1. The displacement of $\alpha$ from $\alpha_0$ is less than
$1.05$ (the crude denominator bounds are $83+85$ for $z\le165$,
$88+89.5$ otherwise), so the rectangles are still uniformly interior:
in particular their maximum $\beta$ is bounded by
$\max(165/83.95,176.1/88.45)<2$. Along the walk let $g=c\ell$;
then $u_t=u_1+gt$ and $v_t=v_1+9t$. In the exact split formula the
only variable term modulo 9 is now $-82\operatorname{rev}_{p-k}(u_t)$:
$y_t$ is constant modulo 9, and all other terms modulo 9 not involving
the reversal are constant. Thus imposing $H\equiv r$ amounts to
imposing one residue of $\operatorname{rev}_{p-k}(u_t)$.

To see bounded gaps in legal times, factor $g=2^a g_0$ with $g_0$
odd and prime to 9. In every aligned interval of $u$ of length
$2^L$ inside the word width, the high prefix is fixed. Set its low
$a$ bits to $u_1\bmod2^a$. In the low block reserve at least $9g_0-1$
positions spaced by $\operatorname{lcm}(6,\operatorname{ord}_{g_0}(2))$
starting at position $a$, using order 1 if $g_0=1$. Take $L$ sufficiently
large to fit (and fix all other low bits to zero). Each reserved bit has
the same invertible weight modulo $g_0$ for $u$ and the same invertible
weight modulo 9 for $\operatorname{rev}_{p-k}(u)$. The number of set
reserved bits can be any integer from 0 to $9g_0-1$; the Chinese
remainder theorem selects it to meet both residue conditions. Since
only finitely many $g$ occur, the acceptable times have a uniform
bounded-gap constant $W$, away from the grid endpoints.

Choose $w=(N-K+M+1)/S^2+C_0$ instead, with constant $C_0$ exceeding
the worst possible product loss from $A_1B_1$ at all $0\le t\le W$.
This constant exists by the quadratic formula and bounded rounding.
The target chart still has $z=N/D+O(D^{-1/2})$; there is an acceptable
time at or before $W$ with $H>N$. There is another acceptable time
within $W$ of $\lceil C\sqrt R\rceil$ with $H<N$, by choosing $C$
sufficiently large as before. Restrict the crossing argument to
acceptable times: consecutive ones differ by at most $W$, so their
actual count jump is still $O(S^2\sqrt R+M)$. This proves the refinement.


### Concrete modular phase bound

The elementary reserved-bit proof establishes existence; the following
much smaller finite computation makes the modular construction practical.
For all odd $d\le83$ prime to 3, the 14-bit pair
$(w\bmod d,\operatorname{rev}_{14}(w)\bmod9)$ covers all $9d$ pairs.
This is checked by starting at $\{(0,0)\}$ and adding MSB $e$ at stage
$j$ via $(a,r)\mapsto(a+e2^{j-1}\bmod d,\,2r+e\bmod9)$.
Our chart yields $\ell\le83$, so the odd part $g_0\le83$ and $a\le7$.
Fix the low $a$ bits of $u$ first and solve the two residues by this DP
using the next 14 bits. Thus every aligned block of width $2^{14+a}$
contains a permissible $u$, and from any admissible search threshold
the next such $t$ can be found within $W=2^{16}$ since $g\ge2^a$, ignoring the
arithmetic grid endpoints where the proof already has ample margin.
The separate exact checker decodes a legal time from any threshold
using the lookup DP and tests the constructive bracket.

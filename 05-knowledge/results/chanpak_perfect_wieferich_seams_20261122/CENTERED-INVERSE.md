# A square-root inverse for the full upper perfect-number neighbor

**PROVED.** For the symmetric apex construction below and each odd exponent $p=2d+1\ge3$, there is just one possible set of noncentral marks. Its binary address is the first $d$ fractional binary digits of $1/\sqrt2$. For $d\ge3$, the next digit determines the central mark as well. No two consecutive $d$ can yield solutions. A bounded number of binary runs cannot realize arbitrarily large exponents in this construction.

**FINITE-EXACT.** The resulting integer test fails for all $1\le d\le100{,}000$. This excludes this specific symmetric family through odd exponent 200001. The unrestricted tournament problem and the all-exponent existence question within this family remain open.

## Which marks are admissible?

Use the transitive spine, isolated highest mark and weighted inversion formula of the [three-contact construction](PROOF.md). Let $E_p=2^{p-1}(2^p-1)$. Start with the mark at binary exponent $x=2p-1$, which gives $1+2^{2p-1}$. Seek $E_p+1$ by adding a set $S$ symmetric about $c=p-1$. We require only that all new positions are positive and have distance at least two from $x$. The marks **inside** $S$ may be adjacent. These boundary conditions describe the symmetric construction being classified.

For any set of positive marks, the individual term of a mark is $2^t$ at position $t$, and the correction for any pair a distance $h$ apart is
\[
                 w(h)=2^{h-1}+[h=1].
\]
Here “individual” means the contribution with every other bit unmarked. Subtraction of pair corrections gives the exact count even for adjacent marks, by counting every marked-to-unmarked pair in the weighted inversion formula. For $u\in S$, the interaction with $x$ costs $2^{x-u-1}=2^{2c-u}$. Reflection $u\mapsto2c-u$ pairs all of these terms with the individual contributions from $S$, including the central fixed point when present. The sole remaining deficit is the sum of $w(h)$ over unordered pairs inside $S$. Thus the requested count is equivalent to that sum being $2^{p-1}$.

Empty sets and singleton sets have deficit zero. For any other symmetric set let $D=\max_{u\in S}|u-c|\ge1$. The two extreme marks contribute $2^{2D-1}$. If every point from $c-D$ through $c+D$ were marked, the total would be $2^{2D+1}-2$. Indeed an interval of length $n$ has $\sum_{h=1}^{n-1}(n-h)2^{h-1}=2^n-n-1$ and the adjacency correction adds $n-1$. An arbitrary subset lies between these bounds. Since $p-1$ is even, the unique possible even power exponent in this interval is $p-1=2D$. Hence $D=d$. In particular the requested radius satisfies both boundary conditions even for $p=3$.

## One binary address and one remaining equality

Let $R\subseteq\{1,\ldots,d\}$ be the positive offsets, so $d\in R$, and write
\[
 v=\sum_{a\in R}2^{a-1},\qquad
 W(v)=\sum_{a>b,\ a,b\in R} w(a-b).
\]
Let $e\in\{0,1\}$ say whether the center is marked. Pairs across the center have distances $a+b\ge2$ and total weight $2v^2$. The two sides contribute $2W(v)$. The center contributes $2e(v+[v\text{ odd}])$. Put $X=2^{2d-1}$. The target equation is exactly
\[
                    X=v^2+W(v)+e(v+[v\text{ odd}]).
\]
For fixed $a\ge2$ the contribution to $W$ from lower offsets is at most $\sum_{b=1}^{a-1}w(a-b)=2^{a-1}$. The smallest offset contributes zero, so $0\le W(v)\le v-1$. The correction to $v^2$ is consequently at most $(v-1)+(v+1)=2v$, strictly below the gap to $(v+1)^2$. Thus
\[
                       v=\lfloor\sqrt X\rfloor
                        =\lfloor 2^d/\sqrt2\rfloor .
\]
This address has $d$ binary digits, as required by $d\in R$.

Write $t=X-v^2$, which is positive. The next binary digit of $1/\sqrt2$ is $b=[t>v]$ since $\sqrt{4X}\ge2v+1$ iff $4t\ge4v+1$. For $d\ge3$, the word starts with 101, hence $R$ contains at least two elements and $W>0$. When $e=0$ the equality forces $t=W<v$; when $e=1$ it forces $t=W+v+[v\text{ odd}]>v$. Therefore $e=b$. Appending a binary digit shifts all old offsets upward by one and adds an offset 1 iff that digit is 1. Its interactions give
\[
 W(2v+b)=W(v)+b(v+[v\text{ odd}]).
\]
The exact necessary and sufficient test for $d\ge3$ reduces to **one equality**, $t=W(2v+b)$. At $d=1$ the values $v=1,t=1,W=0$ satisfy neither center choice; at $d=2$ the values $v=2,t=4,W=0$ likewise satisfy neither.

## Two structural obstructions

**Neighboring exponents.** All pair weights $w(h)$ are even, as is $v+[v\text{ odd}]$. Since $X$ is even, any solution has even $v$. Suppose $d$ and $d+1$ both work. The cases $d=1,2$ were excluded. At $d+1$ the address $2v+b$ must be even, so $b=0$ and the first solution has $e=0$, $t=W(v)>0$. At the next radius the remainder is $4t$, while the pair cost remains $W(2v)=t$. The second center must be present and $4t=t+2v$. Thus $v=3u,\ t=2u$ with $u>0$, giving $X=u(9u+2)$. Both factors would have to be powers of two. For $u=1,2$ the second factor is 11 or 20. For any larger dyadic $u$, it is $2\bmod4$ and exceeds 2. This is impossible. In particular the centered construction cannot work at both odd exponents $p,p+2$.

**Number of runs.** Let $\alpha=1/\sqrt2$. If a run after the first $l\ge1$ binary digits continues through digit $m$, then $\alpha$ is at distance at most $2^{-m}$ from one of the two enclosing dyadic endpoints $a/2^l$. For a run of zeros use the lower endpoint, and for ones the upper endpoint. For either endpoint,
\[
 |\alpha-a/2^l|
 =\frac{|2^{2l-1}-a^2|}{2^{2l}(\alpha+a/2^l)}
 >2^{-2l-1}.
\]
The numerator is a nonzero integer and the denominator's last factor is strictly less than 2. Hence $m\le2l$. The first run has length one because the expansion starts with 10. Each later run can at most double the cumulative length, including the last partial run of a prefix. A $d$-digit address must therefore have at least $1+\lceil\log_2 d\rceil$ runs. In particular any fixed budget on one-sided marked intervals eventually fails. A single filled centered interval is excluded already by its deficit $2^{2d+1}-2$, which has an odd factor greater than one for $d\ge1$.

## Reproduction

    python3 -B 05-knowledge/results/chanpak_perfect_wieferich_seams_20261122/centered_check.py

The check compares the formula to all admissible symmetric masks at $p=3,5,7$, with ordinary tournament path DP also run at $p=3,5$. It then performs the exact base-two square-root recurrence at every $d=1,\ldots,100000$. No floating-point logarithms, unproved assumptions on irrational binary digits, or primality filters enter the computation.

# Both-one closure: simultaneous even and odd congruences with literal +2 choices

The QRH residual packet gives a physical three-center bank:
$$
F_j(x,y)=432xy+129x+258y+15[x\text{ odd}]+30[y\text{ odd}]
 +45(p(x)+p(y))+135+2j,\quad j=0,1,2 .
$$
Here $p(2u)=p(u)+u+[u\text{ odd}], p(2u+1)=p(u)$, and $p(0)=0$.
No assertion about ordinary Collatz or prime distribution is needed.

**Proposition.** For every fixed positive integer $M$ and all sufficiently large odd targets $N$, there are two actual bank triples, one strictly below $N$ and the other strictly above $N$, within $O_M(N^{1/3})$. In each triple the same designated branch has the residue $N\bmod M$, and the two +2 differences are literal. We can in addition take the designated branch to agree modulo 3 even if $3\nmid M$. This is a fixed-modulus result, not exact placement.

## Paired-digit invariant and the obstruction on the old rail

For pre-addresses $X,Y$ with odd sum $S$, set $W=p(X)+p(Y)$, $r=[X\text{ odd}]$, $K=W-X+r$. Append complementary digits $e,1-e$ to $X,Y$. Then
$$
 S'=2S+1,\qquad K'=K+e(S+1)-2e(X+r).
$$
In fact the increment of $W$ is $X+r$ if $e=0$ and $Y+1-r$ if $e=1$; substituting proves the formula. Thus on $S\equiv-1\pmod4$ the state $K\bmod4$ is frozen. Also $p$, hence $W$ and $K$, are even. On the original odd-sum rail a direct substitution gives $F_j\equiv K+3+2j\pmod4$ for either parity orientation. If $j$ is already chosen modulo 3, this rail cannot repair the locked mismatch by another sufficiently deep complementary suffix.

**Escape via a noncomplementary ending.** Append a 1 to *both* pre-addresses: $x=2X+1,y=2Y+1,R=x+y=2S+2$. The checksum remains $W$, and
$$
 F_j=432x(R-x)-129x+258R+180+2j+45W .
$$
For pre-address $X$ even this reduces to $F_j\equiv K+3+3X+2j\pmod4$. Therefore selecting $X=0$ or $2\bmod4$ repairs either odd residue with $j$ held fixed. This is an exact change of terminal rule, not a relabeling of the old odd-sum rail.

## Independent clocks

Write $M=2^k m$ with $m$ odd; we may increase $k$ to at least 2 and replace $m$ by $\operatorname{lcm}(m,3)$. After enough complementary digits to force $S\equiv-1\pmod{2^k}$, reset $X\equiv0\pmod{2^k}$ with $k$ zero digits. For each $h=2,\ldots,k-1$ choose $\ell_h\ge k+h+1$ with $2^{\ell_h}\equiv1\pmod m$. Use alternative equal-length blocks
$$
 U_h=0^{\ell_h},\qquad V_h=1\,0^{h-1}\,1\,0^{\ell_h-h-1}.
$$
Both restore the reset. At the first 1 the dyadic charge is zero, and at the second 1 the current $X\equiv2^{h-1}\pmod{2^k}$ is even; hence $V_h$ adds $-2^h\bmod2^k$ to $K$, hence $W$, while $U_h$ adds zero. The $k$ zeros at the end are essential.

The two block maps on $(X,Y,W,1)\bmod m$, on the even-$X$ input and output parity sheet, are affine unipotent with nilpotent part cubing to zero: the $X,Y$ rows are translations, and the $W$ row adds a linear function of $X,Y$ and a constant. Consequently their $m$-th powers both act trivially modulo $m$, since $m$ is odd and $m\mid\binom m2$. The dyadic difference of $V_h^m$ and $U_h^m$ is $-m2^h$. Choosing one word in each slot for $2\le h<k$ therefore gives every allowable dyadic correction (the sum is $-4m$ times an arbitrary $(k-2)$-bit number) without changing the true sum or the odd-modulus residues. This stronger choice of slots costs $O_m(k^2)$ bits. Append the complementary two-bit word 00 or 10 to fix $X\bmod4$; its effect on $W\bmod2^k$ does not depend on which clock words were used, since their dyadic endpoints were identical. Then apply the both-one ending. Multiplication by 45 remains invertible on the dyadic clock.

## Odd congruence and length bookkeeping

Fix $j\in\{0,1,2\}$ by $2j\equiv N\pmod3$. Force $3R\equiv N-2j\pmod5$ if necessary; this is achieved while varying $S\equiv-1\pmod g$ in an arithmetic progression of step at most $5g$, where $g$ is dyadic. If $5\mid m$, choose $x\equiv0\pmod5$. If $9\mid m$, first choose $x\bmod3$ via $F_j\equiv-3x+6R+2j\pmod9$. The chosen congruences make $N-[F_j-45W]$ divisible by $d=\gcd(45,m)$; solve for $W\bmod(m/d)$ and lift arbitrarily modulo $m$.

For completeness, on the odd-$X$, even-$Y$ sheet the physical words $A=1$, $B=01$ (with complementary $Y$-word) induce
$$
 A(X,Y,W)=(2X+1,2Y,W+Y),\quad
 B(X,Y,W)=(4X+1,4Y+2,W+X+2Y+3).
$$
All arithmetic in the following paragraph is modulo $m$. The maps are permutations, and $A$ has a finite period $L$ dividing $\operatorname{ord}_m(2)$ since
$A^t(X,Y,W)=(2^t(X+1)-1,2^tY,W+(2^t-1)Y)$. In composition order (rightmost first),
$$
 C=A^{-2}B:(X,Y,W)\mapsto (X-1/2,Y+1/2,W+X-Y+3/2).
$$
It has $C^m=1$, since the accumulated linear term in $W$ over $m$ iterations is a multiple of $m$, using oddness of $m$. Then $E=C^{-2}ACA^{-1}$ fixes $X,Y$ and changes $W$ by $-3(X-Y)/2-1/2$. Use $C$ to reach $X=Y$, powers of $E$ for a freely chosen $W$-shift, and $C$ again to reach the desired $X$. These moves preserve the modular sum. All inverses are literal positive words: take $L=\operatorname{ord}_m(2)$ and $C=A^{L-2}B$, and use $C^{m-2}$ and $A^{L-1}$ in $E$. Both words have length a multiple of $L$. The exponents can be chosen from $0,\ldots,m-1$, so their total length is bounded solely in terms of $m$, and one can append whole $A^L$ periods to achieve one common length. The new equality of $X,Y$ is modular only; neither address is forced equal as an integer.

To attain this pair $(x,W)\bmod m$, use the prior odd paired-affine steering result (QRH packet, QRH_TRANSFER.md) on an odd-$X$ complementary prefix, then the reset zeros, chosen two-bit ending and both-one ending. All maps after the steerer and before reading the residues are bijections on the appropriate parity sheets modulo odd $m$: the address multipliers are powers of 2 and the checksum adds known linear or constant terms. So backsolve the desired preimage. The clock is the identity on these odd residues. For identical physical suffix length: steerer words built from $C=A^{-2}B$ and its conjugates have length multiple of the period of $A$, with the inverse represented by positive powers of $A$; use whole periods of $A$ for padding. The steerer starts and ends at odd $X$. First append $k$ complementary zeros and then a complementary 1 to enter that parity and ensure $S\equiv-1\pmod{2^k}$ during the entire programmable part. Thus $K\bmod4$ is available before selecting the 00/10 ending, even though the odd steerer must be backsolved relative to that ending.

## Absolute bracket with the same modulus

Set $g\asymp N^{1/6}$ dyadic. Choose $S\equiv-1\pmod g$ just above $2\sqrt{N/1728}$ (plus a fixed additive constant), with the stated extra mod-5 restriction if required. The prior checksum identity gives, uniformly on the near-diagonal corridor $X=S/2+d$,
$$
 W= (g-1)S/(2g) + O(|d|+g+S/g+1).
$$
Indeed writing $X=gu+b,Y=gv+(g-1-b)$, the reversed low-digit coefficients sum to $g-1$; their weighted error is $O(g|u-v|)=O(|d|+g)$, and the prefix checksum and remaining low terms cost $O(S/g+g+1)$.

On both-one closure the product term gives the main quadratic $432(S+1)^2-1728d^2$, plus an affine function of $S,d$; the checksum's common term is independent of $d$. Near $d=0$ all three branches are above $N$ by order at least $\sqrt N$. At $d=C\sqrt{Sg}$ with fixed large $C$, all branches are below (the quadratic loss dominates the $O(Sg)$ overshoot and all other terms). Partition the pre-address range into aligned blocks of the fixed power-of-two size corresponding to the steering/reset/clock suffix length. Each block contains one physical pair with the designated residue and the same true $S$: fixing $X$'s prefix fixes $Y$'s prefix because $S\equiv-1\pmod g$, and $g$ eventually contains every programmable suffix. Selected addresses from consecutive blocks differ by a bounded constant depending on $M$. Their designated counts differ by $O_M(\sqrt{Sg}+g+S/g)=O_M(N^{1/3})$. Apply sign crossing to lowest and highest branch separately to get strictly below and strictly above windows. Monotonicity is not assumed. The step size of the three choices is untouched.

## Meaning and limits

The Collatz analogy is a valuation gate, not an implication: the complementary rail forces a low-bit invariant, while a terminal change of branch accesses the missing low bits. This resembles tracking actual residue-transition data rather than treating all numbers with a given valuation as interchangeable. This does not prove a local-to-global principle: the programmable suffix length is not bounded independently of the modulus, and congruence agreement with a fixed modulus does not eliminate the absolute residual.

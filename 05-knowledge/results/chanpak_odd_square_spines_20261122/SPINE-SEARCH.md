# Weighted sums exclude twelve full tubes around rational directions at exponent 61

**FINITE-EXACT.** Put $N=2^{60}(2^{61}-1)-1$. Consider either remaining canonical Farey chart, with arbitrary nonnegative physical addresses $x,y$. If its marked arc goes from $N$ to $N+2$, then for **each** ordered pair of distinct $a,b$ in $\{2,3,11\}$,
$$
 |ax-by|>2^{40}\quad\text{on spine 4},\qquad
 |ax-by|>2^{39}\quad\text{on spine 8}.
$$
The six inequalities per spine hold simultaneously. These are tubes around six different slopes $y/x=a/b$. Every binary word on each side is allowed inside each tube. The regions between the tubes and the unrestricted exponent-61 question remain open.

The earlier determinant search selected individual offsets with few binary runs, and the balanced square-root search treated ratios near one. The observation here is that a whole interval of determinant offsets compresses **quadratically** when projected onto its transverse weighted sum. This permits trillions of consecutive signed offsets at each direction. The prime triple supplies directions for which the projection can be computed exactly. Its odd-square extremal property is not assumed in the exclusion. [Prior Farey directions and sparse determinants](../chanpak_decimal_kernel_20261122/FAREY-HYPERBOLA.md) · [Centered ordinary sums](../chanpak_farey_complementary_exact_20261122/BALANCED-SECTORS.md)

## The weighted inverse

**PROVED.** For any fixed coprime positive weights $a,b$, an absolute determinant bound $D$, and either fixed chart, the number of weighted sum slices needed for an exact test at a growing odd height $N$ is
$$
 O_{a,b,\mathrm{chart}}(1+D^2/\sqrt N) \qquad (0\le D\le\sqrt N).
$$
This bounds the slices enumerated, not the work needed to read every lattice candidate in them.

Write $x=2X+i, y=2Y+j$, $0\le i,j\le1$. The physical lower height is
$$
 H=K+QXY+TX+CY+h(f(2X+i)+f(2Y+j)),
$$
where $Q=4q, T=2(t+qj), C=2(c+qi), K=k+qij+ti+cj+ri+ej$. The read is $f(0)=0, f(2u+1)=f(u), f(2u)=f(u)+u+[u\text{ odd}]$. Use $0\le f(z)\le z$. Put $S=aX+bY, V=aX-bY$. At a solution, $|V|\le d=\lfloor(D+|ai-bj|)/2\rfloor$ and $4abXY=S^2-V^2$. Consequently
$$
 S\le U=\left\lfloor\sqrt{4ab(N-K)/Q+d^2}\right\rfloor.
$$
When $N<K$ the phase is empty. Let $L=\max(\lceil(T+2h)/a\rceil,\lceil(C+2h)/b\rceil)$. The upper height bound gives
$$
 N\le K+QS^2/(4ab)+LS+h(i+j).
$$
It suffices to start at the floor of the square root of $\max(0,4ab(N-K-LU-h(i+j))/Q)$. The upper and lower squared bounds differ by $O(d^2+U+1)$, and under the stated growth condition $U\asymp\sqrt N$ and the lower square root is also $\asymp\sqrt N$. Taking their difference proves the slice count.

At each $S$, choose $X_0\equiv a^{-1}S\pmod b$ nearest to $(QS+bT-aC)/(2aQ)$, and take $Y_0=(S-aX_0)/b$. For $b=1$ the congruence imposes no restriction. An arbitrary lattice point on the slice is
$$
 (X,Y)=(X_0+bz,Y_0-az).
$$
Define $A=Qab, B=bT-aC+Q(bY_0-aX_0)$ and $J=K+QX_0Y_0+TX_0+CY_0-N$. Nearest-lattice centering gives $|B|\le A$. The possible $z$ satisfy
$$
 J\le Az^2-Bz \le J+h(2\lfloor S/\min(a,b)\rfloor+i+j).
$$
On each sign of $z$ the quadratic is nondecreasing in $|z|$, so two outward integer roots give the entire interval. Intersect with $X,Y\ge0$ and the physical determinant bound, then check the full height.

## How the odd-prime projection selects sums

The code imposes $K+QXY+TX+CY\equiv N\pmod h$ before the literal word reads. For these directions $\gcd(b,h)=1$. Precompute the admissible $X\bmod h$ by solving the following **linear** congruence for $S\bmod h$ at each $X\bmod h$:
$$
 (QX+C)S\equiv aQX^2-(bT-aC)X-b(K-N)\pmod h.
$$
The gcd of $QX+C$ and $h$ decides solubility and the reduced linear equation lists all answers. CRT then combines $X\bmod h$ with $X\equiv a^{-1}S\bmod b$.

**PROVED.** Let $p$ be an odd prime dividing $h$, with $p\nmid abQ$, and $\Delta=Q(N-K)+CT\not\equiv0\pmod p$. As the chart solution modulo $p$ varies, the exact number of distinct $S\bmod p$ is
$$
 \frac{p+\left(\frac{ab\Delta}{p}\right)}2 .
$$
Indeed $u=QX+C, v=QY+T$ satisfy $uv=\Delta$ and $au+bv=QS+aC+bT$. Calling the latter quantity $w$, the discriminant $w^2-4ab\Delta$ must be a square, zero included, and this condition suffices over the odd field. The character sum of $w^2-4ab\Delta$ is $-1$ because $(w-z)(w+z)=4ab\Delta$ has $p-1$ ordered factor pairs. Counting the zero discriminants gives the formula. If $\Delta=0$, every $S$ occurs. For our charts, modulo 5 the coefficient $Q$ vanishes and $T\ne0=C$; fixing $X$ and varying $Y$ admits every $S$. The dyadic information enters explicitly through the physical parities $i,j$.

This predicts the fraction of slices passing the sieve before any integer search. On spine 4 the two other primes are 17 and 29, on spine 8 they are 29 and 257. In every parity the character signs of $\Delta$ are respectively $(-,+)$ and $(+,-)$. The following numbers of occupied classes are exact in each parity. Exchanging $a,b$ changes individual absolute candidates while preserving these local support sizes.

| Unordered weights | Occupied mod 2465 at spine 4 | Occupied mod 37265 at spine 8 |
|---|---:|---:|
| 2, 3 | 675 | 9675 |
| 2, 11 | 675 | 9600 |
| 3, 11 | 600 | 9675 |

## Finite search

The first count below includes the four parity branches. A slice passes if its modular table is nonempty. The last numeric column counts full integer readings after the root and congruence tests and after clipping to the tube. No full reading equaled $N$. These are counts for this replay, and the exclusions can intersect older searches in other coordinates.

| Spine | $a,b$ | Slices | Passing slices | Full reads |
|---:|---|---:|---:|---:|
|4|2, 3|229,220,420|62,768,266|436,695|
|4|3, 2|229,220,423|62,768,272|435,645|
|4|2, 11|119,706,490|32,779,667|435,562|
|4|11, 2|119,706,502|32,779,672|435,815|
|4|3, 11|97,739,952|23,790,655|290,404|
|4|11, 3|97,739,966|23,790,663|290,189|
|8|2, 3|455,360,492|118,223,805|3,683|
|8|3, 2|455,360,491|118,223,888|3,578|
|8|2, 11|237,804,326|61,261,804|3,576|
|8|11, 2|237,804,334|61,261,803|3,601|
|8|3, 11|194,166,430|50,410,845|2,384|
|8|11, 3|194,166,438|50,410,828|2,409|

The positive chart data come from the [physical inverse](../chanpak_decimal_kernel_20261122/FAREY-HYPERBOLA.md). For the recorded bounds, $S<2^{53}$, $Qab<2^{35}$, and the centered product $QX_0Y_0$ differs from $N$ by less than $2^{110}$. All chart arithmetic fits signed 128-bit integers. The floor of $4abM/Q$ is evaluated by splitting $M$ at division by $Q$, before multiplying by $4ab$. Floating roots only seed candidates. Exact integer comparisons certify both root endpoints. The independent brute-force controls include axes and just-inside and just-outside determinant boundaries. Large planted pairs test those boundaries in all twelve directions/spines.

From this directory run the following commands. The half-width can also be reduced to reproduce shorter subsets.

~~~sh
c++ -std=c++17 -O3 -Wall -Wextra -Werror spines-weighted.cpp -o /tmp/chanpak-weighted
/tmp/chanpak-weighted --check
for n in 4 8; do
  D=1099511627776
  if [ "$n" = 8 ]; then D=549755813888; fi
  for a in 2 3 11; do
    for b in 2 3 11; do
      if [ "$a" != "$b" ]; then /tmp/chanpak-weighted "$n" "$a" "$b" "$D"; fi
    done
  done
done
~~~

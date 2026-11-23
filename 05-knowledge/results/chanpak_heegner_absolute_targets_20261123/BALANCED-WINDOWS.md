# A relative-height selector for prescribed modular Hamilton and path-slot readings

Use an apex over a transitive chain with contacts $O,I$. The [inherited guarded-tile construction](../chanpak_fork_curvature_20261123/ALL-PRIME-WORDS.md) realizes the first $(p-1)/2$ readings $Q_s$ at each odd prime power, simultaneously at finitely many such components and at chosen parity. Its longer-guard variant also controls $Q_2$ at $3^e$ provided $Q_2\equiv(-1)^{|T|}Q_1\bmod3$. The goal here is to prescribe the approximate **integer** height as well as those exact modular readings.

**Square-root-window theorem.** Fix a finite odd modulus (or the finite list of its prime-power observations). There is a constant $C_M$ such that for every sufficiently large real $X$, every choice of order parity, and every compatible modular target, some strong singleton-apex tournament has exactly the modular readings and
$$
|H(T)-X|\le C_M\sqrt X.
$$
In particular each fixed relative window $[(1-\epsilon)X,(1+\epsilon)X]$ works for all sufficiently large $X$, with threshold uniform over the finite target states and both parities. More precisely, for every sufficiently large fixed chain length $n$ in the chosen parity class, the scaled heights $2^{-n}H(T)$ for the target form an $O_M(2^{-n/2})$-net of $[1/4,1]$. No eventual exact integer coverage is claimed.

## Modular read corrections with a prescribed short suffix

For the given modulus let $L$ be the lcm of the unit exponents $\varphi(p^e)$ of its prime-power components. Use the long guarded macros $C=O^gI^g,P=OC,S=CI$, where $g$ is at least the exponent needed at $3^e$ if present, and otherwise $g=1$. Every nonempty guarded macro word starts with $O^g$ and ends with $I^g$. Write $\mathcal O=PC^{-1}$, $\mathcal I=C^{-1}S$ in the finite group for the unit channels. These have precisely the usual one-letter matrices
$$
O_\lambda=\begin{pmatrix}1&1&0\\0&\lambda&0\\0&0&1\end{pmatrix},
\qquad
I_\lambda=\begin{pmatrix}1&0&0\\0&\lambda&1\\0&0&1\end{pmatrix},
$$
length one apiece, and zero increments in nilpotent scalar channels. All the relevant unit channels satisfy that $1-\lambda$ is a unit. Hence $\mathcal O^L=\mathcal I^L=1$ by summing geometric series in each channel, including the length coordinate.

Choose an arbitrary fixed finite suffix $T$ ending $I$, not necessarily a concatenation of macro tiles. Place it after a nonempty guard word $W_0$. Before $W_0$ use one of two long prefixes, $O^a$ or $OI^a$, with $a$ a sufficiently large multiple of $L$. At any unit channel the prefix product is respectively the identity or $O_\lambda$, independent of $a$. In a nilpotent channel, the first prefix has no internal $OI$ pairs; the second has $g_\lambda=1+\lambda+\cdots+\lambda^{a-1}=(1-\lambda)^{-1}$, independent of large $a$. Crossings from $O^a$ into $W_0$ vanish there because $W_0$ begins with $O^g$, and crossings from the single initial $O$ of $OI^a$ into later pieces vanish because $a\ge g$. Crossings out of $W_0$ into $T$ likewise vanish because $W_0$ ends with $I^g$. Thus nilpotent-channel contributions are the sum of the three independent pieces. These observations include $\lambda=0$ for every prime and $\lambda=3/2$ when it is used at $3^e$.

Choose any desired residue class for $|W_0|\bmod L$. Start with the element $\mathcal O^r$ in this class. Append the inherited central corrections to this element before converting to positive macro words. After that conversion, all physical boundaries used above are present. To see the correction equation, the unit-channel matrix for a complete word is the prefix matrix times the $W_0$ matrix times that for the suffix $T$. A central top-right change in $W_0$ adds the same change to the full word's top-right read because it commutes and is a translation with bounding diagonal entries one. The nilpotent changes add by the vanishing of crossings just proved. The total order class remains fixed, so the factors $s^{n+1}$ that convert the normalized observations into the specified residues remain fixed and are units. The inherited central span therefore solves the target independently at all components. At $3^e$ the required correction lies in the inherited congruence sublattice: the baseline physical word has the same order parity and hence satisfies the necessary mod-three relation, just as the target does. The construction uses only the residue class, prefix mode and actual suffix $T$, not the large value of $a$.

Crucially, the length bound $B$ on $W_0$ is **independent of suffix length**. The corrected word is obtained from one of the $L$ bases and a tuple of central correction coefficients in finitely many fixed residue rings. Choose a preferred positive macro realization in advance for each such tuple; there are only finitely many, so their maximum length is a common $B=B(M)$, independent also of the target and both prefix modes. One can bound it explicitly using the signed-word construction and finite group exponent. This is stronger than taking a maximum over the literal suffixes, whose lengths will grow with $n$.

Consequently both modes work at all sufficiently large chain lengths $n$: choose the guard length class $n-k$ or $n-k-1$ mod $L$, and let $a$ fill the rest as a nonnegative multiple of $L$. With $|W_0|\le B$, in either mode $a\ge n-k-B-1$. The nonempty guard word ensures strength even if the suffix begins with $I$: the complete contact word always begins with $O$ and ends with $I$.

## Exact real height change and dyadic balancing

Write $W=W_0T$, $m=|W|$, and put $b_W=\sum_{j:W_j=I}2^{j-1}$. This is an integer whose high digits are the last contacts, in reverse order. The literal singleton-apex formula
$$
H(V)=1+\#\operatorname{adj}_V(OI)+
\sum_{i<j:V_i=O,V_j=I}2^{j-i-1}
$$
holds for all guarded contact words. Applying it at the above prefix boundaries yields exactly
$$
\begin{aligned}
H(O^aW)&=H(W)+(2^a-1)b_W,\\
H(OI^aW)&=H(W)+2^a(1+b_W).
\end{aligned}
$$
The first equation counts all new pairs between a prefix $O$ and a suffix $I$, and there is no new adjacent $OI$ on the boundary. In the second, the first $O$ contributes $2^a-1$ to the internal positive pairs, one adjacent $OI$, and $2^ab_W$ to the pairs ending in $W$. The boundary ends in $I$ and begins with $O$, so again no other adjacent $OI$ is created. All quantities are exact integers.

In the first mode $n=a+m$, so scaled height is within $2^{-a}$ of $b_W/2^m$, since $0\le H(W),b_W\le2^m$ (a Hamilton history chooses its two sides at each old vertex). In the second mode $n=a+1+m$, scaled height is within $2^{-a-1}$ of $(1+b_W)/2^{m+1}$.

Now vary the $k$-letter suffix $T$ ending $I$. The leading $k$ bits of $b_W$ are prescribed by $T$; its arbitrary earlier part changes $b_W/2^m$ by less than $2^{-k}$ inside the corresponding dyadic interval. Thus model ratios in the first mode approximate any target in $[1/2,1]$ to within $2^{-k}$, independently of the correction word $W_0$. The second mode approximates any ratio in $[1/4,1/2]$ to within $2^{-k-1}$: the extra $1$ in $1+b_W$ stays within the upper endpoint of the same half-scaled dyadic interval. At endpoints choose the adjacent available interval. Take
$$k=\left\lfloor\frac{n-B-1}{2}\right\rfloor.$$
Then $a\ge k$, for every possible correction. For sufficiently large $n$, both $k\ge1$ and $a\ge g$, so the modular arguments still hold. Adding suffix quantization error to the prefix error bounds the scaled height error by $2\cdot2^{-k}=O_M(2^{-n/2})$.

Finally the intervals $[2^{n-2},2^n]$ cover all sufficiently large $X$ even when $n$ is restricted to one parity. Choose such an $n$ with $X/2^n\in[1/4,1]$. Multiplying the scaled error by $2^n$ gives $O_M(2^{n/2})=O_M(\sqrt X)$, since $2^n\le4X$. This completes the operational bridge from arbitrary modular targets to shrinking relative height windows: the finite correction occupies a bounded middle, while the high binary digits of the **integer** height are chosen by the terminal tape and a long physical prefix.

## A separate exact-length padding recipe within the original macros

If arbitrary suffixes are not needed, any fixed group word can be padded within the original guarded alphabet to every sufficiently large chain length in its residue class modulo $L$. Let $c=|C|=2g$ and $E$ be the finite group exponent bound. The two positive identity words $C^E$ and $(PC^{E-1})^L$ have chain lengths $A=cE$ and $B=L(cE+1)$. Their gcd is exactly $L$, because $L\mid E$. The conductor for their two-generated length semigroup, in physical length units, is
$$
L(A/L-1)(B/L-1)=cE(cE-L)<c^2E^2.
$$
Thus exact-length padding is deterministic by the elementary two-coin formula. The balancing proof uses the longer uniform prefix instead, as it permits control of the leading real height digits.

## Native factor-seven fork height: a second, base-$3/2$ window

The same modular correction gives a quantitative window for the **physical fork height** $h=280Q_2+7H$, not just the filler height. For any fixed compatible modular target and either order parity, every sufficiently large $X$ has a strong singleton-apex filler whose fork parent is at distance
$$
O_M(X^\theta),\qquad
\theta=1-\tfrac12\log_3(3/2)<1
$$
from $X$, with the specified exact filler residues. The same holds for the center and either exit, with their respective rescaling and additive multiple of $H$.

To see it, work now over the reals and set $\lambda=3/2,\ \mu=1/2,\ r=\lambda^{-1}=2/3$. For a word $V$, let $g_\eta(V)=\sum_{i<j:O_i,I_j}\eta^{j-i-1}$ and $b_\lambda(W)=\sum_{j:W_j=I}\lambda^{j-1}$. The two-channel formula gives
$$
\frac{Q_2(V)}{3^n}=2r^n+\tfrac12r^n\bigl(g_\lambda(V)+g_\mu(V)\bigr),
\qquad n=|V|.
$$
For the two physical prefixes the exact positive-channel formulas are
$$
\begin{aligned}
g_\lambda(O^aW)&=g_\lambda(W)+2(\lambda^a-1)b_\lambda(W),\\
g_\lambda(OI^aW)&=g_\lambda(W)+2(\lambda^a-1)+\lambda^a b_\lambda(W).
\end{aligned}
$$
The first follows from summing over the possible prefix $O$ locations. The second includes the contribution from all $I$'s inside $OI^a$ and the one prefix $O$ continued to $W$.

Now $g_\lambda(W)=O(\lambda^m)$, $b_\lambda(W)=O(\lambda^m)$ for $m=|W|$, by summing geometric series (for the former use the bound $\sum_{\ell\ge1}\ell\lambda^{-\ell-1}$ after division by $\lambda^m$). Also $g_\mu(V)\le 2n$. Therefore the modeled $Q_2/3^n$ ratios in the two modes are
$$
\lambda^{-m}b_\lambda(W),\qquad
\tfrac13\lambda^{-m}(2+b_\lambda(W)),
$$
respectively, with errors $O(r^a+nr^n)$ uniformly.

The quantity $\lambda^{-m}b_\lambda(W)$ reads the reversed contact bits in weights $r,r^2,\ldots$, with top digit one. Because $r=2/3>1/2$, infinite strings with top digit one fill the entire interval $[r,r/(1-r)]=[2/3,2]$: peel successive digits using the two overlapping intervals $[0,r^2/(1-r)]$ and $[r,r/(1-r)]$. Prescribing $k$ leading digits leaves error at most $2r^k$, uniformly over earlier contacts of the correction word. The second mode scales this range by $1/3$ and adds only $O(r^m)$, so it covers $[2/9,2/3]$ at the same precision. Use the same $k,a$ split as for the binary selector. The relative scaled error is $O_M(r^{n/2})$.

Finally the term $7H/3^n$ is $O(r^n)$, so $h/(280\cdot3^n)$ has an $O_M(r^{n/2})$-net over $[2/9,2]$. These normalized bands tile all sufficiently large scales even at fixed parity (neighboring allowed $3^n$ differ by a factor nine). Multiplying the scaled error by the physical scale gives the displayed $\theta$. The center $3h$ and exits $3h\pm2H$ have the same dominant scale after division by three, since $H=O(2^n)$. This result neither forces $h$ to equal a specified integer nor forces its cofactor to be prime.

Replay the two exact real prefix formulas with `python3 chanpak-height-balanced-apex-check.py`.

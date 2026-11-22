# A fixed strong lower source for each parity, and what a single port can't do

**PROVED.** Outgoing masks $(2,16,11,19,133,31,319,111,191)$ and ordered ports $(8,5,1)$ give a strong source for every even $p\ge6$. The usual return-arc test certifies at most two appearances of every port in a compressed spanning path: each return uses an arc from an outneighbor of the port to an inneighbor and no two possible return arcs are disjoint after identifying singleton vertices but treating marked vertices as replicable. The numbers of compressed spanning paths in lexicographic run order $111,112,\dots,222$ are $(93,63,33,24,31,21,11,8)$. A one-run slot has one filling and a two-run slot of width $w$ has $2^w-2$ fillings, by selecting a nonempty proper subset for its first increasing run. The resulting polynomial in the slot powers $X,Y,Z$ is
$$
  -1-X-5Y-5XY+5Z+5XZ+8YZ+8XYZ.
$$
With widths $(p,p/2-2,p/2-2)$ we have $Z=Y$, $X=16Y^2$, giving $128Y^4-8Y^2-1=E_p-1$. The exact source connectivity, cap, signature and an independently evaluated literal inflation at $p=6$ are in [check_even.py](check_even.py). Strong connectivity persists under module inflation: concatenate the base paths between blocks or use a base cycle to connect two vertices of the same block. The order is $2p+2$. With the [existing odd rail](../chanpak_support_birth_grammar_20261122/ODD-PERFECT-PORT.md), this gives two fixed constructions for all $p\ge5$; it does not construct the upper endpoint or give a common edge.

## The holes are actual boundary conditions for one growing port

Fix a source with one transitive growing slot of width $w=p-1$, other vertices fixed. The path count has a finite expansion
$$
    \sum_r c_r \sum_{j=1}^r(-1)^{r-j}\binom{r}{j}j^w,
$$
where $c_r\ge0$ counts compressed words with $r$ visits and each run must receive a nonempty subset. If this agrees with $2\cdot4^w-2^w+d$ on infinitely many $w$, uniqueness of a finite exponential expansion (divide by the largest surviving base) gives equality for **every** $w\ge1$. In particular the source at $w=1$ would have $6+d$ paths. This excludes $d=1$ and $d=15$ by the actual missing 7 and 21, rather than treating their arithmetic images as an abstract phase collision. It also excludes $d=-1$: the largest nonzero $c_r$ must be $c_4=2$, but every return between visits to the marked port needs at least two distinct ordinary vertices in the port's strong component. Four visits require six such vertices, so the marked strong component would have order at least 7, impossible with total count 5 by Busch's floor and strong-component factorization. No implication is claimed for changing the width slope/offset or for several growing slots.

## Two restrictions on extended module grammars

For a fixed source, let transitive growing slots have widths $a_i m+b_i$, with every $a_i>0$ an integer. Suppose each slot has at most three runs in every spanning compressed word. If some word has a three-run slot, maximize the leading base $\prod_i r_i^{a_i}$ among the words of that kind. It is divisible by 3 and its coefficient is positive. All lesser bases of each word are strictly smaller than its leading base, and words with no three-run slot produce only powers of two. Thus the non-dyadic mode cannot cancel, even if one asks for equality with a fixed dyadic polynomial only on an infinite subsequence. This argument does not extend unaltered to four-run slots: their principal mode may be dyadic and their 3-mode negative.

For any module $M$, its $r$-run filling count is $B_r(M)=r!\sum_{\mathcal P}\prod_{U\in\mathcal P}H(M[U])$ over unordered partitions into $r$ nonempty blocks. For equal-order $M,M'$ this implies $2r!\mid B_r(M')-B_r(M)$ by Rédei oddness. Writing the inflated count as a linear combination of the $B_r$'s with coefficient $H(Q)$ on $B_1$ gives
$$
   \Delta H_{\rm inflated}\equiv H(Q)\Delta H_{\rm module}\pmod4.
$$
A dormant internal flip or mod-4-zero flip therefore cannot become a +2 move after homogeneous embedding, independent of the external quotient. If both modules agree through index $r-1$, the entire difference is divisible by (2r!).

For a growing slot of the form $U_i\oplus \operatorname{Trans}_{a_i m+b_i}$ with fixed nonempty $U_i$, the same run count is
$$
  B_r=\sum_{j=1}^r(-1)^{r-j}\binom rj Z_{U_i}(j)j^{a_i m+b_i},
$$
where $Z_{U_i}(j)$ counts ordered covers by $j$ possibly empty paths. This is inclusion–exclusion on empty slots and multiplication across an ordinal cut. Assume a fixed finite run cap. In the exponential expansion the constant term is an integer multiple of $\prod_iH(U_i)$, because base 1 requires $j=1$ for every slot. A perfect rail has constant $+1$ or $-1$, so each $H(U_i)=1$ and each $U_i$ is transitive. Bounded ordinal padding by a fixed nontransitive module cannot repair these polynomial profiles.

## A strong upper rail for every even exponent from 6 onward

**PROVED, separate scalar rail.** The eight-vertex source with outgoing masks
$(2,16,11,19,133,31,63,111)$ is strong. Its ordered ports are $(1,3,5)$, with the same two-run return certificate. Its compressed signature in order $111,112,\dots,222$ is $(53,22,27,12,37,16,18,8)$, giving
$$
   1-Y+(X-2Z)(1+2Y)+8XYZ.
$$
At even $p\ge6$, use widths $(p/2-1,p-1,p/2-2)$. Then $X=2Z$ and $8XZ=2^p$, leaving $1-Y+2^pY=E_p+1$, on $2p+1$ vertices. [The exact checker](check_even.py) independently counts the literal order-13 inflation at $p=6$. The short slots now cancel at ratio two, not at equality as in the lower rail. The lower and upper sources are not a literal single-arc pair. Even exponents beyond two do not give even perfect numbers, but the identity applies to the Euclid-shaped integers as requested.

For three cap-two slots with positive slopes and a positive triple tensor coefficient, isolate the normalized variable $t=2^{(p-\epsilon)/2}$. If the chart is a perfect rail, the slopes must sum to 4 (the top monomial is unique) and some proper subset sums to 2 (or there would be no quadratic term). By complementation one of the slots therefore has slope 2 and the other two sum to 2. If their positive slopes differ, each of their singleton and long-slot mixed terms has a distinct off-target degree and must vanish. Such two shorter slots interact only through their product. Algebraically setting both their slopes to 1 preserves the polynomial along the same parity progression; for rational slopes rebase at an exponent where the original widths are integral and use the old initial widths as the new offsets. Thus varying affine slopes alone cannot repair the failed odd-upper case of this positive triple, cap-two, three-port architecture. This does not exhaust other source sizes or architectures.

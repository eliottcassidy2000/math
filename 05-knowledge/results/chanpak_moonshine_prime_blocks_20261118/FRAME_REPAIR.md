# Same-lattice height-neutral switches and an exact growing-width compiler

This follows the [Gaussian lattice construction](LATTICE.md). The main distinction is no longer just that an unmarked lattice can have different Hamilton counts. **One and the same lattice can have distinct frames at exactly the same two-slot height, and even a single-vertex cut can switch between them.** The smallest nonzero change permitted by the congruences is attained. A common exterior construction then amplifies the switch at a power-law scale. This is a partial prime selector, not all-prime or Chan–Pak coverage.

## 1. Exact cut and divergence formula

Fix an exterior tournament $U$ and add $v$, writing $b_x=1$ if $v\to x$. Put $T_b=U+v$. The single-vertex cut replaces $b$ by $1-b$ and leaves all arcs of $U$ unchanged. For $x\ne y$ let
$$
 M_q(x,y)=\sum_{\substack{S\subseteq U,\;|S|>0\ {\rm even}\\x,y\in S}}
       \#\{\text{Hamilton paths in }U[S]\text{ from }x\text{ to }y\}\,
       f_{U\setminus S}(q).
$$
Every odd directed cycle containing $v$ corresponds to exactly one such internal path from an out-neighbor of $v$ to an in-neighbor. All remaining cycles avoid the internal support. The odd-packing identity gives
$$
 f_{T_b}(q)=q f_U(q)+2q\sum_{x,y}b_x(1-b_y)M_q(x,y).
$$
Define $W_x(q)=\sum_y(M_q(y,x)-M_q(x,y))$. The two pieces of the response are therefore
$$
\begin{aligned}
 f_{T_{1-b}}(q)-f_{T_b}(q)&=2q\sum_x b_xW_x(q),\\
 \tfrac12(f_{T_{1-b}}(q)+f_{T_b}(q))&=
 q f_U(q)+q\sum_{x<y}|b_x-b_y|\big(M_q(x,y)+M_q(y,x)\big).
\end{aligned}
$$
In the second line the symmetric edge weights are positive for positive integral $q$, since the two-vertex support already contributes. In the first line the weights sum to zero. Thus the **center is a weighted cut, but the signed change is a linear divergence read**. To hold the two-slot height fixed, choose $b$ with $\sum b_x W_x(2)=0$; the remaining signed shift in ordinary count is determined by $W(1)$. This makes the target problem a coupled cut-and-subset-sum problem, not simply a search for a point in the integer lattice.

## 2. Sharp elementary step and actual prime pairs at fixed height

For two frames on the same order, the mod-3 identity $B\equiv(-1)^{|T|}H\pmod3$ and Rédei oddness give $\Delta H\equiv0\pmod6$ if $\Delta B=0$. The following literal order-nine example attains 6 using just one cut:

| Exterior outgoing masks | Out-contact mask | Common $B$ | $H$ before and after |
| --- | ---: | ---: | ---: |
| $(200,29,57,176,65,19,46,118)$ | 101 | 38344 | 2399, 2393 |
| $(110,84,152,242,353,454,388,275,15)$ | 65 | 181928 | 9623, 9599 |

The first is a new search witness; the second comes from the background pass. For the second, the marked $W(2)$ entries at 0 and 6 are $-7526,+7526$, while the $W(1)$ entries are $-1069,+1057$, giving $\Delta H=-24$. The subset DP checks both full orientations independently of the identity. The sharpness claim concerns the **size of the step**, not its least possible vertex order.

From that same background pass, the nine-vertex lattice represented by $(18,100,89,51,162,197,153,15,255)$ has switch masks 185, 251 and 164 giving $(H,B)=(1127,22156),(1139,22156),(1151,22156)$. These are three frames of **one** lattice, unlike selecting separate equal-$B$ cores. Apply the existing $\alpha=6$ positive paired reader after 23 ordinal singletons. The common center is $6\cdot22156\cdot2^{23}=1115147993088$, and masks 251 and 164 select the two primes $1115147994227$ and $1115147994239$, along with their sevenfold partners. Primality is checked by complete trial division. Similarly the second elementary move in the table, with reader $\alpha=20$ and five ordinal singletons, gives the two prime quotients $116443543$ and $116443519$. The claims of isometry concern the incoming frames, not the whole composed tournaments.

## 3. A repeatable growing-width compiler without a prime bank

Let the selected neutral pair have $(H,B)=(h_-,B_0),(h_+,B_0)$. The four-vertex count-five tournament $(2,4,9,3)$ has $(H,B)=(5,32)$, so joining $j$ copies to either frame, with $k$ ordinal singletons, gives
$$
 p_\pm=6B_0\,32^j2^k+h_\pm 5^j,\qquad j,k\ge0.
$$
The two signed-height readers with slopes 6 and 42 and intercepts 1 and 7 give **all four actual outputs** $p_-,p_+,7p_-,7p_+$, respectively. The common exterior choices do not depend on which incoming frame is selected. The offset width is $(h_+-h_-)5^j$. With $k$ fixed and $j$ growing it is of order $p^{\log_{32}5}$, an exponent about $0.4644$, larger than the reader-bank lower threshold proved in the main note. This does not evade the sparse-support issue: the two indices only give $O((\log X)^2)$ distinct output requests below $X$.

The target-based inverse loops over $j$ with $6B_0 32^j<p$, subtracts each candidate $h_\pm5^j$, and tests whether the quotient by $6B_0 32^j$ is a positive power of two. This is $O(\log p)$ such integer tests, and it builds literal tournaments when it succeeds. The sharp-step pair at $j=1,k=0$ yields prime quotients $7{,}374{,}013$ and $7{,}374{,}043$. The 24-step pair at $j=2,k=18$ yields prime quotients $293{,}015{,}554{,}074{,}983$ and $293{,}015{,}554{,}075{,}583$, a certified separation of 600. No claim is made that primes continue indefinitely in this family, or that the enlarged source frames stay in the same switching class. The lattice choice is made in the marked inner block.

## 4. Why the simplest way of growing a single lattice switch fails

A different growth idea *does* preserve the switching class within each larger pair: put an equal transitive block of order $a$ on each side of $U$ in its ordinal order, make $v$ lose to both external blocks, then switch $v$. But the old $q=2$ zero need not persist. Define $E_q$ to be the weighted even-internal-support response $\sum b_x W_x(q)$ before adding the padding. Define $O_q$ by the same endpoint difference but summed over **odd** internal supports. Allowing one exterior end in a cycle gives the exact formula
$$
 \Delta f(q)=q^{a+1}\big( ((q+1)^a+(q-1)^a)E_q
                         +((q+1)^a-(q-1)^a)O_q\big).
$$
Indeed internal paths have forward block order. A left-only exterior part adds a prefix to the path ending internally, a right-only exterior part adds a suffix to the path starting internally. A path using both exterior ends has equal contact bits and contributes zero to the difference. For an exterior subset of size $\ell$, the binomial coefficient counts its unique transitive path, and the internal parity required is $\ell\bmod2$. The $\ell=0$ contribution is exactly the original switch, multiplied by the two unused blocks. Taking even and odd binomial sums yields the formula.

In the sharp example the responses are $(E_1,O_1)=(-3,1)$ and $(E_2,O_2)=(0,2)$. The padding gives $\Delta H=-2^{a+1}$, but $\Delta B=2^{a+2}(3^a-1)$ at every $a\ge1$, so it ruins neutrality. The second example likewise has $(E_2,O_2)=(0,15)$, not a double zero. If a primitive source could instead make $E_2=O_2=0$ but $E_1+O_1\ne0$, the same padding would give a true enlarged **same-lattice** unbounded-height-neutral switch with $|\Delta H|=2^a|E_1+O_1|$. The two known examples fail the newly exposed odd-support condition. Sampled searches in small orders did not furnish that stronger source; no nonexistence claim is made.

**Scope and next obligation.** The power-law width result uses ordinal composition *after* selecting the inner frame; the fully enlarged same-lattice route needs the odd-support condition or a different interface. Even with growth, neither selector fills every target height. The exact weighted-cut center and linear divergence above provide a way to test that joint requirement, and prevent mistaking a zero signed response for a feasible baseline.

Replay: run **python3 frame_replay.py** alongside **replay.py**. For instance **python3 frame_replay.py 2051108878529081 --scale 7** materializes a partner from the growing-width family.

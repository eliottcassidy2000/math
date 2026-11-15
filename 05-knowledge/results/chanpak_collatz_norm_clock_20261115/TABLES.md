# The signed mod-eight clock: tables, Collatz, and a literal tournament residue repair

**Scope.** Theorems below are elementary: an exact table encoding, a cycle-word cancellation condition, and a physical mod-24 refinement of the earlier QRH bracket. They are not a proof of Collatz or Chan–Pak. The encoding uses only four distinguished table entries, not a graph isomorphism involving all rows.

## 1. What the middle four multiplication-table entries encode

In the modulo-$m$ table indexed by $1,\ldots,m-1$, take the central rows and columns $i=(m-1)/2, i+1=(m+1)/2$. Since $i+1=-i$ and $2i=-1$ modulo $m$, the diagonal entries are the least positive residue of $4^{-1}$, and the off-diagonal entries are its negative. For $m=4x+\varepsilon>1, \varepsilon\in\{+1,-1\}$, the two representatives are therefore $x$ and $y=3x+\varepsilon$, each twice, with $x+y=m$. Indeed $4x=m-\varepsilon$, and $4y=3m+\varepsilon$; the inverse is $y$ for $\varepsilon=+1$ and $x$ for $\varepsilon=-1$. The “two copies of m” are the sums of the two complementary pairs, not extra entries.

* If $x$ is odd, $x\to y$ is the raw odd step $3x+\varepsilon$, and $y$ is even. The labels $m\equiv5,3\pmod8$ encode every odd input of the $+$ and $-$ maps respectively.
* If $x=2d$, $y=6d+\varepsilon$ is odd, and $3m+\varepsilon=4y$ is an exact two-halving step from the odd source $m=8d+\varepsilon$. The labels are $1,7\pmod8$ respectively. This sheet re-encodes a subfamily of the odd-input sheet, not a second disjoint half of the odd steps.
* The full cross-sign bridge on this second sheet is
  $$
(4d+\varepsilon)\xrightarrow{\ 3t-\varepsilon,\ v_2=1\ }
6d+\varepsilon
\xleftarrow{\ 3t+\varepsilon,\ v_2=2\ }(8d+\varepsilon).
  $$
  Both arrows include the indicated halvings. Thus the center of the table contains the midpoint $y$ and the common step size $2d$. Every positive odd (y>1) coprime to 3 has exactly one such bridge. The fixed point $y=1$ is the degenerate $d=0,\varepsilon=+1$ case (with table modulus 1 excluded).

There is also a complete difference dictionary. For the odd-to-odd map $F_\varepsilon(n)=(3n+\varepsilon)/2^k$ and signed half-difference $d=(n-F_\varepsilon(n))/2$,
$$
6d+\varepsilon=(2^k-3)F_\varepsilon(n).
$$
The shallow exponents $k=1,2$ are precisely the two unit multipliers $2^k-3=-1,1$; thus an all-gaps shallow construction is automatic, but its absolute midpoint is rigid. This identity for $+$ and the significance of the unit multipliers already appear in [THM-4527's result note](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_label_map_paley_snarks_20261001.md). The inverse braid (n\mapsto 4n+1) and its 3-adic regular fibers are inherited from the [reverse-tree lane](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_mod6_20260917_reverse_tree_pieces.md); the signed version is $F_\varepsilon(4n+\varepsilon)=F_\varepsilon(n)$ because $3(4n+\varepsilon)+\varepsilon=4(3n+\varepsilon)$.

## 2. The quadratic character and the cycle-word cancellation condition

On odd labels define the real characters by $\chi_{-1}(m)=(-1)^{(m-1)/2}$, $\chi_2(m)=(-1)^{(m^2-1)/8}$, and $\chi_{-2}=\chi_{-1}\chi_2$. From $m=4x+\varepsilon$,
$$
\chi_{-1}(m)=\varepsilon,\qquad \chi_2(m)=(-1)^x,\qquad
\chi_{-2}(m)=\varepsilon(-1)^x.
$$
In particular the two table sheets are exactly the two values of $\chi_2$, while the signs are exactly the two values of $\chi_{-1}$. The $x^2+2y^2$ prime filter $1,3\bmod8$ uses $\chi_{-2}$, not $\chi_2$: for an odd value of that form $x$ is odd and reduction modulo 8 gives $1$ or $3$. The indefinite norm $x^2-2y^2$ instead gives $1$ or $7$. The homogeneous equation $x^2=2y^2$ has no nonzero integer solution by parity descent and does not give the positive norm's prime filter. Any tournament origin for the positive form must establish it on physical path-pair coordinates, not identify these distinct forms by their approximate ratios.

Take a cyclic valuation word $(k_0,\ldots,k_{r-1})$ for fixed sign $\varepsilon$, $k_i\ge1$, with total $K\ge3$. Write its affine fixed-point equation as $Dn_0=\varepsilon C_0$, where
$$
 D=2^K-3^r,\qquad
 C_0=\sum_{j=0}^{r-1}3^{r-1-j}2^{\sum_{i<j}k_i}.
$$
Length $r$ odd gives $D\equiv5\pmod8$; length $r$ even gives $D\equiv7\pmod8$. Thus $\chi_{-2}(|D|)=-1$ is **not** sign-independent: $\chi_{-2}(-1)=-1$. For the **positive $D$ branch** ((3n+1) on positive odd cycles), $D$ has a prime divisor $p\equiv5,7\pmod8$ to odd exponent, and integral (n_0) forces $v_p(C_0)\ge v_p(D)$. Example for the usual fixed point $1$: take its length-two repetition ((2,2)); $D=7=C_0$, showing cancellation does occur and primitivity matters. On the **negative $D$ branch** ((3n-1) on positive odd cycles), $\chi_{-2}(|D|)=+1$, so the same inference is false. Here a different character yields an inert receipt when available: $|D|\equiv3\pmod8$ for odd $r$ and $1\pmod8$ for even $r$; odd $r$ therefore forces a $\chi_2$-inert prime (class (3) or (5)). In either case a rotated word cannot evade the relevant prime check: $2^{k_i}C_{i+1}=3C_i+D$ implies (gcd(D,C_i)) is independent of $i$, since $D$ is odd and not divisible by 3. A primitive $x^2+2y^2$ representation of (C_i) precludes all inert primes, so would rule out that word on the positive $D$ branch. A nonprimitive norm representation by itself does not suffice.

**Transfer principle:** a character projects a return equation to a lane, but an actual integer witness needs a compatible divisibility receipt and an absolute address. The old [Gold/Pell/Collatz half-step](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/gold_pell_collatz_halfstep_20261019.md) supplies an exact finite analogy: the full modular Pell curve has ten points, but only five lift in the positive admissible orbit.

## 3. A noncomplementary mod-eight repair in the literal three-center bank

Use the literal bank from [target address](https://github.com/eliottcassidy2000/math/blob/2aa5cd874debced736c3461a65bd544f41d17f95/05-knowledge/results/chanpak_target_address_20261115/README.md) and the low-block identity from the [QRH residual session](https://github.com/eliottcassidy2000/math/blob/955d1ed6fe8fdd6e278449d2208a186e9cb81f1e/05-knowledge/results/chanpak_qrh_residual_20261115/README.md):
$$
H_j(x,y)=432xy+129x+258y+15[x\text{ odd}]+30[y\text{ odd}]
 +45(p(x)+p(y))+135+2j,qquad j=0,1,2.
$$
Here (p(2x)=p(x)+x+[x\text{ odd}]), $p(2x+1)=p(x)$, $p(0)=0$, $0\le p(x)\le x$. For $g=2^s$, $0\le b<g$,
$$
p(gu+b)=p(u)+A_b u+p(b)+z_b[u\text{ odd}],\quad
A_b=\operatorname{rev}_s(g-1-b),\quad z_b=[b<g/2].
$$
The existing full complement at the low end need not steer all odd residues modulo 8. For complementary appends (x'=2x+b, y'=2y+1-b) with $s=x+y\equiv-1\pmod4$, let $w=p(x)+p(y)$ and (e=[x\text{ odd}]). Both choices update $w$ by $x+e\bmod4$: their two increments differ by $s+1-2x-2e\equiv0\pmod4$. Thus $w-x+e\pmod4$ is invariant. With final $x$ odd the bank gives $H_j\equiv 2(x+y)-x+6+5w+2j\pmod8$, hence one class modulo 4 is forced for fixed prefix and fixed $j$. Waiting for a longer completely complementary suffix is not a solution.

**Four-choice bypass.** Given $u+v$ odd, take $x=32u+b, y=32v+c$, choosing among
$$
(b,c)=(17,6),(19,4),(21,2),(23,0).
$$
Each pair has $b+c=23$, $A_b+A_c=33$, and opposite top bits, of which $b$'s is 1. The checksum constants (p(b)+p(c)) are (12,8,10,2); the two coefficient rows are (14,6,10,2) and (19,27,23,31). Removing the common term (5(p(u)+p(v)+[v\text{ odd}])+2j), the constant residues are (7,1,1,7). The first two linear phases coincide mod 8, and the last two differ from that phase by $5(4v-4u)=4\bmod8$, giving all four odd residues. The common term plus linear phase is even, either by direct parity using $p(n)\equiv0\pmod2$ (each recursion increment is even), or by inspection. This achieves mod 8 **without changing the selected $j$**; the same pair of positive addresses still produces the full three tournaments with spacing +2.

**All-large bracketing with a mod-24 clock.** For every sufficiently large odd $N$, there are two such physical windows, whose selected members bracket $N$ and are each $N\pmod{24}$, at distance $O(N^{1/3})$. The windows can be taken entirely below and entirely above by bracketing their upper and lower members separately. All selected addresses can share the same (x+y).

Proof: Put $U=u+v, d=u-U/2$, $L=\sqrt{N/110592}$, choose dyadic $g\asymp L^{1/3}$, and let $U$ be the smallest integer $\ge L+1$ with $U\equiv-1\pmod g$. First take $j=2N\pmod3$ from $\{0,1,2\}$; this selects the correct residue since all other terms are multiples of 3. Then choose the suffix pair for $N\pmod8$. For $u=g a+r$, $v=g k+(g-1-r)$, $a+k=(U+1)/g-1$. Thus the high-block expansion gives, uniformly in (u,v),
$$
p(u)+p(v)=\frac{g-1}{2}\left(\frac{U+1}{g}-1\right)+O(|d|+g+U/g).
$$
Indeed $A_r+A_{g-1-r}=g-1$, the deviation is bounded by $g|a-k|/2=O(|d|+g)$, and the low residues and high checksums have total bound (O(g+U/g)). Substituting the fixed five-bit suffixes now gives
$$
H_j=F(U)-442368d^2+O(|d|+g+U/g+1),
\quad
F(U)=110592U^2+\tfrac{331821}{2}U+
45\tfrac{g-1}{2}\left(\tfrac{U+1}{g}-1\right).
$$
The coefficient of $U$ is the same for every suffix because both its spatial sum and its reverse-complement sum are fixed. The omitted linear coefficients of $d$ are bounded absolute constants. At (u=\lceil U/2\rceil), this is above (N+4) for large $L$: the leading overshoot is at least (221184L) and tends to dominate the error. At $d=C\sqrt{Lg}+O(1)$ with sufficiently large constant $C$ the count is below (N-4), as $F(U)-N=O(Lg)$. Both endpoints are in the positive (u,v) domain. Along this integer $u$-interval, each adjacent jump, even when the chosen suffix changes, is $O(\sqrt{Lg}+g+U/g)=O(N^{1/3})$, by comparing two expressions with the same (F(U)). First sign crossing of the selected coordinate brackets $N$; the same argument applied separately to the upper and lower members places whole windows on opposite sides. Finally $x+y=32U+23$ holds identically. Proof complete.

**What remains.** Mod 24 removes one local phase obstruction, not the $O(N^{1/3})$ absolute residue. The proof preserves two independent conditions on the same physical choices: prescribed residue and signed proximity. Any future exact repair must preserve the continuation interface as well. The table bridge identifies the failure of an all-gaps proof to do this; the digit bypass supplies a concrete example where changing the representation, while keeping both its spatial sum and its reverse coefficient sum, repairs a phase with no change in the asymptotic absolute cost.

Arithmetic replay: run `python3 mod8_clock_check.py` (also passes with `-O`). It checks 4,999 central blocks, 9,600 fixed-j local clocks, and 192 complete mod-24 bracket controls through 512 bits. Literal tournament realizability uses the inherited bank theorem, not the arithmetic replay.

# Mirror-framed repairs and a quadratic target selector

Two older experiments suggested different ways to retain what a local switch needs. The [Gray chart work](../gray_selector_recursion_20261018e.md) distinguishes independent physical components from individually legal flips, and the [toothpick recovery](../tournament_joint_parity_20261102/toothpick.md) retains exact incidence at the exposed boundary. Applying both ideas to a painted transitive spine gives an exact two-switch construction. With a different exterior, the same construction supplies arithmetic selection rules for prime and $7p$ targets.

## 1. Contiguous switches: the exact cross term

For a binary transitive spine $w$ with a common apex beating the 1-sites, write $H=1+2I(w)$. All odd cycles use the apex. A marked site $i$ and later unmarked site $j$ specify the two extremes of the spine part of a cycle; the intervening sites in it have even cardinality. Hence
$$
 I(w)=\sum_{i<j}2^{\max(0,j-i-2)}w_i(1-w_j).
$$
This is also the familiar $p$-checksum divided by two when the initial zero digits are ignored.

Take two switches supported in disjoint intervals $L,R$ with $L$ preceding $R$. Give their actual digit differences $\delta,\eta$, each in $\{-1,0,1\}$, and let $s$ be the distance from the last site of $L$ to the first site of $R$. Write $A=\sum_{i\in L}\delta_i2^{\max L-i}$, $B=\sum_{j\in R}\eta_j2^{j-\min R}$. Their mixed difference is
$$
 \Delta_L\Delta_R I=-\tfrac12\big(2^{s-1}AB+[s=1]\delta_{\max L}\eta_{\min R}\big).
$$
Indeed the geometric weight factors across the split; the nearest-neighbor weight is 1 rather than $1/2$. The nonzero binary differences $A,B$ cannot vanish, since each comes from two different binary words. So **separated intervals never give independent additive switches**. For adjacent intervals the mixed difference vanishes exactly when $|A|=|B|=1$ and exactly one switch uses a carry or borrow at its least significant site (the least significant digit change has opposite sign to its one-unit address change). In particular there is no three-dimensional additive cube of disjoint contiguous interval switches. This does not forbid interleaved toggles or a marked nontransitive block.

## 2. Interleaving five sites gives a unit move in every mirror frame

The central transitive word is $(u,1-v,0,v,u)$ with $u,v\in\{0,1\}$. Thus $u$ changes the two outer sites simultaneously while $v$ exchanges the two separated inner marks. They interleave and are not governed by the separated-interval obstruction. Let $L$ be any nonempty painted tournament, with ordinary count $g$, its own apex count $n_L$, and weighted start coordinates $S,d$ as in the [ordinal cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md). The frame is $L\to C\to R$ with one apex, where $R$ is the converse of $L$ with complemented marks. The mirror has end coordinates $S,d$. Assume the first SCC of $L$ contains a marked vertex: then the apex reaches that component, and the mirror gives a return from the last component, so every output is strong.

For an arbitrary transitive word $C$ of length $m$, use its apex count $n_C$, start/end values $s_C,t_C$, and sign differences $e_C,f_C$. The transfer gives exactly
$$
 H=2g n_L+(n_C-2)g^2+gS(s_C+t_C)+gd(e_C+f_C)+2^mS^2.
$$
For this five-site family,
$$
n_C=9+2u-6v,\qquad s_C+t_C=37-12v,\qquad e_C+f_C=1.
$$
These follow from the binary addresses: $s_C$ reads the marks from the left in ordinary order, $t_C$ reads the reverse complement, $e_C$ is the last mark, and $f_C$ is one minus the first. Thus the exact two-switch family is
$$
 H_{u,v}=2g n_L+7g^2+37gS+gd+32S^2+2g^2u-(6g^2+12gS)v .
$$
The first switch is a true unit half-step when $g=1$, and an exact square-scaled step for larger $g$. The other switch has a distinct controlled step. At arbitrary odd $g$, $H_{u,v}\equiv32S^2\bmod g$; the resulting counts are coprime to $g$ whenever the start response is.

For a transitive mirror envelope with binary address $X\ge1$, $g=1,S=X,d=X\bmod2,n_L=1+p(X)$. This gives
$$
 H_{u,v}(X)=32X^2+(37-12v)X+9+2u-6v+(X\bmod2)+2p(X).
$$
Every $X$ yields two adjacent odd pairs on $2\lfloor\log_2 X\rfloor+8$ vertices. The two pairs are disjoint from one another and all other $X$'s: their separation at fixed $X$ is $12X+6$, and $0\le p(X)\le X$ gives $H_{0,1}(X+1)>H_{1,0}(X)$. Each branch has an exact monotone target inverse (binary search and literal test), giving $\Theta(\sqrt T)$ distinct counts through $T$ in this fixed family. The step survives **an arbitrary complementary envelope**, not just repetitions of one fixed word. It is not a common forward-profile switch, and its many narrow pairs do not cover an interval.

## 3. Glue the mirror repair to the independent positional bank

The reflection construction also works if the left and right reservoirs are independent members of the [one-sided positional bank](ONESIDED.md), with the right member conversed and complemented. Replace $2g n_L$ by $g(n_L+n_R)$, since their exposed $g,S,d$ are identical; all other terms survive unchanged. With $m$ cells use $g=37^{2m+1}$, $K=3315840$, $W_j=37^{2(m-1-j)}704^{2j}$ and $n_0$ the bank baseline. The exact family is
$$
 H=B_m+gK\sum_{j=0}^{m-1}\zeta_j W_j+2g^2u-(6g^2+12gS)v,\qquad
 \zeta_j\in\{0,1,2\},
$$
where $B_m=2gn_0+7g^2+37gS+gd+32S^2$. The two actual binary choices on each side give any ternary sum digit $\zeta_j$. The successive ratio $(704/37)^2>3$ makes the digits greedy-decodable; test the four central branches and exact integer division by $gK$. This gives $3^m$ distinct counts within each fixed central branch. The gcd of differences across the full family is **exactly $2g$** for $m\ge1$: the positional differences have gcd $gK$, the first central switch contributes $2g^2$, and all other differences are multiples of $2g$. As $37\nmid K$, $\gcd(K,2g)=2$. The bank begins with a marked singleton; its recurrence has $S\equiv16\bmod37$ after the final guard (before it the value is $822\equiv8$, unchanged across subsequent 37-copies since $704\equiv1$). Thus every mirrored output is coprime to $37$, despite the common difference factor $g$. This combines physical placement and repair modules, but still leaves both a large shared residue modulus and a sparse digit set.

## 4. All-marked envelopes expose quadratic norm sheets

Mark every vertex of $L$, write $Q=P_L(2)$ for its optional ordered two-slot count, and let $g=H(L)$. Then $n_L=g,\ S=Q-g,\ d=g$. The first is ordinal insertion of a source. For the second take all optional covers except the one with empty first slot. For the last the optional first-path signed sum at slot weights $(-1,1)$ vanishes: every odd packing block has total weight $(-1)^j+1=0$, since all its block sizes $j$ are odd. Removing the empty first slot gives the stated sign. The general mirrored output is
$$
 H=R Q^2+B gQ+C g^2,\quad
 R=2^{|C_{\rm word}|},\quad A=s_C+t_C,\quad
 B=A-2R,\quad C=n_C+e_C+f_C-A+R,
$$
where $C_{\rm word}$ is the central mask and $C$ in the form is its constant coefficient. Its discriminant is $\mathcal D=A^2-4R(n_C+e_C+f_C)$.

The five-site four-way switch gives $R=32,\ B=-(27+12v),\ C=5+6v+2u$, hence discriminants $89,-167,113,-143$. In particular the sheet $(u,v)=(0,1)$ obeys
$$
(64Q-39g)^2-113g^2=128H,\qquad
x^2-9xg-8g^2=32H,\quad x=32Q-15g.
$$
The second equation is the very same quadratic norm as the characteristic polynomial of the [alternating two-unit transfer](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-5212-the-bottleneck-free-growth-bound-floor-units-compose-by-exact-hankel-laws-and-every-bottleneck-free-unit-word-grows-at-least-1.7699-per-vertex.md), using its root of $x^2-9x-8$. Here the equality is literal algebra; it does not assert the alternating transfer preserves physical $(g,Q)$.

This also shows why merely walking a Pell orbit will not find more hosts indefinitely. For the 113 sheet $Q\ge2g$, so with $a=64Q-39g$ one has $a/g\ge89$. A positive solution of norm $128N$ with physical ratio has
$$
1<\frac{a+\sqrt{113}g}{a-\sqrt{113}g}\le\frac{89+\sqrt{113}}{89-\sqrt{113}}<2.
$$
Multiplication by a positive norm-one unit $\eta>1$ of $\mathbf Z[\sqrt{113}]$ multiplies this ratio by $\eta^2>2$ (its integer trace is at least three). Thus at most one member of an orbit under these positive units can lie in the physical cone. Norm arithmetic is a candidate generator, not a free physical recontact operation.

For an odd prime $p\nmid 2\mathcal D$, $H=p$, or more generally a squarefree target divisible by $p$, requires $(\mathcal D/p)=1$: complete the square,
$$
(2RQ+Bg)^2-\mathcal Dg^2=4RH .
$$
If the discriminant is a nonsquare, both coordinates vanish modulo $p$, forcing $p^2\mid H$. At 37 all four five-site discriminants are nonsquares; by reciprocity the same is true at every prime $p\equiv37\bmod 480342434$. There are infinitely many by [Dirichlet's theorem](https://arxiv.org/abs/0808.1408). This is an obstruction to this four-way mirror family, not to the tournament spectrum.

## 5. Select the central word to clear both the prime gate and mod 3

A useful three-sheet illustration is
$$
8z^2-g^2,\qquad16z^2-3g^2,\qquad32z^2-3g^2,\qquad z=4Q-3g,
$$
with central words 0001001, 00011001 and 001000101 respectively. Their discriminant squareclasses are $2,3,6$; at least one splits at every $p\ne2,3$ because the product of the three characters is $+1$. Individually these give literal prime counts 199, 397 and 797 with $L$ a marked singleton. But the three-sheet test ignores the compulsory physical relation $Q\equiv(-1)^{|L|}g\bmod3$, so it need not clear that constraint at the same target.

The following four small menus do clear the two named gates simultaneously. A row lists three actual central words and the associated signed squareclasses of their discriminants.

| Target lane | Target mod 3 | Central words | Discriminant classes |
|---|---:|---|---|
| $p$ | 1 | 0011, 011001, 001100 | $-39,-23,897$ |
| $p$ | 2 | 01101001, 0, 0000000 | $-1007,-15,15105$ |
| $7p$ | 1 | 01100001, 0110001, 001100 | $-663,-391,897$ |
| $7p$ | 2 | 01001101, 0100110, 001110100 | $-983,-87,85521$ |

In each row the product of the three classes is square; consequently some sheet splits at $p$ as long as $p$ does not divide that row's discriminants. In the $7p$ rows **every** sheet splits at 7. For every word in the indicated row, one of $R+C+B$ and $R+C-B$ equals the indicated residue mod 3. This lets us choose the allowable parity of $|L|$ at the arithmetic-profile stage. These are finite identities in the displayed words, checked by the replay. In the last row $85521=983\cdot87$; in other rows any overlap of prime factors is removed as a square. Some forms are positive definite and some indefinite; the physical source cone is always $Q\ge2g$.

In particular for the stated unramified targets one can choose a sheet and a **formal** primitive integer pair $(g,Q)$ satisfying the target congruence to arbitrary powers of each prime dividing the target, together with the physical mod-three sign and $Q\equiv0\bmod4,\ g$ odd. To see this take $g=1$ for the local calculation, lift a nonzero square root of $\mathcal D+4RN$ at each odd prime dividing $N$, solve for $Q$, and combine with the mod-three sign and parity conditions by CRT. Positive representatives can be chosen. This does **not** say that $g=1$ with an arbitrary such $Q$ is physical: a transitive spine has $Q$ a pure power of two. It does not assert solvability at all other ramified primes either. What has been removed is the quadratic-character obstruction at the target's prime factors together with the mod-three sign; equality and the supply of one physical source carrying the resulting $(g,Q)$ are open.

More generally, a finite set of discriminants removes an infinite *common nonsplit-prime* obstruction exactly if it contains an odd squareclass dependence. The forward direction is the product of characters. If there is no odd dependence, the classes admit a linear functional over $\mathbf F_2$ equal to 1 on all of them. Quadratic reciprocity and CRT realize its values by one coprime progression (including the sign and the prime 2 via mod 8), and Dirichlet gives infinitely many primes. This is the useful piece of the old Gray-chart idea: a compatible set of physical-looking forms is a linear relation in the correct coordinates, not a count of parameters.

Run `python verify_mirror_patches.py` for literal Hamilton counts on small arbitrary mirror envelopes, independent transitive-mask sums at all 1,024 test addresses, and the exact word and character tables.

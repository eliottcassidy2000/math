# What a two-coordinate edit can hide: Pollock shells and Chan–Pak rungs

This extends [Pollock's one-back shell descent](../../hypotheses/HYP-2491-pollock-tetrahedral-defect-pair-descent.md), the [centered moment filter](../pollock_joint_local_moments_and_selector_20261126/README.md), and [one-apex rungs at prescribed gap and residue](../chanpak_apex_flip_local_global_20261122/README.md). In each problem, impose a very small edit and read a second observable rather than just the integer target.

* **Pollock:** adjacent anchors admit representations modulo every modulus with three unchanged residual summands. Nevertheless almost every shell forbids this three-shared-summand edit **anywhere in the exact shell**, even before solving for the shared terms. An exact divisor-and-square test takes the place of the modular lift.
* **Chan–Pak:** on a strong one-apex rung of height difference +2, there is a complete classification of the change in every odd-cycle count. If triangles do not change, **exactly one net 5-cycle is gained**, and all longer odd-cycle counts agree. In particular no such rung is invisible to both triangle and 5-cycle counts. Every permissible profile still realizes every allowable residue, at every sufficiently large even order.

Neither statement establishes the unrestricted conjecture.

## 1. Pollock: a locally free edit fails at almost every exact shell

Set $F(x)=(x^3-x)/6$; for $x\ge1$ these are precisely the nonnegative tetrahedral numbers. Use the canonical shell
$$
N=F(k+1)+r,\qquad 0\le r<\binom{k+2}{2}.
$$
Using the previous anchor $F(k)$ instead raises the remainder by $\Delta_k=\binom{k+1}{2}$. If the two four-term remainders share three indices, then for some $u,v\ge1$ and $S$ a sum of three tetrahedral numbers,
$$
r=F(u)+S,\qquad r+\Delta_k=F(v)+S.
$$
The proposed edit changes the anchor and one other index. Since $r+\Delta_k<(k+1)^2$, it imposes $v=O(k^{2/3})$, and in fact $v-u\gg k^{2/3}$: on that range a single unit increment of $F$ is $O(k^{4/3})$.

### Every finite modulus accepts the shared-three edit

First, two tetrahedral summands already suffice modulo **every** positive integer, with positive indices. For $p=2,3$, $F(6t)=36t^3-t$ is a bijection mod $p^e$ for every $e\ge1$. It is negation mod $p$, and its derivative is a unit, allowing any next digit. For $p>3$, $f(x)=x^3-x$ has
$$
|f(\mathbb F_p)|=\frac{2p+\chi_p(-3)}{3}>p/2.
$$
Indeed $f(x)=f(y)$ with $x\ne y$ iff $x^2+xy+y^2=1$. This conic has $p-\chi_p(-3)$ points, by completing the square and summing the quadratic character (the identity $\sum_t\chi_p(t^2-a)=-1$ for $a\ne0$ follows by factoring $t^2-s^2=a$). If $c$ diagonal points satisfy $3x^2=1$, there are $c$ double-value fibers (one double root and one simple root), and $p-\chi_p(-3)-c$ ordered off-diagonal collisions. A double-value fiber supplies two collisions, a three-point fiber six. Counting fibers gives the displayed size. Every attained value has a simple preimage: a cubic with a repeated root has a residual simple root unless it is a triple root, impossible here. The image of $F$ is equally large, so a translate of its negative intersects it at every desired sum. One simple coordinate lifts the sum to every $p^e$. The same argument covers differences since $F$ is odd. Combine the indices by CRT, using periods $p^{e+1}$ at 2 and 3 and $p^e$ elsewhere, and take positive representatives.

Now choose $u,v$ with $F(v)-F(u)\equiv\Delta_k$, and use two shared summands to supply $r-F(u)$ after subtracting any prescribed third summand. This works for every $r,k,M$, simultaneously across the prime powers of $M$. It says nothing about the size of these indices.

### Exact selector and a density-zero theorem

Put $d=v-u>0, m=u+v, z=2k+1$. The exact edit is equivalent to
$$
 d(3m^2+d^2-4)=3(z^2-1),\qquad
 z^2-dm^2=C_d:=1+(d^3-4d)/3.
$$
The first identity can also be read before multiplying by four:
$d(3u^2+3ud+d^2-1)=3k(k+1)$.
For a particular shell only test divisors $d\mid 3k(k+1)$ within the cube-root scale. Recover $m$ by an exact square test on
$$
\frac{12k(k+1)/d-d^2+4}{3}.
$$
Require $m\equiv d\pmod2, m>d$, and then require $r-F((m-d)/2)$ to be a nonnegative sum of three tetrahedral numbers. This is necessary and sufficient; for existence somewhere in a shell, the shared sum may be zero, so just check $F((m-d)/2)<\binom{k+2}{2}$.

**Theorem.** On $K\le k<2K$, the number of shell indices at which a shared-three edit is possible is $O_\varepsilon(K^{2/3+\varepsilon})$, for every $\varepsilon>0$. In fact the same bound counts eligible pairs $u,v$.

* Here $d,m=O(K^{2/3})$. The solution $d=1$ has $m=z$, hence $u=k$ and is too large for the shell once $K$ is large.
* For square $d>1$, factor $(z-\sqrt d\,m)(z+\sqrt d\,m)=C_d>0$. There are at most $\tau(C_d)$ choices.
* For nonsquare $d>1$, the algebraic integer $\alpha=z+m\sqrt d$ generates an integral ideal of norm $C_d$ in the maximal order of the real quadratic field. There are at most $\tau(C_d)$ ideals of this norm, by unique factorization of ideals and the split/inert/ramified alternatives. Two positive generators of the same ideal and norm differ by a totally positive norm-one unit. Any such unit $\eta>1$ has integer trace $\eta+\eta^{-1}\ge3$, so it scales the positive embedding by at least $(3+\sqrt5)/2$. In the range here $1\le\alpha=O(K)$; at most $O(\log K)$ generators of any ideal fit. This bound is uniform in $d$.
* Now $C_d=O(K^2)$; the divisor bound $\tau(n)=O_\delta(n^\delta)$ with arbitrary positive $\delta$, and $O(K^{2/3})$ values of $d$, give the theorem.

The scan finds just 171 eligible shells out of the 16,384 shells in $[16384,32768$), containing 172 eligible pairs, even with the shared sum allowed to be zero. This is **not** a defect in the earlier one-back strategy: that strategy needs one of the remainders to work, not for both to work using three shared summands. It rules out a tempting way to couple the two certificates. To couple them, the next genuinely different route is two shared summands and the four-coordinate difference equation.

## 2. Chan–Pak: the first nonzero cycle response is at length 3 or 5

Take an oriented transitive spine and one apex. In an apex flip $U0V\to U1V$, the bit records an outgoing apex arc. Let $u$ be the binary value of $U$ and $w$ the binary value of the reverse complement of $V$. For a strong pair the word starts at 1 and ends at 0. Write $c_j$ for the number of directed $j$-cycles up to cyclic rotation, and $\Delta$ for the change in the indicated direction.

A changed cycle must use the flipped arc. Its other apex contact, if it is at distance $i+1$ from that arc's spine vertex, can be completed with any $j-3$ intervening vertices in their unique increasing order. Thus, for $j=3,5,7,\dots$,
$$
\Delta c_j=\sum_{i\ge0}(\operatorname{bit}_i w-\operatorname{bit}_i u)\binom i{j-3}.
$$
All cycles run through the apex, so by [Grinberg–Stanley](https://arxiv.org/abs/2307.05569), $H=1+2\sum_{j\text{ odd},j\ge3}c_j$. Equivalently the inherited one-site rule is $\Delta H=2(\lceil w/2\rceil-\lceil u/2\rceil)$. When $\Delta H=2$, write
$$
u=2s-\epsilon,\qquad w=2(s+1)-\eta,\qquad \epsilon,\eta\in\{0,1\},\quad s\ge1 .
$$

* **No carry:** $s$ even and $u$ even, or $s$ odd and $w$ odd. Only bits 0 and 1 change: $\Delta c_3=1$ and all higher odd-cycle changes vanish.
* **Carry:** if $s$ is even and $u$ odd, take $r=\nu_2(s)$; if $s$ is odd and $w$ even, take $r=\nu_2(s+1)$. In the high bits the increment loses bits $2,\ldots,r$ and gains bit $r+1$. A direct popcount including the low two bits gives $\Delta c_3=1-r$: in the even case $\operatorname{popcount}(s-1)=\operatorname{popcount}(s)-1+r$; in the odd case $\operatorname{popcount}(s+1)=\operatorname{popcount}(s)+1-r$. For $j\ge5$ the two low bits do not contribute, and the sum of binomial coefficients gives
$$
\Delta c_j=\binom{r+1}{j-3}-\binom{r+1}{j-2}.
$$

In particular the triangle-neutral move is exactly $r=1$. It has $\Delta c_5=1$ and no longer odd-cycle changes. The next examples are $(\Delta c_3,\Delta c_5,\Delta c_7)=(-1,2,0),(-2,2,1),(-3,0,4)$ at $r=2,3,4$. The same scalar increment is made by very different changes in the cycle lengths. A broader two-port construction that preserves both short-cycle counts must leave the single-apex model; it will need to use longer cycles or disjoint cycle packs.

### Fix the profile without losing any residue

Use the earlier rung lower-height law
$$
L=2uw+u+[u\text{ odd}]+h(u)+h(w)-1,quad
h(2q+1)=h(q),\quad h(2q)=h(q)+q+[q\text{ odd}],\quad h(0)=1 .
$$
For the no-carry profile choose even $s$, $u=2s,w=2s+1$; then $L\equiv1+3s\pmod4$, giving both odd classes by $s\equiv0,2\pmod4$. For any specified carry $r\ge1$ take $s=2^r(2z+1)$ even, $u=2s-1,w=2s+1$. Here
$$
L\equiv3+h(s)-h(s-1)\pmod4,quad
h(s)-h(s-1)=(2^r-1)(2z+1)+1-z-[z\text{ odd}].
$$
The difference is $z+2-[z\text{ odd}]$ mod 4 for $r=1$, and $z-[z\text{ odd}]$ for $r>1$. Thus the fixed low bits of $s$ can set either odd class while retaining the entire profile.

For target modulus $M$, put $Q=\operatorname{lcm}(4,M)=2^a m$ with $m$ odd. Choose $\ell\ge a$ with $2^\ell=1\pmod m$. Build a shared high prefix $q$ beginning with an arbitrary nonempty word starting in 1 and ending a guard of $\ell$ ones. Follow it by exactly $Q/4-1$ independent equal-length blocks, each $1^{2\ell}$ or $0^\ell1^\ell$. Both blocks preserve $q\bmod m$ and return to $-1\bmod2^a$. The second block adds an extra $D=(2^\ell-1)q+1$ to $h$, with $D\equiv1\bmod m$, $D\equiv2\bmod2^a$, while the first adds zero. Then form $s=2^Kq+s_0$ by a fixed low suffix from the mod-four choice above, using enough bits that both $u,w$ have $q$ as common high prefix. The suffix increments of $h$ depend only on the prefix residue mod $Q$; thus each active block changes $L\bmod Q$ by $2D$, which generates the whole additive subgroup $4\mathbb Z/Q\mathbb Z$. This retains the cycle profile exactly while covering all desired residues. Since the initial prefix can have any positive length, and both flanks have the same bit length, it works at every sufficiently large even total order.

## 3. What transfers

The Pollock sum and the Hamiltonian height both allow rich modular steering of very short edits, but the exact edit reads additional geometry. For Pollock, positivity places the two changing indices below a cube-root ceiling and then a norm equation controls the edit. For a tournament, the ongoing binary carry quantifies how a +2 can trade one cycle length for another, and rules out hiding the edit from both of the first two cycle counters. The next proof targets are different for that reason: Pollock needs a four-coordinate positive inverse (or a defect-pair estimate not based on shared terms); Chan–Pak needs a physical construction that can use the 5-cycle response or go beyond apex-local cancellation **at an exact absolute height**.

Replay: `python3 replay.py`; see `replay.out` for the counts. The independent Hamiltonian subset DP checks both endpoints on 64 small strong pairs; the remaining tests check the derived arithmetic, not the conjectures.

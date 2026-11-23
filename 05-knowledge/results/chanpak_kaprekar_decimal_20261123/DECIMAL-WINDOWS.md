# Even-modulus height balancing and exact two-adic boundary phases

Consider the strong singleton-apex family: the apex sees a transitive chain of length $n$ via a contact word in $O,I$; $O$ means the apex beats that vertex, $I$ means it beats the apex. A word beginning $O$ and ending $I$ is strong. Put $D=H(T),Q=Q_2(T)$.

The [earlier guarded central construction](../chanpak_fork_curvature_20261123/ALL-PRIME-WORDS.md) controls the compatible odd-prime-power observations, with both order parities. The following adds an independent two-adic control for every odd-coefficient read $AQ+BD$ (in particular $D$, the factor-seven fork parent, center and exits). Separately, it characterizes which *pairs* $(D,Q)$ can be achieved at 2-powers by arbitrarily large members of this family.

## A physically additive two-adic phase

Fix $2^k$ and a fixed odd part of a modulus. Use long guarded physical macros $C=O^gI^g,P=OC,S=CI$, with $g\ge k$ and also at least the exponent of the odd component at three when needed. In every macro, the directed-pair expression
$$
H(V)=1+g_2(V)+g_0(V),\qquad
g_\lambda(V)=\sum_{i<j:V_i=O,V_j=I}\lambda^{j-i-1},
$$
contributes
$$
g_2+g_0=(2^a-1)(2^b-1)+1\equiv 2\bmod 2^k
$$
for its $O^aI^b$ pattern ($a,b\ge g$). Across any guarded macro junction every eligible pair has at least $2g$ intervening contacts, so all new terms for $\lambda=2,0$ vanish modulo $2^k$. Thus a guarded word has $H\equiv1+2N_{\rm tile}\bmod2^k$, where $N_{\rm tile}$ now counts **all** its macro tiles, not just copies of $C$. Track this as an extra additive coordinate modulo $2^{k-1}$ in the earlier finite macro group.

The augmented $\mathcal O=PC^{-1},\mathcal I=C^{-1}S$ each have zero macro-count coordinate. The pure center word $Z=C(\mathcal O^g\mathcal I^g)^{-1}$ has length zero and macro-count coordinate one. Its odd-prime-power coordinates lie in the (additive, central) odd nilpotent-channel group. Raise it to an odd common multiple of their orders: this kills the odd effect while translating the two-adic count phase by a unit. Conversely, any desired odd central shift may be isolated from this extra phase by multiplying its exponent by a suitable number congruent to one at the odd prime power and zero modulo $2^{k-1}$. Hence the extra phase and the previously controlled odd reads are independent. When $k=1$ the extra coordinate is trivial, as it should be.

If an arbitrary suffix follows a long guarded correction word, two-adic crossings again vanish: the final $I^g$ of that word separates every eligible pair by at least $g$. In the two balancing modes, long prefixes $O^a$ and $OI^a$ similarly have no two-adic cross effect, by the first $O^g$ guard or the long $I^a$. Their own $g_2+g_0$ reads are respectively zero and $2^a$, both zero modulo $2^k$ for $a\ge k$. The remaining suffix contribution is known and independent of the correction. The finite correction word length has a uniform bound depending only on the modulus: it is selected from the length class and the fixed number of finite central coordinates, not from a bound on suffix length.

## The two-slot phase is carried by the two ends

From the exact two-channel expression,
$$
Q=2^{n+1}+\sum_{i<j:V_i=O,V_j=I}2^{n-j+i}\bigl(3^{j-i-1}+1\bigr).
$$
For chain length sufficiently large in terms of $k$, only contacts within the first and last $k$ positions matter modulo $2^k$, since all other terms have $n-j+i\ge k$. In particular the long-prefix mode and chosen suffix in the balancing proof determine $Q\bmod2^k$ **before** choosing the bounded correction word. For boundary polynomials
$$
A(x)=\sum_{i=0}^{k-1}[V_{i+1}=O]x^i,\qquad
B(x)=\sum_{j=0}^{k-1}[V_{n-j}=I]x^j
$$
one can factor the formula:
$$
Q\equiv 2\bigl(3^{n-2}A(2/3)B(2/3)+A(2)B(2)\bigr)\pmod{2^k}.
$$
Here division by three is understood in the ring, and $n$ can be taken large enough that the indicated end blocks are disjoint. The formula follows by setting the pair gap to $n-2-i-j$ and discarding powers $2^{i+j+1}$ outside the modulus.

For sufficiently large order, the complete possible $Q$ residue classes at modulus $2^k$, with prescribed parity, are exactly
$$
\begin{cases}
Q\equiv4\pmod{\gcd(16,2^k)},& |T|\text{ odd},\\
Q\equiv0\pmod{\gcd(8,2^k)},& |T|\text{ even}.
\end{cases}
$$
Indeed for $n$ even, the mandatory first/last pair contributes $4\bmod16$. All pairs with $r=n-j+i=2$ have odd gap and a factor 4 in $3^{j-i-1}+1$; pairs with $r\ge3$ are also multiples of 16. For $n$ odd, the $r=1$ term has odd gap and is a multiple of 8; all remaining terms are as well. The large-order qualification is needed in the odd-order clause: e.g. the three-cycle has $Q=12$.

To realize every class, take fixed boundary polynomials $A=B=1$ (initial $O$, then only $I$'s in the first block; final $I$, preceded by only $O$'s in the last). For fixed parity, $3^{n-2}$ runs through one coset of $1+8\mathbb Z$ modulo $2^{k-1}$ as $n$ varies in that parity: $9$ generates this subgroup, by $v_2(9^j-1)=3+v_2(j)$ for $j>0$. For very small $k$ the group is understood as trivial. With $n$ even these boundaries therefore give all $4\bmod16$ residues; with $n$ odd they give all $8\bmod16$ residues. In the latter case also use $A=1+x,B=1$: here $A(2/3)\equiv7,A(2)\equiv3\bmod 8$, and $3^{n-2}\equiv3\bmod8$, yielding $Q\equiv0\bmod16$. In each case the affine factorization has odd coefficient after dividing by 2, so it gives the entire stated coset and not just its lowest digit.

Choose such a length phase $n$, place the fixed left boundary followed by $I^a$ for large $a$ a multiple of the odd length modulus, then a guarded correction word and the fixed right boundary. At every odd unit channel the long $I^a$ is the identity, and at every odd nilpotent channel all crossings through it vanish; the bounded prefix affects only a known baseline. The same is true for the two-adic $H$ phase. Choose correction length class to allow the long $a$, solve the odd and two-adic Hamilton reads by their independent central maps, and take $n$ large in the requisite phase. Thus the stated $Q$ image holds jointly with any odd $D\bmod2^k$ and any compatible odd-prime-power targets. In particular there is a complete joint residue description at $10^k$, with arbitrarily large actual strong witnesses of either prescribed parity.

## Balanced scalar-height extension to every modulus

In the [square-root balancing construction](../chanpak_heegner_absolute_targets_20261123/BALANCED-WINDOWS.md), choose the prefix mode, full length $n$, and binary/weighted suffix digits as before. The factorized end formula fixes an explicit $Q\bmod2^k$ independently of the yet-unknown guarded correction. We can now target any feasible scalar read $AQ+BD\bmod2^k$ with $B$ odd simply by solving for the odd value of $D$ and using the independent two-adic phase. The odd part is controlled as before. No change to the height error is needed: the combined finite correction state still has a uniform maximum word length depending only on the modulus.

Thus every fixed odd Hamilton residue modulo $10^k$ (indeed modulo any integer modulus), together with compatible odd continuation targets, has an actual strong witness within $O_M(\sqrt X)$ of every sufficiently large real $X$ at either chosen parity. For the native factor-seven fork parent, center and both exits, the same extension of their base-$3/2$ balancing proof gives error $O_M(X^{1-\frac12\log_3(3/2)})$ at every sufficiently large scale, with any feasible scalar residue. These statements prescribe a scalar even-part residue, not arbitrary independent even-part $D,Q$ while maintaining the same window.

## Why the two-adic pair phase cannot be balanced freely in this family

Fix a parity and take $X=2^n$, along chain lengths $n$ of that parity in any fixed sufficiently long two-adic length phase. In a sufficiently small fixed relative window about $X$, only chain lengths $n$ and $n+2$ can occur among guarded words of the same parity: the universal height bounds are $2^{\ell-2}+1\le H\le2^\ell$ for chain length $\ell$.

At length $n$, $H/2^n$ close enough to 1 forces the first $k$ contacts to be $O$ and the last $k$ to be $I$. An $I$ at initial index $i\le k$, for example, invalidates all $2^{n-i}$ sidecars whose first right-side vertex is $i$; the terminal assertion is its reverse-converse. At length $n+2$, a height close enough to $2^{n}$ instead forces the only $O$ in the first $k$ contacts to be the initial one and the only $I$ in the last $k$ contacts to be the final one. Otherwise there is an extra $O,I$ pair involving an endpoint contributing at least $2^{-k-1}$ in scaled height beyond the mandatory first/last pair. All contributions are nonnegative. One may take a window small enough in terms of $k$ to force both statements.

The first case has $A(2)=B(2)\equiv-1$ and $A(2/3)=B(2/3)\equiv3$ modulo $2^{k-1}$; its forced $Q$ residue is $2(3^n+1)$. The second case has $A=B=1$ and length $n+2$, giving the *same* residue. Whenever the exact image above contains two residues for this parity, prescribe a different $Q$ residue along an infinite fixed phase of $n$. No such word then hits these relative windows with both readings prescribed, even though the exact pair residues do occur at other heights. This obstruction is specific to the guarded singleton-apex family and does not rule out different tournament hosts.

For example, take even tournament order and prescribe $Q\equiv0\pmod{16}$. This is attainable at arbitrarily large orders, but near every sufficiently large $X=2^n$ with $n$ odd the forced boundary residue is $8\pmod{16}$. At odd tournament order the first such obstruction appears modulo $32$, where the admissible residues are $4$ and $20$ and either is ruled out on one of the two phases of even $n\bmod4$.

Replay: `python3 chanpak-height-balanced-decimal-check.py`.

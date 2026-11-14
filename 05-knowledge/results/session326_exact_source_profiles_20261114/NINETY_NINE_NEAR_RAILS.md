# A literal 99/100 near pair with a reversal of exact continuation cost

**Proved:** both strong rails, the complete all-scale intercept screen for transitive sides at odd slope part 99, and opposite cost order after selective-quotient composition. The result supplies exact targets on two infinite families, not a uniform odd-target selector.

## 1. Source lineage and exact meaning of the numbers

[Least resistance and Paley run sums](../paley_least_resistance_runs_one_modular_form_20261014.md) calls a reduction “least resistance” when series/parallel eliminations retain the exact downstream character response. A small coarse answer is not the certificate; the preserved response is. The adjacent [six-vertex terminal](../paley_k6p6_adjoint_weight3_level32_20261014.md) records 3,243 primes through 30,000 for its modular evaluation, typed finite-exact, separately from its proved Legendre reduction. The [128 atlas](../one_twenty_eight_atlas_balance_gold_20261014.md) extends a related scan beyond that boundary. These explain the recovered 30,000 reference, not a count threshold or a proof that physical 99-based targets are prime. The [prime-power bias packet](../prime_power_bias_20261014d.md) uses another 30,000 cutoff in a different test.

This note uses the literal \(99\cdot2^m+23\) target from the earlier prompt and compares \(100\cdot2^m+23\) at the same tournament order. Ratios \(99/100=198/200\) are exact coefficient ratios; replacing the slope pair 99,100 by 198,200 only shifts \(m\) by one and does not double the offset.

## 2. Two literal strong rails at exactly the same order

Use the [one-return-star transfer](../session317_exact_anchor_atlas_20261113/ONE_RETURN_STAR.md) for \(I\to U\to W\) plus apex \(v\), with \(U\to v\) and \(v\) beating a specified nonempty subset \(S\subseteq I\). Its exact law when \(U=TT_\ell\) is
\[
H=\alpha\beta\,2^\ell +h_I d-h_W\lambda,
\]
where \(\alpha,\lambda\) are the guarded two-piece counts on \(I\), \(\beta\) is the prefix-guard count on \(W\), and \(d=H(W+v)\).

| Chart | \(I\), with apex out-cut | \(W\), with vertices beating \(v\) | \(h_I,\alpha,\lambda\) | \(h_W,\beta,d\) |
|---|---|---|---|---|
| 99 | \((10,0,3,6)\), mask 3 | \(TT_4\), second and fourth | \(3,9,4\) | \(1,11,9\) |
| 25 | singleton beaten by \(v\) | \(TT_5\), fourth and fifth | \(1,1,0\) | \(1,25,23\) |

Rows on \(I\) are outgoing bitmasks in zero-based vertex order. The 99 block \(I\) is a directed triangle followed by a singleton; its out-cut from \(v\) includes one triangle vertex and that singleton. The table is checked by direct permutations of the small sides, not by fitting the output rail. In both cases \(v\) can reach all of \(I\), \(I\) reaches \(U,W\), and every \(W\)-vertex can reach one that beats \(v\) via the chain; thus the entire tournament is strong for every nonempty \(U\).

For the 99 chart take \(\ell=m\); for the 25 chart take \(\ell=m+2\). For each \(m\ge1\) we obtain two strong tournaments of order \(m+9\) with
\[
H_{99}=99\cdot2^m+23,\qquad H_{100}=100\cdot2^m+23.
\]
Their ratio exceeds \(99/100\) by exactly \(23/[100(100\cdot2^m+23)]\). The absolute gap is exactly \(2^m\). The 25 chart by itself works also at padded exponents below the aligned range.

A small modular guard comparison shows why nearby slopes are not interchangeable for the earlier prime sieve. For \(m\ge1\), the 99 chart gives a multiple of 5 exactly at \(m\equiv3\bmod4\), while the 100 chart gives none. The 99 chart gives no multiple of 11, while the 100 chart gives one exactly at \(m\equiv5\bmod10\). These follow from the exact orders of 2 modulo 5 and 11. They say nothing about primality at surviving exponents.

## 3. All-scale obstruction before paying for a cyclic side

If both \(I\) and \(W\) are transitive, the complete possible intercepts at a specified odd part of the slope are explicitly finite. Let \(L\) denote [the exact binary defect](../session317_exact_anchor_atlas_20261113/ONE_RETURN_STAR.md), with
\[
L(2x)=L(x),\qquad L(2x+1)=L(x)+x+(x\bmod2),\qquad L(0)=0.
\]
On \(I\), its membership mask supplies any \(\alpha>0\) and forces \(\lambda=L(\alpha)\).

On \(W\), write \(A\) for the positions beating the apex, using 0-based positions in chain order, and \(c=\sum_{i\in A}2^i\). Then \(\beta=1+c\): a valid nonempty prefix ending at such a position \(i\) has \(2^i\) choices. To count \(H(W+v)\), discard from those \(\beta\) covers the cases in which the nonempty suffix starts in \(A\). For a nonempty prefix those failures are exactly \(L(c)\), by reversing the chain and exchanging the pieces; for an empty prefix there is the additional failure \(c\bmod2\). Hence
\[
d(\beta)=\beta-L(\beta-1)-((\beta-1)\bmod2).
\]
The displayed recurrence gives \(d(2\beta)=d(\beta)\). Thus the intercept depends only on the odd parts of \(\alpha,\beta\). Every division of the odd part \(P=\operatorname{oddpart}(\alpha\beta)\) into ordered positive factors \(x,y\) is allowed; the associated intercept is \(d(y)-L(x)\).

For \(P=99\) these intercepts are \(73,31,5,1,-13,-75\), for factor pairs \((1,99),(3,33),(9,11),(11,9),(33,3),(99,1)\) respectively. None equals 23 modulo 99. Therefore no chart **with odd slope part 99 and both side blocks transitive** hits any target \(99\cdot2^m+23\), even at an isolated exponent or with asynchronous padding. No assertion is made about other odd slope parts hitting individual such targets. For \(P=25\) the intercept list is \(23,3,-13\). Thus the successful 99 construction genuinely needs something beyond the all-transitive side mechanism; the single directed triangle in \(I\) suffices here.

## 4. Exact \(Q(2),Q(3)\) and native cost inversion

Use empty-slot counts \(Q_T(k)\). For a fixed outer star and variable transitive \(U=TT_\ell\), compress every visit to \(U\) in each slot to one token. There can be at most \(k+1\) such tokens across all slots: no slot revisits \(U\) without using \(v\), and \(v\) belongs to only one slot and occurs only once. The finite coefficient \(c_j\) counts compressed covers with \(j\) occupied \(U\) runs. Each admits exactly \(j!{\ell\brace j}\) fillings, by onto choices of the \(U\)-vertices and their forced forward order inside a run. Thus
\[
Q_T(k)=\sum_{j=1}^{k+1}c_j\,j!{\ell\brace j}\qquad(\ell\ge1).
\]
The first \(k+1\) sizes recover the full coefficient vector by an invertible triangular system with diagonal \(j!\). This is an exact finite basis argument, not extrapolation from small graphs. The replay evaluates those sizes by direct path subset DP and ordered slot convolution. In the original exponent \(\ell\) the results are:

| chart | \(Q(2)\) | \(Q(3)\) |
|---|---|---|
| 99 | \(3630\cdot3^\ell+1696\cdot2^\ell-330\) | \(46155\cdot4^\ell+35829\cdot3^\ell-7755\cdot2^\ell\) |
| 25 | \(334\cdot3^\ell+592\cdot2^\ell-58\) | \(2073\cdot4^\ell+5373\cdot3^\ell-609\cdot2^\ell\) |

In the same-order comparison \(\ell=m\) versus \(m+2\), both differences \(Q_{99}-Q_{100}\) are strictly positive:
\[
\begin{aligned}
\Delta Q(2)&=624\cdot3^m-672\cdot2^m-272,\\
\Delta Q(3)&=12987\cdot4^m-12528\cdot3^m-5319\cdot2^m.
\end{aligned}
\]
They are positive at \(m=1\), and after dividing by the largest power the negative terms strictly decrease in magnitude. Meanwhile \(\Delta H=-2^m<0\).

This is also an exact **composition** reversal. In either of the [selective quotients](../session324_busch_inverse_and_gap_structure_20261113/COMMUTATIVE_INVERSE.md), let \(A=-3H+6Q(2)+\gamma Q(3)\), with \(\gamma=2\) or 6. Then \(\Delta A>0\) term by term. Plug the two aligned charts in as sources with any identical nonempty reservoir \(V\). The total difference
\[
\Delta A\,Q_V(2)-\Delta H\,H(V)
\]
is strictly positive, although the uncomposed 99 sources had smaller \(H\). The source vertices are also equal in number, so this is not a padding artifact. Consequently choosing the locally cheaper Hamiltonian count is provably the *more expensive* choice after this native continuation. A least-resistance construction needs a vector order depending on the next operation, rather than a greedy order of isolated counts.

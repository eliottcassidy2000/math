# Exact prime-power clocks on two marked Euclid lower witnesses

This extends the two fixed contacts of [BACKCHAIN-LOWER](../chanpak_decimal_kernel_20261122/BACKCHAIN-LOWER.md) without asserting an absolute upper neighbor. The lower tournament in each branch still has exactly \(E_{p_j}-1\), where \(p_j=1+60\cdot2^j\). This note classifies the previously reported small simultaneous phases exactly, lifts both branches to every depth at 13, and gives clock incompatibilities, including with 11. It does not assume that \(2^{p_j}-1\) is prime.

## 1. The two marked contacts have exact affine width laws

Retain the same prefix registers \(F=35183298379775,\ G=10504481267175\) at exponent 61. For the first contact use suffix widths \(15,p_j\), and flip from its initial boundary to ordinary vertex 53, measured from that boundary. Its growing tail beyond that contacted vertex has width \(w=p_j-37\). For the second use suffix widths \(p_j,15\) and contact 56, with growing width \(w=p_j-56\). The exact changes of the two suffix reads \(k,t\) entering the inherited marked-join formula are
\[
\begin{array}{c|cc}
 &\Delta k&\Delta t\\
53&4503530907893760+9836602018824134393855\,2^w&
4503599627370496+9836602018824134393856\,2^w\\
56&1180627649514430267392-32769\,2^w&
1180627649514430267392.
\end{array}
\]
Here \(k\) is the Hamilton count of the suffix including its initial boundary and \(t\) is the inherited separator count with its designated endpoint restriction.

For the exact width argument use the three fixed boundary vertices as special positions in a spanning path and order them. The intervening ordinary vertices are transitive and distribute into at most four ascending slots. The growing block is entirely beyond the fixed contacted vertex and lies between two neighboring boundary positions; its vertices have equal external incidence to each fixed vertex. Once all fixed assignments are chosen, permitting or forbidding each slot to receive some of this growing block is a count of onto assignments into a fixed subset of the four slots. Inclusion–exclusion expresses each read as a linear combination of \(1^w,2^w,3^w,4^w\) for \(w>0\), where the \(0^w\) term is zero. The endpoint restriction on \(t\) changes only which slots may be occupied, since all vertices of the growing block have the same tag. The companion independent four-width check, which allocates ordinary vertices by endpoint masks, verifies the displayed coefficients; the distinct four bases give an invertible Vandermonde system. In particular no asymptotic assumption about discarded higher bases is used.

Set \(A=2^{60},N=2^j\), and \(P_j=(A^N-1)/(A-1)\). Along the Euclid transport, the prefix registers are multiplied by \(P_j\) and the growing width is \(60N-s\), where \(s=36,55\), respectively. The upper-minus-lower change \(\Delta_j\) from the same physical flip therefore satisfies the single formula
\[
 \Delta_j=P_j\bigl(C+D\,2^{-s}A^N\bigr).
\]
The integer coefficients, obtained as \(F\Delta t+G(\Delta k-\Delta t)\), are:

| contact | \(s\) | \(C\) | \(D\) |
|---|---:|---:|---:|
| 53 | 36 | 158450767610363602280225177600 | 346084103871386661670010152033395225 |
| 56 | 55 | 41538374868278620992215173614796800 | -344221346644057575 |

The formula uses rational notation with dyadic denominator only; at odd primes it is simply multiplication by a unit. Each \(\Delta_j\) is an integer keeping the same literal reversal and the same exact lower.

## 2. Exact simultaneous phases and every 13-adic depth

The full phase matches for the previously marked modulus occur exactly at \(j=0\bmod132\) for contact 53, and exactly at \(j=0\bmod84\) for contact 56. The individual rules before taking intersections are:

| contact | prime | allowed \(j\) |
|---|---:|---|
| 53 | 3 | \(0\bmod2\) |
| 53 | 13 | \(0\bmod12\) |
| 53 | 19 | \(0\bmod2\) |
| 53 | 47 | \(0\bmod11\) |
| 56 | 13 | \(0\bmod12\) |
| 56 | 59 | \(0,3\bmod28\) |
| 56 | 71 | \(0\bmod3\) |

These exact small orbits can be replayed by the elementary update \((x,P)\mapsto(x^2,P(x+1))\), starting at \((A,1)\). The condition is \(P(C+D2^{-s}x)=2\). For contact 56 the nonzero residue 3 at 59 is incompatible with the class \(j=0\bmod12\), giving the stated intersection. Prime 2 is automatic for every tournament count.

On each simultaneous phase progression \(j=d k\), \(d=132\) or 84, there is a unique class of \(k\bmod13^{r-1}\) with \(\Delta_j=2\bmod13^r\) for every \(r\ge1\). In fact \(v_{13}(A-1)=1\) since \(2^{12}-1=13\cdot315\) and \(60=5\cdot12\). As \(d/12\) is prime to 13, \(k\mapsto 2^{dk}\) ranges bijectively over the class \(N=1\bmod13\) modulo \(13^r\); this follows by the order (or LTE) formula. If \(N\) changes by \(13^\ell u\), the geometric factor \((A^N-1)/(A-1)\) changes by \(13^\ell u\bmod13^{\ell+1}\), while the affine factor \(C+D2^{-s}A^N\) is unchanged at this precision. The affine factor is \(2\bmod13\) on the desired class, so each candidate has precisely one next digit. The first lifts in \(j\)-coordinates are
\[
\begin{array}{c|cc}
 &r=2&r=3\\
53&528\bmod1716&3960\bmod22308\\
56&1008\bmod1092&4284\bmod14196.
\end{array}
\]
Thus **there are exact lower tournaments with arbitrarily deep 13-adic upper agreement while retaining every other prime in the marked modulus**. Precision requires varying the exponent and the host. It does not identify one host whose upper count is exactly \(E_{p_j}+1\).

## 3. Formal roots that the exponent clock cannot execute

When \(A-1\) is a unit, finding an upper match modulo \(\ell^2\) is equivalent to solving the quadratic
\[
 (x-1)(C+D2^{-s}x)=2(A-1)\pmod{\ell^2},
 \qquad x=A^{2^j}.
\]
The first equation by itself does have roots at all the four remaining odd primes listed below. The obstruction is that those roots cannot be reached by the compositional exponent clock. The middle column lists all residues \(N=2^j\bmod \operatorname{ord}_{\ell^2}(A)\) whose \(A\)-powers are roots of the quadratic. The last column gives the respective clock constraints on \(j\) from the two prime factors \(m,\ell\) of that order; a dash means the power is not in the subgroup generated by 2.

| contact, \(\ell\) | allowed \(N\) in the formal \(A\)-clock | \(j\) clocks modulo the orders at \(m,\ell\) |
|---|---|---|
| 53, 19 | \(34\bmod57\) | \(0\bmod2,\ 11\bmod18\) |
| 53, 47 | \(868,1013\bmod1081\) | \((-, -),\ (0\bmod11,-)\) |
| 56, 59 | \(182,929\bmod1711\) | \((3\bmod28,6\bmod58),\ (0\bmod28,27\bmod58)\) |
| 56, 71 | \(183,444\bmod497\) | \((0\bmod3,-),\ (-,33\bmod35)\) |

Every row fails either a subgroup test or a common-clock parity test. Thus these four primes cannot lift to their square at any \(j\), not just at the known simultaneous progressions. At 3 for contact 53 the square does match at \(j=2\bmod6\), but 13 requires \(j=0\bmod12\), so they cannot be combined. On the phase progressions the other original primes therefore all have precisely valuation one on \(\Delta_j-2\).

There is a sharp 11/13 comparison already at first depth. Contact 53 matches 11 exactly at \(j=1\bmod10\); it matches 13 exactly at \(j=0\bmod12\). The classes have opposite parity, so no one Euclid exponent in this contact family matches both. Separately the 11 branch lifts to every 11-adic depth: \(v_{11}(A-1)=1\) and \(2^{10}\) generates the principal units, and in the class \(N=2\bmod11\) the affine factor in the normalized geometric formula is \(1\bmod11\), so every next digit lifts uniquely. At \(11^2\) the class is \(j=1\bmod110\), and at \(11^3\) it is \(j=881\bmod1210\). Contact 56 never matches 11. This arithmetic is about the upper contact error, distinct from divisibility of the exact lower by 11 or 13.

The comparison to the [PRIMEGAME arithmetic-vs-execution distinction](../metagraph_primegame_20261102/PRIMEGAME_FANO.md) is precise but limited: a formal solution of the quadratic does not imply an exponent \(j\) respecting both clocks, just as a formal return multiset need not admit a priority-legal execution. No equivalence between the two systems or new PRIMEGAME orbit is asserted.

# Which signed Collatz steps are square-sum edges?

For positive odd integers use the two accelerated maps
\[
 U_e(x)=(3x+e)/2^j,\qquad e\in\{1,-1\},\quad j=v_2(3x+e).
\]
Call an edge *square-sum* if its actual endpoint values add to an integer square. The [earlier positive-sign calculation](../chanpak_middle_out_square_sum_20261122/SQUARE-SUM.md) leaves a natural possibility open. A positive-sign square edge goes from split to inert modulo four. A negative-sign edge can go from inert to split. Can alternating signs concatenate them?

## Classification and a mod-three obstruction

**PROVED.** For a fixed sign and division exponent, the following table gives exactly the cases admitting square-sum edges. Each admissible case occurs infinitely often.

| Sign | Admissible exponents | Classes of input and output modulo 4 |
|---|---|---|
| \(+1\) | Every even \(j\ge2\) | \(1\to3\) |
| \(-1\) | Odd \(j\ge3\) for which every prime factor of \(2^j+3\) is \(1\) or \(3\bmod8\) | \(3\to1\) |

Write \(A=2^j+3\). If \(x+y=m^2\), then \(m\) is even and
\[
 x=\frac{2^jm^2-e}{A},\qquad y=\frac{3m^2+e}{A}.
\]
The congruence \(2^jm^2\equiv e\pmod A\) suffices for both quotients to be integral. Taking sufficiently large even roots gives positive odd endpoints, so \(3x+e=2^jy\) has exact valuation \(j\).

For the positive sign and even \(j\), take \(m\equiv 2^{-j/2}\pmod A\). At odd \(j\) the required Jacobi symbol \((2/A)\) is negative since \(A\equiv5\bmod8\) at \(j=1\) and \(A\equiv3\bmod8\) afterward. For the negative sign at even \(j\), the congruence would make \(-1\) a square modulo \(A\equiv3\bmod4\). At odd \(j\) it is equivalent to \(-2\) being a square modulo \(A\), after multiplying by the square of a unit or its inverse. The case \(A=5\) fails. For odd \(j\ge3\), quadratic reciprocity gives the displayed condition at each prime. A nonzero root lifts through every power of its odd prime, and CRT combines the roots and the prescribed even parity. For instance guards 3 and 7 work for the negative sign, while guard 5 fails because \(2^5+3=35\).

In every successful case \(j\ge2\), the formulas reduced modulo four give \(x\equiv e\), \(y\equiv-e\). More importantly, the sign and exponent satisfy \(e=(-1)^j\) in both rows. Since \(2^j y\equiv e\pmod3\), **every square-sum edge ends at \(y\equiv1\pmod3\)**. Any possible next square-sum edge would likewise end at 1 modulo three, making its endpoint sum 2 modulo three. Hence **no two such steps can occur consecutively under any sign schedule**, including signs chosen from the history.

There is also a structural description. At any odd input, exactly one sign has valuation at least two, so the directed graph of all these edges has outdegree at most one. Sources and sinks are disjoint, by the mod-three argument. Every nontrivial weak component is therefore a star directed toward its center. All incoming signs at any given center coincide, determined by the center modulo four.

## Infinitely many centers have two distinct incoming edges

**PROVED.** There are infinitely many centers with simultaneous positive-sign square edges of guards 4 and 8. The first example needed here has center 427 and leaves 2277 and 36437. The square roots of their endpoint sums are respectively 52 and 192.

More generally suppose distinct guards have \(A=2^j+3\), \(B=2^k+3\), and a common center with sign \(e\) and positive even roots \(u,v\). The common-center condition is
\[
     A v^2-Bu^2=e(B-A)/3 .
\]
The identity follows by cross multiplying the formulas for the center. If \(AB\) is nonsquare, a single example gives infinitely many at those same guards and sign. Choose a positive nontrivial Pell unit \(t+w\sqrt{AB}\) of norm one. The integer matrix taking
\[
       (v,u)\longmapsto(tv+Bwu,\;Awv+tu)
\]
has determinant one and preserves the displayed form. Some positive power acts as the identity modulo \(2AB\), since this group of invertible matrices is finite. Its iterates preserve positivity, evenness and the integrality of the common center, and give unbounded roots. The center stays odd and the two guards remain exact. For the explicit example \(A=19\), \(B=259\) and \(AB=7\cdot19\cdot37\) is nonsquare.

## A gap of two, the forbidden height 7, and perfect midpoints

**PROVED.** The only positive odd pairs that are both square-sum signed Collatz edges and two units apart are \(3\to1\) and \(9\to7\). Indeed every square edge has \(j\ge2\) and strictly decreases (the only possible equality here, \(1\to1\), has sum 2). Writing the input as \(y+2\) gives \((2^j-3)y=6+e\). At \(j=2\), the candidates are \(9\to7\) and \(7\to5\), of which only the first sums to a square. At \(j=3\) only \(3\to1\) is integral. Larger denominators exceed 7. Direct enumeration on the first 36 labels gives exactly these same two edges even if the gap condition is removed.

Only the pair of counts 1 and 3 survives as a physical tournament reversal: flip the long arc in the transitive triangle. The global impossibility of count 7 removes the other arithmetic edge. For 21, its input sign must be positive, giving \(21\to1\) of sum 22. An output of 21 is already impossible modulo three. The two exceptional tournament counts thus play different roles in this graph. Their global exclusion has its own tournament proof. [THM-1370, all-order spectrum exclusions](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md)

At a general even midpoint \(M\), the pair \(M-1,M+1\) is square-sum exactly when \(M=2r^2\). For every Euclid midpoint \(E_p=2^{p-1}(2^p-1)\), \(p\ge2\), the odd part of \(2E_p\) is \(3\bmod4\). Thus even a successful one-arc crossing of \(E_{61}\) must fall outside this direct square-sum encoding. Neighboring midpoint pairs can still be linked horizontally through square sums. The finite construction on 36 midpoints does exactly this. [Tournament lift](LIFT.md) · [Current exponent-61 near witnesses](../chanpak_odd_square_spines_20261122/NEAR-WITNESSES.md)

# Deeper on the same terrace: a tribonacci corridor with an exact inverse

On the [first strong triangle terrace](TERRACE.md), take two positive gap lengths $a,b$ with $q$ zero gaps between them. A zero gap means that neighboring reverse-chain anchors are adjacent in rank. The resulting whole count is
$$
 H=e_q+f_q(2^a+2^b)+g_q2^{a+b},
$$
where the three coefficient sequences all obey the tribonacci recurrence. Their first three rows $(e_q,f_q,g_q)$ are $(1,1,1),(1,2,2),(3,3,4)$, at $q=0,1,2$. In particular the product and shifted product are the first two cases of one interval compiler.

To prove it, after the first positive gap the four transfer states $(E,A,O,V)$ (closed, open with positive, open pure-zero odd, open pure-zero even) are $(1,2^a,0,0)$. Across a zero gap,
$$
 (E,A,O,V)\mapsto(E+A+2V,\ A,\ E+V,\ O).
$$
The positive component $A$ remains fixed through the corridor. Put $\widetilde O=O+A/2$ and $\widetilde V=V+A/2$. Then $(E,\widetilde O,\widetilde V)$ follows the homogeneous map
$$
      (E,\widetilde O,\widetilde V)\mapsto
      (E+2\widetilde V,\ E+\widetilde V,\ \widetilde O),
$$
whose characteristic polynomial is $x^3-x^2-x-1$. Appending the final gap gives
$H=(E+2\widetilde V)+2^b(E+\widetilde O+\widetilde V)$. Thus all coefficients obey the recurrence; reflection of the gaps makes the linear coefficients in $2^a,2^b$ equal. The first three displayed rows follow directly from the transfer.

**Exact prescribed-target inverse.** For target $N$, iterate $q$ only while $e_q\le N$ (the sequence grows at a tribonacci rate). Iterate $s=a+b\ge2$ only while $2^s<N$. At this $q,s$ the prospective sum of the two powers is
$$
       D=\frac{N-e_q-g_q2^s}{f_q}.
$$
Test integrality and whether $t^2-Dt+2^s$ has two power-of-two roots at least 2. This uses $O((\log N)^2)$ exact arithmetic tests and recovers every word in the family, up to reflection. Neither the size nor the zero run needs to be specified in advance.

For example $q=3$ gives $(e,f,g)=(5,6,7)$; $a=2,b=4$ gives $H=5+6(4+16)+7(64)=573$, this time on the minimum-triangle terrace at order **12**, not the older order-ten hub-free intercyclic class or the more deeply overlapping word of the same count. The total 53 has two solutions even inside this restricted corridor. The decoder rejects 385 in this restricted family; this is not a nonrealization statement for the full terrace or for tournaments generally.

Run **python zero_corridor.py** to check the recurrence and the literal interval inverses.

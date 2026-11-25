# A single fixed arc with the exact Euclid midpoint at every exponent from four onward

**PROVED.** Put \(E_p=2^{p-1}(2^p-1)\) for any integer \(p\ge4\), and \(X=2^{p-3}\). There is a single arc reversal between two strong tournaments on \(2p\) vertices whose counts are exactly
\[
 E_p-D_p,\quad E_p+D_p,\qquad D_p=6X^2+2X+1.
\]
In particular \(E_{61}\equiv1\bmod3\) is the **exact** midpoint of one physical flip on 122 vertices. Its lower and upper counts are respectively
\[
 2159995493150488291739304216045289471,\qquad
 3156916489989175197570081015862394881 .
\]
The half-gap here is \(498460498419343452915388399908552705\). The desired half-gap of one is a further simultaneous integer constraint.

Start with outgoing masks \((246,252,152,241,96,196,132,16)\) on vertices zero through seven. Replace 0 and 7 by increasing chains of widths \(m,k\ge1\), respectively, keeping every cross edge in the direction of its quotient edge. Reverse \(1\to5\) to obtain the second tournament. Ignoring the chain-internal edges, write \(U=2^m,V=2^k\). The two counts and their midpoint are
\[
  26UV-U-5V-1,\qquad 38UV-U-V+1,\qquad 32UV-U-3V .
\]
To see this directly, consider all spanning words in the quotient; two distinct successive visits to a port use distinct runs of its chain. Every directed triangle at vertex 0 contains a common ordinary vertex, as does every directed triangle at 7. In either orientation a return to the port must use such a vertex, since somewhere it must go from an outneighbor to an inneighbor of the port. Thus at most two visits are possible, even if both chains grow. The counts of the four legal word types in visit order \((1,1),(2,1),(1,2),(2,2)\) are respectively \((91,51,47,26)\) and \((149,75,75,38)\). A single visit has one filling, and two visits have \(U-2\) or \(V-2\) fillings, by assigning the chain labels to two nonempty ordered slots. Expansion proves the formulas. Both quotient orientations are strong, and a cycle through each port preserves strong connectivity under substitution. A joint small-width literal check appears in the [certificate replay](center_check.py).

Take the single width choice \(m=k=p-3\). The midpoint is \(32X^2-4X=E_p\) and the half-gap is \(6X^2+2X+1\). The two equal chains cancel all additional powers in the midpoint. Their difference retains a nonzero quadratic coefficient. The [three-port search](CENTER-THREE.md) quantifies how often exactly this distinction persists in larger descriptions.

The [two-port ray search](center_two.cpp) uses the earlier return-cap test, subset recurrence and the exact leading-coefficient selector. For a sum of endpoint profiles \(AUV+BU+CV+D\), put \(m=p-1+i,\ k=p-1+j\). The equations for an exact Euclid midpoint for every such exponent are \(D=0,\ A2^{i+j}=4,\ B2^i+C2^j=-2\). Valuation of the last two summands gives at most three candidate shifts. The search supplies positive certificates and makes no global claim about all tournaments.

# The five Heegner registers, nine vertices, and exact Ramanujan targets

**PROVED scoped / FINITE-EXACT.** These are separate maps. The CM approximation, the Ramanujan–Nagell equation, and the order-nine arc representation have different inputs and different lost information.

## Five prime-discriminant pairs and 36 arc coordinates

For the odd Euler lucky primes \(p=3,5,11,17,41\), the principal norm form is
\[
N\!\left(x+\frac{1+\sqrt{-(4p-1)}}2\right)=x^2+x+p.
\]
Their class-number-one radicands are \(11,19,43,67,163\). The remaining Euler lucky prime \(p=2\) maps to 7; the Heegner radicand 3 corresponds to the nonprime parameter 1. The prior [Euler horizon](../../../07-reflections/heegner-prime-horizon-tournament-s650.md) records the exact forced square at \(x=p-1\), and the [Legendre tournament theorem](../../../01-canon/theorems/THM-5125-residue-tournament-model-heegner-domination-and-cycle-gas-identities.md) gives a separate, literal seven-prime QR tournament. Neither theorem identifies class number with a tournament's ordinary count. See [Stark's class-number-one theorem](https://doi.org/10.1307/mmj/1028999653).

The same primes give the continuation register lengths \(d_p=(p-1)/2=\{1,2,5,8,20\}\) from [base-prime slot observations](../chanpak_fork_curvature_20261123/CONTINUATION-ORBITS.md), equivalently the quotient sizes \(|\mathbb F_p^\times/\{\pm1\}|\). They sum to 36, the edge count of a nine-vertex tournament. Here is an actual marked linear decomposition with the same sizes. Partition the labels into cyclically marked blocks of size four and five, and work in the 36-dimensional signed arc space over \(\mathbb Q\).
- The six internal arcs of the 4-block give one rotation-fixed line (oriented rim sum) plus a five-dimensional complement. Opposite chords reverse under the half-turn and give no fixed line.
- The ten internal arcs of the 5-block give two rotation-fixed lines (distance-one and distance-two rim sums) plus an eight-dimensional complement.
- The twenty cross arcs form the remaining block, a single length-20 orbit under simultaneous rotation.

The integer direct-sum lattice has index \(4\cdot5\cdot5=100\): on each untwisted rim orbit of length \(m\), constants and integer zero-sum vectors together leave just a sum residue mod \(m\); all opposite chords already lie in the complement. These cyclic actions are labelled relabellings, not tournament automorphisms. [The companion](NINE-REGISTERS.md) proves that the actual prime-register doubling permutation is not a native relabelling of these 36 edges, and tests exact order-nine constructions.

There is also an information bound. The [joint synthesis theorem](../chanpak_fork_curvature_20261123/ALL-PRIME-WORDS.md) produces all \(3^1 5^2 11^5 17^8 41^{20}\) continuation vectors, about \(2^{163.38}\), at either parity. Every unlabelled order-\(n\) tournament has a Hamilton path, so pinning its vertices along one such path leaves only \((n-1)(n-2)/2\) free bits. Summing that bound through order 19 gives less than \(2^{154}\). Thus a bounded-order universal selector for this joint register must allow order at least 20 (or 21 if odd order is required). There is no contradiction with a scalar CRT check at order nine: the [split bound](NINE-REGISTERS.md#4-prime-field-capacity-and-a-genuine-exact-lift-on-order-nine) gives \(H\le9450\), and \(2\cdot11\cdot17\cdot41=15334\) exceeds that bound. At order nine, three matching residues and oddness force the exact scalar *if a witness already exists there*.

The [five ternary exits](../chanpak_six_five_boundary_atlas_20261123/README.md) are another instance of a five-part continuation, but their weights \(1/3,1/9,1/9,1/9,1/3\) are prefix-code masses, not these vector ranks. Likewise the [five guarded readings](../chanpak_guarded_ternary_packing_20261123/README.md) are measurements of one apex context, not independent copies of the five prime registers. The new literal bridge is the marked 4+5 edge split and its representation obstruction, not a bijection among all five-sets. The earlier [36-dimensional two-apex recurrence](../balanced_ternary_returns_20261103/interior/check.json) is a third object, a Cayley–Hamilton state bound rather than an arc lattice; equality of dimensions is not an intertwiner.

## An exact family through the nine-vertex square \(11^2\)

Put a single apex over a transitive chain; \(O\) means it beats that vertex and \(I\) means it loses. For \(a,b,c\ge1\), take contact word \(OI^aOI^bOI^c\). Every path segment avoiding the apex follows increasing chain order, and the [ordinary apex count](BALANCED-WINDOWS.md) weights an \(O_i,I_j\) pair by \(2^{j-i-1}\), with one more for adjacent \(OI\)'s and a baseline 1. Writing \(x=2^a,y=2^b,z=2^c\) and summing the three geometric runs gives
\[
 H=4xyz+(2y+1)(z-x)+1-y .
\]
Thus equal flank lengths cancel the entire linear correction, leaving the literal strong family
\[
 H=2^{\,2a+b+2}-(2^b-1),\qquad a,b\ge1 .
\]
This gives an absolute inverse on the stated family: for input \(N>1\), take \(b=v_2(N-1)\). Accept exactly when \(1+(N-1)/2^b\) is \(4^r\) for \(r\ge2\), and take \(a=r-1\). The number of accepted heights through \(X\) is only \(O((\log X)^2)\); this is not general height coverage.

At \(b=3\) these are exactly \(2^q-7\) for odd \(q\ge7\). The [Ramanujan–Nagell theorem](https://projecteuclid.org/journals/arkiv-for-matematik/volume-4/issue-5/The-diophantine-equation-2nx27/10.1007/BF02591625.pdf) then gives exactly two square outputs of this rail: \(11^2\) at \(q=7\) on **nine vertices** (\(OIOIIIOI\), eight contacts plus apex), and \(181^2\) at \(q=15\) on 17. The [seven-core signed shell](MERSENNE-RAIL.md) removes the odd-exponent restriction: its exact read is \(512Q_2-7H\), giving \(2^q-7\) on \(q-1\) vertices for every \(q\ge10\) by a transitive filler. Small cases complete the range: every \(q\ge4\) is attained strongly. In particular the shell also gives \(181^2\) on 14 vertices. The other square values \(1,9,25\) occur separately, and the replay uses an ordinary singleton or the apex words \(OIII,OIIOOI\). Of Heegner radicands only 7 is \(7\bmod8\), which is necessary for \(x^2+d=2^q\) with odd \(x\), \(q\ge3\). This is a genuine dyadic norm connection, distinct from the CM near-integer below.

## Exact counts straddling Ramanujan's constant

For \(\tau=(1+i\sqrt{163})/2\), the exact value \(j(\tau)=-640320^3\) and \(j=q^{-1}+744+196884q+\cdots\), \(q=-e^{-\pi\sqrt{163}}\), put \(e^{\pi\sqrt{163}}\) just **below** \(J=640320^3+744=262537412640768744\). See [Zagier's modular-forms notes](https://www.its.caltech.edu/~matilde/Zagier123ModularForms.pdf). The even integer \(J\) is not a tournament count. The nearest odd integer is \(J-1\), and we realize it **exactly**:
\[
 262537412640768743=3677\cdot4787\cdot91433\cdot163129 .
\]
The companion lists adjacency masks and verifies the four respective counts by independent subset DP, on \(10+10+12+12=44\) vertices. Their ordinal join multiplies the counts. On the upper side the same check constructs \(J+3=3^2\cdot47\cdot211\cdot229\cdot431\cdot2143\cdot13907\) on 62 vertices. This is a numerical straddling pair, **not** an arc flip or a claim about \(J+1\); the latter still has the unhandled prime factor \(509781383768483\) in this join route. No assertion is made that the exhibited orders are minimal.

Replay: `python3 heegner-witness.py` for both joins and the two Nagell square cores. The proof uses literal witnesses, not numerical closeness as a proxy for a count.

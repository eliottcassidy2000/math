# The split/inert sign, the square-sum graph, and the signed Collatz clock

**PROVED** for the maps and congruences below; **FINITE-EXACT** for the graph/tournament censuses. Neither fixed-sign Collatz nor the full tournament spectrum is resolved here.

Inherit three different objects from the repository: the [first birth-color ambiguity at \(Q_{36}\)](../session40_marked_ports_square_sum_20261023/square-sum/COLOR_FIRST_FAILURE.md), the [8–9 reflection block and its \(9/8\) scale](../chanpak_middle_out_square_sum_20261122/SQUARE-SUM.md), and the [signed ternary descent with hole bypass](../session304_collatz_peeling_20261113/ARITHMETIC.md). The new connection is a phase map that transports actual square edges, and a sharper criterion for when such an edge is literally a signed Collatz step.

## 1. Two different parity tables

Call odd \(n\) **split** for \(n\equiv1\pmod4\) and **inert** for \(n\equiv3\pmod4\). These are congruence names, not a primality test. With \(\chi(n)=(-1)^{(n-1)/2}\), one has \(\chi(mn)=\chi(m)\chi(n)\), giving
\[
S S=S,\quad SI=I,\quad II=S.
\]
This is XOR written as multiplication. In the odd coordinate \(n=2u+1\), multiplication transports to \(u\star v=u+v+2uv\), which reduces mod 2 to \(u+v\). Ordinary even/odd tags, by comparison, use addition and multiplication in \(\mathbb F_2\): adding odd toggles the tag, and any even factor makes a product even. They are not the same multiplication table.

For an odd *modulus* \(m\), its split/inert tag is the Jacobi symbol \(\left(\frac{-1}{m}\right)\). For an odd prime this detects Gaussian splitting; a composite with tag \(S\) need not have a square root of \(-1\bmod m\): \(21\) is split as a residue class, but \(-1\) has no square root modulo 3 or 7. Also an odd residue modulo an odd \(m\) does not determine the sign tag: \(n\) and \(n+2m\) give the same residue but opposite signs. CRT permits keeping the modulo-four tag as a separate coordinate. On any odd modulus, negation conjugates the plus and minus affine branches with the same halving exponent; there is no intrinsic preference for either sign from the modulus tag alone.

For \(D_s(n)=\operatorname{oddpart}(3n+s)\), \(D_+(-n)=-D_-(n)\). Thus positive \(3n-1\) and negative \(3n+1\) are reflected copies, not two independently discovered families. The three known positive minus cycles in odd coordinates are \((1)\), \((5,7)\), and \((17,25,37,55,41,61,91)\); their negative copies follow by reflection. These are the classical observed cycles ([Lagarias, §2.6](https://www.cecm.sfu.ca/organics/papers/lagarias/paper/html/paper.html)), not a proof that the list is exhaustive. More directly, choosing the numerator sign to agree with the integer's sign gives \(D_+\) on magnitudes; choosing the opposite sign gives \(D_-\) on magnitudes. The positive plus problem is a different fixed-sign dynamics, still open.

Picking the sign from the state is another map again:
\[
s=\chi(n) \Longleftrightarrow 4\mid 3n+s,\qquad v_2(3n-\chi(n))=1.
\]
The fully accelerated first policy strictly decreases positive odd \(n>1\), reaching 1, while the complementary one-halving policy strictly increases \(n>1\). Neither is the fixed-sign Collatz map. For fixed \(+\), the residue \(I\) is the one-halving branch and \(S\) the deeper one; for fixed \(-\) those roles reverse. [Odd modular control and the one-split cycle test](SPLIT-SKEW.md) make this distinction quantitative.

## 2. The exact four-down chart

**Theorem.** The selected-sign chart
\[
\rho(n)=\frac{3n+\chi(n)}4,\qquad n>0\text{ odd},
\]
is a bijection to the positive integers coprime to 3. It acts by \(4u+1\mapsto3u+1,\ 4u+3\mapsto3u+2\). In this chart the split/inert character becomes the ordinary sign of the residue modulo 3, denoted \(\chi_3\): \(\chi(n)=\chi_3(\rho(n))\). Unlike the full odd part, \(\rho\) divides by exactly four and retains the higher binary clock.

For two distinct odd vertices, the following predicates are equivalent:
\[
n+m\ \text{is a square}\quad\Longleftrightarrow\quad
\rho(n)+\rho(m)=3t^2\text{ for some }t.
\]
In the forward direction, \(n+m=(2t)^2\) implies opposite \(\chi\) signs, so the two signs cancel in the sum. In the reverse direction a sum divisible by 3 forces opposite \(\chi_3\) coordinates in the target; the same equation read backwards gives \(n+m=4t^2\). In particular the odd–odd subgraph of \(Q_{36}\) is **isomorphic**, not just homomorphic, to the graph on \([27]\setminus3\mathbb Z\) whose edges sum to three times a square.

Odd multiplication transports exactly as well. With \(c(t)=(\chi(t)-t)/4\in\mathbb Z\), the action of an odd multiplier \(t>0\) in the chart is
\[
(a,e)\longmapsto (ta+c(t)e,\ \chi(t)e),\qquad
a=\rho(n),\ e=\chi_3(a).
\]
Indeed this is \(\rho(tn)\), and the affine carry satisfies \(c(st)=s\,c(t)+c(s)\chi(t)\). For \(t=3\), multiplication is \((a,e)\mapsto(3a-e,-e)\). The same signed affine law reduces mod any odd modulus if the sign is retained separately. In an ordinal tournament join the count multiplies, so adding a triangle realizes the \(t=3\) row physically. No such argument identifies the *Collatz* edge with a tournament reversal.

For the full graph \(Q_N\), put \(\sigma(n)=(-1)^{\lfloor n/2\rfloor}\), extending the odd signs to the even vertices. On a square edge the signs disagree **iff both endpoints are odd**: odd sums are \(1\bmod4\), and even square sums are \(0\bmod4\). Equivalently the odd–odd edge bit is the change of the vertex sign, and its total parity on a path is fixed by the endpoint signs. The mod-four quotient has chain \(0\!-\!1\!-\!3\!-\!2\) with loops only at \(0,2\); not every edge allowed by this quotient is a square edge (e.g. \(1+7=8\)).

The valuation sidecar is visible too. On an odd–odd square edge \(n+m=(2t)^2\) set \(A(n)=v_2(3n+\chi(n))\), \(D(n)=D_{\chi(n)}(n)\). Then
\[
2^{A(n)}D(n)+2^{A(m)}D(m)=12t^2.
\]
With unequal clocks the smaller is \(v_2(12t^2)\), necessarily even. If the clocks agree, their common value is strictly less than \(v_2(12t^2)\), and \(D(n)+D(m)\) is three times a square for an even clock and six times a square for an odd clock. In \(Q_{36}\) the 16 odd–odd edges have 5 equal clocks: \((1,15),(7,9),(31,33)\) have clock 2 and \((3,13),(29,35)\) have clock 3. The edge \((5,11)\) has both clocks at least 3 but has actual clocks 4 and 5. This is one reason a single phase bit cannot stand in for the exact clock.

## 3. Contact theorem: a square edge can be a signed Collatz step, but cannot be followed by another

Call \(x\to y\) a contact for distinct positive odd vertices when \(y=(3x+s)/2^a\), \(a=v_2(3x+s)\), and \(x+y=m^2\). The root \(m\) is even. With \(q=2^a+3\), the condition is exactly
\[
2^a m^2\equiv s\pmod q,\qquad x=(2^a m^2-s)/q,\qquad y=m^2-x.
\]
A root class can always be represented by a sufficiently large even \(m\); then both displayed integers are odd, positive, distinct, and the clock is exact.

- \(a=1\): neither sign works modulo 5.
- Even \(a\ge2\): \(+\) works, since \(2^a\) is a square modulo \(q\); \(-\) would make \(-1\) a square modulo \(q\equiv3\pmod4\).
- Odd \(a\ge3\): \(+\) is impossible since \(\left(\frac2q\right)=-1\) for \(q\equiv3\pmod8\). The minus case works exactly when \(-2\) is a square mod \(q\), i.e. every prime divisor of \(q\) is \(1\) or \(3\bmod8\) by the supplementary laws, odd-prime-power lifting, and CRT. For example clocks 3 and 7 work, while 5 does not.

Thus every contact uses \(s=(-1)^a\), and has \(y\equiv1\pmod3\). The tail cannot be \(1\bmod3\), since \(x+y\equiv2\bmod3\) cannot be a square. Also each tail has at most one exit: exactly one sign has clock at least two, and clock one is forbidden. Therefore the entire signed contact graph is a disjoint union of inward stars and isolated vertices. **No two consecutive accelerated steps, even with freely changing signs, can both be literal square-sum edges.**

This strengthens the earlier [fixed-plus obstruction](../chanpak_middle_out_square_sum_20261122/SQUARE-SUM.md#1-raw-square-sum-edges-are-isolated-on-every-accelerated-collatz-trajectory). It explains why the earlier \(8\)--\(9\) construction used whole square paths and retained a centered clock. If the signs alternate and the two exponents are \(1,2\) in either order, the affine branch is
\[
(9x+\delta)/8,\quad \delta=s(3-2^a)\in\{\pm1\},\qquad x'+\delta=(9/8)(x+\delta).
\]
If the first exponent is \(a\) and first sign \(s\), the allowed residue is \(x\equiv15\bmod16\) for \(\delta=+1\), and \(x\equiv1\bmod16\) for \(\delta=-1\) (substitute to verify both exact clocks). Reduction modulo any odd \(m\) makes the centered step a multiplication by \(9/8\), invertible iff \(3\nmid m\); exact integer legality still uses the mod-16 sidecar. Repetition of these two-clock maps increases positive \(x>1\), so their modular returns alone do not give positive cycles. This is separate from the [fixed-sign unit-gap mirror](../chanpak_mirimanoff_unit_lifts_20261122/MIRROR-UNIT.md), which uses carries \(5,7\) instead of \(\pm1\).

## 4. The first 36 vertices, and a first strong tournament seam

There are exactly two contacts in \(Q_{36}\): \(3\xrightarrow{-,3}1\) and \(9\xrightarrow{+,2}7\). The following are unoriented Hamiltonian cycle counts, grouping by membership of the two underlying undirected edges.

| Order | Neither | \(1\!-\!3\) only | \(7\!-\!9\) only | Both |
|---:|---:|---:|---:|---:|
| 32 | 1 | 0 | 0 | 0 |
| 33 | 1 | 0 | 0 | 0 |
| 34 | 1 | 10 | 0 | 0 |
| 35 | 29 | 10 | 15 | 3 |
| 36 | 15 | 6 | 9 | 1 |

Orders 3 through 30 have a vertex of degree below two, and 31 has zero cycles by the exact scan. So both contacts first coexist on a Hamiltonian cycle at 35, and the first unique both-contact cycle is at 36. Its order beginning at 18 is
\[
18,7,9,27,22,3,1,8,17,32,4,5,11,25,24,12,13,36,28,21,15,34,30,19,6,10,26,23,2,14,35,29,20,16,33,31.
\]
Close this list back to 18. This is a new incidence-constrained cycle statistic, separate from the inherited three doubleton *path* fibres at 36.


### One signed bit repairs the inherited \(Q_{36}\) ambiguity

The birth-color observer keeps the actual endpoints and the ordered word of three birth colors, with each path oriented from its smaller endpoint. The inherited census gave 12960 paths, 12957 labels and three doubletons. The new observation is that the **single bit \(\sigma\) at the 19th vertex** distinguishes all three doubletons simultaneously:

| Endpoint pair | 19th vertices in the two paths | Their signs |
|---|---|---|
| \(6,32\) | \(31,29\) | \(-,+\) |
| \(12,21\) | \(2,13\) | \(-,+\) |
| \(21,32\) | \(35,33\) | \(-,+\) |

Position 19 is the **only** single vertex position where this sign separates all three fibres. The odd–odd edge count alone does not help: paired totals are \(11/11,\ 10/10,\ 10/10\). A single ordered odd–odd edge bit at edge position 17 or 23 (one-based) would also suffice; this follows by testing the same three fibres, and the full ordered bit is equivalent to the signed vertex word given the starting endpoint. The finite [decoder census](q36-sign-repair.cpp) regenerates the three pairs with anchored orientation and compares only the additional observations. This is a one-bit repair of this concrete failure, not a claim that any fixed bit decodes every larger square-sum graph.

The triangle gives the tournament flip \(3\leftrightarrow1\). The arithmetic contact \(9\to7\) is **not** realizable by tournament counts. The next contact by increasing source is \(57\xrightarrow{+,2}43\), from the square \(100\). Outgoing rows \((32,57,43,113,101,64,7)\) and reversal \((3,4)\) give exactly this transition, with both orientations strong. Its order seven is minimal since order six has maximum 45. The additional exact witness table in the replay gives other contacts, not a census of all such tournament seams. This is a direct physical connection from the arithmetic into the metagraph, but it does not supply a renewable Chan–Pak move: the contact theorem already rules out traversing two such square-edges consecutively in their arithmetic directions.

The forbidden heights are in different phases, \(7\in I,\,21\in S\), exactly as \(21=3\cdot7\) predicts. Their mod-four tags cannot explain their impossibility. A proof that all other heights occur still needs an absolute tournament construction. The contact \(9\to7\) shows plainly why an arithmetic descent or its modular inverse is not by itself a physical theorem.

## Replay

Run `node signed-square-contacts.js` for the new cycle census and tournament masks. Run `c++ -O3 -std=c++17 q36-sign-repair.cpp -o /tmp/q36-sign-repair && /tmp/q36-sign-repair` for the finite path-fibre test. Run `python3 split-skew-check.py` for literal modular trajectories and one-split signatures. These distinguish numeric residues, signed clocks, and actual graph incidence.

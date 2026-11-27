# Which tournaments do odd-modulus Collatz multiplication tables select?

This continues the [four-state binary matrix](../mixed_moore_four_state_tiles_20261127/README.md) and [two-odd clock](../cullen_woodall_signed_collatz_clocks_20261127/README.md) work. There are two different interpretations of a modular Collatz edge; only the second gives a tournament on the whole residue set.

## 1. From a sparse routing matrix to a tournament

For odd \(m>1\), write \(a=3/2\) and \(g=a^2=9/4\) mod \(m\). In the coordinate \(x+1\), the halved odd Collatz branch multiplies by \(a\), and two odd steps multiply by \(g\). The matrix of arrows \(u\to gu\) has one entry per row before removing loops. It is a functional digraph, not generally a tournament. On an odd orbit of length \(h\), one can orient \(u\to g^iu\) for \(1\le i\le(h-1)/2\), giving the cyclic carousel on that orbit, but separate orbits remain disconnected.

The whole-residue construction instead uses the **additive difference** \(y-x\). Choose \(D\subset\mathbb Z/m\setminus\{0\}\) with exactly one from each pair \(r,-r\), and orient
\[
x\to y\quad\Longleftrightarrow\quad y-x\in D.
\]
This is a regular circulant tournament. The multiplication table says how to group the difference labels: require \(gD=D\), so the centered two-odd map acts by an automorphism. This does not assert that every modular step is executable at an arbitrary positive integer representative; odd-run parity guards still have to be met independently.

**Existence and count.** If \(3\mid m\), the map \(g\) collapses distinct vertices, so cannot act as a tournament automorphism. Otherwise an invariant tournament of this circulant kind exists if and only if \(\operatorname{ord}_m(g)\) is odd, equivalently if each \(\operatorname{ord}_p(g)\) is odd for \(p\mid m\). If the order is odd, no orbit of a nonzero \(r\) meets its negative: an odd-order cyclic action cannot take \(r\) to \(-r\) without fixing \(r\). If the order is even, some prime factor has even local order, since the prime-power kernel has odd order. On the additive order-\(p\) layer, the halfway power of \(g\) is \(-1\), so it reverses a required edge. In fact an even-order permutation cannot preserve *any* tournament, because its generated subgroup contains an involution that swaps a pair.

At an admissible modulus each orbit must be wholly in \(D\) or wholly out. The exact number of choices on the fixed residue labels is
\[
2^{\frac12\sum_{q\mid m,\ q>1}\varphi(q)/\operatorname{ord}_q(g)} .
\]
The summand counts orbits of differences of additive order \(q\), not merely units. Each orbit is paired with its negative. For instance \(m=25\) has four orbits of length one and four of length five, giving 16 choices.

## 2. One odd step can preserve, reverse, or mix the choices

Whenever \(g\) is admissible, the order of \(a\) at each prime-power factor is either odd or twice odd. Within a difference layer of additive order \(q\):

* If all prime powers in \(q\) have odd \(a\)-order, then \(a\in\langle g\rangle\) mod \(q\); it acts within each orbit and preserves the layer.
* If all have twice-odd \(a\)-order, an odd common multiple of their odd parts gives \(a^j=-1\bmod q\). Thus \(a\) swaps each \(g\)-orbit with its negative and reverses the layer.
* If both kinds occur, \(a\) is neither in \(\langle g\rangle\) nor in \(-\langle g\rangle\) mod \(q\); it moves an orbit to another, non-antipodal orbit. On paired-orbit variables it couples distinct choices.

At \(m=161=7\cdot23\), \(a\) has order 6 mod 7 and 11 mod 23. Thus one step reverses the order-7 layer, preserves the order-23 layer, and exchanges distinct orbit pairs in the unit layer. Two steps preserve the whole tournament. This is an exact reason that the coarse sign \(m\bmod4\) does not by itself predict the matrix: \(5\) and \(29\) admit examples; \(13\) and \(17\) do not.

## 3. Exact isomorphism classes at prime orders

For prime \(p\ne3\) with \(h=\operatorname{ord}_p(g)\) odd, set \(2r=(p-1)/h\). Number the \(g\)-cosets by the successive powers of a primitive root. A connection set is described by a binary string of length \(2r\) with opposite bits \(r\) positions apart. Multipliers give cyclic shifts, and **these shifts classify up to all graph isomorphisms**, not just affine ones. Indeed, in the automorphism group of a prime-order circulant the translation group is Sylow \(p\). An isomorphism transports that group to another Sylow \(p\) of the destination; compose with an automorphism making the two groups agree. The resulting isomorphism normalizes the translation group and is therefore affine.

Burnside gives the number of classes
\[
\frac1{2r}\sum_{\substack{0\le j<2r\\ \gcd(2r,j)\nmid r}}
             2^{\gcd(2r,j)/2}.
\]
A shift \(j\) fixes strings constant on \(\gcd(2r,j)\) cycles. If the antipode is within a shift-cycle, the required opposite bits forbid a fixed string; otherwise shift-cycles are paired freely. In particular if \(h=(p-1)/2\), the subgroup is the quadratic residues, and the tournament is forced to be [Paley](https://arxiv.org/abs/1905.08568), up to converse (which is itself isomorphic). This happens at \(7,11,23,31\) in the replay. At \(5\) the two-odd map is the identity, and the regular 5-carousel is the unique circulant class.

### First split: 19

At 19 the two-odd multiplier has order 3 and the eight choices form exactly two isomorphism classes. The table gives exact Hamiltonian path counts and the histogram in one row of \(A+A^2\) away from its diagonal.

| Order | Class | One-plus-two-step routes per destination | Hamiltonian paths |
|---:|---|---|---:|
| 5 | carousel | \(1:2,\ 2:2\) | 15 |
| 7 | Paley | \(2:6\) | 189 |
| 11 | Paley | \(3:10\) | 95,095 |
| 19 | non-Paley | \(4:6,\ 5:6,\ 6:6\) | 1,167,595,581,285 |
| 19 | Paley | \(5:18\) | 1,172,695,746,915 |

The representative connection sets and exact subset DP are in the replay. Thus one modular multiplication table can give two macroscopic route profiles and two absolute heights even after the same translation and Collatz symmetry are imposed. For a doubly regular tournament, \(A+A^2=(p+1)(J-I)/4\). This is the balanced multiroute cousin of the earlier exact-one Moore identity, not the same equation except for the directed triangle.

Finally the translations and \(g\) generate a group of order \(m\operatorname{ord}_m(g)\) acting freely on ordered Hamiltonian paths, so every count here is divisible by that number. For \(m\ge5\), \(D\) contains at least two different units (one per opposite unit pair), each giving a different Hamiltonian cyclic order from every start by repeatedly adding the chosen unit. Consequently these counts are composite; the high symmetry gives useful test families for Chan–Pak, but cannot by itself cover prime targets.

Run `node replay.js` to reproduce the table and check difference-orbit gating for odd moduli up to 401 prime to 3.

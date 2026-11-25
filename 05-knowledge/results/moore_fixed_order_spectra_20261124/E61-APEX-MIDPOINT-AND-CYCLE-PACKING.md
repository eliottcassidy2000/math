# Midpoint 1 mod 3 at exponent 61, and the phase packing test

This follows the [actual 133-site exact lower and modular upper exit](../chanpak_decimal_kernel_20261122/BACKCHAIN-LOWER.md), the [compact 91/93-site exact lowers](../chanpak_e61_compact_volchkov_20261122/STRONG-93.md), and the [prime shell tests](RADIUS-TWO-PRIME-BALANCE.md). Below, \(E_p=2^{p-1}(2^p-1)\), whether or not \(2^p-1\) is prime. The exact straddle \(E_{61}-1\leftrightarrow E_{61}+1\) remains open outside the tested grammars.

## 1. The complete binary inverse for a shared transitive core

Take an apex \(z\) over a transitive sequence of ordinary vertices. Let \(s_i=1\) when \(z\) beats vertex \(i\). Every directed odd cycle contains \(z\). If its first and last ordinary sites are \(i<j\), then \(s_i=1\) and \(s_j=0\); the internal subset of sites between \(i,j\) has even cardinality. Thus [the odd-cycle packing identity](https://arxiv.org/abs/2307.05569) gives
\[
      H=1+\sum_{i<j:s_i=1,\ s_j=0} f(j-i),\qquad
      f(1)=2,\quad f(d)=2^{d-1}\ (d\ge2).
\]
All the cycles here intersect at \(z\), so there is no further packing correction.

Mark the apex contact to the ordinary pivot at position zero. Write \(u\) for the binary integer formed by the left out-bits \(s_{-i}\) with \(i=1\) the low bit; \(v\) is the integer of right in-bits \(1-s_i\) in the same convention. Trailing far zero bits do not affect the read. Write \(G(w)=\sum_{i>j} f(i-j)\) over positions \(i\) of 1 and \(j\) of 0 in \(w\), with positions starting at 0. The nonpivot count is \(1+G(u)+G(v)+2uv\); this is because every left-to-right pair has distance at least two and its exponential term separates. Changing the pivot from in to out replaces the contribution \(u+(u\bmod2)\) by \(v+(v\bmod2)\). Hence this edit gives \(+2\) exactly when
\[
        u\in\{2q-1,2q\},\qquad v\in\{2q+1,2q+2\}
\]
for some \(q\ge0\), with negative values omitted. Its midpoint is
\[
                     M=2uv+G(u)+G(v)+2q+2.
\]
This makes the useful Collatz analogy precise: divide by two **after recording the low bit**, with two inverse branches on each flank. There is no assertion about actual Collatz dynamics. Midpoint 1 mod 3 is very easy in this grammar: \(u=2^k,\ v=2^k+1\) for every \(k\ge2\) gives \(M=2u^2+9u/2+2\equiv1\pmod3\). These are physical unit flips, not just modularly near-unit jumps. The old [two-mark family](FOURIER-OVERLAPS-AND-THREE-DELETIONS.md) also gives this class via \(54\cdot2^{m-1}-2\) as midpoint.

### Why none hits a perfect midpoint

For each 1-bit at index \(i\ge1\), the sum of all possible weights to lower positions is \(2^i\); hence \(0\le G(w)\le w\). It follows that
\[
                     8q^2+2q\le M\le8q^2+14q+4
\]
(the \(q=0\) cases obey the same bounds). Put \(x=2^{p-2}\), so \(E_p=8x^2-2x\). For \(p\ge5\), \(q\le x-2\) implies \(M-E_p\le-16x+8\), and \(q\ge x+1\) implies \(M-E_p\ge20x+10\). At the only two remaining values, insert \(G(2x-3),\ldots,G(2x+2)=x/2,x,0,2x,x,3x/2+2\). This gives the exact signed offsets:

| \(q\) | \(u=2q-1,v=2q+1\) | \(u=2q-1,v=2q+2\) | \(u=2q,v=2q+1\) | \(u=2q,v=2q+2\) |
|---|---:|---:|---:|---:|
| \(x-1\) | \(-23x/2+6\) | \(-11x/2\) | \(-7x+4\) | \(-x\) |
| \(x\) | \(5x\) | \(19x/2\) | \(11x+2\) | \(31x/2+4\) |

Thus a common transitive core with an apex flip cannot hit any \(E_p\) exactly from \(p\ge5\). The closest apex-flip midpoint is \(E_p-x\). For odd \(p\ge5\), the closest one **also in the target residue class 1 mod 3** is \(E_p+11x+2\): it is the only such offset among the eight closer candidates, and the outside bounds are larger in absolute value. At \(p=61\) these signed errors are \(-576460752303423488\) and \(+6341068275337658370\). Thus attaining the right midpoint class is not an exact repair.

There is one other possibility if the *same* apex deletion is transitive on both ends: reverse an adjacent ordinary pair. In the ordinary order this exchanges two consecutive bits. If they differ, choose the sign so the change is from 01 to 10. The pair contributes \(f(1)=2\). Each earlier 1 and later 0 contributes a nonnegative increment \(f(d+1)-f(d)\), which is zero just for \(d=1\). A unit change therefore requires all sufficiently earlier bits to be 0 and all sufficiently later bits to be 1. The unaffected outside prefix and suffix are ordinal singletons; the active graph has at most the apex and four ordinary vertices, hence height at most \(5!<E_5-1\). A nonadjacent ordinary reversal cannot preserve the common transitive deletion. This exhausts the **shared-apex transitive-deletion** grammar, not all graphs that individually admit different transitive deletions.

### One extra ordinary reversal does not repair the eight nearest parameter forms

For all eight seeds in the table at \(p=61\), I also tried a common extra ordinary-arc reversal at both endpoints, leaving the marked apex contact unchanged. This exceeds the transitive-deletion grammar. In **59,050** such repairs, neither respective endpoint reaches its target even modulo \(2^{128}\). This is a certified exclusion for those seeds only; it does not restrict arbitrary farther two-defect charts.

The replay computes a change by the difference of two forced-adjacency counts, one for each direction of the ordinary arc. Contract that forced adjacent pair to a block, remembering its incoming and outgoing endpoint. There are two orders for the block and the apex, with three increasing ordinary runs. Assign vertices in ordinary order, recording each run as empty, with legal exit, or with illegal exit. First entrances are tested at insertion, so the 27 states give exact totals in any modulus; the finite scan uses unsigned arithmetic mod \(2^{128}\). Literal permutation replay independently checks every ordinary reversal in the small controls. An exact target would in particular have survived the modular test.

## 2. Reciprocal routes determine phase conflicts, not just cycle counts

In a cyclic voltage chart, a lifted cycle is specified by its projected orbit word and its successive phase differences, together with a starting phase. Closure requires total phase zero, and a valid simple lift must not repeat a pair (orbit, phase). For each valid cycle \(C\), retain the **footprint** \(P_C(i)\) of its phases in each orbit \(i\). To translate two cycles \(C,D\) by \(t_C,t_D\), the exact collision rule is
\[
    t_D-t_C\ \notin\ \bigcup_i (P_C(i)-P_D(i)).
\]
Fixed vertices, if present, have their own ordinary occupancy check. The rule follows by equating the phase at a shared orbit. It works for repeated projected visits as long as the individual lifted cycle is simple. Projecting and retaining only the orbit *set* loses this information.

The reciprocal row transports supply the successive phase moves for the known graph edges. In the [fixed-pentagon atlas](PENTAGON-PORT-ATLAS.md), the port slots determine which physical destinations are available at each step; in the [order-13 voltage branch](SINGER-ORDER13-MODULES.md) they can be represented by residues. These are **separate conditional symmetry branches**, not assumed actions on the same graph. Either way, the actual compatibility test is the footprint disjointness rule after lifting. A useful sufficient condition for an ordered list of cycles is that when placing each new cycle, the sum of the sizes of its forbidden difference sets against the earlier cycles is less than the phase modulus. Greedy choice then always works. When that bound fails, solve the finite simultaneous difference constraints; pairwise nonemptiness is not a certificate.

There is also an exact per-flip arithmetic target for that information. For an arc \(e\), take all directed odd cycles that use \(e\) before or after reversal, assigning signs minus or plus, respectively. In a disjoint collection at most one cycle can use the flipped arc. Every cycle not using it is unchanged. Therefore
\[
  \frac{\Delta H}{2} = \sum_C \epsilon_C\,H(T-V(C)).
\]
Each term counts compatible odd-cycle packs on the complement, not just the number of possible words through the edit. To hit \(E_{61}\) one needs an exact lower at \(E_{61}-1\) and this signed sum **equal to 1**. Row reciprocity can certify which lifted cycles occupy the remaining sites, but gives neither their tournament orientations nor that integer balance for free.

### A concrete Hoffman–Singleton control

Use the order-five action \(A_{i,x}\mapsto A_{i,x-i}\), \(B_{j,y}\mapsto B_{j+1,y}\) on the [Robertson–Hafner model](https://www.math.auckland.ac.nz/~hafner/his.pdf), which fixes a pentagon. The [previous cycle atlas](THIRTEEN-ODD-ALPHABET.md) counts 126,000 induced 13-cycles, the longest odd induced length. Enumerate them now by undirected support and test their own translates under this fixed action:

| Disjoint from translate by \(g\)? | Disjoint from translate by \(g^2\)? | Number of 13-cycles |
|---|---|---:|
| no | no | 125,900 |
| yes | no | 50 |
| no | yes | 50 |
| yes | yes | 0 |

The tests at \(g^{-1},g^{-2}\) coincide with those at \(g,g^2\). In particular at most **two** translates of a single such cycle can coexist in these five phases, even though the trivial volume bound is three. In the hypothetical order-13 Moore graph, the internally cyclic orbits are different: the symmetry rotates an induced 13-cycle and its footprint is already **the entire orbit**. Such an internal cycle excludes any other cycle that visits that orbit, while distinct internally cyclic orbits are automatically vertex-disjoint. The prior census allows 30, 120 or 210 such orbits (only 30 in the cyclic crossing branch). Other lifted cycles still compete for phase slots. The numerical length 13 therefore means different compatibility constraints in the two graphs.

### Why the degree-two infinite family is a good control

A directed chordless odd cycle, viewed as a digraph of prescribed arcs, has odd-pack polynomial \(1+z\) regardless of its length; a disjoint union of \(r\) such cycles has \((1+z)^r\). Its large diameter is allowed at degree two; only the pentagon is diameter two. A legal tournament completion illustrates why the length and the completion nevertheless matter: take transitive order on \(\ell\) vertices and turn the last–first edge backwards. All nontrivial cycles use that backward arc, and the odd ones are counted by choosing an odd subset of the \(\ell-2\) interior vertices. Thus the packing polynomial is \(1+2^{\ell-3}z\) and \(H=1+2^{\ell-2}\). The planted spanning odd cycle is only one of the letters on the completed board. In an ordinal join, such polynomials multiply; introducing mixed arcs requires the compatibility calculation again. This is a constructive control for the distinction between allowed letters and allowed simultaneous letters, not a conversion of undirected cycle lengths into prescribed heights.

Replays: `python3 apex-midpoint-check.py`; compile `apex-one-defect-check.cpp` and `hs-cycle-phase-check.cpp` with GNU C++17 and run both. They check the arithmetic offsets, the bounded exact-exclusion scan, and the full Hoffman–Singleton translate census.

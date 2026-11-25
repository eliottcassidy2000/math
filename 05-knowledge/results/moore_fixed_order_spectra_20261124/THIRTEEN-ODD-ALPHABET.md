# Thirteen as a cycle letter, an orbit length, and an exact arithmetic selector

The [Chan–Pak conjecture](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf) asks for all positive odd Hamilton-path counts other than 7 and 21. An induced cycle length in an undirected Moore graph, a tournament Hamilton-path count, a path-orbit count, an eigenvalue multiplicity, and a gap in one construction bank are different observables. This session compares them by tests that preserve the observable, rather than by equating them.

## 1. The odd alphabets of the small Moore graphs

The exact replay constructs Petersen as a pentagon–pentagram matching and Hoffman–Singleton in the [Robertson–Hafner affine model](https://www.math.auckland.ac.nz/~hafner/his.pdf). It checks the Moore equation physically and visits every starting edge.

| Graph | Odd induced cycle lengths | Other exact observations |
|---|---|---|
| Pentagon | 5 | The whole graph is one letter. |
| Petersen | 5 | Ordinary cycle lengths are 5, 6, 8, 9; **7 is absent even without the induced restriction**. |
| Hoffman–Singleton | 5, 7, 9, 13 | **21 is absent as an induced length, but present as an ordinary cycle**. Also its \(-3\)-eigenspace has dimension 21. |

For Hoffman–Singleton the entire induced-length list is 5, 6, 7, 8, 9, 10, 12, 13, 16, 18. Every edge has 9,360 induced 13-cycles through it, oriented with that edge first. Counting by edges, the graph has **126,000 induced 13-cycles**. One is
\[
(0,1,2,3,28,8,7,6,44,15,16,17,30).
\]
Give its vertices distinct colours and all others an initial common colour. Refining by the current colour and the multiset of neighbour colours gives 44 classes, then 50. Hence the **pointwise stabilizer is trivial**. By [Hafner's group-order computation](https://www.math.auckland.ac.nz/~hafner/his.pdf) the whole group has order 252,000, with no factor 13. A setwise stabilizer injects into the dihedral group of a 13-cycle, so it has order at most two. The orbit therefore has at least 126,000 elements, exactly the total just counted. This proves all induced 13-cycles form **one orbit**, with stabilizer a reflection. Keeping the cyclic orientation, but not a start, removes that stabilizer and gives a free transitive action of the full group.

This is a stronger role for 13 than an isolated coincidence: it is the longest odd induced letter in Hoffman–Singleton, and its oriented copies give a regular indexing set for the group. It does not make 13 forbidden as a tournament height. Nor can missing induced lengths alone explain 7 and 21: Hoffman–Singleton has induced 7-cycles, and also omits induced lengths 11 and 15, both realized by tournaments.

## 2. A taxicab test of the “next odd gap” analogy

The adjacency identity \(A^2+A=56I+J\) implies that the hypothetical degree-57 graph has nonconstant eigenvalues 7 and −8 of multiplicities **1729** and 1520. The odd multiplicity \(1729=7\cdot13\cdot19\) is also \(1^3+12^3=9^3+10^3\), returning to the old [March Ramanujan dossier](../ramanujan_1729_bridge_s133.out). Its [path-orbit observation](../../../04-computation/bridge_1728_1729_s122.py) is exact: the Paley 11-tournament has 95,095 paths, with 55 automorphisms, giving 1,729 path orbits. Paley 7 has \(H=189\), and **21 nonclosing** Hamilton paths, not 21 total. The cyclic 13-carousel (differences 1 through 6) has \(H=3,711,175\), giving 285,475 path orbits. Its group is the translations: the count of common out-neighbours distinguishes step 1, so every automorphism preserves that directed 13-cycle. A quadratic-residue Paley **tournament** on 13 is not available (\(-1\) is a square).

**1729 is not a new odd gap.** There is a literal one-arc step **1729 → 1731** on nine vertices. Encode an arc \(i\to j\) for \(i<j\) by bit 1 in lexicographic pair order. Reverse \(\{1,5\}\) in mask 815318702, giving 815320750. The [Grinberg–Stanley odd-cycle collection law](https://arxiv.org/abs/2307.05569) sets \(H=F_T(2)\) for \(F_T(z)=\sum_{\mathcal C}z^{|\mathcal C|}\), with the sum over collections of pairwise vertex-disjoint directed odd cycles of length at least three. The actual alphabets have weight enumerators
\[
  1+478z+181z^2+6z^3,\qquad
  1+485z+180z^2+5z^3.
\]
Both agree with independent endpoint subset DP. The change is \(7z-z^2-z^3\) and reads 2 at \(z=2\). The seven extra single letters do **not** give displacement 14, because they interact with the other letters. This is a concrete counterexample to reading either a single coefficient or a graph multiplicity as a raw tournament gap.

## 3. What order 13 really predicts for the odd-cycle alphabet

Suppose a tournament admits an automorphism of odd prime order \(p\), with \(m>0\) moving orbits and fixed induced tournament \(T_0\). Modulo \(p\), for every \(a\in\mathbb F_p\),
\[
                 F_T(a)=F_{T_0}(a)(1-a/2)^m .
\]
This is an identity of **values** (or polynomials modulo \(z^p-z\)), not a formal factorization in \(\mathbb F_p[z]\).

Proof: noninvariant collections cancel in orbits of size \(p\). An invariant collection is a permutation on its occupied vertices commuting with the prime action. Fixed sites can occur only in fixed cycles, entirely inside \(T_0\). On moving sites the occupied orbits project to permutation cycles. An odd cycle of \(l\) distinct orbit labels lifts, depending on summed shift, to one cycle of length \(pl\) or \(p\) cycles of length \(l\). Evaluated at \(a\), both have weight \(a\) modulo \(p\). Even projected lengths do not contribute. For two distinct orbit labels \(i,j\), the numbers \(b_{ij},b_{ji}\) of shift options sum to \(p\). Reversing the first projected cycle of odd length at least three gives its negative weight, so these contributions cancel. A lone moving orbit is either unused or carries one of \((p-1)/2\) nonzero directed step cycles, giving \(1+(p-1)a/2=1-a/2\). This proves the law.

For the Moore graph itself the [order-13 branch](FIFTH-POWER-HALVES.md) has 30, 120 or 210 internally cyclic orbits (only 30 in the [cyclic-crossing subbranch](SQUARE-REFLECTION-CROSSING.md)). Each is an **induced 13-cycle**, with five, twenty or thirty-five per step class. The action rotates those cycles, unlike Hoffman–Singleton, where the cycle has a reflection stabilizer and no rotation of order 13. An equivariant *tournament* built over that action necessarily has \(13\mid H\). That restriction belongs to the chosen symmetry and disappears when we permit a nonequivariant edit; it says nothing about global unattainability of heights.

## 4. Breaking the symmetry: a literal +2 in every 13-power residue

Two physical families of single-arc edits have the lower/upper readings
\[
12\,2^{m-1}-1\ \leftrightarrow\ 12\,2^{m-1}+1,\qquad
54\,2^{m-1}-3\ \leftrightarrow\ 54\,2^{m-1}-1,\quad m\ge1.
\]
The second is the [previous two-deletion-neutral family](FOURIER-OVERLAPS-AND-THREE-DELETIONS.md). The first uses the six-site mask 12, editing \(\{2,4\}\), and inflating vertex 1 to a transitive chain. That vertex has a unique outside entrance or exit in both tournaments. Every Hamilton path uses at most two runs in the chain, so the height is \(u+v(2^m-2)\). At \(m=1,2\) the readings are **11→13** and **23→25**, respectively; these determine the formula for all widths. Thus the prime itself is also a small exact height, not just a cyclic-orbit label.

The [November 24:13 ledger](../marked_24_13_descent_wheel_20261116/ARCHIVE.md) used the fact that 2 has order 12 modulo 13. Here it extends to a selector: \(2^{12}\equiv40\not\equiv1\pmod {169}\), so 2 generates the units modulo \(13^k\) at every \(k\ge1\) by binomial lifting. Since 12 and 54 are units, the first family selects any prescribed lower residue r with \(13\nmid r+1\); the second selects it whenever \(13\nmid r+3\). These conditions cover **every residue at every precision**, using actual single-arc +2 moves. More generally the same argument works for any prime other than 2 and 3 at whose powers 2 generates the units (e.g. 5 and 11). It does not assert simultaneous CRT selection across all such primes.

This precise local result still leaves absolute targets: the two lower families are sparse, and both omit all lower targets \(1\bmod3\). The old [133-site exponent-61 witnesses](../chanpak_decimal_kernel_20261122/BACKCHAIN-LOWER.md) exhibited the opposite tension (exact absolute lower and an upper phase containing 13, but no exact +2 within the audited family). To attack additional Chan–Pak gaps, the missing ingredient is a height-changing control that preserves the arc displacement while closing the absolute, not just 13-adic, error.

## 5. Historical carry-forward ledger

| Repo lane | What 13 was doing there | What transfers here |
|---|---|---|
| [March Ramanujan/taxicab](../ramanujan_1729_bridge_s133.out) and [1728/1729 path-orbit audit](../../../04-computation/bridge_1728_1729_s122.py) | 1729 factors through 13, and occurs as a path-orbit count at 11. | New exact raw-height seam, and a different equality of graph eigenvalue multiplicity. |
| [Paley/interval race](../../../04-computation/paley_race_p19.py) and [Paley 7 nonclosing](../euler_torus_hamilton_transport_20261020/PUBLICATION_AUDIT.md) | 13-label carousel and the 21 nonclosing paths at 7. | Keep cycle/closure statistics separate from total H. |
| [September LRC(14)/Collatz](../procgen_bridges_20260923_lrc_amm_collatz.md) and [July essential-region theorem](../../../01-canon/theorems/THM-1142-essential-region-criterion-proved.md) | 13 speeds, and exact which-pieces-cover-what tests. | Carry the overlap data; not a conversion from runner clocks to tournament heights. |
| [Marked 24:13 and Riesel](../marked_24_13_descent_wheel_20261116/README.md) and [historical ledger](../marked_24_13_descent_wheel_20261116/ARCHIVE.md) | Order-12 phase clock, two signed carries, Bernoulli denominator and the \(A_{12}^2\) Niemeier Coxeter label. | The order of 2 genuinely lifts to all 13-powers here. The lattice and Monster comparisons alone supply no physical edit. |
| [Signed height / pure-A Leech labels](../chanpak_absolute_liouville_leech_20261122/LEECH-LABELS.md) | 13 among the eight height labels, with observed shadow factors and independent flanks. | A genuine small 13 witness and literal upward seam; retain arithmetic versus geometric labels. |
| [Square-sum middle-out](../chanpak_middle_out_square_sum_20261122/SQUARE-SUM.md) and [our square-reflection crossing test](SQUARE-REFLECTION-CROSSING.md) | 13 beside the centre of the first Pythagorean path; later a 13-phase lift of the two 125-label panels. | Keep position and reflection domain; literal square shifts cannot implement the 25-transversal. |
| [Heegner target selector](../chanpak_heegner_selector_20261122/HEEGNER.md) | A real filled-port \(-13P_1\) term on one branch, not a universal height correction. | Exact width reads and independent exit guards, rather than naming equal coefficients. |
| [Keith–Ulam–Kaprekar](../chanpak_lrc_two_edit_ports_20261123/FILLER-KEITH-ULAM.md) | \(6+13\) is one witness rejecting 19 in its Ulam rule; 19 is a Keith number. | Multiplicity of representations can decide a different grammar. 13 itself is not a two-digit Keith hit. |
| [133-site exact lower at exponent 61](../chanpak_decimal_kernel_20261122/BACKCHAIN-LOWER.md) | Factors of 13 in both upper moduli despite an exact lower height. | Demonstrates why the local selector still needs an absolute inverse. |

Run the two exact companions for small graph alphabets, the marked 13-cycle and its pointwise stabilizer, the taxicab seam, the invariant-packing checksum, and all displayed finite seeds. No degree-57 graph is assumed to have been constructed.

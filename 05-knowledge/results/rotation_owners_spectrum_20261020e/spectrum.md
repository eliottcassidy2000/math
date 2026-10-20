# Ordinary Hamiltonian-path spectrum: two proved omissions, open completeness

**CURRENT STATUS: OPEN** whether the ordinary Hamiltonian-path counts of
finite tournaments are exactly the positive odd integers other than 7 and
21. **PROVED / CITED + FINITE-EXACT:** 7 and 21 never occur. **PROVED:**
the inherited two-prime-lane reduction and the reduced-witness size bound
below. **FINITE-EXACT:** the small independent controls and the displayed
regular cycle-reversal hostile. This note supplies no new infinite
coverage theorem.

## 1. Pin the object before interpreting a spectrum

For a tournament on a fixed finite vertex set, let
\[
 H(T)=\#\{(v_1,\ldots,v_n):v_i\to v_{i+1}\text{ for every }i\}.
\]
Paths are distinct vertex words. Relabelling the entire tournament does
not change this count. We neither identify paths under automorphisms nor
multiply the count by the number of labellings of an isomorphism class.

This is different from:

* \(H(T)/|\operatorname{Aut}T|\), the number of automorphism orbits of
  Hamiltonian paths, since this action is free;
* any blue-incidence count or black-degree statistic of the tiling
  metagraph;
* the spectrum restricted to regular tournaments;
* attainable values at one fixed order, rather than over all orders.

For example, the regular five-vertex tournament has \(H=15\), automorphism
order five, and quotient three. Paley's seven-vertex tournament has
\(H=189\), automorphism order 21, and quotient nine. The checker verifies
these directly. Thus a metagraph count must retain its path-orbit
multiplicity if it is to certify the ordinary spectrum.

The inheritance board is deliberately small:

| Role | Mechanism or hostile | Exact use |
|---|---|---|
| Anchor | THM-4094, two prime lanes | The unclosed infinite existence obligations |
| Niche | Moon--Busch strong minimum | Global exclusions and a finite witness-size cutoff |
| Niche | SCC multiplication | Exact arithmetic construction, with an indivisible seven-adic carry |
| Wildcard | Regular cycle reversal | Degree-preserving motion, but no unit interpolation of \(H\) |
| Hostile | THM-4118 ear image with gcd two but a hole | Lattice compatibility does not imply coverage |

The closest corrected route is
[THM-1370 / h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md),
followed by
[THM-4094 / hamiltonian-matching-deficit-and-two-prime-lane-completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md).
The least-used sidecar for a constructive extension is the full
endpoint/insertion response, including orphan paths, not just its parity.
THM-1221 / seven-wall-strict-spectrum-hunter-floor concerns LRC comb
overlaps; its word “spectrum” does not make it a tournament-value result.

## 2. Primary-source status and a finite-endpoint correction

The primary sources opened on this pass are:

1. Chan--Pak,
   [Computational complexity of counting coincidences](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf),
   *Theoretical Computer Science* 1015 (2024), 114776;
   [arXiv:2308.10214v2](https://arxiv.org/abs/2308.10214v2).
   Section 6.5, Conjecture 6.2 states the exact proposed range as a
   conjecture, attributing its origin to bof and Gordon Royle in 2016.
2. The original
   [MathOverflow question and Royle answer](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament).
   Royle reports coverage of allowed odds **strictly below 80,557**
   using tournaments on at most eleven vertices. Thus the safe inclusive
   endpoint is 80,555. Chan--Pak's wording “up to 80557” should not silently
   change that endpoint.
3. Grinberg,
   [Lecture 12, version April 27, 2026](https://darijgrinberg.gitlab.io/t/22s/lec12.pdf),
   Remark 1.1.11, still describes the higher-value question as open and
   states coverage through 80,555.

This is affirmative documentation of an open problem, not an inference
from a failed search. No later closure was found in this targeted pass.
The original 2016 discussion also already explains how a strong-order
bound and the small enumeration exclude 21. Therefore the “appears new
as a stated theorem” sentence in the historical narrative of THM-1370
is not a defensible current priority claim.

For locally replayable provenance,
[THM-4104 / selected-order-eleven-strong-ear-solid-interval](../../../01-canon/theorems/THM-4104-selected-order-eleven-strong-ear-solid-interval.md)
certifies a selected construction bank and arithmetic coverage through
80,405. It is explicitly not a full order-eleven census. Its first
unforced values are relative to that selected bank, not asserted global
gaps and not the larger literature prefix.

## 3. The two omissions: an exact finite proof with a published cutoff

Busch,
[A Note on the Number of Hamiltonian Paths in Strong Tournaments](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf),
*Electronic Journal of Combinatorics* 13 (2006), N3, Corollary 1,
proves the exact minimum \(f(n)\) among strong tournaments:
\[
f(n)=
\begin{cases}
3\cdot5^{(n-3)/3},&n\equiv0\pmod3,\\
5^{(n-1)/3},&n\equiv1\pmod3,\\
9\cdot5^{(n-5)/3},&n\equiv2\pmod3,
\end{cases}\qquad n\geq3.                              \tag{1}
\]
This is a minimum theorem, not an attainable-value completeness theorem.

The first values, for orders three onward, are
\[
3,5,9,15,25,45,75,125,225,375,625,1125,\ldots.
\]
They satisfy \(f(n+3)=5f(n)\). The exponential base is \(5^{1/3}\),
not \(\sqrt3\). The false Leonardo recurrence gives 41 at order eight
instead of 45.

For completeness, the standalone checker enumerates every labelled
tournament through order six and independently checks path counts by
permutation enumeration through order five. The strong spectra are:

| Order | Strong ordinary \(H\) values |
|---:|---|
| 3 | 3 |
| 4 | 5 |
| 5 | 9, 11, 13, 15 |
| 6 | 15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45 |

By (1), a strong witness to \(H=7\) would have order at most four;
a strong witness to \(H=21\) would have order at most six. The table
excludes both.

Every tournament's strong components are totally ordered by dominance,
and every Hamiltonian path concatenates a path of each component in
that order. Thus
\[
                    H(T)=\prod_{C\text{ SCC of }T}H(C). \tag{2}
\]
A product equal to prime seven requires a factor seven. A product equal
to 21 requires either a factor 21 or a factor seven. Both possibilities
have just been excluded for strong components. This proves the two
omissions at every order.

### Correction lineage that must stay visible

* MISTAKE-559 explicitly retracts the claim that the repository already
  proved *exactly* these two omissions; it routes completeness to the open
  THM-4094.
* MISTAKE-055 corrected the spurious recurrence giving 41 and 67, but its
  replacement sentence beginning “Empirically the data satisfies” contained
  another stale fit: \(f(n)=3f(n-2)\) except at seven, followed by growth
  \((\sqrt3)^n\). At order ten this predicts \(135\), while (1) gives
  \(125\). This session repairs that sentence to (1), preserves the
  correction lineage, and removes the adjacent unsupported completeness claim.
* [THM-200 / h7-impossibility](../../../01-canon/theorems/THM-200-h7-impossibility.md)
  retains its explicit retraction of the proposed \(7\cdot3^k\) gap tower.
  A strong \(H=63\) witness exists; adjoining a source preserves 63 and
  makes it non-strong. Divisibility by seven does not force an SCC
  whose *whole* path count is seven.

The parent session applied these corrections to the maintained mistakes
ledger and retracted THM-1370's unsupported novelty assertion. The proof of
the two exclusions and the open completeness status are unchanged.

## 4. A finite decision cutoff, without an existence theorem

Here is a useful corollary of (1) and (2). For \(h>1\), if any tournament
has \(H=h\), there is such a tournament on at most
\[
                  \left\lfloor3\log_3h\right\rfloor     \tag{3}
\]
vertices. More precisely, deleting all singleton SCCs from any witness
already achieves this bound.

Proof. For every nontrivial strong component of order \(r\geq3\),
\[
                         H(C)^3\geq3^r.                \tag{4}
\]
Check (4) for \(r=3,4,5\) using the minima \(3,5,9\); equality occurs
only at three. Increasing \(r\) by three multiplies the cube of the
minimum by 125 and \(3^r\) by 27, so the inequality is strict thereafter.
Delete the singleton SCCs, whose factors in (2) equal one. Multiply
(4) over the remaining components to obtain \(h^3\geq3^n\), proving
(3). Equality holds exactly for an ordered chain of directed triangles.
This sharpness is for the deletion-only bound on the reduced witness;
it does not assert that this chain has smallest possible order for its
value. For instance \(H=9\) also has a strong order-five witness. \(\square\)

The code uses \(3^n\leq h^3\) with integers, not rounded logarithms.
For \(h=1\), a singleton suffices. Consequently membership in the
ordinary spectrum is decidable by a finite search over
\(2^{O((\log h)^2)}\) labelled candidates, with exact subset-DP path
counting. This is an existence-of-a-finite-algorithm observation, not
a suggestion to run that search at large \(h\), nor a polynomial-time
claim in the bit length of \(h\).

For a value required to be strong, use the sharper direct cutoff
\(\max\{n:f(n)\leq h\}\). For example 80,557 has reduced-witness cutoff
30 and strong-witness cutoff 22. These finite bounds do not assert that
80,557 is absent or present.

## 5. Exactly what remains to close the ordinary spectrum

The proved THM-4094 reduction says that, using (2), the two omissions,
and known strong carry witnesses \(49,63,343\), completeness is
equivalent to both infinite statements:
\[
\begin{split}
&\text{every odd prime }p\ne7\text{ is strongly attainable};\\
&\text{every }7p,\ p\text{ an odd prime},\ p\ne3,
  \text{ is strongly attainable}.                       \tag{5}
\end{split}
\]
These are necessary because factoring a prime SCC product gives its
prime factor, while any proper factorization of \(7p\) requires the
forbidden whole factor seven.

For sufficiency, write an allowed odd value \(h=7^am\), \((m,7)=1\).
When \(a=0\), use prime witnesses. Positive even \(a\) uses powers of
49; odd \(a\geq3\) uses 343 and powers of 49. When \(a=1\), either
\(m=3^k\) with \(k\geq2\), giving \(63\cdot3^{k-2}\), or \(m\) has
a prime divisor \(p\ne3,7\), giving \(7p\cdot(m/p)\). Order joins
realize these products exactly. The checker verifies this conditional
arithmetic recipe through 10,001; it does not supply the unproved
strong witnesses in (5).

A finite generating bank closed under order join and source/sink
padding cannot settle (5): those operations only multiply its existing
values, so a new prime outside the bank cannot appear. This is the
actual finite-bank hostile of THM-4094, not a conclusion that arbitrary
vertex insertion is powerless.

For insertion of a vertex \(v\), the exact count increment is
\[
H(T)-H(T-v)=\sum_{P\in\operatorname{Ham}(T-v)}(a_v(P)-1)+O_v, \tag{6}
\]
where \(a_v(P)\) counts legal positions in \(P\), and \(O_v\) counts
new Hamiltonian paths whose deletion of \(v\) fails the old shortcut.
A construction for a prescribed prime must control this exact
increment, including the orphan term. A selected insertion matching
proves monotonicity but omits precisely these extra witnesses.

## 6. Regular flips preserve a different predicate

Reversing one arc of a regular tournament decreases one endpoint's
outdegree by one and increases the other's by one. It leaves the regular
class immediately. Reversing a directed cycle preserves every outdegree,
since each affected vertex loses one outgoing arc and gains one.

There is an exact positive completeness statement for these moves:
**every fixed labelled score-vector fibre is connected under directed-cycle
reversals.** For two tournaments with the same scores, form the digraph of
arcs in the first that point oppositely in the second. At every vertex
its indegree equals its outdegree. It therefore decomposes into directed
cycles. Reverse these edge-disjoint cycles one at a time; each reversal
preserves scores and removes precisely its edges from the disagreement
set, eventually reaching the second tournament. The checker independently
constructs these paths for all 576 ordered pairs of regular five-vertex
tournaments. This statement permits directed cycles of arbitrary length.

But the latter does not give unit steps in the odd \(H\) spectrum.
Take the regular carousel on seven vertices with arcs of differences
\(\{1,2,3\}\) modulo seven. Its bit rows are
\[
                     (14,28,56,112,97,67,7)
\]
and \(H=175\). Reverse the directed triangle \(0\to2\to4\to0\).
The rows become
\[
                     (26,28,41,112,100,67,7),
\]
still all of outdegree three, but \(H=171\). In fact the fourteen
directed triangles of this source give seven changes of zero and
seven changes of minus four. Thus even this natural regular motion
does not supply the increments of two needed for an intermediate-value
argument.

The regular restriction also loses attainable ordinary values altogether.
Every nontrivial regular tournament is strong: a source SCC of size
\(s<n\) would have average outdegree
\(n-s+(s-1)/2>(n-1)/2\). Its order is odd. Consequently a regular
tournament cannot have \(H=5\): orders one and three give one and
three, while orders at least five have strong minimum at least nine.
The strong order-four tournament does have \(H=5\).

The appropriate next construction is therefore an exact response
mechanism with a specified target class, not regularity alone. For
unit-step interval methods the relevant inherited boundary is
[THM-4118 / ear-response-lattice-and-stateful-unit-component-intervals](../../../01-canon/theorems/THM-4118-ear-response-lattice-and-stateful-unit-component-intervals.md):
even an ear-response gcd of two can coexist with missing values. A
connected state component whose responses change by zero or two
certifies a solid odd interval, and overlapping/touching certified
intervals can be combined. What is still absent is a uniform
unbounded construction of enough such intervals to cover both lanes (5).

## Reproduction

~~~text
python 05-knowledge/results/rotation_owners_spectrum_20261020e/spectrum_check.py
python -O 05-knowledge/results/rotation_owners_spectrum_20261020e/spectrum_check.py
~~~

The saved receipt is spectrum_check.out. Explicit exceptions enforce the
checks in both modes. Universe: all labelled tournaments through order six
(33,867 including orders one and two); independent literal permutations
through five; all fourteen directed triangles of the displayed regular
seven-tournament; all 576 ordered pairs of labelled regular five-tournaments
for the score-fibre transport; exact integer cutoff and conditional
arithmetic controls.
No order-seven-or-higher census is claimed.

# Placing a physical +2 window at 30011 and 30013

**Exact target result.** The two primes 30011 and 30013, fixed before the
search, are Hamiltonian-path counts of two tournaments on the same
eleven vertices differing by two old arc reversals. This does not
give a rule that succeeds for every prime.

## 1. The prescribed target and the actual cancellation

Number vertices 0 through 10, with bit \(j\) of row \(i\) indicating
\(i\to j\). The rows are
\[
\begin{array}{c|l}
30011 & 578,\ 1164,\ 1449,\ 673,\ 1039,\ 659,\ 1086,\ 1105,\ 1787,\ 1238,\ 41\\
30013 & 578,\ 1164,\ 1449,\ 689,\ 1063,\ 643,\ 1086,\ 1105,\ 1787,\ 1238,\ 41 .
\end{array}
\]
Reverse \(4\to3\) and \(5\to4\) in the first. Reversing only the first
gives 29895 (change \(-116\)); reversing only the second gives 27637
(change \(-2374\)). The net change after both is +2, not the sum of
those first-order effects. The mixed finite difference is +2492:
\[
-116-2374+1007+1485=2.
\]

This is a general identity, not a fit to the two counts. In the
mixed difference \(H_{11}-H_{10}-H_{01}+H_{00}\), all ordered paths
testing at most one changed edge cancel. A Hamiltonian path testing
both edges at their common vertex must contain the consecutive
chain \(5\to4\to3\) in the original graph or \(3\to4\to5\) in
the final graph. In each intermediate, both changed edges have
the same direction relative to vertex 4, so neither chain is possible.
The two chain contributions are 1007 and 1485. To verify them
independently, contract a three-vertex chain to a block with the
incoming arcs of its first vertex and outgoing arcs of its last,
and count Hamiltonian paths by subset DP in the resulting digraph.

The generator takes the requested lower count as input. It computes
subset endpoint tables for a fixed-seed source bank. From a first
candidate source (or after a first flip), the exact response of a
reversal \(i\to j\) equals the number of paths using \(j\to i\) after
the flip minus those using \(i\to j\) before. Each is
\(\sum_S E_S(i)S_{\bar S}(j)\) with the appropriate endpoint order,
where \(S\) contains the first endpoint but not the second. The
induced orientations inside the two subsets do not depend on that
cross edge. This searches the second edit without guessing its
response. The 30011/30013 pair occurs at source sample 50071, and
its counts are independently checked from the literal rows.
All 55 single reversals of this particular source miss 30013.
No global minimum edit distance or minimal order is claimed.

## 2. Why a new host was necessary, and an exact arithmetical inverse

Neither 30011 nor 30013 has a one-apex realization at **any** order.
After stripping neutral source/sink padding, a strong one-apex spine
of length \(n\ge3\) has \(H\ge1+2^{n-2}\) by the odd cycles through
the first and last spine vertices. Thus all possible spines have
length at most 16. Exhaustive checking rules out both counts,
and also 29989, the closest prime below 30000. The exhibited pair
needs five vertex deletions per graph to reach a transitive
tournament, by covering all directed triangles.

For one apex, there is now a complete [much narrower two-swap
inverse](TWO_SWAPS.md). Even permitting unequal exteriors does
not lead to a dense family; the only unequal-exterior case is
when the two pairs touch, with \(z=x+1\) for odd \(x\). It gives
the additional twin-prime window 12539/12541. Summed over all
allowed lengths, the classified bank has only \(O(\sqrt B)\)
windows up to \(B\).

There is a separate arithmetic route to the two *individual*
prescribed targets that does not need the seed search, although
it does not put both counts on the same graph. Take an ordinal
join \(C_3\to T\) and reverse one arc \(x\to y\) with \(x\) in
the triangle. For a subset \(S\) of \(T\), let \(E_S(y)\) count
spanning paths of \(T[S]\) ending at \(y\). Put
\[
R_y=\sum_{S\ni y} E_S(y)H(T[V\setminus S]),\qquad
H_{\rm new}=3H(T)+3R_y-S_y(T),
\]
where \(S_y\) counts full paths of \(T\) starting at \(y\).
The lost paths are exactly the ordinal paths whose exposed edge
was \(x\to y\). For each new path using \(y\to x\), split \(T\) into
its before and after runs at the boundary, and choose one of the
three possible nonempty final triangle runs starting at \(x\);
its complement has a unique linear path (including the empty
option). This gives \(3R_y\). Splitting each full path of \(T\)
right after \(y\) injects it into \(R_y\), so \(H_{\rm new}\ge5H(T)\).
Consequently this is not a small way to escape an **already
scaled** ordinal-triangle window; the source must be chosen
from the target first.

Let \(f(z)=H(\operatorname{bin} z)\) be the one-apex word recurrence
from the parent packet, with \(f(0)=1\). There are two readily
invertible slots:

* Take \(T\) with word \(\operatorname{RC}(\operatorname{bin} z)\) for
  positive even \(z\), and take \(y\) to be its apex. Here
  \(R_y=1+z\) (each admissible last selected spine vertex of a run
  can choose any subset before it), and \(S_y=1\). The output is
  \(3(z+f(z))+2\).
* Take the word \(1\,\operatorname{RC}(\operatorname{bin} z)\) for
  \(z\ge1\) and take \(y\) to be its first spine vertex. If \(b\)
  is the first bit of the remaining word, the only incoming
  edge of \(y\) is from the apex. Then \(S_y=f(z)-b\) and
  \(R_y=H(T)+b\): the only split that cannot be concatenated
  uses \(y\) as a singleton followed by the apex-start path,
  if it exists. Also \(H(T)=f(z)+z+1-b\) and
  \(b=1-(z\bmod2)\). Thus the output is
  \(6z+5f(z)+4+2(z\bmod2)\).

The literal inverse works upward from a **requested number** using
\(F_{a,b,c}(z)=az+bf(z)+c(z\bmod2)\). Strip the last bit:
\[
F_{a,b,c}(2u)=F_{2a+b,b,b}(u),\qquad
F_{a,b,c}(2u+1)=a+c+F_{2a,b,0}(u).
\]
At each branch parity and the positive linear lower bound prune
impossible choices. This is an exact inverse for those slots,
not a promise of a short runtime or a preimage for every given
prime. It selects \(z=6246,6558\) for 30011 in the apex slot and
\(z=4727\) for 30013 in the first-spine slot. Independent subset
DP replays all three output tournaments.

## 3. Why 10000 and 30000 actually arise, and what 24 does not prove

In the [previous mirrored +2 construction](../collatz_rows_prime_descent_20261116/MIRROR.md),
replace the middle word \(00\) between the two contact pairs by
\(0\,1^k\,0\), \(k\ge0\). For reflected exteriors with ordinary
binary value \(x\), the lower window endpoint is
\[
N=2^k\big((8x+4)^2-1\big)+2f(x)-x+(x\bmod2)-1-[k=0].
\]
The middle words are \(01\,0\,1^k\,0\,10\) and
\(10\,0\,1^k\,0\,01\). They are palindromes of the same length
\(k+6\), their \(X+Z\) profiles are \(2^{k+6}-1\), and their
counts are \(15\cdot2^k+1-[k=0]\) and two more by the binary
append recurrence. The two exterior effects therefore cancel
exactly. With \(k=0,\ x=12\), the squared address is
\(8x+4=100\), and \(f(12)=11\) gives endpoints \(100^2+7\)
and \(100^2+9\). This explains why the first witness sits near 10000.

More generally take \(k=2j-4,\ x=(5^j-1)/2\), \(j\ge2\).
Then \(2^k(8x+4)^2=10^{2j}\) **exactly**, because the address
is \(4\cdot5^j\). This gives physical +2 windows within \(O(5^j)\)
of the decimal-square anchors: \(1\le f(x)\le x+1\) bounds
the explicit remainder and \(2^k=O(4^j)\). For \(j=3\) the
endpoints are 999999 and 1000001, both composite. It does not
give a prime-selection or a twin-prime infinitude theorem.

Ordinally gluing a directed triangle to the old twin-prime pair
produces \(30021,30027\), with center
\[
3\cdot10008=3\cdot100^2+24=30000+24,\qquad
30000=3\cdot2^4\cdot5^4.
\]
Of course both scaled outputs are composite. The *different*
physical repair at 30011/30013 has center 30012, one half-step
of 24 below the scaled center. Their residues mod 24 are
11 and 13; the old pair had residues 23 and 1. This is wheel
arithmetic, not primality or lattice geometry.

The older [tournament–Leech bridge](../lattices_codes_tournaments_20261021.md#5-score-vectors-permutohedra-and-ii_251-session27--s52-lane-my-cross-checks)
is different and exact in its own domain: transitive 25-vertex
score lifts in \(II_{25,1}\) are null Weyl lifts whose rank-24
quotient is Leech. The [deep-hole
swap](../leech_deep_hole_swap_theorem_20261022.md) exchanges
coset contributions, leaving a controlled theta difference.
That suggests the useful design principle of canceling the
large paired contributions while retaining a small local one.
It does not turn a lattice shell, or the dimension 24, into a
Hamiltonian-path count or a prime selector. Here the two exact
cancellations are the reflected binary profile at 10007/10009
and the shared-vertex two-edge interaction at 30011/30013.

There is also a *literal obstruction in the archived lattice map*,
not merely an analogy. The following two different, score-labeled
order-11 witnesses have counts 30011 and 30013:
\[
\begin{split}
&(160,145,259,519,205,94,1295,876,1083,375,703),\\
&(12,145,1282,102,173,327,1175,813,603,127,955).
\end{split}
\]
Both have the same labeled outdegrees
\((2,3,3,4,5,5,6,6,6,7,8)\) and 38 triangles.
Append an ordinal transitive suffix of order 14 to each. This
preserves \(H\), shifts each old outdegree up by 14, and adds
the same scores \(0,\dots,13\). After the same permutation to
sort the coordinates, their score lifts \(v=(s\mid70)\) are
**identical**, with Lorentzian norm \(-76=-2\cdot38\). No
function of this score-lattice vector alone can distinguish
these neighboring prime counts. This point is not a null Weyl
vector; the assertion is about the common ambient score lift,
not that a Leech lattice quotient has been produced from it.
These witnesses are separate from the two-reversal pair of
section 1. The replay also verifies a score-preserving switch between
these exact rows: reverse the three edge-disjoint directed cycles
(0,5,3), (6,3,9,8,5,4), and (0,7,6,8,10,2).

## Replay

    python3 verify.py
    python3 inverse.py
    python3 score.py
    python3 two_swaps.py
    c++ -O3 -std=c++17 search.cpp -o /tmp/prime-window
    /tmp/prime-window 11 30011 50071

The mathematically useful next question is whether a target-dependent
family can force the mixed interaction to differ from the two
first-order responses by exactly 2, without relying on a finite
trial of source graphs.

Follow-up: [literal Leech octad baselines and a second prescribed prime pair](OCTADS.md).

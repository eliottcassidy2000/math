# A symmetry-preserving route to \(7p\), and where the Paley module fails

**Proved.** For prime \(p\ge11\), a 7-times-prime count \(7p\)
can have an automorphism only if its full automorphism group has
order 7. Its moving orbits must form a *connected* nonuniform
coupling graph. If there is only one moving orbit, the carousel
is its only possible circulant core: a homogeneous Paley
7-module forces an extra factor 9 in the cofactor.

There is a literal escape from that Paley obstruction: couple two
Paley orbits by two offsets. It yields \(H=7\cdot2163467\) with
prime cofactor, and a strong infinite family with an exact padding
inverse. This is independent of the 13-gap determinant repair and
is not claimed to cover all \(7p\).

## 1. The prime and 7-prime symmetry sieve

An automorphism group acts freely on the *ordered* Hamiltonian
paths: fixing one path fixes every vertex. Its order divides \(H\).
Also strong components are canonically ordered in a tournament,
and each Hamiltonian path uses their paths in that order.
An ordinary prime count therefore has just one nonsingleton strong
component. A \(7p\) count with \(p\ne7\) also has just one,
because otherwise one factor would be the forbidden count 7.
Singleton components are individually fixed under automorphisms.

For a prime count \(p\ge5\), the strong core has fewer than
\(p\) vertices, since [Busch's bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
\(H\ge5^{(n-1)/3}\) exceeds \(p\) at \(n=p\).
No permutation of fewer than \(p\) vertices has order divisible
by \(p\), so such a witness is rigid. For \(7p\), \(p\ge11\),
the same bound already exceeds \(7p\) at \(n=p\).
Thus the only possible nontrivial automorphism group has order
7 (and is cyclic). This is a useful asymmetry:
nontrivial symmetry is a viable construction lane for \(7p\),
but cannot be the universal mechanism at ordinary prime targets.

## 2. The coupling graph of moving orbits

Let a generator \(\sigma\) of order 7 act on a witness. Label each
moving orbit by \(\mathbb F_7\), increasing under \(\sigma\).
Between two distinct orbits \(i,j\), the directions are determined
by an offset subset \(S_{ij}\subseteq\mathbb F_7\):
\(i_x\to j_y\) exactly when \(y-x\in S_{ij}\).
A fixed vertex has a uniform relation to each orbit.

Join two moving orbits in the coupling graph exactly when
\(S_{ij}\) is a nonempty proper subset. An independent shift
on each orbit preserves its internal arcs, and preserves
\(S_{ij}\) precisely if the difference of shifts stabilizes
that subset under translation. A proper nonempty subset of
\(\mathbb F_7\) has no nonzero period. Thus the subgroup of
independent shift automorphisms is exactly \(C_7^c\), where
\(c\) is the number of coupling components. In particular
\(7^c\mid H\). A primitive \(7p\) witness with a nontrivial
automorphism must have one coupling component, not several
independently rotatable substructures.

## 3. The two moving-orbit types separate a single slot from all others

If there is only one moving orbit, it is a module: every other
vertex is fixed and sees it uniformly. There are exactly two
circulant 7 vertex tournaments up to isomorphism. Their difference
sets may be taken as \(R=\{1,2,3\}\) (carousel) and
\(P=\{1,2,4\}\) (Paley). Indeed there are eight skew choices
of one of each \(\pm1,\pm2,\pm3\); multiplication by
\(\mathbb F_7^\times\) yields the two quadratic-residue choices
in one orbit and the other six in one orbit.

Let \(F_m(U)\) sum products of ordinary Hamiltonian counts over
\(m\)-colored partitions, with empty color slots allowed. The
[odd-cycle expansion](https://arxiv.org/abs/2307.05569) gives
\[
   F_m(U)=m^7+28m^5+c m^3+d m .
\]
The two regular graphs share 14 directed triangles (use their
outdegrees 3), and all other nonempty packings have one 5-cycle,
two triangles, or a 7-cycle. Thus only \(c,d\) remain. The
literal path subset DP gives \((F_1,F_2)=(175,1988),(189,2016)\),
respectively, so \(c=112\) for both and \(d=34,48\).
In particular
\[
                        F_m(P)-F_m(R)=14m
\]
as a polynomial, not just on checked values of \(m\).

Write \(B_r(U)\) for the number of ways to distribute the orbit
into \(r\) ordered nonempty directed paths. Inclusion-exclusion
from \(F_m\) gives
\(B_r(U)=\sum_{j=0}^r(-1)^{r-j}\binom rj F_j(U)\).
The first two templates therefore have the exact normalized
vectors
\[
\begin{array}{c|rrrrrrr}
 r&1&2&3&4&5&6&7\\\hline
 B_r(R)/7&25&234&954&2160&2880&2160&720\\
 B_r(P)/7&27&234&954&2160&2880&2160&720 .
\end{array}
\]
Every \(B_r\) is divisible by 7 directly, by free rotation of
the ordered slot paths. The vectors agree at every multi-slot
channel, which makes the change predictable in any context.

For any quotient tournament \(T\) with a distinguished slot
vertex \(v\), replace \(v\) by either \(U\) and retain all outside
arcs uniformly. Fix any order of the outside vertices and the
locations of the \(r\) nonempty runs through \(U\).
Legality of the transitions depends on those locations and
the outside order, not on the internal slot paths. Thus for
nonnegative integers \(G_r(T,v)\) independent of the chosen
7-module,
\[
 H(T[v\gets U])=\sum_{r=1}^7G_r(T,v) B_r(U).
\]
For \(r=1\) the single run contracts to \(v\), so \(G_1=H(T)\).
Consequently the Paley replacement increases the normalized
cofactor by *exactly* \(2H(T)\) with all other contributions
unchanged. The Paley normalized cofactor is always divisible
by 9; all multi-run coefficients are multiples of 18 and
the single-run coefficient is 27. For the carousel the
normalized cofactor is \(7H(T)\pmod{18}\), so it at least
does not require divisibility by 3. A primitive \(7p\) witness
cannot hide **any** homogeneous Paley 7-module, not just a
specified single moving orbit.

For example, use a cyclic triangle with modules \(U,TT_r,\{x\}\)
in that cyclic order. The inherited two-module formula
\(H=(F_2(U)-H(U))(F_2(TT_r)-H(TT_r))+2H(U)H(TT_r)\)
follows by whether either module is split into two visits around
\(x\). The carousel cofactor becomes \(259\cdot2^r-209\) and
the Paley one \(261\cdot2^r-207\); at \(r=2\) the carousel
gives \(7\cdot827\). The phase is explicit, but not every prime
belongs to this rail.

## 4. A genuinely nonuniform Paley pair

Take two Paley orbits \(i,j\). Between them use
\(i_x\to j_y\) if \(y-x\in\{0,1\}\), otherwise \(j_y\to i_x\).
This tournament is strong: each internal orbit is strong
and the cut has arcs in both directions. The shared rotation
is an automorphism but neither orbit is a module. The exact
literal counts are
\[
 H=15144269=7\cdot2163467,\qquad F_2=400589504.
\]
The cofactor is prime by trial division. By the sieve the
automorphism group is therefore exactly the simultaneous
\(C_7\). This is not a disguised portal-matching or terrace
example: its degree multiset is seven 5's and seven 8's, hence
no \(n-2=12\) vertex as every outer matched-port template has;
it has 98 rather than 12 directed triangles. It also is not
an ordinal sum.

Use this entire new core as \(U\) in the cyclic triangle with
\(TT_r\) and the singleton. The strong cofactor rail is
\[
                          55063605\cdot2^r-50736671,\quad r\ge1 .
\]
For a prescribed cofactor \(p\), the padding inverse is exact:
\((p+50736671)/55063605\) must be an integer power of two at
least 2. At \(r=2\) the cofactor is the prime 169517749.
This is a concrete nonuniform \(7p\) construction with repeatable
change of size and a sharp inverse. It does not prove infinitely
many cofactor primes on the rail or target selection in general.
The next positive design question is how to change the connected
orbit coupling at a prescribed cofactor while retaining the unit
phase of the carousel, or the equivalent nonuniform phase of
this two-orbit core.

Run `python orbital.py` for the literal subset DP, coefficients,
and trial-division checks.

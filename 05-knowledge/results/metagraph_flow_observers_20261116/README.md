# Information carried by a metagraph line and by a tournament gadget

**What connects the two problems.** A neutral ordinal append on one
side of a framed metagraph line becomes a marked one-arc reader
on the other side. We can recover exactly the rooted profile it
reads, and prove why this operation does not keep a fixed small
nonzero repair as the tail grows. On the gadget side, a bipartite
matching measures precisely how many path-cover coordinates an
unchanged module exposes. A family of matching readers detects
those coordinates over every modulus. This gives a physical
prime output from one rigid \(H49\) birth class and proves another
rigid \(H49\) birth class cannot do the same under *any* homogeneous
reader. These are structural and arithmetic results, not a proof
of Chan–Pak.

## 1. The prior information-flow lenses

The archived [integer-clock observer](../clock_channels_20261021/FINITE_OBSERVERS_INTEGER_CLOCKS.md)
gives an exact coarse-graining test: identical observed states
may still have different successors or different costs. The
[Session33 fixed-path observer](../session33_continuation_states_20261023/metagraph-memory/NOTE.md)
distinguishes the stored frame, its node, the node plus its
partner, and the signed linear observer. Its fixed-size
two-node state works for the complement involution alone,
not for new growth moves. The older
[Klein dictionary](../black_detours_closed_blue_lines_klein_dictionary_20261020.md)
treats the Hamilton path as the clock and the complement as
the charge flip. The [flow and cycle appendix](../flow_decorated_tournament_series_20261021/session23b-flows/PROOF.md)
separates binary cycle supports, nonnegative integer
circulations and actual Hamiltonian paths. In particular,
a group that permutes the support need not act on the directed
graph.

On the growth side the
[birth-fibre packet](../birth_fibers_binary_growth_20261116/README.md)
already distinguishes \(H\), a minimal witness, and the exact
rooted information for a continuation word. Our point is to
make those distinctions operational across the metagraph and
to prove precisely what an exterior gadget can read.

## 2. Ordinal join becomes a one-arc reader across the line

A **frame** lists a directed Hamilton path. Write \(\tau\)
for the involution that complements every nonconsecutive arc
in that frame and leaves the path arcs fixed. A metagraph line
is \(\{T,\tau T\}\) with framed presentations projected to nodes.
This is not full converse and not one single arc flip.

Take frames on \(A,B\) with last vertex \(a\) and first vertex
\(b\), respectively. Concatenate them in the ordinal join
\(A\to B\). The literal result of \(\tau\) is
\(\tau B\to\tau A\) except that the seam \(a\to b\) remains
forward. All other cross edges reverse. For marked \(U=\tau A\),
let \(S\) count paths starting at \(a\) and \(R\) its end-split
response: all weighted prefix paths ending at \(a\), with an
arbitrary possibly empty path on their complement.
For \(V=\tau B\), let \(E_V\) be its end count at \(b\) and
\(R_V\) its start-split response. The [one-arc join law](../../../01-canon/theorems/THM-4866-single-cycle-floor-bottom-band-uniquely-hamiltonian-one-arc-lorentzian-joins-and-crevasses.md)
in these coordinates is
\[
H(\tau(A\to B))=H(U)H(V)-S E_V+R R_V.
\]
Indeed a path avoiding the marked seam is the old ordinal path,
unless it needed \(b\to a\). A new path has four consecutive
possible pieces \(V_{\rm free}|U_{\rm end}(a)\to
V_{\rm start}(b)|U_{\rm free}\). The two outside free pieces
can be empty. Independent choices give the product of the two
split responses. The transfer law is inherited; the exact
series of observations under neutral framed growth is used
below.

### Three neutral padding probes recover three rooted coordinates

Take \(B=TT_m\), a transitive suffix of order \(m\) in its
forward path. This never changes \(H(A)\). Write \(F_m=\tau TT_m\),
the [tribonacci frame of THM-4731](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md),
and \(h_m=H(F_m)\), with \(h_0=h_1=h_2=1\) and \(h_m\) the
tribonacci recurrence. Put \(q_m=\sum_{j<m}h_j\).
A path starting at the first vertex of \(F_m\) is a consecutive
initial segment, so \(R_V=q_m\). Taking \(A\) a singleton
in the join formula gives \(E_V=q_{m-2}\) for \(m\ge2\);
for \(m=1\), \(E_V=1\).

Let \(h=H(\tau A)\), \(P_m=H(\tau(A\to TT_m))\). Then
\[
P_1=h-S+R,\qquad P_m=h h_m-S q_{m-2}+R q_m\quad(m\ge2).
\]
In particular
\[
S=(P_3-3P_1)/2,\qquad R=P_2-P_1-S,\qquad
h=2P_1+2S-P_2 .
\]
Thus two marked coordinates are exactly the extra information
in *all* neutral suffix counts if the current partner count
is already known. Three counts recover all three coordinates
even without that prior count. No amount of these probes
separately recovers the deletion count in general.

For a literal witness, the two order-four frames with free bits
010 and 101 in positions (02,03,13) are both blue, literal
complements of one another, and in the same ordinary H5 node.
Their endpoints have the same classes and color. Adding the
identical neutral sink to each leaves both primary counts 5
but gives partner counts **11 and 13**. The endpoint-class pair
and color do not determine the new readout.

### The missing coordinate is exposed by endpoint deletion

For the rooted converse of \(\tau A\) marked at \(a\), the
previous \(S\) is its full marked end count \(e\), and \(R\)
is its total even plus odd rooted split response. The earlier
growth identity \(R-e=2r\) gives \(r=(R-S)/2\), \(s=(R+S)/2\).
Deleting the last frame vertex commutes with \(\tau\), so the
partner count of the shortened frame gives \(g\). Thus three
padding probes plus one endpoint deletion reconstruct the
**full** \((h,g,r,s)\) needed to predict
[the binary marked growth program](../birth_fibers_binary_growth_20261116/README.md)
on this complementary rooted core. This is a physical bridge
between the two workflows; it is not a catalogue of every
possible marked operation.

For example the primary order-five frames with rows
\((26,4,25,18,2)\), \((22,4,24,19,2)\) both have \(H=9\).
Both have \((h,S,R)=(13,2,18)\), hence they have the same
*entire* neutral suffix count sequence, beginning 29,49,91.
The deletion probe returns 5 and 3. In the rooted companion
the first identical marked append gives 29 on both, but the
second gives 55 and 59. Its formula is
\(h''=h+6r+2e-2g\). Repeating a fixed-root append is not the
same operation as regenerating a neutral suffix on the
primary frame. This is precisely the omitted memory coordinate.

### A small gap cannot keep its size by this neutral-padding scheme

Set \(D(z)=1-z-z^2-z^3\). The series for \(h_m,q_m\) are
\((1-z^2)/D(z)\) and \(z(1+z)/D(z)\). For two fixed marked
cores the partner-count difference therefore has, up to a
polynomial from the small boundary lengths, series
\[
\frac{(1+z)\,[\,\Delta h(1-z)+\Delta R\,z-\Delta S\,z^3\,]}{D(z)}.
\]
The denominator is irreducible by the rational root test,
and its simple positive root \(\rho\) is the unique one of
smallest modulus (triangle equality in
\(z+z^2+z^3=1\) at equality of moduli). If the differences
stay bounded along even an infinite subsequence, that pole
must cancel. The degree-three numerator in brackets has no
\(z^2\) term, so divisibility by \(D\) forces it to vanish.
Thus the marked observables agree and all differences are
zero, or else the difference grows in magnitude like
a nonzero constant times \(\rho^{-m}\).
The order-four +2 is real, but **cannot** be kept
microscopic by repeatedly lengthening that neutral suffix.
Other port choices or changed cores are outside this result.

A complementary control starts on the gadget side: from the
same marked \(H13\) host with \((g,r,s)=(3,6,9)\), the binary
words 100 and 000 both give 97 and even the same degree
multiset, but their \(r\)'s are 54 and 48. With the same
extension of the original frame their line partners have
counts 199 and 197. One case loses track of the rate under
scalar observation; the other loses the frame under node
observation. The two losses are not interchangeable. Moreover the archived
[capacity-transport lemma](../capacity_flow_metagraph_transport_20261021/session23b-capacity/PROOF.md)
gives these equal-score graphs bijective bounded-circulation sets
at every uniform capacity (match the labelled scores, then
complement capacities on the differing arcs). Thus even that
entire exact flow-count sequence cannot recover their distinct
rooted growth and frame-partner counts.

## 3. Exactly how much a homogeneous exterior can read

A homogeneous module \(T\) has all its vertices sending to,
or all receiving from, each exterior vertex separately.
Let \(Q_j(T)\) count **ordered** covers of \(T\) by \(j\)
nonempty directed paths. Set \(D_j(T)=(j-1)!Q_j(T)\).
The sufficiency of the ordinary path-cover profile is
inherited from [THM-3121 and THM-3134](../../../01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md).
Here the exact **integer** readout and visibility cost
follow from a matching.

Partition exterior vertices as \(A\) beaten by the module,
and \(B\) beating it. Make a bipartite graph with edge
\(a b\) iff \(a\to b\) in the host. If its matching number
is \(\nu\), the maximum readable number of module runs
is \(\min(|T|,\nu+1)\), and every smaller positive
run number has a positive coefficient. Indeed each
nonempty internal exterior gap begins in \(A\), ends in
\(B\), and contains an adjacent \(A\to B\) crossing.
A collection of disjoint gaps gives a matching. Conversely
prescribe \(j-1\) matching edges and place their pairs
between \(j\) ordered anchors in a transitive inflation of
the module slot. Insert all remaining exterior vertices,
preserving the prescribed path as a subsequence by the
usual tournament path-insertion lemma. The old gaps still
separate the \(j\) anchors. Thus the coefficient is positive.

Furthermore the coefficient of \(Q_j\) is divisible by
\((j-1)!\): freely permute the \(j-1\) disjoint internal gaps.
Each starts in \(A\) and ends in \(B\), so all orders are
legal, independent of the beginning and end tails.
Every actual host therefore gives a nonnegative integer
linear combination of \(D_j\), not an arbitrary real
combination of ordinary cover counts.

### Matching probes recover the scaled profile over any modulus

The \(m\)-th probe uses \(m\) vertices \(a_i\in A\) and
\(m\) vertices \(b_i\in B\) with \(a_i\to b_j\) iff \(i=j\).
Orient within each side transitively. Its highest module
run count is \(m+1\). To reach that run count use all
exterior vertices in matching length-two gaps, with no
exterior tails. The \(m!\) possible orders are exactly the
coefficient of \(Q_{m+1}\). The probe matrix in the \(D_j\)
coordinates, \(m=0,\dots,|T|-1\), is lower unitriangular
over the integers. Consequently:

- The integer span of actual count outputs is exactly the ideal
  generated by the \(D_j\); the largest common divisor of all
  outputs is their gcd. This does not say that every multiple
  in the ideal is itself an actual count.
- Two modules are indistinguishable modulo any integer in all
  homogeneous readers iff their scaled coordinates agree
  componentwise (pad absent coordinates by zeros).
- With at most \(b\) exterior vertices, only coordinates up to
  \(\lfloor b/2\rfloor+1\) are available, and the matching probes
  attain that bound. If the first detectable disagreement or
  prime-release coordinate is \(j\), the exact minimum budget
  is \(2(j-1)\).

There is an **absolute** cost as well. If \(\nu=0\) the
module slot is a singleton SCC of the quotient, so the
result is a scalar multiple of \(H(T)\). Otherwise it is
in a nontrivial SCC, with at least three once-visit skeletons
and a two-visit skeleton. Thus the least activated count is
exactly \(3H(T)+Q_2(T)\), attained by the two-vertex matching
probe (a directed triangle quotient). More generally if
\(k=\min(|T|,\nu+1)\ge2\), the count is at least
\(3H(T)+\sum_{j=2}^kD_j(T)\). Cutting a full path at
\(j-1\) places gives \(Q_j\ge H(T)\binom{|T|-1}{j-1}\).
These inequalities show why an **unchanged homogeneous**
near-target module cannot supply a small positive
nonmultiplicative correction: it needs selective contacts
or a different source. The exact modulo freedom does not
remove this absolute-price constraint.

## 4. What the group and the prime actually encode

The automorphism group of \(T\) acts freely on every ordered
cover: fixing it fixes each vertex in the ordered paths.
Its order divides every \(Q_j\), hence the readout gcd. But
the gcd is not a group order. The order-seven H49 birth
fibre contains rigid witnesses of both types:

| H49 masks | Common divisor over all readers | First three probe counts |
|---|---:|---|
| (2,4,73,19,39,15,59) | 7 | 49, 861, 15925 |
| (2,12,17,37,11,87,31) | 1 | 49, 873, 16069 |

The first can never become a prime in any homogeneous
reader: every output is a multiple of 7 and is at least 49
by induced monotonicity. The second **literally gives
prime 16069** with four exterior vertices (verified by
trial division and direct path DP). A two-vertex probe
already breaks 7, but gives the composite 873. Equal
birth count and equal trivial automorphism group do not
determine either release or primality. One H343 witness in
the replay first releases 7 only at budget four. Both
rigid H63 birth classes have readout gcd 9, not 63:
a two-vertex reader can remove 7 but not the shared 9.

For a prime \(p\), a useful number-theoretic chart is the
polynomial \(P_T(t)\) which colors all vertices with
\(t\) colors and sums the product of the counts on each
color class (empty class count 1). In the archived
[packing expansion](../../../01-canon/theorems/THM-4756-forward-step-polynomial-mod-3-is-h-times-eulerian.md)
each disjoint odd-cycle pack with \(j\) cycles and
\(\ell\) vertices contributes \(2^j t^{n-\ell+j}\).
It follows also that \(P_T(t)=\sum_j\binom tj Q_j(T)\)
at nonnegative integers. Since \(j!\mid Q_j(T)\),
evaluation at \(0,\dots,p-1\) is triangular and
shows that \(p\) is permanent under homogeneous readers
iff the polynomial vanishes at every element of
\(\mathbf F_p\). On the multiplicative group
\(\mathbf F_p^\times\), exponents can be reduced modulo
\(p-1\), and they all have the same parity as \(n\).
Only \((p-1)/2\) channels need be read for odd \(p\).
In particular at 3 this reduces exactly to \(3\mid H(T)\),
not merely to whether \(T\) has an order-three symmetry.

Also \(P_{A\to B}(t)=P_A(t)P_B(t)\) because
each colored class factors over the ordinal cut.
Nonzero support modulo \(p\) intersects. In the concrete
rigid H35 and H119 graphs of the replay, their
nonzero supports modulo 7 are \(\{\pm2\}\) and
\(\{\pm3\}\), respectively. Each individually has
a reader that escapes divisibility by 7. Their ordinal
join has no such homogeneous reader, since the supports
are disjoint. This is a genuine arithmetic obstruction
created by composition **without** creating any
automorphism. It is the readout profile, not the abstract
group order or the bare product count, that encodes it.

## 5. The limiting missing set: what changes in a proof search

Ordinal padding makes the spectra nested. The assertion
that their limiting missing set equals \(\{7,21\}\)
is exactly the open coverage direction, not something
to assume from the finite census. The existing
[32-bit certificate](../../../session203/coverage-recertification/README.md)
and the [all-moduli constructions](../one_apex_all_residues_20261116/README.md)
have different scopes. For each fixed height \(B\ge3\), the
missing set below \(B\) has stabilized by
\(\lfloor3\log_3 B\rfloor\), since every attainable count
there has a reduced core of that bounded order. That does
not identify the stabilized set at all heights.

The concrete design rule is now two-sided. In the
metagraph, choose and retain the marked coordinates that
control the next read; a +2 at one length can diverge
exponentially under neutral growth. In the coverage
construction, **first test whether the chosen core can
leave the unwanted residue through the permitted reader**,
then test whether activation already overshoots the
absolute budget. For example a prime target near a
multiple of 3 cannot come from a homogeneous reader of
a \(3\)-divisible core, and a target only two above a
large unchanged module cannot come from any activated
homogeneous reader. The tests show exactly where a proof
needs to change the port or the physical core. They do not
rule out selective seams, coordinated editing, or choosing
different birth witnesses from the target.

## Replay

    python3 verify.py
    python3 atlas_table.py ../birth_fibers_binary_growth_20261116/births.tsv

The first runs direct path-count subset DP on the graph
witnesses and covering probes; the second recomputes
the finite release-cost census at the three selected
birth counts.

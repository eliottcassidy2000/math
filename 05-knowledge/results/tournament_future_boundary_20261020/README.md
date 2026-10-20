# Future responses compute tournament counts; run ports compress tribonacci line partners

2026-10-20. **PROVED:** the finite-horizon response quotient and its exact
recurrences, the marked three-letter run carrier, and the two-run collapse.
**FINITE-EXACT:** the attached computation. **CITED/INHERITED:** the Burnside
formula, canonical transitive-chain decomposition, Eisenstein identities and
earlier graceful/square-sum results. **OPEN:** a bounded-state or
polynomial-complexity all-order recurrence, a small grammar for all tournament
cores, and the transfers to fixed-shape graceful existence or Collatz coverage.
No historical-priority claim is made.

The principal advance over the inherited weighted-core packet is an implemented
recurrence that calculates **ordinary tournament and pair-free core counts
without constructing cores or visiting every cycle partition**. It sums
histories having identical responses to every still-possible cycle insertion.
It is an exact compression of Burnside's arithmetic, not a replacement for its
proof. A second advance gives a concrete sparse-defect language that preserves
the chosen Hamiltonian path and its ports under line complementation.

**Concurrent integration.** Mid-session commit `602286dc6f` independently
supplied the same future-observability principle using live divisor
coordinates, together with much stronger fixed-feedback classification
results. The overlap is credited in [incoming.md](incoming.md).
The extensions here are degree-dependent pruning, direct ordinary core and
boundary computations, the order-200 run, and the sharp marked-run carrier.
Combining the incoming feedback viewpoint with that carrier then gives the
exact feedback formula in Proposition 6 below.

## 1. Inheritance and the concept board

The closest proved mechanisms are
[THM-514, unit-tail core kernel and prime-tail ladder](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md),
[THM-516, divisor-profile quadratic form](../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md),
and the incoming
[weighted-core cycle-index hierarchy](../weighted_core_transport_20261020/inheritance.md), §8.
The new work uses all three; rediscovering the odd-cycle derivative hierarchy
alone would not advance that frontier.
The unchanged [all-weight cubic bound](../continuing15_20260909_universal_cubic_bound.md)
and [marked allocation transport](../weighted_core_transport_20261020/README.md)
remain separate parts of the structural toolkit.

The canonical hostile is two histories with the same mass and number of cycles
but different response to a future three-cycle. A second hostile, already in
the weighted-core packet, is two rigid cores with identical cycle index but
different child-automorphism distributions. Our recurrence counts globally;
it does not identify or transport individual child objects.

The corrected near miss is the March tribonacci note's scalar
clique-separator quotient. The least-used relevant sidecars are its actual
interval-packing boundary, July's cancellation-aware graceful monomial
selector, and August's ordered boundary-extension vectors.

The live board is:

1. histories versus their responses to remaining operations;
2. fixed-point-free symmetry mass versus actual pair-free cores;
3. run lengths, endpoint ports, and transitive clone weights;
4. boundary extension vectors and graphical-matroid contractions;
5. odd-cycle edge-sum inversion and observer-dependent zeta identities.

The assigned anchor is A000568 counting and sparse defects. Graceful/square-sum
continuation interfaces are the niche; the analytic/arithmetic E-atlas and
small complement graphs are the wildcard. No LRC(14) consequence is inferred.

## 2. The exact finite-horizon response quotient

Let \(N\) be a requested maximum degree. Process odd cycle lengths in a fixed
order; the implementation normally uses descending order. For the processed
cycle multiset \(\mu\), write

\[
w(\mu)=2^{e(\mu)}/z_\mu,\quad
e(\mu)=\sum_{a\in\mu}(a-1)/2+\sum_{\{a,b\}\subset\mu}\gcd(a,b).
\]

The centralizer denominator is \(z_\mu=\prod_\ell \ell^{m_\ell}m_\ell!\).
Only odd cycles can preserve a tournament: an even cycle has a power that
reverses an unordered pair. Conversely the edge orbits of an all-odd
permutation can be oriented independently. These are the inherited Burnside
facts, giving \(a_n=\sum_{\mu\vdash n,\ \mathrm{odd}}w(\mu)\).

At a given stage, let \(d\) be the already committed degree. Initially
\(d=|\mu|\); the core substitution below introduces a separate excess.
For every **unprocessed** odd \(j\le N-d\), retain only

\[
R_j(\mu)=\sum_{a\in\mu}\gcd(a,j).                 \tag{1}
\]

Histories with the same \(d\) and the same vector \(R\) share one state;
their rational weights are added. No cycle list is retained.

**Proposition 1 (future-response equivalence).** At fixed stage and degree,
equality of (1) is sufficient for every possible continuation to multiply
the two historical weights by the same factor. It is also necessary if
every *individual* permitted continuation weight is to be preserved.

**Proof.** Adding \(k\) cycles of the current length \(\ell\) multiplies by

\[
\frac{2^{k(\ell-1)/2+\ell\binom{k}{2}+kR_\ell}}
     {\ell^k k!}.                                      \tag{2}
\]

The first term counts internal edge orbits, the second counts interactions
between new cycles, and the third counts their interactions with the past.
The remaining profile updates by

\[
R'_j=R_j+k\gcd(\ell,j),\qquad d'=d+k\ell.          \tag{3}
\]

Drop coordinates for processed lengths and those greater than \(N-d'\).
They cannot be queried by any continuation. Equations (2)–(3) inductively
prove sufficiency. If two retained \(R_j\)'s differ, inserting exactly one
\(j\)-cycle and no other future cycles distinguishes their normalized
weights by a nontrivial power of two. This proves necessity for this
wordwise notion, not optimality among arbitrary symbolic summation methods.
QED.

The implementation initializes weight one at \((0,0)\), applies (2)–(3)
for every allowed \(k\), and sums equal successor states. After the last
length, the weight at degree \(n\) is exactly \(a_n\). This is a closed
finite recurrence with no unknown core bank.

**Why this improves the inherited profile.** THM-516 retains a divisor
profile sufficient for *all* possible future lengths. Here a coordinate
disappears as soon as no remaining continuation can observe it. For
example, \((11,5)\) and \((9,7)\) both have mass 16 and two cycles.
Their \(R_3\)'s are 2 and 4, so they must stay separate while length three
is available. Once only length one remains, both have \(R_1=2\) and can
merge. Equal mass and cycle count were insufficient; equal remaining
responses are sufficient.

Descending elimination is a choice, not a theorem of optimality. The
checker also runs ascending elimination through degree 24. The two
answers agree even though the intermediate state counts differ.

## 3. The same recurrence computes pair-free cores and a boundary observer

Let \(Z_T\) be the tournament cycle index and \(Z_C\) the cycle index of
pair-free tournament cores, each including the empty object. The inherited
canonical-chain theorem gives

\[
Z_C(p_1,p_2,\ldots)=
Z_T\left(\frac{p_1}{1+p_1},\frac{p_2}{1+p_2},\ldots\right). \tag{4}
\]

The full odd-cycle differential hierarchy was already proved in the incoming
packet; it determines the missing boundary but had no demonstrated
computational compression there.

Specialize \(p_\ell=x^\ell\) in (4). When choosing \(k>0\) actual
\(\ell\)-cycles, additionally choose substitution excess \(s\ge0\):

\[
\left(\frac{x^\ell}{1+x^\ell}\right)^k
 =\sum_{s\ge0}(-1)^s\binom{k+s-1}{s}x^{\ell(k+s)}. \tag{5}
\]

Multiply (2) by the signed binomial coefficient, set
\(d'=d+\ell(k+s)\), and update the response by **\(k\)**, not \(k+s\).
The excess is not another automorphism cycle. For \(k=0\), only \(s=0\)
is allowed.

**Proposition 2.** This signed recurrence returns the number \(c_n\) of
pair-free unlabelled tournament cores. Omitting length one returns

\[
B(x)=Z_C(0,x^2,x^3,\ldots)
 =\sum_{[Q]}\frac{\#\{g\in\operatorname{Aut}Q:
                     g\text{ fixes no vertex}\}}
                    {|\operatorname{Aut}Q|}x^{|Q|}.       \tag{6}
\]

**Proof.** Expand the finite truncation of (4) using (5). Its interaction
weights depend only on actual cycle multiplicities, so Proposition 1 applies
unchanged to the excess-graded states. The standard ordinary cycle-index
specialization proves the assertion about \(c_n\). Setting \(p_1=0\)
retains exactly fixed-point-free automorphism terms, proving (6). QED.

This computes the required ordinary counts directly. It does **not** recover
the full multivariate boundary from the one-variable series \(B(x)\), and
does not treat its rational coefficients as counts of a subclass of cores.
For example,

\[
B(x)=1+\tfrac23x^3+\tfrac45x^5+\tfrac{10}{9}x^6+
\tfrac87x^7+\tfrac{16}{15}x^8+\tfrac{1958}{81}x^9+\cdots.
\]

The core counts start

`1, 1, 0, 1, 2, 5, 26, 255, 4602, 147095, 8228360,
814376271, 144675514036` (orders zero through twelve).

## 4. Exact computation and its boundary

The [checker](verify.py) and [receipt](receipt.json) use rational arithmetic.
The final normal and optimized runs agree byte for byte and execute
34,705 always-active checks. [Validation scope](VALIDATION.md).

| Quantity | Degree bound | Maximum stored states | Transition contributions |
|---|---:|---:|---:|
| A000568 | 100 | 2,953 | 90,219 |
| pair-free core classes | 100 | 5,480 | 1,697,975 |
| fixed-point-free core mass \(B\) | 100 | 4,784 | 280,479 |
| A000568 extension | 200 | 96,836 | 2,994,814 |

Controls:

- Direct, unmerged odd-partition Burnside summation through 60.
- Independent direct polynomial substitution through 18 for all three series.
- Opposite elimination order through 24.
- Direct enumeration of invariant orientations and homogeneous pairs for
  every odd permutation cycle type through six vertices.
- The previous actual core bank through eight vertices.
- All 77 entries of Keith Briggs's
  [A000568 table, orders 0–76](https://oeis.org/A000568/b000568.txt);
  the exact retrieved table is [stored here](a000568-reference.txt).
- Increasing the bound to 200 preserves every earlier coefficient.

The checked answer through 200 is generated by a proved recurrence; the
external-table control stops at 76, not 200. Normal and optimized Python
replays are compared byte for byte. No independent-agent audit is claimed.
The receipt includes all coefficients and the reference hash, not merely
integrality checks or final decimal digits.

This establishes a practical compression at the tested sizes. It does not
establish polynomial asymptotic complexity or a scalar recurrence on the
previous \(a_n\)'s alone. Burnside already avoided enumerating tournament
classes; the further gain here is merging many distinct permutation
histories before their tails are processed.
For scale, the numbers of odd cycle partitions at masses 100 and 200 alone
are 444,793 and 487,067,746. These are partition counts, not memory bounds;
the table above reports the different quantities actually used by the DP.

## 5. A marked sparse-defect grammar for the tribonacci family

Let \(T^*_n\) have \(i\to i+1\), and \(j\to i\) for \(j\ge i+2\).
The March proof and
[THM-4731, statement 6](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md)
already establish its tribonacci Hamiltonian-path count, rigidity for
\(n\ge4\), and blue multiplicity. Those conclusions are inherited.

**Proposition 3 (run/port form).** Hamiltonian paths correspond bijectively
to compositions \((a_1,\ldots,a_r)\) of \(n\) with no consecutive parts
both equal to one. Decode by taking consecutive label intervals in
descending interval order, but traversing each interval increasingly.
The initial and final vertices are
\[
(\text{source},\text{target})=(n-a_1,a_r-1).          \tag{7}
\]

**Proof.** The only ascending arcs increase by one, so every maximal
ascending run is an interval. Disjoint intervals must appear in descending
order: an ascending jump is either forbidden or merges the runs. Across
two adjacent runs, the downward jump has size \(a_i+a_{i+1}-1\);
it is allowed exactly when the two sizes are not both one. The inverse
construction and (7) follow immediately. QED.

Retaining \(u\) for the number of runs gives the exact series
\[
\sum_a u^{r(a)}z^{|a|}
 =\frac{1+uz}{1-u\,\frac{z^2}{1-z}(1+uz)},           \tag{8}
\]
including the empty composition. For \(u=1\) this is
\((1-z^2)/(1-z-z^2-z^3)\).
More finely, with \(r-s>0\) and \(n-2r+s\ge0\), the number with \(r\)
runs, \(s\) singleton runs, is
\[
\binom{r-s+1}{s}\binom{n-r-1}{r-s-1}.               \tag{9}
\]
Otherwise it is zero, apart from the one-singleton composition of one.
Equation (9) chooses the separated singletons and then distributes the
excess above two among the other runs.

Reflection/reversal sends the composition to its reverse. Thus the
anti-palindromic paths are exactly palindromic compositions; the checker
recovers the already-proved half-size blue law. Endpoint rows, rather than
only (8), are retained when a subsequent operation uses the marked path.

### The line-complement carrier

Let \(B_\uparrow\) be the increasing transitive tournament, \(I\) the set
of adjacent-label pairs, and \(E(P)\) the unordered edges of the chosen
Hamiltonian path. The line partner preserving that path is

\[
L_P(T^*_n)=(T^*_n)^{op}\mathbin{\triangle}E(P)
         =B_\uparrow\mathbin{\triangle}(I\triangle E(P)). \tag{10}
\]

There are exactly \(2(r-1)\) defects: one adjacent cut edge and one long
joining edge at each run boundary. Long joining edges have length at least
two, so these two collections are disjoint.

**Proposition 4 (three-letter weighted carrier).** Replace every run
length \(a\ge3\) by three, recording the middle vertex's weight \(a-2\).
Keep lengths one and two unchanged. The line partner of the original path
is exactly the transitive blowup of the line partner of this capped path.
The blowup restores the original labelled graph and marked path.
In particular, its canonical pair-free core has order at most \(3r\),
independently of the total \(n\).

**Proof.** Every defect in (10) touches endpoints of runs only. Strict
interior vertices of one run have identical adjacency to everything outside
that interior, and form an increasing transitive module. Replace that module
by one vertex; the endpoint defect pairs become exactly those of the capped
path. Expanding that vertex into its recorded increasing chain restores
every edge and every position of the chosen path. Canonical chain
contraction then gives the same pair-free core before and after expansion.
QED.

The normalized *syntax* has alphabet \(\{1,2,L\}\), forbids adjacent `11`,
and attaches a positive integer to each `L`. It is not a canonical
isomorphism classification: different words may still have isomorphic
cores. There are \(R_r\) permitted type words, where
\[
R_0=1,\quad R_1=3,\quad R_r=2R_{r-1}+2R_{r-2}.
\]
### Sharp core size: exactly \(3r-2\)

**Proposition 5.** Among all permitted \(r\)-run codes, the largest
canonical core order of the line partner is exactly \(3r-2\), for
every \(r\ge1\).

**Proof of the upper bound.** Work with the capped word. Its vertex
count is \(3r-\delta\), where \(\delta=2\,\#1+\#2\).
If \(\delta\ge2\), the bound is immediate.

If \(\delta=0\), all blocks have length three. For \(r\ge3\), number
blocks from bottom to top by \(k=0,\ldots,r-1\), with vertices
\(L_k=3k,M_k=3k+1,R_k=3k+2\). The pairs
\(\{R_0,M_1\}\) and \(\{M_{r-2},L_{r-1}\}\) are homogeneous.
For \(r=3\) they form one three-vertex chain; for larger \(r\) they
are disjoint. They reduce the core by at least two.
For \(r=1\) the child is transitive; for \(r=2\) the two-run law
below gives a core of order three.

If \(\delta=1\), there is exactly one length-two block. For \(r\ge4\),
at least one of the second block from the bottom and the second block
from the top still has length three. At that end, the extreme block's
inner endpoint and the neighboring block's middle vertex are a
homogeneous pair. This follows directly from (10): the endpoint has
one defect, to the intervening vertex, while the middle has none.
Thus at least one vertex is removed from \(3r-1\). For \(r=3\), the
same argument works unless the middle block has length two; the
remaining code `(3,2,3)` has homogeneous pair \(\{2,5\}\).
For \(r=2\) use the two-run law; \(r=1\) is transitive.

**Attainment.** For \(r=1\), a transitive child has core order one.
For \(r=2\), code `(1,3)` gives a directed triangle plus a sink,
which is pair-free of order four. For \(r\ge3\), take every part
equal to three. We show that the two displayed pairs are the only
homogeneous contacts, so exactly two vertices are removed.

For a transitive tournament with defect graph \(D\), vertices \(i<j\)
are homogeneous exactly when
\[
\bigl(N_D(i)\triangle N_D(j)\bigr)\setminus\{i,j\}
   =\{i+1,\ldots,j-1\}.                            \tag{12}
\]
This compares adjacency inside and outside the interval. Here the
defect edges are \(\{R_k,L_{k+1}\}\) and
\(\{L_k,R_{k+1}\}\). Middle vertices have no defects.
Two middle vertices cannot satisfy (12). Two `L` vertices have
symmetric difference supported on `R` vertices only, while their
intervening interval contains a middle; the same applies to two `R`s.
An `L` and an `R` also have an intervening middle unless they are
the adjacent endpoints across a block cut. At that cut their
remaining defect neighborhoods cannot both be empty when \(r\ge3\),
so they are not homogeneous.

Finally, an interior-block endpoint has two defect neighbors six
positions apart, which cannot be the consecutive interval in (12)
when paired with a defect-free middle. Among the four outer-block
endpoints, exactly \(R_0\) and \(L_{r-1}\) have the one intervening
defect neighbor needed for the two stated contacts. This exhausts
all pairs. QED.

The complete three-letter census for \(r=1,\ldots,8\) independently
gives `1,4,7,10,13,16,19,22`, agreeing with the all-order proof.

**Two-run law and the missing split.** If \(a,b\ge2\), every code \((a,b)\)
has line partner \(C_3[1,a+b-2,1]\). Indeed the only defects from
\(B_\uparrow\) are \(\{0,n-1\}\) and \(\{b-1,b\}\).
The middle vertices are transitive in the order obtained by interchanging
the adjacent \(b-1,b\); the two extreme vertices and the middle chain
form a directed triangle. The target core and weights therefore forget
the split \(a,b\), whereas the ports are \((b,b-1)\).
That split is essential for reversible marked transport.

### Exact feedback when every run has length at least two

**Proposition 6.** If every \(a_i\ge2\), then the minimum feedback arc
number of \(L_P(T^*_n)\) is exactly \(r-1\), independently of the
run lengths.

**Proof.** Write the blocks in increasing numeric order as
\([L_i,R_i]\), \(i=0,\ldots,r-1\). The defect pairs are
\(\{R_i,L_{i+1}\}\) and \(\{L_i,R_{i+1}\}\).
Swap each adjacent pair \(R_i,L_{i+1}\) in the increasing vertex order.
Because every block has at least two vertices, these swaps are disjoint.
They remove precisely the \(r-1\) short backward arcs without changing
the relative order of any other pair. The long defects remain backward,
so this order has cost \(r-1\).

For each \(i<r-1\), the vertices
\((L_i,R_i,R_{i+1})\) form a directed triangle. The three vertices are
distinct, the first two arcs are forward, and the last is the long
defect. These triangles have pairwise disjoint edge sets, even where
successive triangles share a vertex. Every total order has a backward
arc in each triangle. Thus its cost is at least \(r-1\). QED.

This is a concrete normalization of the incoming
[fixed-feedback grammar](../defect_boundary_20261020/sparse_theory.md).
The displayed reference order had \(2r-2\) backward arcs; its true
minimum is half that. For all-three words with \(r\ge3\), the core
has order \(3r-2=3f+1\). The general incoming sharp bound is \(4f+1\);
these are different families, not contradictory constants.

**Coverage limit.** This grammar covers all line partners of all marked
tribonacci paths. It is not a grammar for all tournaments. The sources
\(T^*_n\), \(n\ge4\), are rigid, so they themselves contribute no
nonidentity fixed-point-free automorphisms. The two constructions in this
note solve different tasks: global symmetry-weighted counting and exact
transport in a restricted object family.

## 6. The three small graphs locate the missing information

There are precisely three unlabelled simple graphs with four vertices and
three edges: \(P_4\), \(K_{1,3}\), and \(K_3\sqcup K_1\).
Their labelled multiplicities are 12, 4, 4. Complementation fixes the
isomorphism type of \(P_4\) and exchanges the latter two.

For the edge-sum map \(x\mapsto(x_u+x_v)_{uv\in E}\), all three have rank
three over a characteristic-zero field. Yet their one-dimensional kernels,
with suitable labelled representatives, are

\[
(1,-1,1,-1),\qquad(1,-1,-1,-1),\qquad(1,0,0,0).
\]

Thus the same rank and nullity describe three different locations of
uncertainty. On the star, changing one root value changes every label.
On the triangle-plus-isolate, the triangle is completely determined and
only the isolate is free. This extends the organizing use of the
[E-atlas signless-incidence lemma](../e_atlas_20261019.md), not its
existence conclusions.

It explains a concrete historical repair. The odd-cycle conflict graph
of \(T^*_9\) contains an induced claw: the interval `[0,8]` meets the
three disjoint triangle intervals `[0,2]`, `[3,5]`, `[6,8]`.
It is an interval graph but is not claw-free. Moreover, gluing two edges
along a vertex gives \(P_3\), and
\[
I(P_3,t)=1+3t+t^2\ne(1+2t)^2/(1+t).
\]
The March scalar clique-cutset quotient is false.
For a clique separator \(K\), the correct carrier retains its state:
empty, or the identity of its one occupied vertex. Sum the products of
extension counts over those compatible states. The number of states
cannot generally be replaced by a quotient of total partition functions.

## 7. Graceful, square-sum, quadratic, and zeta transfers

The precise reusable move is to compare **continuations**, not merely
statistics of the objects already built.

For a matroid with future ground set \(R\), partial independent selections
\(I,J\) have identical admissible future subsets exactly when
\[
(M/I)|R\quad\hbox{and}\quad(M/J)|R
\]
have identical independent sets. This follows from the definition:
\(X\subset R\) is independent in \(M/I\) iff \(I\cup X\) is independent
in \(M\). In a graphic matroid, retaining the connected-component
partition on future attachment vertices is sufficient. Components having
no future incident edge can be forgotten for this continuation test.
This is the exact graphical counterpart of (1), with a different predicate.

That observation connects to yesterday's stronger
[quadratic upper-tail reconfiguration diameter](../graceful_polynomial_diameter_20261019/POLYNOMIAL_STAIRCASE_DIAMETER.md):
tight suffixes are contracted, and completed inner flats can be fixed.
It supersedes the weaker exponential-diameter route. The
[factorial-fiber construction](../graceful_factorial_fibers_20261019/GRACEFUL_FACTORIAL_FIBERS.md)
shows why one-color moves are insufficient even when two-color buffers work.
None of these operations preserves an arbitrary requested tree shape.
That extra constraint must remain a sidecar in a graceful-existence program.

For square-sum Hamiltonian constructions, connectivity of partial pieces
must be supplemented by endpoint pairing, degrees, and actual labels.
The [exchange study](../graceful_square_exchange_20261019.md) already
distinguishes endpoint rotations from two-end splices; its endpoint-blind
degree-two forcing was previously corrected. A component quotient alone
does not certify that a spanning path, or an arithmetic contact, extends.

The Petersen/\(K_5\)/\(K_{3,3}\) distinction also becomes exact through
edge-sum continuation: Petersen and \(K_5\) are connected nonbipartite
graphs and have injective signless incidence; \(K_{3,3}\) has one
alternating gauge. The
[integral E4/Petersen lattice identification](../e_atlas_20261019.md), §3,
is stronger than a coincidence of vertex counts and is retained.

For a quadratic three-cycle, the triangle's edge sums
\(s_i=x_i+x_{i+1}\) give \(x_i=(s_i-s_{i+1}+s_{i+2})/2\).
The secant product is one, but the derivative multiplier is
\(8x_0x_1x_2\). At \(c=-29/16\), the cycle
\((-7/4,5/4,-1/4)\) has multiplier \(35/8\); at \(c=-7/4\), the
algebraic period-three cycle has multiplier one. These inherited exact
identities distinguish a cycle's combinatorial return from its perturbation
response. They do not identify a return map with a tournament automorphism.

Finally the direct analytic \(E_3,E_4,E_5\) connection is already stronger
than a numerical analogy. With
\[
E(z,s)=\tfrac12\sum_{\gcd(c,d)=1}y^s/|cz+d|^{2s},
\]
its constant Fourier coefficient is
\[
y^s+\frac{\sqrt\pi\,\Gamma(s-\tfrac12)\zeta(2s-1)}
                 {\Gamma(s)\zeta(2s)}y^{1-s}.             \tag{11}
\]
Thus \(s=3,4,5\) intrinsically selects \(\zeta(5),\zeta(7),\zeta(9)\).
The normalization was rechecked against
[Bell, Theorems 1 and 3](https://jordanbell.info/LaTeX/mathematics/eisenstein/).
The [AP square-sum height series](../ap_square_odd_zeta_20261019.md) is
another exact observer of the same Eisenstein family at \(z=i\), but uses
parameters \(5,7,9\) for those zeta values. Its angular positivity filter
must accompany the height when only positive graph labels are allowed.
Neither observer supplies a new irrationality result.

## 8. Next questions with decisive tests

1. **Profile width.** Bound the number of distinct feasible future-response
   profiles, or choose an elimination order provably better than descending.
   The observed table is evidence, not a complexity theorem.
2. **Run normalization.** Classify precisely when different three-letter
   words give the same weighted core, now that the sharp \(3r-2\)
   maximum is proved. Ports must remain recoverable.
3. **Symmetric sparse families.** Replace the rigid tribonacci source by
   an explicitly equivariant defect carrier. Test its fixed-point-free
   cycle-index terms against (6), preserving the action rather than only
   the underlying graph.
4. **Arithmetic continuation interfaces.** For a bounded set of free
   graceful colors or square-sum attachment ports, derive the exact
   contraction/endpoint response and test whether it stays bounded under
   refinement. A rank, total count, or graph minor alone fails the test.

Reproduce from this directory:

```text
PYTHONINTMAXSTRDIGITS=0 python3 verify.py --max-order 100 --large-order 200 > receipt.json
PYTHONINTMAXSTRDIGITS=0 python3 -O verify.py --max-order 100 --large-order 200 > /tmp/tournament-replay.json
cmp receipt.json /tmp/tournament-replay.json
```

Use Python 3.10 or later. The environment setting permits decimal output
of the large computed integers on Python versions with a digit limit.
The more detailed [inheritance map](inheritance.md) records the historical
side themes and precisely what was retained or rejected.

# Sparse-defect grammars and the fixed-point-free boundary

20 October 2026. Continuation of
[weighted cores and allocation transport](../weighted_core_transport_20261020/README.md).
**PROVED:** the scoped constructions and counting theorems linked below.
**FINITE-EXACT:** their declared compilations and independent controls.
An efficient structural recurrence covering unrestricted A000568 remains
**OPEN**. This packet does not claim a Collatz proof.

## The strongest advance: finite templates count every fixed-defect layer

Let `a_(n,<=k)` count tournament isomorphism classes whose minimum feedback
arc number is at most `k`. The [proof and exact compiler](sparse_theory.md)
establish:

\[
\begin{gathered}
\text{every automorphism moves at most }3k\text{ vertices},\\
\text{the canonical transitive-clone core has at most }4k+1\text{ vertices},\\
A_k(x)=\sum_n a_{n,\leq k}x^n
       =1+\sum_{\text{finitely many capped template orbits }s}
           \frac{x^{v(s)}}{(1-x)^{h(s)}}.
\end{gathered}
\]

Only the defect endpoints, the intervening chains, and weights capped at `k+1`
are needed. Automorphisms fix every unbounded weight coordinate individually.
There are at most `2k+1` such coordinates. Therefore, for each fixed `k`,
the counts are eventually polynomial of **degree exactly `2k`**, with leading
coefficient `1/(2^k k!)`. Perfect matchings on `2k` ordered endpoints classify
the leading terms.

[The next layer is also classified](next_coefficient.md): shorten one gap,
omit one gap, or let two defects share an endpoint. These give the universal
second coefficient:

\[
a_{n,\leq k}=
\frac{n^{2k}-3k\,n^{2k-1}}{2^k k!}+O_k(n^{2k-2}).
\]

The first two nontrivial all-order outputs are

\[
a_{n,\leq1}=1+\binom{n-1}{2}\quad(n\geq1),\qquad
\boxed{a_{n,\leq2}=
\frac{n^4-6n^3-n^2+62n-72}{8}\quad(n\geq4).}
\]

The budget-two formula comes from 4,032 ordered defect templates, normalized
to 2,691 capped weighted-core orbits. It is not an interpolation. An independent
complete tournament census and feedback subset DP agree through order seven.
This construction avoids generating all cores at each target order. Its
cost as the **budget** grows is the remaining difficulty.

## Three complementary advances

1. **A000568 through 100 without constructing its cores.**
   [The live divisor transfer](boundary.md) retains only divisibility
   coordinates that a future cycle length can observe. At bound 100 it has
   3,663 peak states and 66,404 transitions, versus 444,793 unmerged
   multiplicity states. Independent partition calculations through 40 and
   the previously recorded 1,333-digit value at 100 agree. The same note
   gives signed ordered-list clone removal and compatible differential
   equations for every cycle coordinate. This still uses the classical
   fixed-permutation edge-orbit weights; it is not Burnside-free.

2. **The tribonacci grammar survives one nonlocal defect.**
   [The two-tail chord family](grammar.md) reverses one distance-two chord in
   the inherited tribonacci tournament. Its Hamiltonian-path counts satisfy
   tribonacci in each tail, with an exact rational bivariate generating
   function and a necessary three-dimensional transfer. Its simple
   complement partner is `TT_m -> C3 -> TT_r`, connecting this grammar to
   actual metagraph leaves. Centered defects preserve the blue count even
   while changing the full path count.

3. **A forgotten local-interaction question has an all-order answer.**
   [The signed Hessian proof](web.md) shows that mixed Hamiltonian-path
   differences on two incident edges never vanish. Their signs recover the
   product of the two edge orientations. A spanning tree of the line graph
   therefore reconstructs the entire tournament up to converse; one edge
   direction restores the last bit. This proves the previously empirical
   June endpoint-star property. It changes the representation, not the number
   of required bits or the cost of computing Hamiltonian-path counts.

## The four-vertex observation, made exact

The user's three graphs are precisely `P4`, `K_(1,3)` and `K3 + K1`.
Graph complementation exchanges the last two and preserves the type of `P4`.
The ambiguity also appears in their line graphs: both the claw's and triangle's
three edges give `K3`.

[The tetrahedral calculation](web.md) distinguishes graph complement from
opposite-edge duality. Signed local interactions distinguish vertex stars
from triangular faces, restoring the incidence information lost by the
unsigned octahedron. Complement also swaps independent twins with adjacent
twins, so a proposed cancellation of weighted graph defects must transport
that module type. The independent triple blowup of `P4`, of defect +2, is a
concrete hostile: quotient complement preserves it instead of negating it.

## Inheritance, archive excursions, and the live board

The portfolio is **Anchor:** the automorphism boundary and exact class counts;
**Niche:** sparse path grammars and actual line partners; **Wildcard:** the
four-vertex complement/duality observation and signed reconstruction.

The six live concepts are **canonical chains, capped feedback, automorphism
cycles, two-tail path transfer, signed local relations, future observability**.
Their roles are different:

| Source | Map | Preserved target | Lost data / necessary sidecar |
|---|---|---|---|
| Ordered sparse tournament | Endpoint-gap template, chain normalization | Feedback budget and class expansions | Redundant order removed by weighted isomorphism |
| Automorphism cycle history | Live divisor profile | Every future continuation weight | Discard only expired coordinates |
| Clone cycle index | Signed ordered-list transform | Fixed-permutation pair-free count | Repeated cycle lengths and phase powers retained |
| One-chord path family | Three-state tail transfer | Hamiltonian paths | Not the number of tournament classes |
| Tournament edge bits | Signed adjacent Hessian | All arcs up to converse | One orientation bit and incidence addresses |
| Tetrahedral graph | Complement or opposition | Explicit character/parity laws | The operations and twin types cannot be conflated |

Closest inherited mechanisms: THM-4463's homogeneous chains, THM-514/516's
unit-tail and divisor kernels, and the March tribonacci run proof as recovered
by THM-4731. The canonical hostiles remain the six-core examples from the
previous packet: equal full cycle index does not determine individual growth.
The corrected near miss is THM-4730's affine-only anti-automorphism test;
its alphabet consequence remains provisional here.

The archive excursions deliberately recovered side themes:

- [June observer-extension payload](../../../07-reflections/observer-extension-payload-codex-s218.md):
  rootless `[3,3]` coupling and the rule that forgotten data must be irrelevant
  to every future operation. Here that rule becomes an exact variable
  elimination algorithm.
- June [endpoint-star](../../../07-reflections/character-ratio-is-blind-to-dichromatic-the-OCF-governs-both-it-and-the-delta-field-s683c.md)
  and [arc interaction](../../../07-reflections/arc-interaction-graph-refines-H-s686.md):
  support versus signed response. The new block-path argument resolves their
  specific adjacent-edge nonvanishing question.
- August [tetrahedral representation atlas](../../../07-reflections/tetrahedral-k4-representation-atlas-tournaments-xor-berggren-h1-codex-20260816.md):
  incidence, opposition, parity and the extra orientation bit. The new map is
  explicit; no geometric ancestry equivalence is inferred.
- [THM-3079's Laguerre row transform](../../../01-canon/theorems/THM-3079-laguerre-pf-row-transform-and-strict-integer-mesh-terminal-minus-one.md):
  its derivative-weighted coefficient row is exactly the ordered-list array
  needed for clone removal.
- [The March zeta reflection and its correction](../../../07-reflections/riemann-zeta-tournament.md):
  shared prime labels and tribonacci trace values do not transport a counting
  predicate. Only the actual divisor identity and transfer matrix are used.

No historical reflection is treated as proof of its surrounding speculation.
The maps above and their hostile controls specify what has been retained.

[A late incoming commit](incoming.md) supplied boundary-support and metagraph
degree work. Its 105-matching atlas gives an explicit new dictionary: the
105 leading budget-four defect templates correspond to the 105 four-chart
partitions of `H8+I2`, with the same `7+42+56` split after a binary address is
fixed. The tournament's endpoint order remains essential. The incoming
support-versus-weight hostile also explains why path existence alone cannot
replace the weighted counting transfer.

## Boundaries that survived hostile testing

- A defect polynomial for a fixed **ordered** permutation is explicit, but
  conjugate permutations can have different polynomials in that same order.
  [Minimal four-vertex witness](defect_polynomial.md).
- Equal defect count and equal score sequence do not determine the two-defect
  path count: the order-six examples have counts 25 and 29.
- The live divisor bank cannot omit `M5` while length five remains possible:
  histories `(7,15)` and `(9,13)` give continuation increments eight and four.
- An unsigned interaction graph is orientation-blind at order four.
- Fixed-budget polynomial recurrences do not become a fixed scalar recurrence
  for A000568 when the budget grows.

The next strongest target is a structural rule for adding one defect to the
finite endpoint grammar, including shared endpoints, missing gaps and clone
mergers. The top two polynomial degrees are now classified. The lower-degree
collision rules are the remaining work.

See [VALIDATION.md](VALIDATION.md) for exact universes, commands, outputs,
independent audits and Git integration notes.

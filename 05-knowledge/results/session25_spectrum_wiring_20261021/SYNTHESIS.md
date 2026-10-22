# Session 25: exact interfaces, the Hamiltonian spectrum, and prime-power closure

Session started 21 October 2026; research continued into 22 October UTC.

## Summary and scope

This session gives several all-order results, not a claimed solution of the
complete Hamiltonian-count spectrum or merged-metagraph classification.
All principal proofs below have independent mathematical review and exact
controls. Computational universes and remaining limitations are stated in
the linked packets. No external literature-priority claim is made.

The common mechanism is **retaining exactly the information needed by the
next operation**: a path's ports, a parity sheet, an equality partition, a
true clock, or an exponent in a finite group algebra. A numerical resemblance
or a matching group order is not enough.

The strongest closures are:

1. Every parity-admissible residue modulo every integer is attained by a
   strong tournament with feedback vertex number one, even within an
   extremely sparse explicit family.
2. The odd girth of the **labelled** tournament complement graph is determined
   at every order. A separate, exact two-sheet covering explains the
   parity information lost on merging isomorphism/converse classes.
3. HYP-9358(d′), the additive prime-power quasi-periodicity of strong labelled
   tournament counts, is proved at all prime powers. Its CRT consequence
   removes the earlier squarefree/cube-free restrictions.
4. Named port signatures give exact wiring/colouring interfaces, while a
   six-colour two-port assembly cannot manufacture equality from soft
   biases: its hard-inequality skeleton must already be six-chromatic.

## 1. The spectrum has no finite congruence obstruction beyond parity

Let H(T) be the actual number of directed Hamiltonian paths, not H(T)/Aut(T).
For every positive modulus M, every residue compatible with H being odd is
attained by a **strong** tournament that becomes transitive after deleting
one vertex. Each such residue contains infinitely many distinct attained
counts.

The construction inserts one vertex into a transitive tournament with
binary orientation word (1^L 0^L)^t, L,t>=1. Put Q=2^L. Its exact count is

    H = 1 + Q²(Q^(2t)−1)/(Q+1)² + 2t/(Q+1).

Choose L>=max(1,v2(M)) and Q=1 modulo the odd part of M. At block boundaries,
the exact affine count state resets modulo M and each block adds

    D = (Q−1)²+1,    gcd(D,M)=gcd(2,M).

Positive t therefore traverses every allowed residue. First bit1 and last
bit0 guarantee strongness, not merely an arbitrary one-apex presentation.

Crucially this same two-parameter family has only O(log X log log X) distinct
counts at most X. Full congruence coverage can coexist with extreme numerical
sparsity. Thus no finite collection of residue tests, including a combined
CRT test, can distinguish the missing values 7 and 21 from every attainable
count. This is an obstruction to a proof strategy, not a construction of
all other odd integers.

The older THM-4094 two-prime-lane reduction and its much larger finite
constructive ranges remain relevant. This session does not present a small
new census as an improvement on those ranges.

[Proof and exact recurrence](session25-spectrum/PROOF.md) ·
[Independent spectrum review](session25-spectrum-review/REVIEW.md).

## 2. Complete finite wiring certificates, with ports retained

Two tournaments A,B have a complement line exactly when a vertex bijection
makes their agreement digraph one directed Hamiltonian path. Counting the
compatible bijections, dividing by the correct automorphism factors, and
then treating converse reversal gives exact blue/black multiplicities and
loops. This is a necessary-and-sufficient finite realization certificate;
it is not an all-order description of the graph's shape.

For substitution, a path may revisit a block in several runs. Maximal-run
path covers together with actual boundary entry/exit vertices reconstruct
the complemented substituted tournament exactly. Cover counts alone do not.
The old triangle-lift failure, H=2983 versus27, is thereby explained at the
operation level rather than hidden behind matching degrees.

A local directed triangle has ordinary C3 symmetry and signed S3 symmetry.
Signs in different modules must synchronize; arbitrary independent local
reflections do not define a global tournament/converse symmetry.

The two independent wiring engines agree on every class pair through order 6.
At order 6 there are 34 merged nodes, 512 raw lines, 32 blue lines and 240
fundamental black edges.
The distinction between raw black lines and fundamental twin-pair edges is
retained throughout. The order 6 black-bridge calculation finds three bridges,
all leaf-incident; it does not classify bridges at larger orders.

[Wiring and maximal-run transfer](session25-wiring/RESULT.md) ·
[Independent review](session25-wiring-review/REVIEW.md).

## 3. All-order labelled odd girth

Let L_n have all labelled tournaments as vertices. An edge preserves a
Hamiltonian path and reverses every other arc. For n>=3:

- If n=0 or 3 modulo4, L_n is bipartite.
- If n=1 or 2 modulo4 and n>=5, its odd girth is

      2 floor(n/4)+1.

For an odd closed k-walk, XOR of the k underlying path-edge sets must be
all pairs. Thus k(n−1)>=binom(n,2), giving k>=ceil(n/2). The parity of the
number of flipped arcs supplies the bipartite cases.

At orders2 modulo4, an explicit Hamilton-path decomposition attains the
bound. At orders4r+1, a Walecki cycle decomposition followed by deletion and
matching repairs gives2r+1 paths: every pair appears once, except a matching
of r pairs appearing three times. A paired-triples constraint graph is a
functional pseudoforest. Reordering one spare path in each cyclic component
removes its orientation obstruction; path reversals then give the actual
directed cycle. This proves all orders, not just the n5/n9 examples.

Independent checks include every n=5,9,...,241 and orders301,401. Those
checks supplement the construction proof. They do not establish it by
extrapolation.

A quotient walk can close only after relabelling or converse. Its full
holonomy correction is explicitly recorded; labelled cycles can collapse
to short quotient walks or loops. **The formula is not the odd girth of the
merged metagraph.** The newer all-order proof supersedes the still-open
odd-residue question recorded in the earlier wiring packet's section5.

[All-order odd-residue construction](session25-odd-girth/RESULT.md) ·
[Independent construction review](session25-odd-girth-review/REVIEW.md).

### The construction contains an exact modular doubling clock

Index the repaired path pairs by a modulo M=r+1. Their parent relation is
exactly a -> 2a modulo M, with 0 the unpaired base (its formal fixed loop
is omitted from the constraint skeleton). Writing M=2^s u with u
odd, the number of nonzero cycles is

    C(M)=sum_(d|u, d>1) phi(d)/ord_d(2).

For each compatible path order there are exactly 2^(C(M)+1) reversal
assignments. Every path order is compatible precisely when M is a power
of two. For M=21 the number is 64; for M=341=11*31, the divisor cycles
contribute 1+6+30=37, giving 2^38 assignments. These count reversal choices
in the named construction, not Hamiltonian paths or distinct tournaments.
The clock is genuinely conjugate to modular doubling; no conjugacy to
classical Collatz or payment of source-credit debt is asserted.

[Doubling-clock proof and typed examples](session25-doubling/DOUBLING.md).

## 4. A genuine S2 covering of the merged metagraph

For n>=3, let eta(T) record orientation parity relative to fixed labels, and let
Gamma=S_n×C2 act by relabelling and optional converse. The character is

    chi(pi,epsilon)=sgn(pi)(−1)^(epsilon binom(n,2)).

The quotient by ker(chi) is an honest two-sheet covering of the full
fundamental merged coloured multigraph. There are no edge inversions:
an ordinary edge swap would fix its common directed path, and a signed swap
would incorrectly preserve the nonpath arc joining its endpoints.

For n=0,3 modulo4 this is the canonical bipartite double cover. For n=1,2
modulo4 it is two disjoint copies. A base loop in the first case lifts to
TWO parallel edges; discarding one breaks the covering.

This gives an exact bridge-transfer law. In the bipartite-cover case a lift
is a bridge iff its base edge is a bridge and at least one of its two sides
is bipartite. If both sides contain odd cycles, the two lifted edges form
a two-edge cut but neither alone is a bridge. The sides refer to the edge's
original connected component. In the two-copy case bridges are preserved.

This is a classification of bridge behaviour under the covering, conditional
on actual base wiring. It is not a classification of which tournament
metagraphs contain non-leaf bridges.

[Covering proof, loop controls and bridge law](session25-double-cover/PROOF.md).

## 5. Prime-power closure and an unrestricted CRT criterion

Here s(n) counts strong **labelled tournament objects**, while
c(n)=2^((n−1)/2) for odd n and 0 for even n. Write D(n)=s(n)−c(n).
This D is not the older metagraph generating polynomial D(z), and neither
s(n) nor D(n) counts paths inside a fixed tournament.

For every odd prime p, r>=1, and n>=2r−1, put epsilon=(2/p). Then

    s(n+phi(p^r)) = epsilon s(n)  (mod p^r),

and the same holds for D. This closes exactly HYP-9358(d′). The stated
threshold is sufficient; individual minimal thresholds are not classified.

The proof preserves an exact finite group-algebra carrier. A Teichmüller
Fourier expansion places its exponents on one additive line. Integer-valued
Newton expansion puts every p^j correction in divided-power degree at most
2j. That filtered algebra survives inversion. Beyond2r−1, nonunit
exponents vanish at the required precision; unit exponents have the claimed
quasi-period. There is no illicit factorial division or assumption that the
representation into the Hurwitz ring is injective.

A separate orbit-contraction proof establishes the compatible dilation law

    s(mp^r) = epsilon^m s(mp^(r−1))  (mod p^r),

also for D. It contracts internally strong p-blocks to a weighted
semicomplete quotient and controls the p-valuations of digon coefficients.
Normalized sequences have p-adic limits; their values are not determined.

For arbitrary composite n, and each odd p^r exactly dividing n, set

    t = 2r−1 + ((n−(2r−1)) mod phi(p^r)).

Then n divides D(n) iff every such p^r divides D(t). The2-primary factor
is automatic when n is even. This is the full criterion, with no
squarefree or cube-free hypothesis. At prime squares it recovers
p² divides D(p²) iff p² divides D(p).

The same filtered method works for base a coprime to p, with character
(a/p), and for valid Hurwitz sums, products, shifts and unit-constant
inverses within that fixed-character algebra. It does not mean ordinary
pointwise products. For connected graph counts, the logarithmic derivative
gives the sufficient threshold 2r; no claim of optimality is made.

Incoming3596c5d0/MISTAKE-651 already corrected the source's connected-graph
and census statements and proved the square additive law. Those corrections
are preserved and credited. The new all-power proof extends them; stale
proposed corrections are not applied. The aggregate examples21,341,1023
remain distinct from the forbidden individual-H values 7 and 21.

[Additive theorem and full CRT proof](session25-congruence-review/ADDITIVE_PRIME_POWER_PROOF.md) ·
[Independent additive audit](session25-additive-review/REVIEW.md) ·
[Independent contraction proof](session25-congruence-review/PRIME_POWER_PROOF.md).

## 6. The attached paper and the actual S2×S3 construction

The attachment is Dos Santos et al., *At the Limits of Computability:
Steps Towards Autonomous AI Discovery on the Busy Beaver Frontier*,
arXiv:2610.22603v2. It reports certified accelerations for eleven BB6
machines, not their resolved halting or a classical Collatz proof.

Appendix A.9 merges two configuration families by encoding the binary family
label in an integer's parity. It does not discard the label. The naive
phase-forgetting map fails on an explicit pair of states. This session
also derives a second exact boundary-return map, including its number of
A.9 macro rules. That is not raw Turing time or native Collatz time.

The prior six-vector construction has B³=−I and JBJ=B^−1, yielding
S2×S3 on its specified orbit. A fixed J does not preserve an arbitrary
B-orbit; it does so only on six explicit lines. For a general nonzero vector v,
the transported reflection J_v restores the orbitwise action. Reflections
reverse the directed three-cycle, rather than being ordinary tournament
automorphisms. These qualifications survived explicit hostile tests.

[Paper and repaired group-action account](session25-paper/PAPER_CONNECTIONS.md) ·
[Macro-clock acceleration](session25-crossfield/CERTIFICATES.md).

## 7. Palette forcing and actual clocks

For a graph assembled from internally disjoint two-port fragments, each
admitting a separating k-colouring, retain a skeleton edge exactly when
its fragment forbids equal port colours. The assembly's complete boundary
support is exactly proper k-colouring of this hard-inequality skeleton.
Soft equality biases create no new existence constraint under this interface.

Consequently a nonvacuous six-colour equality assembly needs a six-chromatic
hard skeleton, hence treewidth at least 5. The threshold is sharp for
abstract graphs, not asserted for Euclidean unit-distance realizations.

For r named ports, the exact counting signature is indexed by partitions
into at most k colour classes. It is both sufficient and necessary for
all abstract boundary gluings. At three ports there are five states when
k>=3, not six permutations. S3 permutes the three states identifying which
pair is equal; forgetting that distinction fails even for two single-edge
fragments. This is a proved special case of standard connection-matrix
ideas, with Lovász's general theorem cited, not a claimed new general theory.

On the clock side, n_j=(2^(2+12j)−1)/3 has identical mod15/mod35 guards,
identical accelerated odd-Collatz successor1, and identical valuation phase
mod6, but unbounded exact unaccelerated arrival time for j>=1. Thus those
registers cannot replace an actual clock. The native 153N ledger remains
R′<=M−J′ as its outstanding target: no credit was silently reassigned,
recycled or claimed paid by this result.

[Hard-skeleton and clock proofs](session25-crossfield/CERTIFICATES.md) ·
[Named port signatures](session25-ports/PORT_SIGNATURES.md) ·
[Independent port audit](session25-ports-review/AUDIT.md).

## Remaining frontier

- Exact attainment of every odd H other than 7 and 21 remains open.
- Full merged-metagraph shape, general non-leaf bridge occurrence, and
  all-order black-core connectivity are not settled by these certificates.
- The new labelled odd-girth theorem and one-bit cover do not retain full
  permutation holonomy or determine every adjacency.
- HYP-9358's other density, distribution and exact-threshold questions stay
  open; the additive subproblem is the one closed here.
- General Collatz, uniform 153N coverage, and a new Hadwiger–Nelson bound
  remain open. Abstract colouring constraints require real geometric
  realization before they imply a Euclidean result.

The next proof-oriented targets are constructive strong prime lanes,
port-preserving cycle repairs for actual merged bridges, and genuinely
coupled multi-port geometric gadgets. Congruence completeness, scalar
counts, soft two-port bias and compressed clock phases cannot substitute
for those missing mechanisms.

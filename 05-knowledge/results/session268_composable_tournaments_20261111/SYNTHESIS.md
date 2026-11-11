# Classify what composes, then let the hidden states interact

**Session268, 11 November2026. Chan–Pak remains OPEN.**

This session produces two exact infinite-family enumeration methods, a sharp
physical observation theorem, and a literal interacting cellular-automaton
family. The main gain is not a larger unstructured census. It is identifying
which state can safely be merged, how many distinct graphs each state hides,
and which operations can make that hidden information visible again.

H counts directed Hamiltonian paths. Ordinary tournament isomorphism classes,
classes after identifying converse, labelled tournaments, and distinct integer
H values are different objects. All class counts below are ordinary classes
unless explicitly marked labelled. No complete global class census or full
odd-spectrum theorem is asserted.

## 1. More than a million classes without a million graph records

Let q_r(T) count ordered covers of T by r nonempty directed paths. Retain the
entire vector q, with q_1=H. Uniform tournament substitution consumes precisely
this vector, via an exact quotient-walk/separator generating function.

The new inventory compiler groups already classified children by q, but retains
two weights: their actual number of isomorphism classes and the sum of their
reciprocal automorphism orders. A response-weighted Burnside formula then
counts parent classes exactly. Multiplying the second weight by n! gives the
labelled count. A response bucket is never mistaken for a single graph type.

The executable infinite family allows every canonical prime quotient to be C3,
and allows arbitrary ordered strong-component sums. At order16 it contains:

- **1,343,039 ordinary classes**;
- **83,073 full response states**, a16.17-fold reduction;
- **74,149 distinct H values**;
- **10,997,137,501,622,272,000 labelled tournaments**.

The largest response bucket contains252 classes. The whole inventory through16
ran in9.728 seconds in the recorded environment before JSON serialization,
using69,883 distinct cyclic response contractions and103,793 ordinal-sum
contractions. This is a measured scoped workload, not a speed comparison with
nauty or a claim that all tournament enumeration is polynomial-time.

The essential diagonal term is instructive: if m actual child classes have the
same response, three cyclic ports contribute (m^3+2m)/3 ordinary parent classes,
not one response orbit and not m^3/3. The observable symmetry is larger than
the actual graph automorphism group. Their distinction creates the efficiency
gain without sacrificing multiplicity.

The canonical decomposition, old prime cycle-index census, and original
path-cover substitution machinery are inherited inputs. The new step is their
response-refined class and labelled-mass pushforward. Its implementation and
all H histograms through16 are in [the modular packet](modular/README.md).

## 2. The exact state has a sharp physical observation cost

Equal H alone is not safe for composition. For example
C3[TT3,1,1] followed by one sink, and C3 followed by C3, both have H=9.
Their ordered cover vectors begin (9,150,...) and (9,126,...). Applying the
same cyclic context C3[child,1,1] gives177 and153.

The new tomography theorem supplies an explicit converse to the composition
law. The scalar experiment

    P_d(T)=H(C3[TT_d,TT_d,T])

is triangular in q_1,...,q_(d+1), with diagonal coefficient(d!)^2. Thus H(T)
and these physical probes recover the entire response exactly. In particular:

**Two nonempty tournaments have the same full path-cover response if and only
if every uniform substitution context gives them the same H.**

There is also an optimality statement. If the first difference is q_r, a
separating context needs at least2r−2 auxiliary vertices, and the displayed
probe attains that cost. Between consecutive runs through a module, a
Hamiltonian path must use at least two external vertices: one outside vertex
has a uniform orientation to the module and cannot both receive and return.
That is the mechanism behind the lower bound.

For known order n, the top two coordinates q_n=n! and
q_(n−1)=(n−1)! binom(n,2) are universal. Consequently, for n≥4, H and the
probes d=1,...,n−3 suffice; the largest has3n−6 vertices. This is an exact
information and construction theorem, not a polynomial-time algorithm for
counting paths in an arbitrary raw graph. Nor does it say q is a minimal
nonlinear encoding of graph isomorphism.

[Coefficients, inverse, and sharp observation filtration](spectrum/PROOF.md).

## 3. Previously invisible computation becomes genuinely interacting

The repository's Session168 contained two five-vertex bits differing in one
arc but sharing their complete unrooted path-cover profile. The theorem above
shows that **no pure uniform substitution context can reveal their bit by H**.
More elaborate substitution alone cannot fix this particular blindness.

One deliberately nonuniform seam does. Start with A followed by B and reverse
a single cross arc a→b. The exact new law is

    H(A⋈B)=H(A)H(B)−end_A(a) start_B(b)+L_A(a) R_B(b),

where L and R retain rooted two-path covers, allowing one empty segment.
These rooted data were erased by the old unrooted profile. Applied to the two
bits, the four resulting ten-vertex cells are strong and rigid, with counts

    w00=215, w01=179, w10=245, w11=201.

The mixed additive interaction is−8. The multiplicative determinant is−640.
So the response is neither an independent sum nor an independent product of
two unary bit weights. This is a physical interaction through actual arcs.

For a binary word x on an indexed N-cycle, take the ordered strong-component
sum of the N cells G_(x_i,x_(i+1)). The SCC order is intrinsic, so the word is
recoverable from the unmarked tournament. Rotations are not silently identified.
This gives exactly **2^N ordinary classes on10N vertices**, all rigid.
The exact labelled count is therefore **(10N)!2^N**.

Nevertheless, only **floor(N^2/4)+2 H values** occur. For a nonconstant word
with n0 zeroes and k zero-runs, its count and multiplicity are

    h=215^(n0−k) 201^(N−n0−k) (179·245)^k,
    multiplicity=(N/k) binom(n0−1,k−1) binom(N−n0−1,k−1).

Prime valuations recover n0 and k, so these bins really are distinct.
At N=100,2^100 classes compress to exactly2502 bins; the committed program
computes every multiplicity and checks their sum without enumerating2^100
words. At N=1000 the theorem gives250002 bins for2^1000 classes on10000
vertices. The full explicit list of classes still has exponential output size.

[Interaction law, exact class spectrum, rigidity and labelled counts](automata/PROOF.md).

## 4. A literal CA action, with a visible limit to scalar dynamics

Any binary cellular rule acts on the decoded words and hence on these ordinary
tournament classes. In the four-cell alphabet the overlap condition is local.
Changing one bit changes exactly two actual arcs in the pinned graph. A
radius-r rule has its usual finite light cone in ordered-component geometry;
the H update for a single defect depends on just its two adjacent pair factors.

H is a global readout, not a complete dynamical state. Under Rule110,
101000 and100100 have the same present H88902740880625, but their next counts
are76567725723825 and77701668671025. Longer initial block statistics predict
a longer future horizon: length2rt+2 blocks suffice at time t. There is no
claim that a fixed pair histogram closes every rule for all time.

This implements interacting computation, rather than attaching a computer
analogy to a static tournament. Finite rings remain finite systems. No
undecidability or computational-universality theorem is inferred from the
sample Rule110 trace, and no automatic local rule is extracted from an
undirected arc-flip graph.

There is also a literal statistical-mechanics transfer. Summing H over all
classes of spatial size N gives Z_N=trace(W^N), W=(w_ab), hence

    Z_N=416 Z_(N−1)+640 Z_(N−2),
    Z_N=(208+56 sqrt(14))^N+(208−56 sqrt(14))^N.

This recurrence is in spatial size **after summing all classes**, not in the
time evolution of one CA state. [The joint bridge and loss ledger](BRIDGES.md).

## 5. An independent recursive lens: exact parity dynamics

For a regular odd child B, its F2 bit space splits into the global constant
bit and the zero-sum internal modes. The adjacency automaton of Q[B] splits
exactly into Q shifted by the child's degree parity, plus copies of those
internal modes. Iterating this rule classifies the automaton of every word in
C3 and Paley7 using O(word-length) integer operations.

All those automata have only fixed points and eventual three-cycles, with
explicit multiplicities. Yet composition order is visible: C3[Paley7] and
Paley7[C3] both have21 vertices but binary ranks11 and17, and respective
Hamiltonian counts122621733302733 and123156134945235.

The same physical substitution used in Session267's wreath-action example
therefore supports another exact recursive invariant. A simple parity clock
does not provide the integer Hamiltonian response, however. These are different
quotients of the common module-word carrier.

[All-word block theorem and direct matrix controls](dynamics/PROOF.md).

## 6. What this achieves toward Chan–Pak, and what is still missing

The spectrum lane gives an actual strong construction from every seed T:

    H(C3[1,TT_b,T])=(H(T)+q_2(T))2^b+H(T)−q_2(T).

For a fixed seed, testing a target requires one divisibility and power-of-two
test, not a tournament search. Coupling the old invisible bits through a
transitive probe and seam gives additional explicit two-parameter affine
families. These are genuine constructions, but their integer ranges are
specified rather than assumed to be complete.

The basic singleton-transitive atom h(b,c)=2+(2^b−1)(2^c−1) has a unique
unordered parameter decoder, exact chirality and automorphism multiplicities,
and logarithmic vertex cost. Even closing **all** these atoms under arbitrary
products covers only exp(O((log X)^(2/3))) integers below X. A partition
generating-function bound proves this sparsity. The complete finite product
inventories have95,368,1308,4264 values through10^3,10^4,10^5,10^6.

This explains a real obstruction: huge class families and clean multiplicative
composition need not approach the desired scalar spectrum. The missing step
is target-directed, height-controlled interaction on a physical rooted state
domain. The exact seam law now exposes a useful subtraction term, and the
tomography theorem tells us how much hidden state must be retained, but neither
proves that every desired odd count except7 and21 is reached.

The next focused work is to classify rooted-state image constraints for one
seam grammar and find a decreasing target algorithm or finite covering
certificate. Globally strong interacting chains are another explicit gap:
the present cells are strong, but their ordered sum is not strong for N>1.

All new proofs, programs and controls are included. Historical results were
accepted in their declared scope; no inherited Lean rebuild, independent
global audit or external novelty claim is made. The session advances exact
enumeration, compositional state theory and constructive families while
leaving the global Chan–Pak conjecture honestly open.

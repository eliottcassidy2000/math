# Session 26: cuts, memory, and the forced wiring of low-intransitivity tournaments

22 October 2026. This packet records proved structural results, exact
finite controls, explicit counterexamples, and scoped tools. It does not
claim to solve the full Hamiltonian-count spectrum, unrestricted tournament
enumeration, LRC(14), Collatz, or Hadwiger–Nelson.

The exploration sampled **103 distinct archived repository paths** in four
reproducible passes. It repeatedly moved from an analogy to a precise map,
tested the missing coordinate, and then proved or repaired the resulting
statement. See the [exploration and concept board](EXPLORATION.md),
[root sampling manifest](ROOT_SAMPLE.json), and
[sampling totals](SAMPLING_SUMMARY.json).

## 1. A global constraint on every actual metagraph edge

Write tau(T) for the minimum number of vertices whose deletion makes T
transitive. If T' is obtained by retaining a directed Hamiltonian path
and reversing every other arc, then

    tau(T)+tau(T') >= floor(n/3),
    tau(T') <= floor((n+tau(T))/2).

The first bound is sharp at the transitive node. Partition the path into
disjoint consecutive triples: each triple is cyclic in exactly one of
the two tournaments. The upper bound retains a feedback deletion set,
then covers the remaining path-edge forest to leave a reversed transitive
tournament. No sharpness claim is made for that upper bound.

These are isomorphism- and converse-invariant statements, so they apply
to BOTH colors in the merged metagraph, not merely to a labelled
auxiliary graph. Consequently the stratum tau<=k is independent for
n>=6k+3. Every neighbor of that stratum lies between feedback levels
floor(n/3)-k and floor((n+k)/2).

Sources: [lower bound and exact one-apex incidences](session26-blue/RESULT.md),
[general fixed-feedback analysis and upper bound](session26-fixed-apex/proof.md).

## 2. Complete induced wiring of the one-apex stratum

Here “one-apex” includes every tournament with tau<=1, not just a
particular chosen apex presentation.

At n=8, the complete nontrivial induced graph is

    0001010 --- 1000010 --- 0010110
    1000110 --- 0010010 --- 1001010
    0011010 === 1010010

All six edges are black. The last component has TWO parallel edges.
There are 47 additional isolated nodes, for 55 selected merged nodes
in total. No loop occurs. The two H=15 endpoints in the two path
components are different nodes; H alone is not their identifier.

The configuration is forced by a short structural argument. If both
endpoints have a deletion apex, those apices must occupy path positions
3 and 6. The other six vertices split into the agreement blocks
(12),(45),(78). Only two block orders survive the actual entry/exit
constraints, leaving twelve tilings and their twelve companions.
Their six reflection/complement orbits give the displayed graph.

The complete n=3,...,8 node-and-edge catalogue is supplied. For EVERY
n>=9, the induced one-apex graph is edgeless by the general bound above.
Thus this induced stratum, including its final exceptional wiring,
is classified at every order. This is not a classification of edges
from that stratum into the rest of the metagraph, or of bridges in the
full graph. The displayed path edges are bridges only in the induced graph.

[Structural proof and full catalogue](session26-last-wiring/RESULT.md).
[Independent exhaustive audit through all 2^21 pinned-path masks at n=8](session26-last-wiring-review/AUDIT.md).

## 3. Exact counts, tiling mass, and colored boundary budgets

Let A_n be the number of ordinary one-apex classes, C_n the self-converse
number, and V_n the sum of H(T)/|Aut(T)| over ordinary classes. Then, n>=2,

    A_n = 2^(n-1) - binom(n-1,2) - 1,
    V_n = 4^(n-2) - 3*2^(n-2) - (n^2-7n+4)/2.

For r>=1,

    C_(2r)=r,                   C_(2r+1)=2^r+r-1,
    B_(2r)=2^r+r-2,             B_(2r+1)=2^(2r-1)+r-1,

where B_n is the number of blue incidences contributed by the stratum.
Selected merged nodes number (A_n+C_n)/2.

For n>=9 every selected incidence is a boundary edge, so the exact
fundamental boundary counts are B_n blue and (V_n-B_n)/2 black.
At n=8 the six internal black edges consume twelve black incidences.

| n | selected merged nodes | selected tilings V_n | blue boundary | black boundary |
|---|---:|---:|---:|---:|
| 8 | 55 | 3,898 | 18 | 1,928 |
| 9 | 123 | 15,989 | 131 | 7,929 |
| 10 | 240 | 64,751 | 35 | 32,358 |

Within this stratum, points of symmetry are exactly the transitive class
and, at odd n>=3, the centrally padded C3 class. Other strata may contain
additional points of symmetry.

The proof restores root-forgetting multiplicities rather than treating
apex words as unlabelled tournaments. A strong one-apex tournament of
order at least four has either one deletion apex, or is the exceptional
double-apex blow-up C3(1,1,n-2). These tournaments are rigid; C3 is the
separate automorphism-order-three exception.

Both the class-count and tiling-mass ordinary generating functions are
rational and are given explicitly in the [counting and wiring note](session26-blue/RESULT.md).
[Source-bound review](session26-blue-review/AUDIT.md).

## 4. Hamiltonian counts as cuts and Young-diagram weights

For a chosen apex over an ordered transitive base, let b_i indicate the
arc apex->i. Then

    H(b)-1 = sum_(j<i) b_j(1-b_i) w_(i-j),
    w_1=2,       w_d=2^(d-1) for d>=2.

Every odd cycle contains the apex, and its interior is an increasing
even-sized subset. That independently explains the cut weights.

The inversion cells of the binary word form a Ferrers partition.
The hook length of each cell is exactly its endpoint distance, giving

    H = 1 + sum_cells 2^(hook_length-1) + number_of_corners.

This is not a tableau hook-product identity. It gives sharp extrema:
fixed-content strong minima are hook diagrams and maxima are rectangles.
For strong one-apex tournaments of order n>=5,

    H_min = (2^floor((n-3)/2)+1)(2^ceil((n-3)/2)+1).

The additional small-order equality cases are classified in the proof.
For n>=3,

    H_max = 2^(n-1)-2^floor((n-1)/2)-2^ceil((n-1)/2)+3.

In particular, H=13 is impossible throughout the entire one-apex family,
including transitive padding. A concrete five-vertex tournament with
feedback number two and H=13 is supplied. Thus this restriction is real,
not a global tournament gap.

[Cut, extrema and Ferrers proof](session26-tournament/CUT_ADDENDUM.md).
[H=13 witness](session26-tournament/H13_WITNESS.md).
The exploration note also derives an exact electrical cut-energy
representation, with an indispensable linear field; energy alone fails.

## 5. A spectral counting law and a structural gap beyond congruences

For a uniformly chosen strong apex word of length m, the normalized
count H/2^m converges to U*V, where U,V are independent uniform variables
on [1/2,1]. The convergence has a uniform CDF error O(2^(-m/2)),
proved by disjoint leading-bit coupling, not fitted from data.

The product density is

    4 log(4z) on 1/4<z<1/2,
    -4 log(z) on 1/2<z<1.

Its mean is 9/16 and variance 55/2304. Summing across lengths cancels
the dyadic phase and proves

    W(X) = (log 2) X + O(sqrt X),

for the number of strong apex WORDS with H<=X. Ordinary strong one-apex
classes have the same main term; converse-merged strong classes have
half that main term.

Reverse-complement pairing also proves, without needing the class
enumeration, that the number of DISTINCT attained one-apex H values
up to X is at most

    (log 2)X/2 + O(sqrt X).

Hence the missing positive odd values have lower relative density at
least 1-log(2), about 30.7%, within the odd integers. No existence of
a natural density for the attained set is asserted.

This contrasts with the inherited Session25 theorem that the same
one-apex family realizes every parity-compatible residue modulo every
modulus. Congruence completeness is therefore not actual spectrum
completeness.

All class-counting-by-H statements and the associated Dirichlet
continuation are restricted to STRONG CORES. Without that restriction,
padding gives infinitely many classes at fixed H.
[Full proof, exact conventions and Dirichlet consequence](session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md).
[Independent review](session26-cut-review/REVIEW.md).

## 6. Every fixed feedback budget has an explicit asymptotic structure

Fix k>=1, write n=m+k, and let m tend to infinity. Put

    L=2^(binom(k,2)+km).

Ordinary classes with tau<=k satisfy

    N_(n,k) = L/k! * (1+O_k(m^2*2^(-m))).

Their strong subfamily has leading fraction (1-2^(-k))^2.
Almost all presentations have a unique deleting k-set and are rigid;
self-converse classes are negligible, so converse merging halves the
leading class count.

The selected tiling mass satisfies

    V_(n,k) ~ 2^(binom(k,2)-2k) * [2^k(k+1)]^m.

Blue incidence mass is o(V_(n,k)). Once n>=6k+3, every selected edge
leaves the stratum, and black boundary edges number asymptotically
V_(n,k)/2.

These formulas include their error terms and full weighted
multiplicity proof in the [fixed-feedback theorem](session26-fixed-apex/proof.md).
[Independent audit](session26-fixed-apex-review/AUDIT.md).
They are NOT unrestricted A000568 asymptotics: k is fixed, not allowed
to grow with n inside these estimates.

## 7. Exact attachment transforms and an executable path counter

For a fixed k-vertex core Q, let p_r count its ordered spanning covers by
r nonempty directed paths, and P_Q(t)=sum p_r t^r. Attach m transitive
base vertices with all cross arcs independently fair. If E_Q(m) is the
expected H, then, with y=e^x-1,

    sum_m E_Q(m)x^m/m! = ((y+2)^2/y) P_Q(y/4).

The apparent singularity is removable. The first ceil(k/2) expectations
recover the whole ordered path-cover profile, but that profile does not
classify tournaments. The first two exponential channels are universal;
the first possible triangle sensitivity is quantified exactly.

The channel polynomial obeys a centered reflection law and has a
double, or in even core order at least triple, zero at -1. Positive
expectations do NOT force positive channel weights: the transitive
four-core has a negative positive-frequency coefficient and an exact
negative shifted Hankel determinant.

[Transform, tomography, reflection and hostile controls](session26-attachments/proof.md).
[Independent review](session26-attachments-review/review.md).

Without averaging, the same slot decomposition gives an exact
nonnegative weighted automaton with at most 2*k!*3^k states for a
supplied core and attachment word. It retains first/last ports and
checks core arcs across empty gaps. Its arithmetic-operation count
is linear in the base-word length for fixed k; preprocessing and
integer bit complexity are stated separately.

Run, from its directory:

    python3 counter.py witness13.json

The answer is 13.
[Counter and proof](session26-counter/proof.md).
[Independent audit](session26-counter-review/REVIEW.md).
This is not a minimal-state claim, an algorithm for finding the
deletion set, or an unrestricted tournament canonizer.

## 8. Three orthogonal results from the wider archive

### Minimal continuation memory for the one-apex scalar

The exact minimal deterministic residue-output observer has
2M^2 states for odd M>=3, M^2/2 for even M>=4, and one state for M=1,2.
Its rational Hankel rank is only four; in characteristic two it is one.
The integral Smith fingerprint is (1,2,2,4). CRT products must share
one literal last-bit control, rather than multiply independent copies.

[Proof and independent controls](session26-tournament/PROOF.md).

### Quadratic digit ideals and the real role of the prime five

For a finite reachable abelian group G under the binary affine update
s <- Ms + dv, with M invertible and homomorphic output, let K be its
unobservable subgroup. Validating at most r>=1 consecutive ones while
returning every output residue needs exactly (r+1)|G/K|+1 states. If legality is promised,
only |G/K| arithmetic states are needed.

For Fibonacci and Lucas inputs under the no-adjacent-one grammar, the
validated counts are respectively 2M^2+1 and 2M^2/gcd(M,5)+1.
The distinction is a genuine quadratic-ideal calculation:
the Fibonacci seed is the unit phi^3, while the Lucas seed is
phi*sqrt(5). Reachability loss and unobservability are kept separate.

This provides exact modulo15/modulo35/joint105 readers. It does not
simulate Collatz transitions or preserve an integer stopping clock.

[Full theorem and scope](session26-wildcard/RESULT.md).

### Runner clearance as a signed height measure

The clearance curve f(t)=min_i ||v_i t|| has an exact signed height-atom
representation. Its NORMALIZED moments, not its raw moments, have a
finite rational ordinary generating function. The nearest positive pole
recovers the maximum clearance.

For two speeds (1,N) with N even, the minimal normalized-moment recurrence
order is exactly N-1. Thus runner count alone does not bound this linear
recurrence order. This does not rule out nonlinear reconstruction or
a closed formula within a special family. The (1,4) example refutes
positivity of all curvature atoms.

[Proof, exact controls and limitations](session26-runners/RESULT.md).

## 9. Repairs, current boundaries, and the next obligations

The [archival repair packet](session26-corrections/README.md) distinguishes
a truncated counterterm from its reciprocal-odd parent, preserves a
Coxeter numerical analogy without calling it a causal theorem, and
repairs two stale LRC reflections. The exact large-denominator family
in that packet has maximum clearance 1/8 despite evading every
prescribed finite rational census.

Incoming repository corrections and new source work were preserved,
including the S50 prior-art/finite-census repairs and S51 repunit
corrections. They were not silently promoted into independent replays
by this session.

The most useful next steps are now concrete:

1. Use the two-apex counter to seek actual spectrum realization mechanisms,
   not another finite congruence sieve. H=13 already requires leaving
   the one-apex stratum.
2. Classify the cross-stratum destinations and full-graph bridges left
   by the feedback inequalities; induced-stratum bridges are not enough.
3. Seek sharper fixed-k correction terms and smaller exact observers,
   with a benchmark before claiming practical enumeration acceleration.
4. Restore phase/owner/clock information when using height distributions
   or modular digit readers in LRC and Collatz.
5. Pay the native 153N debt with actual distinct sources and clocks.
   The new cut count 153, or a shared state-space size, is not that payment.

The full {7,21} spectrum question, earlier global metagraph D(z)/M_k
questions, all full-metagraph wiring, LRC(14), general Collatz and
Hadwiger–Nelson remain unresolved by this session.


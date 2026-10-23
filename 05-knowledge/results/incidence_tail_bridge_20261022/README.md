# Incidence retained: cycle capacities, delayed clocks, and angular data

October22,2026. **PROVED scoped:** the linked optional-arrival cycle-capacity
bound below. **FINITE-EXACT:** its controls and frozen native-prefix consumers.
**CITED / ANALOGY:** the lattice and Dehn connections are explicitly separated
from the Collatz statement. Uniform153N remains **OPEN**.

## 1. Inheritance and a concrete transfer from the pasted paper

The [SRFP311T1 incidence paper](https://srfp311t1.com/incidence_moments_metric_decks)
uses shared-tuple incidence to recover information lost by marginal counts.
Its [structured-incidence companion](https://srfp311t1.com/leech_niemeier_incidence)
retains a marked incidence operator before taking its Gram matrix. Here the
transfer is methodological: retain the interaction graph and its shared edges,
then ask which scalar maximizers are jointly possible. No lattice reconstruction
theorem is assumed to apply to tournaments or Collatz.

Closest mechanism: the signed-cycle packing of
[the earlier terminal-cube packet](../collatz_cube_context_20261015d.md).
Canonical hostile: three negative edges cannot all be satisfied by vertex
signs. Corrected near miss: a triangle-only bank can miss a longer negative
cycle. Least-used sidecar: the exact capacity of every shared interaction edge.
Board: immutable source / pair interaction / cycle sign / edge capacity /
terminal incidence / angular operator.

## 2. From an indefinite linked kernel to a sound cycle saving

Use the original linked state, delta, kappa and envelope U from
[future tails](../future_tail_blocks_20261021b/README.md). Form a signed
undirected graph on the literal pending nodes with edge coefficient
k_ij=kappa(b_i,b_j), omitting zero edges. It is not a tournament: its observable
is symmetric and ties really are absent.

For each simple negative cycle C choose weight p_C>=0, with

    sum_(C containing e) p_C <= |k_e|  for every edge e.       (1)

Then at every later native cut

    Delta_h <= U - 2 sum_C p_C.                              (2)

**Proof including arbitrary destinations and nonarrival.** Fix any legal
optional signed terminal assignment. Randomly choose an independent sign
for each occupied terminal group and each omitted vertex. Multiply every
included vertex's actual sign by its group's random sign; give an omitted
vertex its own random sign. The expectation of the full signed edge sum is
exactly the actual pending-pair contribution: cross-terminal and omitted
terms average to zero. For any full sign assignment a negative cycle has an
unsatisfied edge. By(1), the weight of unsatisfied edges is at least sum p_C.
Every unsatisfied edge loses2|k_e| relative to its absolute maximum. Take
expectations and retain the old valid local bounds for the delta/current-port
terms. This proves(2), without positivity of the four-dimensional linked form.

The earlier terminal-Gram argument used a positive quadratic cube and factor4
because its matrix counted off-diagonal terms twice. Here kappa already is
the unordered-pair coefficient, hence factor2. This is a classical signed-graph
capacity method specialized to the linked optional-arrival interface; no
novelty claim is made for signed-cycle frustration itself. See also
[Aref et al., frustration-index models](https://arxiv.org/abs/1611.09030).

If exact cluster bounds already consume internal edges, apply(1) only to
remaining cross-cluster edges. This combines the methods without recycling
an internal edge. Do not subtract a cycle saving from an independently
optimized angular bound; both may already have paid the same incompatibility.

## 3. Why longer cycles and mixed channels matter

For the Q-Gram of the five planar vectors

    (-4,-5), (-3,-8), (3,-7), (-7,-1), (-5,-5),

triangle packing has exact optimum48, while adding four-cycles raises the
exact all-cycle optimum to67. The checker retains matching primal and dual
witnesses: three triangles of weights27,1,20 give48; adding four-cycles
(0,2,3,4) and(1,2,3,4) of weights6,13 gives67. A literal sign vector
(-1,-1,-1,1,-1) has unsatisfied-edge weight67. Thus no numerical optimizer is
needed for either claim. A cycle basis or triangle tally is not the capacity
ledger of all interacting cycles.

For linked four-vectors

    (-3,2,2,-3), (2,0,-2,-2), (3,0,-1,3)

and zero current ports, U=57 and a negative triangle saves36, giving21.
The exact optional three-port assignment maximum is18. In contrast every
scalar Lorentz projection of the prior angular family is at least
(sqrt(2701)-1)/2>21. One fixed signed sum has

    4 Q(g_lambda)=37lambda^2-2+73/lambda^2,

and AM-GM proves the lower bound for all lambda>0, not just a sampled grid.
So the cycle method adds a real mixed-channel capability. It is not an exact
optimizer for general linked states.

On the saved pure-new states, the existing exact angular bound remains at
least as good. Greedy triangles plus at most1000 longer negative cycles give:

| N | Old U | Capacity saving | Capacity bound | Existing angular bound |
|---:|---:|---:|---:|---:|
| 4096 | 557164 | 5320 | 551844 | 551844 |
| 8192 | 1138407 | 19526 | 1118881 | 1118881 |
| 32768 | 4241482 | 99200 | 4142282 | 4142214 |

The first two exactly explain the angular saving by owned cycle payments;
the last is a bounded discovery result, not a claimed packing optimum.
N7 is still hostile. All original source, guard and integer-clock commitments
are retained. No new trajectory suffix is used here.

## 4. Dehn invariance and the rank88 angular operator: precise analogy

For Euclidean polyhedra the classical Dehn invariant is

    D(P)=sum_edges length(e) tensor [dihedral_angle(e)]
      in R tensor_Q (R / pi Q).

It retains angle-length coupling that volume loses, and is additive under
dissection; see [MIT's Hilbert-third-problem notes](https://ocw.mit.edu/courses/res-18-011-algebra-i-student-notes-fall-2021/mit18_701f21_lec35.pdf).
It is not merely a scalar angle sum or an element of R/(pi Q).

The incoming [rank88 theorem](../session28_modular_spectral_20261022/session28-rank88/REPORT.md#7-a-fixed-quartic-sidecar-controls-all-gaussian-hessians)
instead gives Hess_t=a(t)I+b(t)A. Scalar theta determines a,b but loses the
fixed quartic shell operator A. Its [audited scalar correction](../session28_modular_spectral_20261022/session28-cm-audit/scalar-followon/ADDENDUM.md)
shows that the reciprocal-symmetric ratio a/b has a strict local maximum at1
and at least two reciprocal off-center global minimizers. It does not prove
there are exactly two, locate them, or compute the minimum angular eigenvalue.

The actual shared pattern is **a scalar observable plus retained coupling
data**, with different operations in each setting. Dehn data are invariant
under dissection. A is a fourth-moment operator on trace-free symmetric
deformations, not a Dehn invariant. Our kappa graph retains shared source
interactions under an optional destination relaxation, not a lattice shell.
There is no claimed map identifying these objects or their spectra. The
cheapest decisive transfer was to keep shared-edge capacities, yielding(2)
and the five-vector longer-cycle witness above.

## 5. A repaired tournament moment hierarchy

The pasted suggestion that bounded small-subtournament counts refine the
score sequence needs a scope condition. A source joined to a directed triangle
and a directed triangle joined to a sink both have one cyclic and three
transitive triples, but scores(1,1,1,3) and(0,2,2,2). Their three-decks agree
while their scores differ. Conversely the cyclic-triangle count is determined
by the full score multiset. Use joint marked invariants, not an assumed total
hierarchy between different quotients.

Recovery also exposed an invalid inference in
[THM-1970, now scope-corrected](../../../01-canon/theorems/THM-1970-H-sits-at-the-formula-sharp-P-edge.md):
counting hardness alone does not prove that a polynomial-time invariant fails
to determine the count; its decoder may be hard. Nor is a spectral trace the
full induced deck. The finite collision witnesses survive; the unsupported
universal reconstruction/complexity deductions do not enter these proofs.

## Reproduction and next test

    python -B 05-knowledge/results/incidence_tail_bridge_20261022/check.py
    python -B -O 05-knowledge/results/incidence_tail_bridge_20261022/check.py

The compressed bank stores literal cycle nodes and integer weights. The consumer
checks each negative sign and each total edge spend; duplicate credit is a
hostile control. The mixed-control enumeration covers every optional assignment
to three terminals. Native receipts pin all prior source coordinates.

The useful next combination is with
[authenticated delayed pulses](../lagged_tail_pulses_20261022/README.md): first
change the admissible future carrier by a proved native relation, then bound
that whole carrier. A better certificate cannot by itself establish that every
scale eventually admits a passing certificate.

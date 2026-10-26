# Incoming work changes the frontier

2026-10-25. Initial HEAD was 1b4e6dcc0f25. A mid-session fetch brought
5d87c3b009 and 9e79465c2c from the concurrent Session64 packet.

**Current update:** the next fetch brought 2a9a9c342e and f265adc2b8.
Session64 now proves the **sharp n>=8** closed endpoint-four threshold,
with a whole-illegal-leaf repair at n8 and the actual order-seven obstruction.
Its [boundary proof](../session64_lca_reverse_valuation_20261026/session64-repairs/boundary/README.md)
and [independent lower audit](../session64_lca_reverse_valuation_20261026/session64-audit-lower/AUDIT.md)
were read in full. The earlier n>=11 threshold below records the first
incoming checkpoint, not the current strongest result.

The [Newton quotient](../session64_lca_reverse_valuation_20261026/session64-counting/NEWTON_QUOTIENT.md)
now compresses the **signed clone-excess weight bank at absolute precision**
using precisions K,K-1,K-3,K-6,... . Its [LCA proof](../session64_lca_reverse_valuation_20261026/session64-root/newton-lca/REPORT.md)
identifies each triangular loss as a sum of named pairwise valuation depths.
These proofs were read, including the sparse-support and negative-division
hostiles. It leaves the already scalar ordinary state graph unchanged.
This is a DVR evaluation-lattice theorem; our non-UFD ideal example explains
why local valuation data cannot automatically diagonalize a global module.

The separate [Mahler/LCA arithmetic note](../lca_ports_reverse_20261025/arithmetic/README.md)
recovers the exact translation z=b+3 to ceil(3z/2), and the reset family
x_R=(2^(R+2)+7)/9, R=5 mod6. Its old depth is one, its next letter is two,
yet its new depth is R. The congruence extension preserves any fixed finite
modulus while R grows. This strengthens the obstruction against a finite
congruence-only rank; our density argument supplies the separate quantitative
cost of retaining a source through refinement. Neither proves failure on
the marked blank-tape orbit or on a Collatz trajectory.

**Third incoming checkpoint (07695d3952 / 9053ebfc58):** read
[THM-4772's conductor source](../fixed_twists_conductors_reverse_potential_20261025/conductors/DEDEKIND_CONDUCTORS.md)
and [THM-4773's lift source](../fixed_twists_conductors_reverse_potential_20261025/twists/FIXED_TWIST_LIFTS.md),
plus their canonical statements. The normalization formula provides a stronger
general description of the character ideals than our examples alone. It is
consistent with I1 after inserting its normalized coefficient ideals.

The audit found and repaired three precise issues in those incoming sources
and their canon statements: N2's ideal inclusion needed reversal and a
reflexive hull in higher dimension; shared valuations invalidate the claimed
reduction of every unreduced product test to single-root tests; whole-conic
nonemptiness implies open-central nonemptiness only over infinite fields.
`check_incoming.py` supplies the exact four-root and F3 counterexamples.
The full conductor bank, normalization formula and Q-conic application
survive. These are targeted corrections, not a full independent audit of
the incoming packet's Dedekind arithmetic or other claims.

The following statements were read in their proofs and scope notes, rather
than inferred from commit titles:

* [All six failures](../session64_lca_reverse_valuation_20261026/session64-repairs/README.md):
  complete affine blue slices have dimensions 3,3,0,0,0,0. A cyclic-shift
  bridge supplies three-black certificates for all twenty exceptional
  sources. Combined with the earlier residual theorem, the closed
  endpoint-four threshold becomes n>=12.
* [Order eleven](../session64_lca_reverse_valuation_20261026/session64-repairs/next11/README.md):
  seven exact exceptional sources close all 224 affine cells, including
  the split endpair cells. The threshold becomes n>=11. One certificate
  has a genuine intermediate ordered clone merger (5,4), then a split.
  Literal vertices and complete incidence are retained.
* [Reverse implementation](../session64_lca_reverse_valuation_20261026/session64-counting/REVERSE_IMPLEMENTATION.md):
  the actual ordinary DAG now computes normalized residues directly, with
  an exact structural valuation-budget filter. Its independent controls
  include intermediate coordinates, unmerged partitions and full small
  class enumeration.
* [Direct skip contraction](../session64_lca_reverse_valuation_20261026/session64-counting/SKIP_CONTRACTION.md)
  and its [audit](../session64_lca_reverse_valuation_20261026/session64-audit-counting/CONTRACTION_AUDIT.md):
  forced unaffordable-stage skips can be removed before graph generation.
  Active stage and defect coordinates stay. Suffix minima descend through
  contraction; prefix minima and merged weights agree only at surviving
  aggregates. The order-six hostile rules out pointwise prefix descent.
* [Antihydra report](../session64_lca_reverse_valuation_20261026/session64-reverse-valuation/REPORT.md)
  and its [audit](../session64_lca_reverse_valuation_20261026/session64-audit-dynamics/REPORT.md):
  full-shift coding requires removal of the singular preimage trees. Its
  same-word conjugacy is a 2-adic isometry, but does not preserve integers.
  Twisted integer trajectories are unbranched rays. The golden harmonic
  ruin probability is exact in Haar measure; the marked source remains open.

## Consequence for this session

Our separately constructed reverse implementation and repair certificates
are retained as cross-checks, not announced as a newer or stronger frontier.
`reverse_ordinary_check.py` uses the old actual transition helpers, compares
independent odd partitions through N=60 and normalized precisions 1..8.
It does not replace Session64's stronger production/audit package.

`check_repairs.py` independently reconstructs all six blue slices with the
pinned earlier Gaussian solver and directly replays twenty literal sources.
Its order12 path-reframing repair uses five moves; Session64's three-move
recipe is stronger. Our four order14 witnesses were found by a separate
SAT search, then checked without a solver dependency. Neither certificate
set is described as proving unrestricted black-core connectivity.

Our n14 searches restricted to three/four forward contiguous blocks failed
in their declared banks. Unrestricted first-word search found certificates.
This is a recorded reason to enlarge the grammar, not an obstruction theorem.
The inherited initial SAT timeouts likewise mean unknown, never impossibility.

The incoming work redirects the anchor to the non-UFD ideal formula and the
fixed-twist Humbert/elliptic comparison. The wildcard now extends the golden
conditional law and quantifies source-concentration cost. This avoids spending
the remaining session re-proving the stronger coding, order11 or contraction
results already present on the shared surface.

## Fourth incoming checkpoint: integer roots and guarded affine rays

Read commits `6a8ad5d3bf`, `f9738023f6` and `10695de67d`, especially
the complete [integer ancestry proof](../lca_ports_reverse_20261025/arithmetic/ANCESTRY.md),
[inverse-depth density proof](../lca_ports_reverse_20261025/arithmetic/DENSITY.md),
and [late integration](../lca_ports_reverse_20261025/LATE_INTEGRATION.md).
These provide a stronger source-specific coordinate than a valuation prefix:
for the twisted integer map, the terminating inverse algorithm gives
`y -> (root(y), chronological word from that root)` exactly. The source27
condition is `root(y)=27`; the word then supplies the integer credit.
Neither coordinate can be dropped. This is an unbranched-ray decomposition
of that expanding integer map, not a claim about the ordinary Collatz graph.

There are now three rigorously different tails to keep separate:

| Quantity | Measure / varied parameter | Exact formula or bound |
|---|---|---|
| Eventual barrier ruin from credit a | Haar-distributed valuation itinerary | q^a |
| At least j inverse ancestors of X(t)=(128*64^t+7)/9 | Natural density in t | 2^(-j) |
| Length-m cylinder of one fixed itinerary | Haar measure on starting points | 2^(-S_m) |

The density proof transports finite ternary congruences using the isometry
`v3(X(t)-X(s))=v3(t-s)`, then bounds the omitted large-valuation tail by
`2*3^(-M)`. It does not interchange an uncontrolled infinite union and
natural density. Each fixed integer root receives a density-zero set of
reset parameters. Thus it separates typical reset examples from the marked
ray while leaving possible exceptional intersections, and the chronological
credit inequality on that ray, OPEN. This sharpens the role of the source
concentration cost in DYNAMICS_EXTENSION rather than overcoming it.

The new A6 [accelerator source](../../../04-computation/bb6_cryptids_s56/a6_accel_20261025.py)
suggests a useful algebraic control that can be stated independently of its
large simulation. Let one **fixed native word**, with every branch and port
guard retained, have domain

    D={x in Z^d : A x <= b} intersect (c+Lambda)

and affine return F(x)=M x+e. Suppose x0,x0+v belong to D and
F(x0)=x0+v, F(x0+v)=x0+2v. Then Mv=v and v belongs to Lambda.
For any integer K>=1 with x0+Kv in D, every x0+kv, 0<=k<=K,
is in D and returns to x0+(k+1)v. The proof is affine linearity,
convexity of each inequality, and closure of the lattice under integer
multiples. The entire infinite ray is certified **iff A v <= 0**,
given the two starting conditions. Necessity follows from any positive
row slope eventually violating its finite bound. This elementary lemma
is PROVED; applicability to a proposed native-word domain must be proved
separately. Strict rational inequalities on integer points can first be
scaled and converted to weak integer inequalities.

The failure boundary matters: three samples in a union of branch cells
do not establish convexity, and even a very large finite K does not replace
the slope test. We read the accelerator's jump logic but do not certify its
infinite-capture cutoff or reproduce its deep orbit. If native step counts
are needed, sum the actual clock function along the retained word; a phase
count is not that clock. For a finite Boolean port cube a nonzero integer
ray cannot stay inside forever. An unbounded application must therefore
retain the ports while moving genuine size/debt coordinates, and justify
that its clone mergers preserve the same guarded domain. Likewise, fixed
valuation inequalities do not certify unit squareclasses or nonreflexive
coefficient ideals. These are the precise conditions for reusing the ray
mechanism across the three lanes, not a claimed universal construction.

## Remaining targets

A later targeted Humbert search also recovered the already stronger
[five elliptic factors over Q](../six_seven_gale_cube_20261025/GENUS_FIVE_ELLIPTIC_FACTORS.md),
the [branch cube theorem](../six_seven_gale_cube_20261025/BRANCH_CUBE_AND_DOUBLE_COVER.md),
and [Session60's six-branch/Fano/Richelot bridge](../session60_six_seven_geometry_20261025/recovered-bridges/THEOREM.md)
with its [normalization audit](../session60_six_seven_geometry_20261025/audits/curves/CURVE_BRIDGE_AUDIT.md).
They distinguish the base-line class 2H from H and retain the Richelot
determinant's arithmetic twist. F6 now extends the full five-factor theorem
to delta and identifies the three factors belonging to our marked quotient;
we do not replace that older result with the weaker statement that one
elliptic quotient exists. The one-root fibre-product reconstruction and
the two-dimensional Prym are explicitly different objects.

1. Compute coefficient ideals and comparison torsors for additional fixed-twist
   norm-seven boundary components, including simultaneous root collisions.
   Use the normalized operator bank, not just the old visible roots.
2. Apply the explicit elliptic auxiliary equation while keeping the same
   t-port of a genus-three point. Its rational points alone do not solve the
   fibre-product problem.
3. For ROOT transport, either pay the proved concentration cost or certify a
   larger event than a shrinking source cylinder. Uniformly integrable
   changes of Haar weights cannot yield a fixed positive singleton mass.
4. Beyond the incoming administrative skip contraction, any further A000568
   state merger must preserve the complete continuation kernel. Neither
   LCA shape nor coincident low-precision residues supplies that equivalence.

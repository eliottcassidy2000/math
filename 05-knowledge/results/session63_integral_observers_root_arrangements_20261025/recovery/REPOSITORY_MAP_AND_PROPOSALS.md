# Session63 source recovery and ranked theory map

Recovery by `map_new_unifying_results`, 2026-10-25. This is a source/navigation
and finite-control handoff, not a claim of historical novelty. Sources were
read at the commits in `SOURCE_PINS.json`; inherited theorem statuses below
are reported, not a substitute for replaying their independent audits.

## 1. The two user-supplied targets are not yet the same source object

At Session62 closure `61a9ae67`, the exact clone source is
`session62_agent_theme_bridges_20261025/clone/README.md`. It proves the odd
action theorem `ord_(t=-1) P_A=min_g #cycles(g)`, versus invariant XOR rank
`#(A-orbits)`. Its genuine examples are 24 vertices (regular C3^2, depths
1 versus 3), 33 vertices (4 versus 6), and 63 vertices (1 versus 9).
It does **not** contain the 50-vertex/mod-5 depth-21/mod-2 depth-25 claim.
That user-supplied target is assigned to the cancellation producer; it
must not be presented as a recovered Session62 theorem.

The exact number `q^2-6q+12` belongs to the branch complement

    U=P^2 \ {E(E+x)(E-x)(E+z)(E-z)(E+x+z)(E-x-z)=0}.

The seven-root surface has degree 64 over P2, but its **open etale part**
is the degree-64 torsor over U. The complement has 12 double and 3 triple
intersection points in characteristic not 2 or 3. Its point count is
`q^2-6q+12`; after writing `t=q-1`, this is `t^2-4t+7`.
The exact count and fresh 64-twist average are independently checked in
the attached finite checker. The full singular projective cover has fresh
twist average `q^2+q+1`, not the complement polynomial. Counting a resolution
requires the node/exceptions sidecar; the producer owns the general proof.

At characteristic 3 the seven lines remain distinct, but there are 6 double
and 5 triple points, giving `q^2-6q+10`. Characteristic 2 collapses pairs of
lines and destroys the quadratic etale-torsor interpretation. Neither is
an ignorable bad-prime convention.

The base has the class `[U]=L^2-6L+12` by elementary arrangement
inclusion-exclusion. This does **not** make an average of a fixed finite
set of twist motives equal to that class. Fixed nonsquare constants in Fp
become squares over Fp2: fresh twists over each extension and base changes
of the same twists have different averages. The twist producer is retaining
this distinction explicitly.

## 2. Exact nearest source map

All short paths below are under `05-knowledge/results/` unless indicated.

| Source | Exact mechanism | Load-bearing boundary |
|---|---|---|
| `session58_hourglass_complement_geometry_20261025/hourglass-fiber-product/fiber_product.md` | Actual shared-root fiber product, finite locally free degrees 4,4,16; C2-even/odd energy formula | Agree on the same root-valued base, not just equal marginal counts or squared labels |
| `session59_f1_geometry_and_information_20261025/toric-shared-roots/RESULT.md` | AP-root base becomes cross-polytope toric after adjoining i; exact blue-algebra pushout retains additive equations | Original Q-form not split; interior branch obstructs toric cover morphism; literal blue F1 assignment count differs from formal q=1 |
| same packet, `toric-shared-roots/TORIFICATION_OBSTRUCTION.md` | Seven-root cover degree64, complete intersection, 24 nodes, canonical p_g=7; six-root K3 threshold | No-fourfold condition for the displayed crepant calculation; torification is not every F1 framework |
| same packet, `graph-magic-degeneration/RESULT.md` | Graph incidence lattice, seam obstruction A^(k-1), smooth G_t-family u+v+tuv, ordinary root pullback | Include isolated common vertices; saturate cycle equations before t=0; ordinary squaring is not group doubling |
| same packet, `graph-magic-degeneration/branch-addendum.md` and `distinctness-addendum.md` | Counts in every finite abelian coefficient group via Smith factors and h_d=|A[d]| | Cardinality alone loses torsion; distinctness introduces 5-torsion and mixed 2/3 behavior |
| `session60_six_seven_geometry_20261025/thresholds/ROOT_THRESHOLDS.md` | General X(d,n;r), canonical coefficient (r-1)n-r(d+1), coordinate-section ladder | Root count alone loses base dimension, incidence and singular discrepancies |
| `six_seven_geometry_20261025/root_deletion/DELETION_AND_CUBE.md` | Branch components as Boolean cube, exact matrix M_t=2^(4-t)((2-t)I+A(Q_t)), genuine half-cube divisor relations | Smooth quotient is a blowup of K3; seven nonconstant modes are not seven independent divisor classes |
| `session62_agent_theme_bridges_20261025/clone/README.md` | Cycle minimum, common invariants, orientation sidecar, exterior Molien model | Rational positivity proof is not modular multiplicity; even subset actions and complex representations are hostiles |
| same packet, `common-action/COMMON_ACTION_BRIDGE.md` | Same permutation module sends action to min individual fixed space and chosen graph H1 | One element versus whole group; graph holonomy is additional data, not an intrinsic tournament move graph |
| same packet, `incidence/cyclic_common_kernel.md` | Physical line sums become step-three path forest plus pc=0 | Geometry characteristic and coefficient characteristic differ; weights (-3,-9,1) degenerate mod3 |
| same packet, `palettes/RESULT.md` and audit | Exact joint palette quotient and weighted response ranks; integer row-sum lattice index q-1 | Shared exterior assignments and aligned palette labels are essential; empty word adds a unimodular coordinate |
| `hadwiger_six_seven_eight_20261025/cube_boundary/RESULT.md` (latest eaa81fd8) | Exact unit-graph boundary functions; full Walsh coefficients glue by convolution | Pairwise feasibility shadows lose correlations; pinned external colors are part of the example |
| same packet, `projective_eight/RESULT.md` and `SYNTHESIS.md` | Marked conic-to-S(3,4,8) correspondence, seven Walsh partitions, 28 frame choices | Affine shears preserve design but not conic intersection; projective extension is not metric folding |

## 3. Incoming source that changes the priority ordering

Commit `4bb68b99306a0eb7aa046a0c3c102ed3bd66d97d` arrived after recovery began.
Its files are exactly:

* `05-knowledge/results/normseven_f1_20261025/KUMMER_AND_F1.md`
  (Git blob `9cb9e9520b2114208de7be69add56ae36920013b`). Over a nonzero
  commutative ring with unit radicands, marked mth-root algebras transform
  isomorphically iff the exponent matrix is invertible mod m. Inverting m
  is unnecessary for the algebra statement, though etaleness is separate.
  The original norm-seven Boolean stack has determinant 3; changing D mask
  4 to 5 gives determinant 1, with an explicit inverse. The base radicand
  correction is essential. Every mod-2 invertible character matrix admits
  a Laurent unimodular lift preserving each marked mod-2 character, by
  lifting elementary matrices. This is not an affine-space automorphism
  across zero coordinates and does not encode additive E+x relations.
* `05-knowledge/results/normseven_f1_20261025/INTEGRAL_WRAP.md`
  (Git blob `fe248b3bf28230411b6a8e0ead2c0cffa88f15c2`). It extends the
  Session62 physical incidence answer to every commutative coefficient ring
  and cyclic N>=5 coprime to 3: kernel `R^d direct_sum R[D]`, with exactly
  one possible nonunit Smith factor D. Depending on path anchors and parity,
  D is 1, N, or gcd(N,beta). It includes D=4 at N=16 and D=49 at full N=49.
  Lengths divisible by 3 remain an explicit frontier. Do not duplicate
  this as new Session63 work or use field nullities to infer torsion exponent.

These incoming files were read and pinned, but their proof/checker package
was not independently replayed by this recovery lane.

## 4. Ranked unifying proposals, with actual operations and hostiles

### 1. Decorated Kummer pullbacks rather than torification

The faithful universal object is the square

    X_open = U x_(Gm^6,[2]) Gm^6,
    U -> Gm^6, (x,z) |-> (1+x,1-x,1+z,1-z,1+x+z,1-x-z).

The root operation is multiplicative; the embedding is additive arrangement
data. `KUMMER_AND_F1` now supplies all valid marked monomial reparameterizations
and their correction factors. **Next useful test:** transform the entire
character-sum vector with the explicit unimodular decoder, retaining base
units, then compare direct root counts pointwise. **Canonical hostile:**
extend the Laurent inverse through zeros, or drop the additive equations;
both change the object. A blue presentation records these equations but is
not by itself a geometric torification.

### 2. Conditional Fourier transfer as the common gluing calculus

For the six quadratic character signs on U, let n_a count base points of
squareclass signature a in F2^6. Let

    A_S = sum_a (-1)^(a dot S) n_a,
    N_e = 64 n_e = sum_S (-1)^(e dot S) A_S.

The mean keeps only A_empty=|U|. Variance is `sum_(S nonempty) A_S^2`,
and shift covariance is `sum_(S nonempty) (-1)^(h dot S) A_S^2`.
This is the exact same finite Fourier operation as Session58's C2
even/odd covariance and the latest cube boundary's Walsh gluing, on
different marked carriers. **Next useful test:** recover selected twist
counts and common-base quotient products from the full vector rather than
mean counts. **Hostile:** q=5 and q=7 hourglass covariance have opposite
signs; arbitrary gluing cannot be reconstructed from marginals.

### 3. Separate intrinsic integral defects from presentation defects

Keep the original integer map/lattice and then specialize coefficients.
The incoming wrap theorem gives actual cyclic torsion R[D]. Palette
response lattices supply (q-1)-torsion; sign-orbit representations supply
orientation obstructions; modular root multiplicities use Taylor
coefficient valuations. These are related by integral specialization, not
one shared numerical rank. **Next useful test:** a table of exact integer
objects, allowed unimodular operations, and coefficient rings, with an
explicit specialization certificate for each claimed rank jump.
**Canonical hostile:** N=49 has one extra F7 dimension but order49 torsion;
the optional determinant3 of a Boolean Kummer lift can be removed with
base-unit correction, while a physical Smith factor cannot be removed by
unimodular coordinates. This contrast is now proved in the incoming sources.

### 4. Seam complexes with local coefficients and affine gains

Graph-potential overlap obstruction, palette holonomy, and the same-action
bridge all come from keeping a concrete boundary differential and its
kernel/cokernel. The operation is restriction to pieces, followed by the
actual compatibility map; not comparison of dimensions. **Next useful
test:** a fully specified exterior with unit-segment seams, where the
retained joint relation is transported by its actual gain matrices.
**Hostiles:** ignored isolated common vertices erase four K3,3 equations;
two loops can each fix a frame but have no common frame. Minimum individual
fixed space is not the invariant space of all holonomy generators.

### 5. Effectivity tests for scalar invariants

Distinguish an element of a group-completed counting ring from an effective
geometric/graded object. The open complement has t^2-4t+7, the regular C3^2
subset polynomial divided by (1+t)^3 has negative coefficients, and a
twisted singular toric surface can have the same complete zeta function as
a split one. **Next useful test:** declare the exact positive cone being
tested (split-torus strata, graded multiplicities, physical responses),
then exhibit its first forbidden coefficient. **Hostile:** negativity
rules out that effective interpretation only, not every F1 framework,
every free action, or every possible change of carrier.

### 6. Marked cube representations, not numerology of seven

The same marked F2^3 Fourier basis organizes K3 intersection modes and
palette boundary modes. The map is explicit after choosing an origin and
three directions; it does not equate the underlying geometric predicates.
The incoming unimodular Kummer lift strengthens the marked-character
view. **Next useful test:** record the subgroup preserving the marked
directions, intersection form, and arithmetic equations simultaneously.
**Hostile:** the 28 sextic labels and 28 projective quadrangles have
different permutation characters; an affine shear preserves S(3,4,8)
but destroys cube adjacency.

## 5. Recent correction ledger to keep visible

* Finite quotient after root deletion is a blowup of the K3, not the
  minimal K3 itself; branch nodes matter.
* Canonical root signs carry the residue determinant twist.
* Six-root K3 polarization has degree 8 and sectional genus 5, not genus 6.
* Torification obstruction is not a universal F1 no-go; marked monoid,
  blue-assignment, Euler-specialization and rational-point counts differ.
* Partial-rank fixed spaces do not supply simultaneous equivariant toggles.
* Whole-word palette response includes an empty-word coordinate omitted by
  a nonempty-prefix rank calculation.
* Physical field characteristic is not coefficient characteristic; literal
  incidence weights can vanish even while the unweighted common kernel
  theorem remains valid.
* The inherited lattice chromatic-cap argument has a published replacement
  in Session61; historical nonarithmetic direction counts are not evidence.

## 6. Finite artifact handoff

`check_twist_histograms.py` is independent standard-library code. It builds
prime fields and quadratic extensions directly, enumerates projective
base points and root counts through square-root tables, compares the
six-bit histograms and Walsh inversion, and checks all 64 shift
covariances. The declared field universe is q=3,5,7,9,11,13,17,25,49.
Normal and optimized outputs are compared byte for byte. Resolution counts
use the ordinary A1 replacement rule and are marked as such, not an
independent proof of surface geometry. See the companion result JSONs.

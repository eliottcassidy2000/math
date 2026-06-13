# HYP-2450 - Diagonal coefficient tilings are a Cohn/magnetization quotient of tournament space

**Status:** OPEN synthesis; finite quotient atlas verified.
**Source:** codex-2026-06-12.
**Companions:** HYP-2447, HYP-2448, HYP-2445, HYP-2444, HYP-2430,
HYP-2425, OPEN-Q-070.
**Computation:** `04-computation/coefficient_tiling_prime_bridge_codex.py`;
stored output `05-knowledge/results/coefficient_tiling_prime_bridge_codex.out`.
**External anchors:** Singh arXiv:2411.18366; Iravanian arXiv:2410.15880;
Church arXiv:2508.14876.

## Statement

The user's triangular coefficient picture is a legitimate quotient of fixed
Hamiltonian-path tournament tilings.

For `N=n+1` vertices, the gap-`d` diagonal contains `N-d` arcs.  Thus the
degree-5 picture has layer sizes

```text
1, 2, 3, 4, 5
```

from apex to base.  A tournament tiling can be projected to a coefficient
profile in at least two useful ways:

```text
count profile:       c_d = # forward arcs on gap d
magnetization:       A_d = #forward - #backward = 2c_d - (N-d)
```

The count profile gives a Cohn digit-polynomial lane.  If `b >= N` and

```text
P_b = a_0 + sum_d c_d b^d
```

is prime, then the digit polynomial `a_0 + sum_d c_d x^d` is irreducible over
`Z`, provided it has positive degree.  This is not a new irreducibility
criterion; it is Cohn's criterion installed directly on the tournament
diagonals.

The magnetization profile gives the sign/magnitude lane.  A fixed magnitude
vector `(|A_1|,...,|A_{N-1}|)` is a slice through the tiling hypercube, and the
sign of each coefficient records the majority orientation of the corresponding
gap layer.

The resulting object is not "a polynomial is a tournament."  It is:

```text
tournament tiling -> coefficient quotient -> irreducibility/prime certificate
```

and the fiber over a coefficient quotient is where the lost Hamiltonian-path,
SCC, support, and code/LRC information lives.

## Exact Finite Evidence

The stored atlas enumerates diagonal profiles for `N=3..7`.

Full count-profile grids:

```text
N=3: profiles=6,    weighted tilings=8,       Cohn prime profiles=2,   mismatches=0
N=4: profiles=24,   weighted tilings=64,      Cohn prime profiles=11,  mismatches=0
N=5: profiles=120,  weighted tilings=1024,    Cohn prime profiles=28,  mismatches=0
N=6: profiles=720,  weighted tilings=32768,   Cohn prime profiles=271, mismatches=0
N=7: profiles=5040, weighted tilings=2097152, Cohn prime profiles=580, mismatches=0
```

The zero mismatch line is expected from Cohn, but it confirms that the
diagonal-count encoding respects the base-digit side condition.

Fixed-base-path grids treat the gap-1 layer as the constant Hamiltonian path:

```text
N=6: profiles=120, weighted tilings=1024,
     count-polynomial irreducible profiles=115,
     positive-degree Cohn prime profiles=57,
     centered irreducible profiles=96,
     weighted centered-irreducible tilings=859.
```

The coefficient quotient is deliberately lossy.  In the fixed-path `N=6`
fiber audit:

```text
tilings = 1024
profiles = 120
profiles with H-variation = 91
max H-spread inside one profile = 34
max-spread profile = (5,1,1,1,1)
```

Thus the coefficient polynomial usually does not determine `H(T)`.  That is a
feature for this program: the polynomial is the scalar shadow, and the fiber is
the support channel.

## Magnitude Slices

For `N=6`, the magnitude slices show how the user's coefficient-sign idea
behaves.

```text
max_layer_magnitudes      (5,4,3,2,1): 32 distinct polynomials, 26 irreducible
minimum_parity_magnitudes (1,0,1,0,1):  8 distinct polynomials,  8 irreducible
alternating_rigid_soft    (1,4,1,2,1): 32 distinct polynomials, 24 irreducible
apex_rigid_base_soft      (1,0,1,0,1):  8 distinct polynomials,  8 irreducible
```

The parity-minimum slice is especially provocative: many sign assignments
collapse to only 8 distinct polynomials, and all 8 are irreducible in the
pilot.  This suggests a "thin slice" proof route: characterize magnitude
vectors whose parity constraints force irreducibility or force a small family
of factor patterns.

## Procedural Hiding Places

The script generated ten possible tournament hiding places:

```text
diagonal_count_cohn_map
centered_magnetization_slice
fiber_entropy_vs_irreducibility
root_argument_tournament
newton_polygon_edge_tournament
finite_field_factor_race
evaluation_time_tournament
derivative_trace_tournament
galois_orbit_recombination
code_support_polynomial_lift
```

The majority tournament on novelty, testability, prime-polynomial bridge,
tiling bridge, support-transfer, and risk was transitive in this run:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles = 0
edge_flips_vs_novelty_only = 15
leader = diagonal_count_cohn_map
runner-up = centered_magnetization_slice
```

The absence of cycles is not a failure.  It says the first useful procedural
ranking is a spine: start with the Cohn diagonal quotient, then add
magnetization slices, then audit the fiber entropy that the polynomial
forgets.

## Transfer To LRC14 And The 72-Code Gate

For LRC14, a denominator/blocker ledger can be read as a coefficient profile:
low denominators, shell-27 classes, 13-clock residues, divisor fibers, and
Bprime/owner exceptions become layers.  The Cohn lane is the warning that a
scalar prime can certify an atom only when its place-value side channel is
retained.  The LRC analogue should therefore not ask merely whether `q` is
blocked, but which diagonal/layer resource created the block.

For the extremal Type II `[72,36,16]` problem, a scalar weight enumerator is a
coefficient polynomial shadow.  The support/matroid/design realization is the
tiling fiber.  HYP-2430's Tutte support gate is exactly the instruction to
work in that fiber rather than in the scalar quotient alone.

The working dictionary becomes:

```text
coefficient profile        <-> scalar quotient / weight enumerator / q-ledger
tiling fiber               <-> support realization / design incidence / owner data
Cohn prime profile          <-> scalar atom with retained place-value address
centered magnitude slice    <-> sign/magnitude support gate
H-variation inside profile  <-> information destroyed by scalar projection
```

## Assumption Challenge

Alternate vertex sets considered: arcs, coefficient layers, layer profiles,
fixed coefficient-magnitude slices, roots, factor subsets, finite primes,
evaluation times, Newton polygon edges, derivative obligations, code support
moves, and LRC denominator resources.

Chosen primary quotient: diagonal profiles of fixed-path tilings.

Preserved: the triangular coefficient-layer geometry, Cohn prime certificates,
sign/magnitude slices, and the fact that many tournaments share one
coefficient profile.

Destroyed: individual arc arrangement, most `H`/SCC data, and actual
asymptotic Bunyakovsky information.  The fiber-H audit measures this loss.

Challenged assumption: the base Hamiltonian path need not merely be a gauge
choice.  In the coefficient picture it can be treated as the constant term,
while higher diagonals become the polynomial's mutable coefficients.

## Next Moves

1. Prove the exact Cohn-diagonal lemma for fixed-path profiles in arbitrary
   base `b >= N`, including the positive-degree guard.
2. Characterize magnetization magnitude vectors whose sign slices are forced
   irreducible, forced reducible, or have bounded factor patterns.
3. Replace profile-level `H` spread by SCC, directed-cycle, and Hamiltonian
   path count distributions inside each coefficient fiber.
4. Build a finite-field factor-race version: vertices are primes `p`, and
   edges compare factorization type of the same coefficient-tiling polynomial
   modulo `p`.
5. Translate the `code_support_polynomial_lift` row into a `[72,36,16]`
   support atlas: coefficient layers are weight/design incidence constraints,
   and the fiber records binary matroid realizability.

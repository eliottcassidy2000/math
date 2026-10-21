# Independent frozen audit of Session 19 ellipsoid packets

**PASS, 2026-10-21 UTC. No blocking mathematical or artifact defect found.**

This audit covers the two exact source packets identified below, not a
later revision. Their complete written proofs and checker sources were
read independently. The retained audit does not rely on earlier informal
review results. No producer file or repository file was modified, and no
Git write or publication was performed by this audit.

## Exact input identity

- `session19-ellipsoid-optimum/MANIFEST.json`:
  `050fc733781c929b8d2eb36a128cec3bc3d0f7b85c95f491941dcb3dbf849dad`.
- `session19-triangle-ellipsoid-repair/MANIFEST.json`:
  `43ec88a26d022a7e10f9968a1a4e84f2634325c42550ab7defe34f9151ab17fe`.

The complete source packets are retained under `frozen_sources/`. Every
file listed by their manifests, its byte length and SHA-256, both manifest
sidecars, and the exact input inventories were verified before and after
execution. The original sibling packets were also rechecked unchanged.

## Retained results

- Base producer: **71,504** exact gates, normal and optimized output
  byte-identical to the frozen `check_exact.out`.
- Triangle producer: **31,344** exact gates, likewise byte-identical.
- New independent checker: **7,914** exact gates, normal and optimized
  output byte-identical. It imports no producer code.
- All six executions returned zero with empty stderr. The outputs and
  machine-readable execution receipt are in `receipts/`.

The independent implementation uses rational ambient projection matrices,
Gauss--Jordan inverses, determinant evaluation after adding the normal
rank-one projector, exact consistency ranks for the diagonal-lift linear
system, and direct facet-basis matrix construction. Its deliberately
different finite universe uses the six rational speeds
`{1/2,2/3,1,3/2,2,3}` in dimensions 3 through 5, additional boundary
controls in dimensions 3 and 4, every labelled four-vertex tree basis on
an explicit rational box, and the nonempty repair family through
dimension 10. It also checks translated covariance equality controls,
the two requested optimum examples, the threshold equality example,
the two-speed lift continuum, and the compression/inverse hostile.

Finite checks are corroboration, not proofs of the universal real-speed
claims. The written proof audit below is the basis of that conclusion.

## Reproduction

Python 3.8 or newer, standard library only. The recorded interpreter was
Python 3.12.14. From any working directory:

```sh
python3 -B /path/to/session19-ellipsoid-audit/replay.py
python3 -B -O /path/to/session19-ellipsoid-audit/replay.py
```

Both commands verify the frozen sources, execute all three checkers in
normal and optimized modes, and compare the retained outputs and receipt.
No network, repository, or original sibling directory is needed. To also
check the live original packets, append
`--originals /path/to/their/common/parent`. The explicit `--record` option
regenerates only this audit's receipts; do not use it for verification.

The independent companion alone is run with:

```sh
python3 -B /path/to/session19-ellipsoid-audit/check_independent.py
python3 -B -O /path/to/session19-ellipsoid-audit/check_independent.py
```

`MANIFEST.json` freezes this audit's own inventory, including its input
snapshots and receipts; `MANIFEST.sha256` pins the audit manifest without
self-reference. `SOURCE-PINS.md` records the independently inspected
external and repository attributions.

## Written mathematical audit: base packet

1. **Carrier and exact constraints.** For positive real `v`, the support
   identity `h_D(y)=||y||_1` holds on `H=v^perp`. The pair section vertices
   are the opposite-sign edge intersections of the cross-polytope with
   `H`. The independent interval-intersection argument gives the same
   primal slab system and removes reliance on a hidden polarity premise.
   A nonnegative diagonal covariance remains positive semidefinite even
   when singular, so the convex-vertex test remains valid at zero entries.
   No integer-speed or primitive hypothesis is used here.

2. **Exact dominant cone and nondegeneracy.** Substitution of the proposed
   singular lift makes every pivot pair tight. The remaining slack has
   the claimed factorization with positive prefactor. Monotonicity in
   the two leaf speeds makes the largest two leaves decisive, with
   equality included. The star chart rows annihilate `v`, are independent,
   and map the candidate covariance to the identity. Consequently the
   ambient zero does not cause a zero-dimensional direction on `H`.

3. **Global optimum, not merely diagonal optimum.** The selected star
   slabs define an enclosing parallelotope, and its chart image is a
   cube. Coordinate support bounds apply to every translated ellipsoid.
   Hadamard's inequality bounds its determinant by the product of its row
   lengths, which is at most the product of `1-|c_i|`, hence at most one.
   Equality forces zero center and orthonormal rows. Transport by the
   invertible chart proves global uniqueness; it is not a relaxation
   that only optimizes diagonal matrices. Singular competitors have zero
   intrinsic volume.

4. **Intrinsic determinant and normalization.** The Cauchy--Binet identity
   for the compressed diagonal covariance is polynomial in the entries
   of `Q`, so it extends to singular lifts. The denominator `||v||^2`
   is necessary and correct. Our separate computation instead forms
   `PQP+vv^T/||v||^2`: the added projector has eigenvalue one in the
   normal direction, so its determinant is the intrinsic covariance
   determinant on `H`. This independently verifies the Euclidean volume
   normalization and the chart Jacobian. A determinant is squared volume
   divided by squared unit-ball volume, not volume itself.

5. **Diagonal cap argument and endpoint.** Coordinatewise monotonicity
   of the determinant numerator allows replacement by pivot caps even
   when that cap vector is not feasible for the remaining facets. The
   polynomial identity `N_cap=C(f-t f')` and derivative `-Ct f''` are
   correct. The second derivative is positive before the first cap zero
   for dimension at least two. Strict loss at the endpoint follows by
   integrating over a nonempty initial interval; no division by a zero
   cap is made. At pivot coordinate zero the product objective forces
   every other cap tight.

6. **Boundary-optimizer equivalence.** Compactness follows from pairwise
   upper bounds on each coordinate. The feasible identity covariance
   has positive determinant, excluding two or more zero coordinates at
   a maximum. With exactly one zero coordinate, the stated direction
   preserves all potentially tight constraints except the chosen slack
   pivot pair and decreases the other constraints. For sufficiently
   small positive step it remains nonnegative. Its determinant-numerator
   derivative is exactly `2 v_i^2 product_(j outside {p,i}) q_j`, positive.
   Thus all pivot constraints must be tight, and feasibility then forces
   the unique largest-speed pivot and dominant cone. The converse was
   already proved. Hence outside the cone every maximizing nonnegative
   diagonal lift is positive, without identifying that optimizer.

7. **Dimension exception and injectivity.** For two speeds the single
   pair constraint permits a continuum of maximizing lifts representing
   the same interval. For at least three speeds, equality of diagonal
   restrictions gives pair equations for `d_i/v_i^2`; any triangle forces
   all these quantities to vanish. Thus the higher-dimensional lift
   uniqueness statement does not accidentally include the two-speed
   exception.

8. **Nonstar facet-tree obstruction.** Incidence rank after the invertible
   coordinate scaling and restriction to `H` is `k-c`. A basis of `k-1`
   pair facets therefore corresponds precisely to a spanning tree.
   Telescoping along its unique paths gives the stated sum-of-squares
   containment criterion. A length-two path forces its middle reciprocal
   speed strictly below both ends. Two adjacent length-two paths along
   a length-three path contradict one another. Every nonstar tree has
   such a path, proving impossibility for all positive real speeds.
   This rules out only facet-basis parallelotope certificates, not
   overcomplete contact configurations.

## Written mathematical audit: triangle addendum

1. **Regime and all-facet feasibility.** Pulling back the pair slabs into
   star coordinates yields the positive coefficients stated in the
   packet. Each is at most one because `M` is a largest speed. Exactly
   one violating pair is equivalent to the two ordered inequalities,
   with ties counted separately. The prescribed positive correlation
   is less than one since the squared difference of the two positive
   coefficients is less than one. Thus the covariance is positive
   definite. The exceptional support becomes one exactly; every other
   leaf support sees no changed cross term. Axis supports remain one.
   No omitted higher-support facets are present.

2. **Positive certificate and translated competitors.** The stated axis
   and exceptional-pair weights are strictly positive, including cases
   with a coefficient equal to one. Their rank-one sum is precisely
   the inverse covariance, and the contact equalities give total weight
   `m`. For any feasible translated covariance `T`, each support bound
   has `|w.h|<=1`, so squaring gives
   `w^T T w <= (1-|w.h|)^2 <= 1`. Summation produces the trace bound.
   Applying scalar AM--GM to the nonnegative eigenvalues of the
   conjugated covariance gives the determinant bound. Equality makes
   every eigenvalue one; positive axis weights then force every center
   coordinate zero. This establishes uniqueness among all ellipsoids,
   not just centered or diagonal ones. Additional tight nonviolating
   facets need no positive certificate weight.

3. **Determinants and complete triple classification.** The difference
   of squares cancels precisely the two affected star factors; the
   resulting formula is homogeneous of degree zero and has the correct
   intrinsic normalization. A triple has only one leaf pair, so the
   base cone and strict repair regime exhaust all positive triples.
   The repaired diagonal lift has strictly positive entries, matches
   the chart covariance, and meets the zero-pivot lift continuously at
   threshold equality. This includes repeated and tied-largest speeds.

4. **Strict higher-dimensional gap from the entire diagonal class.**
   With at least three leaves, a real ambient diagonal lift has all
   off-diagonal chart entries proportional to its single pivot entry
   by nonzero positive factors. The repaired covariance has exactly
   one nonzero unordered off-diagonal pair, so it cannot be represented
   by any real diagonal lift. This is stronger than nonrepresentability
   by a nonnegative one. Compactness of the feasible nonnegative class
   and uniqueness of the global optimum turn nonrepresentation into
   a strict volume gap; merely proving nonrepresentation without
   attainment would not have sufficed. The separate exact checker also
   verifies inconsistency by augmented linear-system rank in ambient
   coordinates, rather than assuming the chart argument.

5. **Named examples and dimensional nonvacuity.** The determinants
   `22/7` at `(1,2,3)` and `1488/325` at `(1,2,3,5)` are correct.
   The latter has exactly the exceptional leaf pair `(2,3)` and
   correlation `3/28`. The family `(1,...,1,2,2,4)` has one violating
   pair in every claimed dimension, with correlation `1/9` and chart
   determinant `80/81`. The overcomplete star-plus-edge contacts contain
   a triangle and do not contradict the base tree obstruction.

## Scope guards and nonblocking evidence note

- The analytic geometry is valid for positive real speeds. Primitive
  integer hypotheses enter only the stated projected-lattice covolume
  interpretation, not the ellipsoid theorems.
- The gain is volume, not set containment or pointwise norm dominance.
  The `(1,2,4)` source chart covariance has both a larger and a smaller
  directional support than the optimum. Neither body contains the other.
- No multiple-violating-pair classification, quantitative diagonal gap,
  uniform LRC improvement, lattice-point incidence result, lonely-time
  witness, computational-prime audit, or literature-priority claim has
  been established or silently added.
- The producer's base grid branch named `grid_unique_diagonal_lift`
  does not fire in its frozen output: its selected grid does not attain
  the candidate. This is a limitation of that finite control, not a
  defect in the theorem. The written uniqueness arguments are complete,
  and the independent translated-covariance test includes the exact
  optimizers and records both equality cases. Gate totals should not be
  interpreted as a count of distinct mathematical theorems.

**Final disposition: PASS for both exact frozen packets, with the stated
scope unchanged.** This is an independent proof review plus reproducible
finite audit, not machine formalization of the universal arguments.

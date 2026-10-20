# Fixed and moving cylinder obstructions

**Status:** proved; independently reviewed.

This continuation contains three precise mechanisms:

1. A fixed rational plane curve with a normalization collision and unique
   ambient embedding has a constant-cylinder period obstruction. Finite
   intrinsic automorphisms suffice; they need not extend ambiently.
2. For
   \(E_g=\mathbb C[x,y^2,y g(x,y^2)]\), with \(g(x,0)\in\mathbb C^*\),
   generic unique embedding forces the ambient Jacobian values at all nodes
   to agree. Gauss descent and the surviving singular-point conormals make
   that scalar constant. Provided every fiber retains a collision, a
   coefficient comparison then excludes any generating pair whose Jacobian
   is a nonzero polynomial in \(x\).
3. A principal-conductor actual source algebra containing a monic quadratic
   normalization coordinate is exactly such an \(E_g\). This is an
   entire-algebra recognition theorem, not a statement about conductor
   multiplicities alone.

The prime-exponent embedding theorem supplies the moving obstruction whenever
\(2\deg_w(g)+1\) is prime. A separate cited one-line argument handles fibers
where all collisions disappear. In particular all nontrivial cubic-or-lower
\(g\) models admit no generating Keller pair.

Under the earlier actual-source hypotheses **and the additional monic
quadratic-element hypothesis**, noninvertibility forces generic conductor
degree \(m\) to have \(m+1\) composite; the first possible \(m\) is eight.
No unconditional degree-eight bound, arbitrary-pair closure, or resolution
of the complex Jacobian conjecture is claimed.

## Files

- [Moving-algebra theorem, recognition, and exact scope](MOVING_HYPERELLIPTIC_GENERATOR_OBSTRUCTION.md)
- [Fixed-cylinder finite-period lemma](FIXED_CYLINDER_FINITE_PERIOD_LEMMA.md)
- [Source and dependency record](SOURCE_AND_SCOPE.md)
- [Exact finite controls](verify_cylinder_obstruction.py)
- [Results](results.json)
- [Completed independent audit](review/AUDIT.md)

Reproduce with Python 3 and SymPy:

```sh
python3 verify_cylinder_obstruction.py --output normal.json > normal.out
python3 -O verify_cylinder_obstruction.py --output results.json > optimized.out
cmp normal.json results.json
cmp normal.out optimized.out
```

Producer result: **205 exact checks across 26 declared families**.
Independent result: **396 checks across 43 controls**. Normal and optimized
execution are byte-identical for both.

Controls include degree drop from prime leading exponent eleven to composite
exponent nine, repeated generic roots, nonconstant units on the node cover,
an empty-fiber pole trap, and a nonprincipal-conductor boundary.

The next missing input is a criterion for the monic quadratic element in
actual source algebras of low conductor degree. It is not inferred from
the conductor divisor or from arbitrary abstract gluing data.

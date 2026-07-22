# Kelvin--Farey gates and the scaled clock/binding closure

This pass deliberately combined three repo niches that had mostly been kept
apart:

```text
circle inversion / polar duality,
toric regular fans / Farey neighbors,
CRT covering clocks / affine binding pairs.
```

The first two give a clean exact certificate for the rank-two determinant
gate. The third is what actually closes an LRC plane.

## 1. Invert the tangent circles

THM-2055 writes the directional error as the support norm

```text
D(d)=h_K(Rd).
```

The gate compares this degree-one function with `||d||^2/91`. Kelvin inversion
`I(d)=d/||d||^2` removes the degree mismatch:

```text
D(d)<=||d||^2/91
iff I(d) in (1/91)R^(-1)K^o.
```

Thus every tangent circle through the origin is just the inverse image of a
facet of a fixed rational polar polygon. For primitive `d=(a,b)`, the inverse
point is the reciprocal Gaussian integer

```text
a/(a^2+b^2)+i b/(a^2+b^2)=1/(a-i b).
```

This is the lawful arithmetic structure that the failed Heegner analogy was
trying to find: visible Gaussian reciprocals tested by polygon facets. It has
no discriminant or class group unless extra structure is proved.

## 2. A Farey cone can carry a defect budget

On one normal-fan sector the gate is

```text
F(d)=||d||^2-91 w dot d>=0.
```

For a unimodular pair `u,v`, the exact expansion of `F(mu+nv)` has a positive
cross term `2mn u dot v`. If the cone is acute, that cross term can pay the
negative endpoint defects. THM-2056 packages the payment as

```text
2 u dot v >= max(0,-F(u))+max(0,-F(v)).
```

One inequality then certifies every positive integral combination. The
acuteness is essential: the safe unimodular endpoints `(91,1),(-90,-1)` hide
the unsafe mediant `(1,0)` when their dot product is negative.

This is toric resolution used as a proof certificate rather than as a
metaphor. A finite regular fan isolates every bad primitive direction as a
ray and discharges the infinite complement cone by cone.

## 3. The negative result that changed carriers

On the positive primitive one-tail plane

```text
a(1,...,13)+b e_12,
a>0, 12a+b>0,
```

the complete THM-2053 determinant residual has `640702` directions, or
`640690` after removing the twelve collision directions. This count uses the
proved radius `||d||<91*13`; it is not a truncated sample.

So the polar polygon and Farey fan are exact, but they are not by themselves
a useful terminal enumeration. The challenge to the earlier assumption is
sharp:

> the best geometric quotient need not be the best arithmetic proof carrier.

The geometry should address a direction; a phase-producing sidecar must erase
large families of addressed directions at once.

## 4. Scaled unit clocks close the whole plane

The row in this plane is

```text
S(a,w)={a,2a,...,11a,13a,w},      w=12a+b.
```

HYP-2896 closed only the `a=1` version using the `12`, `14`, and `84m`
branches. The missing scaled statement is unexpectedly elementary. If
`q=Na` with `N=12` or `14` and `q` does not divide `w`, reduce by
`g=gcd(w,q)` to a nontrivial modulus `h=q/g`. Every `h>=2` has an explicit
unit at circular distance at least `1/14`:

```text
1                         if h<=14,
(h-1)/2                   if h is odd,
h/2-1                     if h=0 mod 4,
h/2-2                     if h=2 mod 4.
```

Lift its preimage to a numerator coprime to `N`. On the `12a` clock the core
has margin `1/12`; on the `14a` clock it has margin `1/14`. Therefore:

```text
12a does not divide w       -> 12a-clock witness;
12a divides w, 14a does not -> 14a-clock witness;
both divide                  -> 84a divides w.
```

In the last case `w=84am`, so scaling the explicit HYP-2896 phase gives margin
`7m/(84m+5)>1/14`. THM-2057 thereby closes the entire plane, including all
`640690` distinct-speed determinant failures, with three symbolic leaves.

The proof actually yields a general missing-clock sieve. If no member of an
integer core `C` is divisible by some `2<=N<=14`, then a counterexample of the
form `aC union {w}` must have `Na|w`. Thus the tail must be divisible by

```text
a*lcm{N<=14:C contains no multiple of N}.
```

For the zeta core the missing clocks are exactly `12` and `14`, producing the
`84a` ray. This lcm tax is the portable lemma for other star cores.

The first transfer already works. For the adjacent core `{1,...,12}`, the
missing clocks are `13` and `14`, so only `182a|w` survives. Writing
`w=182am`, the phase

```text
t=14m/[a(182m+1)]
```

has exact minimum `14m/(182m+1)>1/14`: the tail is `-14m` modulo
`182m+1`, and the twelve AP residues are all at least `14m` from an endpoint.
Thus THM-2057 closes a second complete AP one-tail plane rather than merely
suggesting a transfer.

## 5. Transfer pattern for the remaining atlas

The promising general pattern is now:

```text
rank-two normal fan
  -> Kelvin/Farey address
  -> choose a labelled lower-rank core
  -> complete unit orbit on a safe clock
  -> killed clock forces a divisibility sublattice
  -> repeated killing forces an affine binding ray or an Euler endpoint.
```

This resembles the repo's covering-system work more than its quadratic-form
work. It also says what data cannot be projectivized away: the scale `a` is
needed because the live clocks are `12a` and `14a`, not `12` and `14`.

The open HYP-8871 target is to find such safe-unit clock orbits on every
THM-2052 star, or prove that a bounded chain of killed clocks lands in the
owner-labelled Euler branch. The finite atlas supplies the possible cores;
THM-2056 supplies exact ray addresses; THM-2057 supplies the model descent.

## Tournament analysis

Vertices are proof carriers, not runners:

```text
scaled_clock_binding
resolved_phase_height
kelvin_polar_polygon
farey_defect_cone
raw_tangent_scan
heegner_form_class
```

Pairwise observable: retention of an actual LRC witness, exactness of the
rank-two address, whole-family compression, sidecar debt, and verification
cost. The resulting diagnostic order is the displayed order: a transitive
tournament with score histogram `{0,1,2,3,4,5}`, no directed triangles,
singleton SCCs, and one Hamiltonian path.

The quotient losses explain the ranking. `scaled_clock_binding` preserves an
actual phase and margin. `kelvin_polar_polygon` and `farey_defect_cone`
preserve only the determinant gate but address every direction exactly.
`heegner_form_class` preserves neither predicate for an arbitrary column
polygon and remains last.

Artifacts: THM-2056, THM-2057, HYP-8871, MISTAKE-225, and
`lrc_kelvin_farey_scaled_core_codex_20260721.py` with its frozen output.

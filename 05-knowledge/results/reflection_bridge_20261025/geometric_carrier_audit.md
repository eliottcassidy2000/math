# Independent review of the projective divisor carrier

2026-10-25. **VERIFIED proof review; FINITE-EXACT independent controls.**
Reviewed `geometric_carrier.md` after deriving the independent seven-set
and projective-support model in `seven.md`. No mathematical failure was
found in the proposed carrier or counts. No producer code was imported.

The coordinate formula on `P^k x P^1 x P^1` is an actual projective
automorphism in the stated pointed-monoid model and after base extension:
each factor permutation takes a standard affine chart isomorphically to
another chart. The factor swap is also an automorphism. They commute
and induce exactly `(a,b,c)->(k-a,1-b,1-c)` and `(a,b,c)->(a,c,b)` on
the coordinate labels. The deleted pair is invariant, so both maps
restrict to its open complement.

Two precision edits were requested before publication:

- Over `Z`, the two deleted objects are the two coordinate **sections**,
  each a closed subscheme isomorphic to `Spec Z`; they are not single
  topological points. In the monoid model they are closed points, and
  in every finite-field fiber they are two rational points.
- Spell out specialization direction: `S subset T` means `x_S` lies in
  `closure({x_T})`. Larger supports are more generic. This has no
  bearing on divisibility of the coordinate-point labels, all of which
  are closed in the monoid model.

Independent direct support enumeration at `k=2` gave 63 support tuples
before deleting `(1,1,1)` and `(4,2,2)` in bit-mask notation, and **61**
after deletion. Exactly **10** remaining tuples have singleton support
in each factor. Independent enumeration of their exponent labels gave
fixed counts `(0,4,2)` for `(C,Q,CQ)`. Coordinate reversal on the seven
supports of `P^2` fixes exactly masks `2,5,7`, namely `{1}`, `{0,2}`,
and the full support. Its three fixed points differ from the one fixed
point of the boundary-contracted complement in `seven.md`.

The count polynomial follows by multiplication and removal of the two
rational sections. At `k=2` it is

    (1+Q+Q^2)(1+Q)^2-2
      = -1+3Q+4Q^2+3Q^3+Q^4
      = 10+24h+19h^2+7h^3+h^4,       h=Q-1.

The last expansion records the remaining standard torus orbits by
dimension: only two zero-dimensional orbits have been deleted. The
support expansion and the product-cell polynomial were checked
independently for `2<=k<=8`, `1<=Q<=7`; the value at `Q=1` in those
checks is a formal polynomial evaluation. It is not a field or a
cyclic-unit group of size zero.

The ten-label rank polynomial `3Q+4Q^2+3Q^3` differs by `Q^4-1`, because
removing the maximal coordinate point punctures the top cell rather
than removing that whole cell. This distinction is correctly stated
in the producer note. At `Q=2`, the two numbers are 61 and 46.

The seven-support model in one `P^2` factor and the three coordinate
points in that factor are different functors/strata, even though they
belong to the same explicit projective model. Consequently the positive
carrier does not circumvent the earlier obstruction to realizing
divisor complement as a tetrahedral affine symmetry, nor the obstruction
to preserving Fano incidence under boundary-contracted complement.
The carrier retains the intended finite-set action and requires the
coordinate flags to retain divisor order and rank.

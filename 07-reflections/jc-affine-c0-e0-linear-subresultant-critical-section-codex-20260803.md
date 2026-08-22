# Affine C0-E0 linear-subresultant critical section

**Status:** FINITE-EXACT PARTIAL SCOUT for the two rational controls
`C=1+x`, `E'=1` in the THM-3212 accessory fields.  This note is not canon and
does not generalize the result to arbitrary affine parameters.

## Exact object

For the THM-3289 localized cubics

```text
R1=2L(2y+1)+VA,
R2=V^3+V^2y+L(-V'y+2V^2),
L=y^2+y+(1+x)V,
```

the chosen standard subresultant sequence has degree profile

```text
(3,3,2,1,0).
```

Its final constant is the critical resultant `V^3 K`.  Its penultimate row
has a direct polynomial Bezout representation

```text
S1=U R1+W R2,
```

where

```text
U=-2V'(4V^2+V')y+4V^2(4V^2+2V'),
W=-8(4V^2+V')y-4V'.
```

The displayed row is one standard normalization, not a canonical raw PRS
iterate.  Only its ideal and its primitive coefficient pair up to a common
nonzero field unit are used; this is the unit discipline required by
MISTAKE-360.

## Boundary content and residual graph

In each accessory field the full `x`-content of the two coefficients of `S1`
is, up to a field unit, exactly the degree-44 passport boundary:

```text
S1 = unit * boundary_i * (a(x)y+b(x)),
deg(a)=36,   deg(b)=44,   gcd(a,b)=1.
```

For the degree-53 saturated residual `H`, exact computation gives

```text
gcd(H,H')=gcd(H,boundary_i)=gcd(H,V)=1,
gcd(a,H)=gcd(b,H)=1.
```

Thus both the boundary and `a` are units in
`A_H=K_i[x]/(H)`.  The Bezout row puts `a y+b` in the localized gradient
ideal, while direct substitution of

```text
y_H=-b/a in A_H
```

makes both `R1` and `R2` vanish.  Consequently

```text
(R1,R2)=(a y+b) in A_H[y]
```

up to the chosen unit normalization.  The class `y_H` has a degree-52
polynomial representative in both fields.

## Hostile branches and meaning

All four immediate specialization hazards are absent in these controls:

- `a` has no zero on `H`, so no selector split is needed;
- the constant leading coefficient `4` of `R1` excludes a common projective
  `y`-root at infinity;
- `H` is disjoint from the boundary and from `V=0`;
- `H` is squarefree, so no repeated residual branch is hidden.

The positive conclusion is precise: the reduced degree-53 residual critical
scheme is one global graph over its `x`-scheme, rather than 53 separately
selected roots.  This is stronger structure than the bare resultant count for
these two controls, but it is still only a first-coordinate critical locus.
It supplies no cofactor, second coordinate, inverse cover, or implication for
`JC(2)`.

## Reproduction

Run

```text
python3 04-computation/jc_affine_c0_e0_linear_subresultant_section_partial_scout_20260803.py
python3 -O 04-computation/jc_affine_c0_e0_linear_subresultant_section_partial_scout_20260803.py
```

and compare LF-normalized output with
`05-knowledge/results/jc_affine_c0_e0_linear_subresultant_section_partial_scout_20260803.out`.

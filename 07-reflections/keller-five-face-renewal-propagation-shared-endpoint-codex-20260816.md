# Renewal is the norm of one shared endpoint, not a second state variable

**Status:** research reflection for PROVED + INDEPENDENTLY AUDITED THM-3522.
The theorem is scoped to
the fixed THM-3495 inverse chart and assumes the cleared norm is polynomial.
It does not supply a later finite-sheet unit, image prime, degree-243 gate,
all-level tower, or general Jacobian statement.

## Inheritance pass

- **Closest proved mechanism:** THM-3513's two hybrid Newton limits for the
  fixed `J -> G` step.
- **Canonical hostile:** MISTAKE-415.  An exposed pair was once treated as a
  closed scalar recurrence before all inverse-branch faces were retained.
- **Corrected near miss:** THM-3513 made its leading `J` value root-independent
  because its exponent happened to be divisible by three.  That convenience
  fails already for `L -> H` and `H -> J`; it is not part of the mechanism.
- **Least-used sidecar:** the coefficient compatibility at the intersection
  of two complete faces.  It fixes both the input and output renewal scalars.

The active concept board was:

1. the singleton maximum-`z` face;
2. the binomial minimum-`gamma` face;
3. their common endpoint coefficient;
4. the `delta_6=gamma-k` hybrid valuation;
5. the `delta_8=gamma-3k` hybrid valuation;
6. the nonmonic product-of-roots factor;
7. polynomiality versus the finite old-`L` sheet.

The board collapsed from seven obligations to two: the first six form one
universal renewal mechanism, while polynomiality/finite-sheet control remains
a separate global gate.

## The mechanism

For packet `A(e,m)`, put

```text
r=e-2m/3,
p=2r,
d=e+r.
```

The minimum-gamma face

```text
eta z^e(27x^2z+y^3)^r
```

has a unique largest-`z` endpoint `eta 27^r x^p z^d`.  The complete
maximum-`z` face is `kappa x^p z^d`, so `kappa=eta 27^r`.

For any other source monomial, let `u` be its gamma gap and `v` its deficit
from maximum `z`-degree.  Its two hybrid gaps are exactly

```text
delta_6 gap = u+v,
delta_8 gap = u+3v.
```

Both vanish only at the common endpoint.  This is the entire uniqueness
argument; no expansion of the input polynomial is needed.

Under the two inverse scalings, that endpoint contributes a power

```text
q^rho,       rho=p-3d=-4e+2m/3.
```

The crucial correction is that `q^rho` need not be constant on the three
residual roots.  Instead one takes its norm:

```text
product(q)=2/(27A^2C)        on the maximum-z chart,
product(q)=2/D, D=27A^2C+B^3, on the gamma chart.
```

The resulting output scalars are

```text
eta'   =kappa^3 2^rho,
kappa' =eta' 27^r'.
```

Thus the output faces meet with the right coefficient automatically.  The
nonmonic Vieta denominator is not a nuisance factor; it creates the full
power `D^r'` and thereby the complete gamma face.

## Hostile controls that mattered

- `rho=-4` for `L -> H/64` and `rho=-26` for `H -> J/2^35`.  Neither is
  divisible by three, so a proof demanding sheetwise root-independence is
  false.  The exact norm scalars nevertheless agree.
- The admissible boundary `(e,m)=(3,3)` has a zero exponent in the beta packet
  and `rho=-10`; all propagated exponents remain nonnegative.
- A raw resultant with `q` equals `2`, not `2/D`.  Dividing by the nonmonic
  leading coefficient is required to obtain the function-field norm.  The
  wrong monic treatment predicts exponent `43` rather than `205` at `J -> G`.
- The residual discriminants are generically nonzero, and the constant term
  `-2` makes every negative power of `q` lawful.
- Polynomiality is retained as a hypothesis.  Hybrid asymptotics alone do not
  establish global denominator clearing.

Two concurrent clean-room audits then rebuilt the mechanism without importing
the primary companion.  One exhausts `2,500` admissible endpoint intersections
through `e<=120`; the other checks `15,250` packets through `e<=300`, including
`10,167` non-cube root-label rows, and adds split-prime Vieta controls at
`109,127,163`.  Both recover the fixed supports and scalars exactly.

## Connection contract

```text
source:      complete maximum-z and minimum-gamma faces of P
target:      the same two complete faces of L^e N(P)
map:         unique endpoint -> two hybrid valuations -> cubic norm
preserved:   endpoint coefficient and packet exponents
destroyed:   higher hybrid layers and individual root labels
sidecar:     nonmonic product(q), plus polynomiality of the cleared norm
test:        residual-cubic reduction, multiplication determinant/resultant,
             and the L/H/J scalar calibrations
```

This is a genuine connection rather than a shape analogy: it transports both
the support and coefficient compatibility of the target faces.

## Fixed consequences and remaining frontier

THM-3506/3513 give polynomial `R_5` from full `G`, so THM-3522 proves

```text
R_5 has A(1699,615).
```

THM-3521 gives polynomial `R_6`; applying the same theorem proves

```text
R_6 has A(10663,3867).
```

The renewal obligation is therefore gone at these rungs and is conditional
only on polynomiality at every later rung.  The next decisive fixed-tower
question is not another face expansion: it is the next old-`L` finite-sheet
or denominator-clearing gate for `R_6`, followed separately by image status
and degree-`243` separability.

The session used the repository method cards **Refine and saturate before
transporting a factor or shadow** (complete faces and the full nonmonic norm)
and **Classify response-state growth before naming the closed form** (repairing
MISTAKE-415's premature scalar recurrence).  No new meta-pattern card is
needed; this is a sharp new instance of both.

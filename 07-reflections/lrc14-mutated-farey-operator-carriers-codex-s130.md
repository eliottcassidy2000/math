# LRC14 Mutated Farey Operator Carriers

The useful thing from the four mutated Farey operations is not that one of them
replaces the denominator.  It is that they rank the ways a quotient can forget
the exact fraction address.

The exact gap formula is stubborn:

```text
M(S)-1/14 = (14p-q)/(14q).
```

So `q` is still the binding-scale coordinate.  The ordinary denominator is not
a decoration; it is the part of the address that turns the integer excess
`e=14p-q` into an actual analytic gap.

The additive mutation `p+q` is almost innocent near the apex.  Along the
unit-excess chain

```text
p/(14p-1) = 1/13, 2/27, 3/41, 4/55, ...
```

the denominator walks by `+14`, and the numerator walks by `+1`.  This is the
same shape as the repo's `n+2` recursion, but viewed through the Farey child
of `1/14`.  It explains why the additive payload is a good ledger: it records
the recursion without trying to be the theorem.

The product mutation `p*q` is more interesting and more dangerous.  It feels
like the `n*2` side of the same picture: a two-coordinate area, not a binding
scale.  On the row bank it reverses many true gap comparisons, so it cannot be
the proof denominator.  But it is the right place to keep the multiplicative
coimage information: divisor profiles, exact-period units, product-Mobius
packets, and the `Div(D) x B_r` side of the p0 route.

The two power mutations are not local coordinates.  They blow up magnitude.
That makes them bad proof denominators but good lie detectors.  If a proposed
tournament invariant survives `q^p` or `p^q` without noticing the distortion,
then the invariant is probably only reading a finite class label and has lost
the unbounded scale that HYP-2926 warns about.

The tournament view made this sharper.  Taking the payloads themselves as
vertices, with edges determined by risk-order agreement and Farey-locality,
gives a transitive majority tournament:

```text
q > p+q > p*q > p^q > q^p.
```

That is a small hierarchy of proof roles.  The first vertex is the theorem
coordinate.  The second is the additive recursion ledger.  The third is the
multiplicative side-channel ledger.  The last two are stress tests for
magnitude leakage.

For LRC14, this does not bypass the hard step.  It says how to avoid turning
the hard step into the wrong object.  The proof should retain the fraction as
a labelled address:

```text
(excess e, binding scale q, additive ledger p+q, multiplicative ledger p*q)
```

and only then compare tournament spectra.  A single quotient, even a clever
one, is too easy to fool.

# Borel Constructive Embedded Maximality, S673

The clean synthesis is:

```text
diagonalization is constructive only when the diagonal address survives.
```

That sounds almost too simple, but the finite audit makes it concrete.  Take an
`m x m` binary matrix and build the antiword:

```text
a_i = 1 - M_ii.
```

For `m=4`, row and column weights together still leave `8769` mixed buckets for
the antiword.  The diagonal vector alone has zero mixed buckets.  So a quotient
can remember an impressive amount of scalar information and still forget the
one coordinate needed to compute the witness.

That is the constructive reading of diagonalization.  Classical existence says
"there is a row-escaping object."  Constructive mathematics asks for the
program:

```text
given the code, produce the witness.
```

The program is tiny in Cantor's case because the diagonal bits are named.  In a
Borel setting the analogous address is a Borel code, selector, invariant
embedding, or naturality condition.  If that address is quotiented away, the
witness still feels morally present but cannot be used.

This plugs straight into embedded maximality.  HYP-2242 said a maximum is not a
unary property; it is:

```text
maximal(object, ambient embedding, allowed extensions).
```

Now the same thing is true of antidiagonal witnesses:

```text
witness(object, diagonal address, allowed quotient/extension).
```

A local maximum leaks when a new cut is admitted.  An antidiagonal leaks when
the quotient forgets the coordinate that builds it.  Both are failures of
embedding, not failures of imagination.

The regressive toy was useful precisely because it was too small.  On finite
intervals, the map

```text
f(t) = min(t) - 1
```

avoids every shifted collision in the tested windows.  That predecessor value
is the whole branch.  If we forbid the predecessor boundary, the toy collapses
immediately.  So the endpoint is not incidental; it is the escape route.  This
is the finite shadow of why Friedman's no-endpoint/critical-order setting is
not just a decorative generalization.  The ambient embedding determines which
diagonal escapes remain legal.

For repo work, this turns constructivism into a practical test:

```text
Can this invariant construct the next witness without hidden labels?
```

For A000568, HYP-2246 has a beautiful half-filter purity result.  S673 says the
next theorem should be a selection theorem: given parent deck data and a pure
half-filter bucket, select the child representative without full hidden
canonicalization.

For LRC14, the visible `Res_27` data is row/column weight.  The owner-private
deletion bit is closer to the diagonal address.  The next target is not just
classifying local carries after the fact; it is constructing the strict-tax or
lonely-time witness from the retained owner/cut address.

For unit distance, edge counts and direction counts are also row/column
weights.  The diagonal address is the unit-spine owner or point-deletion
frontier owner.  Without it, the construction can look large enough while still
not telling us where the mandatory spine witness lives.

The phrase "tangible incompleteness" now has a specific repo meaning for me:

```text
a finite-looking object whose witness exists locally,
but whose uniform construction requires a stronger embedded recursion law
than the visible quotient supplies.
```

This is why the incoming HYP-2247 fourth-face split matters.  The fraction face
is the address.  The recursion face is the law saying the address survives
extension and can keep producing witnesses.

Working slogan:

```text
Existence is a shadow; address plus transition is a proof object.
```

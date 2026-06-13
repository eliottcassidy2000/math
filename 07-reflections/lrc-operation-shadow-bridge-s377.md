# LRC Operation Shadow Bridge

codex-2026-05-31 S377

The natural-number operation graph story and the Lonely Runner recursion story
now look like two faces of one projection problem.

The old "incomplete tournament" phrase was almost right, but it has to be
used carefully.  For addition, the simple one-input shadow of

```text
x + y = z
```

is not incomplete at all once the complementary summand is forgotten.  It is
the complete transitive order:

```text
x -> z  iff  x < z.
```

The missing structure is in the fiber label `y=z-x`.  Addition is dense, but
the density hides the operation fiber.

Multiplication behaves differently.  The simple one-input shadow of

```text
x * y = z
```

is already sparse:

```text
x -> z  iff  x divides z.
```

The divisor DAG is the residue left after forgetting the second input.  This
is the real incomplete-tournament analogue.

## Product-Sum As The Interface

The family

```text
x_1 + ... + x_k = x_1 * ... * x_k
```

is the collision interface between dense additive completion and sparse
multiplicative residue.  The S365/S366 normal form says: delete the `1`s from
a product-sum tuple and call the remaining core `C`.  Then

```text
prod(C) - sum(C) = number of deleted 1s.
```

So `1`s are not decorative.  They are additive slack that repairs
multiplicative excess.

That sentence now has an LRC translation: unit endpoints are not decorative
either.  They are boundary mass.  The initial segment `{1,...,n-1}` leaves
exactly the unit residues `a/n` as safe boundary witnesses, and THM-360 says a
unit endpoint can only be strictly protected by an `n`-divisible speed.  In
both problems, units are where additive slack is paid.

## The LRC Quotient Ladder

The new computation `lrc_operation_shadow_bridge_s377.py` reads the stored S376
recursive LRC atlas and compares it with operation-shadow metrics.  The sharp
finite observation is:

```text
For every composite n<=18 with a stored quotient ladder,
the best ladder divisor equals lpd(n), the largest proper divisor of n.
```

The list is:

```text
n=4,6,8,9,10,12,14,15,16,18.
```

This is exactly the divisor-DAG signal one would expect from the multiplication
shadow.  The first quotient channel is not an arbitrary proper divisor; it is
the one closest to the additive completion, namely `n / spf(n)`.

The pressure ranking also clarifies the fourteen-runner case:

```text
endpoint_defect / gap_ratio:
n=18 > n=16 > n=14 > n=15 > n=12 > n=10.
```

So `n=14` is not just weird.  It is the first really loud member of a family
where the largest proper divisor gives a tiny global gap but leaves many
endpoints exposed.

## A New Invariant Vector

For LRC, the tournament-style feature vector should now include an operation
shadow layer:

```text
(
  C(n),
  phi(n)/(n-1),
  multiplication-shadow density on [n],
  largest-proper-divisor ladder gap,
  exposed endpoint count,
  scalar-puncture moat,
  product-sum minimum at arity n
)
```

The multiplication-shadow density decreases as the additive completion grows.
In the small atlas its correlation with `C(n)` is about `-0.673`, which is
only descriptive, but the sign is right: the more complete the additive
ambient space becomes, the smaller the visible multiplicative residue looks.

## Next Proof Shape

The next proof object should not be a scalar monotonicity statement.  It should
be a divisor-indexed endpoint-transfer matrix:

```text
unit endpoints
-> lpd(n) quotient endpoints
-> lower divisor quotient endpoints
-> residual nondivisor endpoints
```

This is the LRC analogue of the tournament program's projection-defect and
endpoint-transfer matrices.  The theorem to try first:

```text
Among one-divisor quotient ladders at composite n,
the lpd(n) ladder minimizes gap ratio.
```

If false, the exceptions should be valuable: they ought to occur exactly when
a multi-factor product-sum packing beats the binary divisor layer.

# Quadratic pi/e Carrier S636

S635 said that `S=e+pi` and `P=e*pi` are trace and norm.  S636 adds the sheet:

```text
D = e - pi,
S^2 - 4P = D^2.
```

That is the missing branch datum.  Trace and norm define the quadratic; the
discriminant tells us how the two roots split.

This gives a small but real theorem:

```text
At least two of e+pi, e*pi, e-pi are transcendental.
```

Proof: any two of them reconstruct `e` and `pi` algebraically.  `S,D` do it
linearly.  `S,P` do it by `T^2-S*T+P`.  `D,P` do it by `T^2-D*T-P`, whose roots
are `e` and `-pi`.  Since `e` and `pi` are transcendental, no pair of shadows
can both be algebraic.

This is not the dream theorem.  It still does not tell us whether `e+pi` is
irrational, whether `e*pi` is irrational, or whether `e-pi` is irrational.  But
it makes the obstruction geometry more rigid: the algebraic possibility can
occupy at most one of the three coordinate shadows.

The Schanuel route is the clean conditional completion.  Schanuel applied to
`1` and `i*pi` implies algebraic independence of `e` and `pi`.  Once that is
known, any nonconstant algebraic-coefficient polynomial in `e` and `pi` is
transcendental.  So Schanuel implies `S`, `P`, and `D` are all transcendental,
and even that the coordinate pairs `(S,P)`, `(S,D)`, `(P,D)` are algebraically
independent.

That pinpoints why the actual problem is hard.  The question is not one-variable
transcendence.  It is whether `(e,pi)` lies on one of a few very simple
algebraic curves:

```text
X + Y = algebraic,
XY = algebraic,
X - Y = algebraic.
```

Known one-variable machines do not see enough.  Lindemann-Weierstrass gives
`exp(alpha)` transcendental for nonzero algebraic `alpha`; Gelfond-Schneider
handles algebraic bases with algebraic irrational exponents.  But if
`S=e+pi` were algebraic, the lifted equation is only

```text
exp(S) = exp(e)*exp(pi),
```

and that asks for algebraic independence among exponentials of transcendental
inputs.  If `P=e*pi` were algebraic, then `pi=P/e`, and again the exponent is
not algebraic.  The easy lifts all ask for Schanuel-like strength.

The PSLQ part of S636 is just a height sieve, but it is a useful sanity check.
At 100 digits it finds no small scalar relation for `S`, `P`, or `D`, and no
small pair relation for `(S,P)`, `(S,D)`, or `(P,D)`.  It does find the intended
degree-2 sheet relation:

```text
S^2 - D^2 - 4P = 0.
```

So the computational witness agrees with the geometry: the visible relation is
not an accidental algebraic dependency among the shadows; it is the structural
discriminant equation.

This fits the older `transcendental_bases_s116h.py` thread.  That thread asks
whether `ln(pi)` is rational.  If it were, then `e` and `pi` would be
multiplicatively commensurable:

```text
e^p = pi^q.
```

S636 is the additive/quadratic sibling.  Instead of asking whether the log
ratio lands in `Q`, it asks whether trace or norm descends to `Qbar`.  Both are
commensurability questions, one multiplicative and one quadratic.

The LRC translation is unexpectedly good.  A finite reset period is a trace-like
clock: it tells us the orbit returns to the same visible state.  A relation
lattice is a norm-like shell: it tells us which hidden resonances are conserved.
But the lonely predicate is observer-coupled, so we also need a discriminant
sheet: who is on which side of the safe box, which endpoint owns a pinch, which
branch is source or sink.

That is the algorithmic lesson I would carry forward:

```text
Do not store only trace/norm quotients.
Store the discriminant sheet that names the branch.
```

For LRC this means reset-period ledgers should carry observer/source/sink labels
and short-circuit branch ownership.  For unit distance it means edge count must
carry spine/bulk and endpoint-compatible ears.  For tournaments it means rooted
perspectives and anti-cosets, not raw complement merges.

The surprising thing is how small the theorem is and how useful the shape is.
`T^2-S*T+P` says: two shadows can reconstruct a hidden pair, but the sheet still
matters.  That is almost exactly the observer-blind versus observer-coupled
dichotomy the repo has been circling.

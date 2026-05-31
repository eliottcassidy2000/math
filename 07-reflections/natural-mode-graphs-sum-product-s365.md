# Natural Mode Graphs and the Product-Sum Interface

codex-2026-05-31 S365

The prompt pointed back to an old object: an oriented graph on natural-number
modes where a decomposition `x+y=z` creates arrows `x -> z` and `y -> z`.
The repo had already called this the summand graph, with one important
convention: binary parents are distinct.  The multiplicative analogue uses
`x*y=z`, again with distinct parents at least `2`.

The useful new move was to put both graphs next to the equations

```text
x+y = x*y
x+y+z = x*y*z
x_1+...+x_k = x_1*...*x_k
```

and ask where additive coalescence and multiplicative coalescence actually
agree.

## Archaeology

The older summand-graph work had two core facts:

1. Additively, `{2,3}` generates every positive integer except `{1,4,6}`.
2. Multiplicatively, the product graph is much sparser, because every prime and
   every prime square is an atom under distinct-factor multiplication.

That already smelled like an additive/multiplicative asymmetry, but the
product-sum equations explain the small exceptional set better.

## The Low Resonances

Binary equality gives

```text
x+y=x*y  iff  (x-1)(y-1)=1,
```

so the only positive solution is `(2,2)`.  Under the old distinct-parent rule,
this is hidden.  Its target is `4`.

Ternary equality has the unique distinct positive solution:

```text
1+2+3 = 1*2*3 = 6.
```

So the additive complement `{1,4,6}` is not arbitrary:

```text
1 = source / multiplicative identity
4 = hidden diagonal binary resonance
6 = first visible distinct ternary resonance
```

This is much cleaner than treating `{1,4,6}` as just a weird reachability
artifact.

## Defect Normal Form

The session added THM-361.  Strip all `1`s from a product-sum tuple and call
the remaining core `C`.  Then

```text
product(C) - sum(C) = number of stripped 1s.
```

Conversely, any core with nonnegative defect becomes a product-sum tuple by
adding that many `1`s.

This turns the whole family into a slack law.  The multiplicative graph creates
excess; `1`s repair that excess additively without changing the product.

The first few exact arity rows are:

```text
k=2: (2,2)
k=3: (1,2,3)
k=4: (1,1,2,4)
k=5: (1,1,1,2,5), (1,1,1,3,3), (1,1,2,2,2)
```

The infinite two-core spine is especially simple:

```text
core=(2,n), defect=n-2
(1^(n-2), 2, n) has sum=product=2n.
```

## Relation To Tournaments

The repo keeps finding `6`, `7`, `21`, and `42`.  This lens separates their
roles:

- `6` is the unique distinct product-sum resonance and also the first ternary
  summand threshold.
- `7` is not a product-sum equality; it is the next Lucas/summand escape and
  the first forbidden `H` value.
- `21` is triangular `T_6`, the six-part additive threshold and the second
  forbidden `H` value.
- `42=2*3*7` belongs to the hyperbolic/self-product thread, where the product
  of angle denominators is one above the pairwise-product sum.

The new suggestion is: do not compare special numbers by value alone.  Compare
their defect layer.  A number can be special because it is a source, a diagonal
collision, a first distinct resonance, a triangular arity threshold, or a
one-off hyperbolic self-product.

## Files

- `04-computation/natural_mode_graph_s365.py`
- `05-knowledge/results/natural_mode_graph_s365.out`
- `01-canon/theorems/THM-361-product-sum-defect-normal-form.md`
- `05-knowledge/hypotheses/HYP-1821-natural-mode-defect-interface.md`

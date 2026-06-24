# Farey, Fibonacci, Goldbach, Polygonal Numbers, And Zeckendorf

The new connection is that the Fibonacci row you named,

```text
1,
1,
1+1,
1+2,
1+3+1,
1+4+3,
1+5+6+1,
...
```

is not just "the Fibonacci sequence with decorations."  It is the sparse-carry
normal form

```text
F_{n+1} = sum_k binom(n-k,k).
```

The entry `binom(n-k,k)` counts supports of size `k` with no adjacent carries.
That is the local automaton behind Zeckendorf: coverage plus a confluent rule
that destroys ambiguity.

That makes the Farey bridge cleaner.  On the golden Stern-Brocot spine

```text
F_i / F_{i+1},
```

the additive payload of a fraction is literally the next Fibonacci number:

```text
p + q = F_{i+2}.
```

So the operation "add numerator and denominator" is the same sparse-carry clock
when the continued fraction is all ones.  The product payload is a different
clock:

```text
p*q = |E(K_{p,q})|.
```

That is the old HYP-2932/HYP-2934 `K_{p,q}` area/incidence carrier.  It sees the
bipartite packet, not the binding denominator.  The power payloads are a third
clock: `p^q` and `q^p` are not useful proof denominators, but they are excellent
stress tests for whether a quotient has forgotten order.  The scout saw the
power winner flip between `2/3` and `3/5`, which is exactly the kind of signal a
sum/product scalar would erase.

This also clarifies the additive-basis triangle:

```text
Goldbach: many representations, local residue/singular-series correction.
Ternary Goldbach: one extra prime creates enough smoothing to prove coverage.
Fermat polygonal: bounded arity is the resource; k summands pay k-gonal debts.
Zeckendorf: no-adjacent normal form is the resource; uniqueness replaces mass.
```

The punchline for the repo is not that all four theories are secretly the same.
They are four different economies for a representation hypergraph.  If a future
LRC proof uses a sequence shadow, it must declare its economy:

```text
smoothing,
bounded arity,
normal form,
Farey address.
```

Then the controlled-kernel rule can ask the right question.  What coordinate is
this quotient allowed to forget?  For Goldbach, not the residue/singular-series
side channel.  For Fermat polygonal, not the arity budget.  For Zeckendorf, not
the carry automaton.  For Farey, not the separation between `q`, `p+q`, `p*q`,
and the power stress clocks.

That is the transfer back to LRC14: `q` remains the theorem-level binding scale,
`p+q` is the Stern-Brocot/Fibonacci recursion ledger, `p*q` is the Kpq incidence
ledger, and powers are anti-scalarization tests.  Any packet classifier that
uses representation counts without these clocks should be treated like a raw
scalar quotient: interesting, but not proof-safe until the missing address is
reconstructed, annihilated, exact as a cocycle, descended, or named as residual
debt.

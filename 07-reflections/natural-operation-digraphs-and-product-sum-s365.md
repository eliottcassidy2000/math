# Natural Operation Digraphs And Product-Sum Critical Pairs

**Session:** codex-2026-05-31-S365  
**Prompt:** Revisit the old incomplete-tournament idea where natural-number
nodes point into their sums, compare the sparser multiplication version, and
connect both to equations such as `X+Y=X*Y` and `X+Y+Z=X*Y*Z`.

## Repo Hooks

Several older threads were pointing at the same object from different angles.

- `natural_numbers_s116i.py` says multiplication is the natural formal-group
  operation in rapidity space, while addition is external and harder.
- `understand_naturals_s90ay.py` makes this sharper: multiplication is
  ordinary addition of rapidities, but addition of naturals is log-sum-exp.
- `sum_product_tournament.py` uses the Paley/interval crossover as a
  multiplicative-versus-additive phase competition.
- `non_tournament_apps_s253.py` keeps partial/incomplete tournaments as the
  natural non-complete comparison layer.
- `audit.py` had already recorded the small fact that `(1,2,3)` is the unique
  positive sorted triple with `x+y+z=xyz` in its finite check.

The new observation is that the proposed additive natural-number arrow graph is
not an incomplete tournament after forgetting the second input.

## The Two Shadows

Let `A_N` be the additive two-input operation system on `[N]`: for each
`x+y=z`, draw both `x -> z` and `y -> z`.

As a simple directed graph this collapses:

```text
x -> z in A_N  iff  exists y>=1 with x+y=z
                 iff x < z.
```

So the one-input shadow of addition is exactly the transitive tournament on
`1<2<...<N`. It is complete, acyclic, and loses almost all information about
which second input made the edge.

Let `M_N` be the multiplication analogue: for each `xy=z`, draw both
`x -> z` and `y -> z`. Its simple shadow is

```text
x -> z in M_N  iff  x divides z and x<z.
```

This is the divisibility DAG, i.e. a genuinely sparse incomplete tournament
inside the additive transitive completion. Its Hasse cover edges are prime
multiplications `x -> xp`.

The computation `natural_operation_digraphs_s365.py` records the densities:
by `N=100`, addition has all `4950` transitive comparisons, while multiplication
has only `382` if the unit `1` is included and `283` if it is removed.

## Keep The Gate Labels

The interesting object is therefore not the simple arrow graph but the labeled
two-input gate hypergraph:

```text
additive gate:        {x,y} -> x+y
multiplicative gate:  {x,y} -> xy
```

The additive one-shadow answers only "can x contribute to z?" and the answer is
yes for every larger `z`. The product one-shadow answers "can x divide z?" and
retains prime/factor structure. Product-sum equations are the critical pairs
where the two gate systems reach the same endpoint.

For the binary equation,

```text
x+y = xy
```

we get `(x-1)(y-1)=1`, hence the only positive solution is `(2,2) -> 4`.
If distinct inputs are required, the binary collision disappears entirely.

For higher arity, the units are the hidden mechanism. Given a nondecreasing
seed `F=(a_1,...,a_m)` with all `a_i>1`, set

```text
P(F) = product a_i
S(F) = sum a_i
D(F) = P(F)-S(F).
```

If `D(F)>=0`, then adding `D(F)` copies of `1` gives

```text
D(F)*1 + S(F) = P(F),
```

so the witness

```text
1, ..., 1, a_1, ..., a_m
```

solves `x_1+...+x_k=x_1...x_k` with arity

```text
k = m + D(F).
```

This is exactly the product-sum-number mechanism, reinterpreted as an
additive/multiplicative gate critical pair.

## Divisor Layer

The two-nonunit stratum has a clean closed form. With `r=k-2` ones,

```text
r+a+b = ab
```

is equivalent to

```text
(a-1)(b-1) = k-1.
```

Thus two-factor product-sum witnesses at arity `k` are classified by the
divisors of `k-1`. This is a nice bridge: the additive arity defect is governed
by a multiplicative divisor lattice.

The first point where a multi-factor seed beats every two-factor seed is
already `k=5`:

```text
1+1+2+2+2 = 8 = 1*1*2*2*2
```

The best two-factor solutions at `k=5` have products `9` and `10`.

That smells like the OCF pattern: a single local obstruction is not always the
best packing. Several small compatible multiplicative atoms can beat one larger
binary atom.

## A Useful Reframing

The natural-number operation graph should be treated as a two-colored operation
complex:

- vertices are natural numbers;
- additive 2-cells are `x,y -> x+y`;
- multiplicative 2-cells are `x,y -> xy`;
- the additive 1-shadow is the transitive tournament;
- the multiplicative 1-shadow is the divisibility DAG;
- product-sum witnesses are two-colored critical pairs.

This gives a concrete version of the repo's repeated "additive versus
multiplicative" theme:

```text
multiplication path:  factors -> product       (prime-generated DAG)
addition path:        summands -> sum          (transitive completion)
critical pair:        same endpoint
defect:               number of unit paddings needed
```

In rapidity coordinates, multiplication is linear but addition is log-sum-exp.
The equation `sum=product` is therefore the place where a log-sum-exp fold lands
exactly on a linear rapidity fold. That is a clean toy model for the general
number-theoretic pain point: additive questions about multiplicative structure.

## Next Directions

1. Build the two-colored operation complex and compute small homology/critical
   pairs after quotienting by commutativity.
2. Treat the product-sum number sequence as a minimal endpoint problem in the
   additive completion of the divisibility DAG.
3. Compare the defect `D(F)=P-S` to OCF fugacity: small compatible factors can
   win over a single larger binary factor exactly as packings can beat raw
   abundance.
4. Ask whether divisor-layer formulas for two factors have higher-factor
   analogues using multiplicative partitions of `k-c`.

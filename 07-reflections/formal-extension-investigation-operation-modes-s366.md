# Formal Extension Investigation: Operation Modes

**Session:** codex-2026-05-31-S366  
**Prompt:** Spend a long overnight-style session cycling between formalization,
extension, and investigation, using the other two modes as inspiration for the
third.

This session continued the S365 natural-operation/product-sum thread and made
the requested mode loop explicit.

## Loop 1: Formalization To Extension

The clean formal target was the ordinary one-input shadow of the two operations
on positive natural numbers.

Lean module:

```text
04-computation/lean/TournamentH7/TournamentH7/NaturalOperationDigraphs.lean
```

The axiom-free core now proves:

- `AddShadow x z <-> x < z`.
- Unit-allowing multiplication shadow is divisibility on positive targets.
- Nonunit multiplication shadow is proper divisibility on positive sources.
- Shifted binary collision:
  `(a+1)+(b+1)=(a+1)(b+1) <-> a*b=1`.
- Two-factor product-sum layer:
  `r+(a+1)+(b+1)=(a+1)(b+1) <-> a*b=r+1`.
- The universal binary witness:
  `r` ones plus factors `2` and `r+2` solve `sum=product`, giving `m(k)<=2k`.
- Seed-defect padding:
  if `sum(F)<=prod(F)`, then `(prod(F)-sum(F))+sum(F)=prod(F)`.

That formalization corrected the level of abstraction. Addition, as an ordinary
edge shadow, is too complete: it is just the transitive tournament. The product
shadow is sparse because it is divisibility. Therefore the extension should not
study ordinary additive edges; it should study critical pairs of the labeled
operation gates.

## Loop 2: Extension To Investigation

Computation:

```text
04-computation/natural_operation_modes_s366.py
05-knowledge/results/natural_operation_modes_s366.out
```

The extension enumerates all potentially minimal product-sum witnesses through
arity `k<=120`. The key pruning comes from the formal binary witness:

```text
m(k) <= 2k,
```

so products above `2*MAX_K` cannot be minimal for the searched range.

The product shadow remains very sparse inside the additive completion. On
`[1000]`, addition has all `499500` transitive comparisons; multiplication has
only `6069` unit-including product-shadow edges, density `0.01215`.

For product-sum witnesses, the binary divisor layer is not the main story after
small arity. Through `k<=120`, strict multi-factor wins over the best binary
witness occur in `96` of the `119` arities.

Early strict wins:

```text
k=5:  seed (2,2,2),      P=8,  beats binary P=9
k=8:  seed (2,2,3),      P=12, beats binary P=16
k=12: seed (2,2,2,2),    P=16, beats binary P=24
k=21: seed (3,3,3),      P=27, beats binary P=30
k=27: seed (2,2,2,2,2),  P=32, beats binary P=42
```

This echoes a repo-native OCF lesson: raw abundance of one large local object is
often beaten by a compatible packing of smaller atoms.

## Loop 3: Investigation Back To Formalization

The computation suggests three formal targets.

First, the binary layer is solved:

```text
(a-1)(b-1)=k-1.
```

Lean already proves this in shifted form. This should become the base case of a
general "factor-packing beats binary" theory.

Second, the bound `m(k)<=2k` is not merely a computational convenience; it is a
formal normalization theorem. It makes finite minimal product-sum search
possible by bounding the endpoint range.

Third, the real proof target is a packing inequality for seeds. For a seed
`F`, define

```text
P = prod(F)
S = sum(F)
D = P-S
k = |F|+D.
```

The optimization problem is:

```text
minimize P subject to |F| + P - S = k.
```

Equivalently:

```text
P - S + |F| = k.
```

The score `S-|F|` is the "additive discount" a multiplicative seed earns. Small
factors have high discount per product cost, which explains why strings of
`2`s and `3`s dominate early minima.

## New Mental Model

The operation complex has three projections:

```text
two-input gates       -> ordinary one-shadow -> endpoint minimization
{x,y}->x+y            -> total order         -> all additive paths available
{x,y}->xy             -> divisibility DAG    -> sparse multiplicative paths
sum=product witness   -> critical pair       -> minimal common endpoint
```

The additive shadow supplies a full transitive completion. The divisibility DAG
supplies the sparse multiplicative skeleton. Product-sum numbers measure the
least endpoint at which a multiplicative seed can be padded by units to land in
the additive completion with the same value.

That makes product-sum numbers a small exact laboratory for the broad repo
theme:

```text
additive completion + sparse multiplicative skeleton + compatibility packing
```

The same grammar appears in Paley/interval, OCF packings, endpoint-protection
graphs, and quotient bucket transport.

## Next Moves

1. Formalize a `Seed` structure with `prod`, `sum`, `defect`, and arity.
2. Prove the finite search bound `m(k)<=2k` in a small abstract definition of
   minimal product-sum endpoint.
3. Search for a dominance theorem: replacing a large factor by several small
   factors improves `P` when it increases the discount `S-|F|` enough.
4. Build the two-colored operation complex and compute its small critical
   cycles after quotienting by commutativity.

# Pair-First Twin Prime Lens S502

The user prompt was blunt in the best possible way: think of everything in
terms of pairs.  Two nodes in a tournament matter, but most importantly the
connection between them.  Twin primes.

That suggests a basic ontology:

```text
pair-cell = endpoints + relation/label + fiber coordinates.
```

The endpoints are carriers.  The connection is where the structure lives.

## Tournament Pair-Cells

A tournament on `n` vertices has exactly one bit for each unordered pair:

```text
{i,j} -> i->j or j->i.
```

The repo often projects those bits to vertex scores, Hamiltonian paths, odd
cycle conflicts, good cuts, pressure DAGs, or TDA features.  The pair-first
warning is that these projections can erase too much.  Before projecting, the
edge already has useful coordinates:

```text
range = j-i
midpoint = i+j
orientation sign = s_e
flip impact
local cycle participation
```

The S502 script repairs the missing `edge-variable.md` note and makes `s_e`
explicit again as the primitive tournament coordinate.

## Prime Pair-Cells

For primes the same idea is:

```text
{a,b} -> sum=a+b, gap=b-a, midpoint=(a+b)/2, residue data.
```

The coordinate change is invertible for same-parity endpoints:

```text
a = (sum-gap)/2
b = (sum+gap)/2
```

So Goldbach and twin primes are not separate mental universes.

```text
twin primes = fixed gap row, gap=2
Goldbach    = fixed sum column, sum=N
```

The S502 computation shows that fixed gap rows already behave like labelled
arithmetic edges.  Up to `100000`, gap `6` has almost exactly twice as many
prime-pair edges as gap `2`, matching the singular-series row boost.  Gaps
`30` and `60` have about `2.7` times as many.  The connection label is carrying
the local obstruction product.

This is the twin-prime reframing I like:

```text
The twin prime conjecture is not primarily about prime nodes.
It is about infinite recurrence of the smallest nontrivial surviving edge row.
```

The nodes are prime because the edge needs two prime endpoints, but the
conjectural persistence is a statement about the edge type `gap=2`.

## Goldbach As The Perpendicular Slice

Goldbach fixes the other coordinate.  For a fixed even `N`, the edge fiber is

```text
{p,N-p}: p and N-p prime.
```

The S502 script records:

```text
N=100      6 pairs
N=1000    28 pairs
N=100000  810 pairs
```

The same edge can be read as belonging to one gap row and one sum column.  This
is exactly the kind of two-coordinate board the repo likes: tournament ranges
and midpoints on one side, prime-pair sums and gaps on the other.

## Zeckendorf Carry Debt Is Pair Data

S502 also tests a repo-native extra label on prime-pair edges:

```text
Zeckendorf carry debt of endpoint supports relative to the sum support.
```

Some Goldbach edges have zero carry debt:

```text
100    = 11 + 89
5000   = 673 + 4327
50000  = 61 + 49939
100000 = 24971 + 75029
```

That does not prove anything about Goldbach, but it is a useful new metric.  It
asks how hard the endpoint normal forms have to work to become the target normal
form.  Again, the quantity belongs to the connection, not to either endpoint.

## SC Blowup As Twin-Gaining

The old SC blowup note becomes more vivid here.  Each vertex becomes a pair
`{v_0,v_1}`.  Each old pair `{u,v}` becomes four edge-cells:

```text
same lane follows T
cross lane follows T^op
```

This duplicates carriers while preserving connection memory.  At the vertex
level, score variation is flattened into the universal sequence

```text
(n-1)^n, n^n.
```

At the pair level, nothing is forgotten: lane/cross labels still encode `T` and
`T^op`.  That is the tournament analogue of the twin-prime lesson: the twin
relation is a structural edge type, not a mere node property.

## LRC Direction

For LRC, the natural next experiment is to lift pressure data one level down:

```text
node = runner-pair, or endpoint-runner pair
edge label = who changes whose safe distance, chord distance, or pressure debt
```

The current pressure DAGs have been useful, but they are probably projections.
If the real structure is pair-first, then a pressure edge should remember the
pair-cell that generated it.  This might split acyclic pressure searches into
smaller peelable layers, or reveal the first genuine pressure cycle as a
cycle of pair-cells rather than a cycle of runners.

## New Questions

1. Can tournament TDA store an edge table first and compute score sequence,
   cycle counts, good cuts, and pressure features as views?
2. Can prime-pair TDA compare the topology of fixed-gap rows and fixed-sum
   columns in the same `(sum,gap)` grid?
3. Do zero-carry Goldbach edges cluster in predictable residue or singular
   factor classes?
4. Does SC blowup give a general "twin operation" on pair-labelled structures:
   duplicate carriers, flatten vertex statistics, retain connection data?
5. Can LRC pressure DAGs be rebuilt on endpoint-runner pair-cells and then
   projected down only after the peel order is known?

HYP-1965 records this as the pair-cell ontology.

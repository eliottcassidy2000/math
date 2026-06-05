# LRC Contact Invariants Spec

This fixes a reproducible convention for the requested `D(A)` table.

Let `A={v_1<...<v_k}` be primitive speeds and `delta=1/(k+1)` unless stated.
Define `M(A)=max_t min_i ||v_i t||` on `t in [0,1]`.

## Exact Candidates

The checker evaluates all rational candidate times

```text
0, 1,
(2m+1)/(2v_i),
m/(v_i+v_j),
m/|v_i-v_j|
```

for legal integer ranges, then scores each by exact rational arithmetic.

## `D(A)`

`D(A)` is the sorted set of reduced denominators of all exact maximin times:

```text
D(A) = { den(t) : min_i ||v_i t|| = M(A) }.
```

This is intentionally a contact denominator invariant, not all denominators
appearing in the candidate scan.

## Endpoint Activity Sets

For a contact time `t` and score `M(A)`, the activity set is

```text
Act(t) = { (i, side_i) : ||v_i t|| = M(A) }.
```

`side_i = -` if `{v_i t}=M(A)` and `side_i = +` if `{v_i t}=1-M(A)`.
These are the active lower/upper endpoints of the forbidden arcs.

## Denominator Classes

For each `q in D(A)`, the denominator class is the set of distinct activity
patterns among contacts with denominator `q`:

```text
C_q(A) = { Act(t) : den(t)=q and t is maximin }.
```

Two candidates are in the same class iff the runner indices and endpoint signs
match exactly.

## F2 Incidence Rank

Build a binary matrix over `F_2` with one row per contact time and one column per
runner:

```text
I[t,i] = 1 iff i is active at t.
```

The F2 incidence rank is `rank_F2(I)`.

## Cyclomatic Rank

Build the bipartite contact graph with runner vertices and contact-time vertices,
and edges `i--t` iff `i` is active at `t`. Its cyclomatic rank is

```text
beta_1 = |E| - |V| + components.
```

This detects independent endpoint-contact cycles left after forgetting signs.

## SCC And Tournament Entropy

At the first maximin witness, build the half-turn phase tournament on runners:
orient `i -> j` iff `{v_j t - v_i t} < 1/2`, with index tie-break at exactly
`1/2`. Report SCC sizes and `log2 H`, where `H` is the exact Hamiltonian-path
count of this tournament.

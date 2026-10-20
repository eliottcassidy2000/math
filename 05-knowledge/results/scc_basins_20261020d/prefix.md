# A four-black SCC interchange with an unrestricted prefix

**PROVED, all orders:** if both endpoints are black-active, the ordered
joins `A -> 1 -> C3` and `A -> C3 -> 1` are joined by four black path
reversals, for any tournament A of order at least two. For empty A the
two endpoints are already converse-equivalent. For a singleton A the
second endpoint is pure blue, so the activity hypothesis fails.
All-order black-core connectivity remains **OPEN**.

## 1. Inherited mechanism and the new choice

The four paths below are inherited from
[the guarded SCC interchange, section 5](../black_descent_20261020c/scc.md).
That theorem required A to begin with two dominating singletons. Its
literal arc ledger did not need the guard; its blackness proof did.
The missing coordinate is a suitable choice of Hamiltonian path inside A.
The canonical hostile is `1 -> C3 -> 1`, a pure-blue target which cannot
be reached by a black edge. The least-used sidecar is the set of vertices
which can be last in a Hamiltonian path of A.

Here `T^P` reverses the arcs of a directed Hamiltonian path P. It represents
the same converse-merged line as the nonpath-complement operation.
A path is blue if its positional reversal is an anti-automorphism.
In particular, a blue reversal pairs all vertex scores to `n-1`.

## 2. The path ledger

Choose a Hamiltonian path `p=(u_0,...,u_(m-1))` of A, with `m>=2`.
Write R for its last vertex. After A put x, followed by the cycle
`a -> b -> c -> a`, with x dominating that cycle. Successively reverse

```text
P1 = (p,       x,c,a,b),
P2 = (rev(p),  c,x,b,a),
P3 = (x,a,b,c, p),
P4 = (c,b,a,x, rev(p)).
```

These are directed paths for every choice of A and p. The internal path
arcs of A reverse four times; other arcs of A never reverse. The edges
`Rx` and `u_0 c` each reverse twice. Among the suffix edges only `ac`
and `xb` reverse an odd number of times. Thus the final tournament is
literally A followed by `x -> a -> b -> x`, followed by the sink c.
All arcs of A and between A and the four-vertex suffix are restored.

The problem is to choose p so that all four presentations are black.

## 3. The second move is always black

Immediately before P2 the vertex a has score one and unique outneighbor
c. If P2 were blue, its positional reversal r would send a to R, forcing
R to have indegree one. Direct endpoint scores give

```text
d_A(R) = m-2.
```

R originally loses to its predecessor on p. At that extremal score this
is its only loss inside A; P1 reverses that edge. Hence R now dominates
every other A vertex, and its unique inneighbor is x. Reversal of the
arc `a -> c` then requires `r(c)=x`.

In the displayed P2, c is at position m, its reversed position is 3,
and x is at position m+1 (zero-based positions). This forces m=2.
For m=2 the second and penultimate vertex scores are 4 and 2, whose sum
6 differs from `n-1=5`. This also rules out blue. Thus P2 is always black.

## 4. The first and fourth moves have exact SCC tests

Ordered strong components must be paired in reverse order by an
anti-automorphism. Since the singleton and C3 each have just one
Hamiltonian-path orbit, this gives the following exact tests:

* P1 can be blue only if `A=C3`, or if `A=C3 -> 1 -> Z`.
  In the second case it is blue exactly when the chosen path of Z is blue.
* P4 can be blue only if `A=1`, or if `A=1 -> C3 -> Z`.
  In the second case it is blue exactly when the chosen path of Z is blue.

For P4 this can be checked in the final tournament: reverse(P4) is its
Hamiltonian path, and reversing a path preserves its blue/black status.
The tests include empty Z. For example, the P1 test pairs the final C3
with the first SCC of A, then pairs the preceding singleton; only the
middle Z can remain. If A consists of just that initial C3, the two
triangles pair across the central singleton.

The exceptional prefixes `A=C3` and `A=1` make respectively the source
and target pure blue. They are excluded by the theorem's hypothesis.
The two remaining forms of A are mutually exclusive: their first SCCs
have different orders.

If one such form occurs, its corresponding endpoint is pure blue
exactly when Z is pure blue. Activity therefore supplies a black path
of Z. Choose it and any paths in the preceding rigid SCCs. This makes
P1 and P4 both black. Moreover its last vertex satisfies

```text
d_A(R) <= |Z|-1 = m-5.                                (1)
```

If neither form occurs, P1 and P4 are black for every path p of A.

## 5. Choosing the final vertex makes the third move black

Immediately before P3 the scores of its two endpoints are

```text
d(x)=3,             d(R)=d_A(R)+3.
```

Consequently P3 is black whenever

```text
d_A(R) != m-3.                                       (2)
```

This already handles (1). In all other cases choose R in the last SCC
of A with score different from m-3, if possible. Every vertex in that
last SCC can end a Hamiltonian path of A: for a nontrivial strong SCC,
cut a directed Hamiltonian cycle just after that vertex, and precede it
by Hamiltonian paths in the earlier SCCs. A singleton is immediate.

It remains to consider when every vertex of that SCC has score m-3.
If its size is t, its vertices have no outneighbors in earlier SCCs,
so averaging their internal scores gives

```text
(t-1)/2 = m-3,       t=2m-5 <= m.
```

A tournament SCC has size one or at least three, and this regular SCC
has odd size. The only possibilities for m>=2 are:

| m | Last SCC | Entire A |
|---:|---|---|
|3|singleton|transitive tournament on three vertices|
|4|C3|`1 -> C3`|
|5|regular five-vertex tournament|R5|

The middle case has pure-blue target `1 -> C3 -> C3 -> 1` and is excluded.
For transitive A of order three, P3's second and penultimate scores are
2 and 5, summing to 7 rather than `n-1=6`. For regular A of order five,
P3's middle vertex u_0 has score 5, whereas a blue reversal on nine
vertices would require its fixed vertex to have score 4.
Thus P3 is black in both remaining cases.

Together these arguments select a single p for which all four moves
are black. They prove the theorem without a bound on |A|.

## 6. Exact tests, boundary cases, and what is still missing

Run from the repository root:

```text
python3 05-knowledge/results/scc_basins_20261020d/prefix_macro.py
python3 -O 05-knowledge/results/scc_basins_20261020d/prefix_macro.py
```

The script first fixes p and enumerates every orientation of its remaining
arcs for prefix orders one through seven (33,868 presentations in total).
For each it checks all path arcs, the literal final tournament, the exact
P1/P4 SCC tests, the universal P2 claim, and the sufficient P3 score test.
It then uses every ordinary tournament class of orders three through eight
from the committed census, searching its Hamiltonian paths for a witness.
Failure is accepted only when the source or target is certified pure blue
by **all** prefix paths; rotating the suffix C3 is an automorphism.

| Prefix order | Ordinary classes | Four-black certificates | Excluded pure endpoint |
|---:|---:|---:|---:|
|3|2|1|1|
|4|4|2|2|
|5|12|10|2|
|6|56|54|2|
|7|456|452|4|
|8|6,880|6,878|2|

The [receipt](prefix_macro.json) retains difficult witnesses and all
excluded classes. The fixed-path audit includes deliberate blue failures;
it does not merely test a favorable path.

Taking converse yields the unrestricted-suffix version at the left end
of an SCC word. The theorem does **not** yet swap arbitrary interior
letters with both a nonempty prefix and suffix. Nor does it eliminate an
arbitrary strong SCC letter. These are genuine remaining normalization
obligations. The [golden word codec](arithmetic.md) maps this certified move
to the terminal rewrite `code(A)121 <-> code(A)211`, subject exactly to
activity of the two endpoints, and retains the four actual tournaments
visited along the route.

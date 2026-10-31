# Canonical weighted-cut law for tournament vertex extensions

**Status: PROVED here, using Rédei's classical oddness theorem; finite exact
verification is recorded separately. No novelty or global spectrum-completeness
claim.** Session 124, 2026-10-31.

## The theorem

Let `T` be a finite tournament on `V`. Count directed Hamiltonian paths as
ordered vertex sequences, and write `H(T)` for their number; set `H(empty)=1`.
For a subset `S` of `V`, form `T+x(S)` by giving the new vertex exactly the
out-neighbors `S`, and put `F_T(S)=H(T+x(S))`.

For distinct `u,v`, let `W_uv` count **ordered pairs of nonempty directed
paths** whose vertex sets partition `V`, the first path ending at `u` and the
second starting at `v`. Define the symmetric pair weight

```text
b_uv = (W_uv + W_vu)/2.
```

Then the following hold for every tournament, with no strong-connectivity or
bounded-order assumption:

1. Every `b_uv` is an integer, and in fact
   `b_uv >= H(T-{u,v}) >= 1`.
2. The canonical response graph is **the original tournament itself**, with
   weight `b_uv` on its arc joining `u,v`:

   ```text
   F_T(S) = H(T) + 2 sum_{u in S, v outside S, u -> v in T} b_uv.       (1)
   ```

Thus no feasibility search, auxiliary orientation, or exponentially many
inequality tests are needed. The original orientation already has exactly the
required weighted outdegrees.

Here "canonical" means specified directly by `T`, not that a directed graph
representing the same cut function must be unique.

The empty and one-vertex cases have no pair weights and are covered by (1).

## Proof

If `V` is empty, both the parent and its one-vertex extension have count one,
and (1) holds with an empty sum. For the remainder of the proof assume `n>=1`.

### 1. Deletion of the new vertex gives the quadratic polynomial

Let `s_v,t_v` be the numbers of Hamiltonian paths of `T` starting and ending
at `v`. Write `c_v=1` when `v` belongs to `S`.

If a Hamiltonian path of the extension starts at `x`, deleting `x` leaves a
Hamiltonian path starting at some `v` with `c_v=1`. If it ends at `x`, it leaves
a path ending at `u` with `c_u=0`. Otherwise deletion of `x` leaves an ordered
two-path cover, with first endpoint `u`, second start `v`, and necessarily
`c_u=0,c_v=1`. These cases are disjoint and exhaustive. Therefore

```text
F_T(c) = sum_v s_v c_v + sum_u t_u(1-c_u)
         + sum_{u != v} W_uv(1-c_u)c_v
       = H + sum_v L_v c_v - sum_{u<v}(W_uv+W_vu)c_uc_v,              (2)

L_v = s_v-t_v+sum_u W_uv.
```

### 2. Endpoint divergence determines the original orientation

For an original arc `u -> v`, concatenating the two paths counted by `W_uv`
gives a Hamiltonian path using that arc. Cutting that arc reverses the map.
Consequently `W_uv` is exactly the number of Hamiltonian paths containing
the original arc `u -> v`.

In each Hamiltonian path, outgoing-arc use minus incoming-arc use at `v` is
`+1` at its start, `-1` at its end, and zero elsewhere. Summing gives

```text
s_v-t_v = sum_{v -> u} W_vu - sum_{u -> v} W_uv.
```

Add `sum_u W_uv` to obtain the decisive cancellation:

```text
L_v = sum_{v -> u}(W_vu+W_uv) = 2 sum_{v -> u} b_vu.                (3)
```

Substituting (3) into (2) gives exactly (1), since
`c_u(1-c_v)=c_u-c_uc_v`. This establishes the orientation convention, the
factor of two, and the absence of any extra constant term.

### 3. Integrality

Reverse only the arc joining `u,v`. Neither `W_uv` nor `W_vu` changes:
their two distinguished endpoints lie in different paths, so the edge
between them is never used inside either path.

If the original arc is `u -> v`, reversal loses exactly `W_uv` Hamiltonian
paths and gains exactly `W_vu`. Thus

```text
H(T with uv reversed) - H(T) = W_vu-W_uv.
```

Rédei's oddness theorem makes the left side even. Therefore `W_uv+W_vu` is
even and `b_uv` is integral. Equivalently, the mixed Boolean second difference
in (2) is even because all four extension counts are odd.

### 4. Strong pair-deletion lower bound

Fix a Hamiltonian path `P=(y_1,...,y_m)` of `T-{u,v}`. Temporarily allow either
orientation between `u,v`, and insert either consecutive block `(u,v)` or
`(v,u)` into `P`, requiring all edges *outside* that block to point forward.
After cutting the edge inside the block, each such insertion yields an
ordered two-path cover counted by `W_uv+W_vu`.

Put `A_i=[y_i -> u]`, `B_i=[y_i -> v]`, with artificial boundary values
`A_0=B_0=1` and `A_{m+1}=B_{m+1}=0`. The number of valid block insertions is

```text
N(P) = sum_{i=0}^m [A_i(1-B_{i+1}) + B_i(1-A_{i+1})].
```

The elementary termwise identity

```text
A_i(1-B_{i+1})+B_i(1-A_{i+1})
 = A_i+B_i-A_{i+1}-B_{i+1}
   + A_{i+1}(1-B_i)+B_{i+1}(1-A_i)
```

telescopes to `N(P)=2+` a nonnegative integer. Hence each `P` yields at least
two covers. They are distinct: concatenating a resulting cover across its
distinguished endpoints and deleting `u,v` recovers `P`, the insertion
position, and the block order. This also proves distinctness across different
`P`. Therefore

```text
W_uv+W_vu >= 2 H(T-{u,v}),
```

as claimed. The same computation works for `m=0`, giving the two singleton
covers. Rédei's theorem gives `H(T-{u,v})>=1` when the remainder is nonempty;
our empty-set convention handles the other case. This completes the proof.

## Consequences with exact scope

### Equality and strong components

All pair weights are strictly positive, so `F_T(S)=H(T)` if and only if there
is no original arc from `S` to `V-S`. Such a set cannot split a strongly
connected component: a path from its inside part to its outside part would
have a crossing arc. The condensation of a tournament is totally ordered,
with all arcs pointing forward, so its zero-outgoing cuts are exactly suffix
unions of strong components.

If `T` has `k` strongly connected components, there are precisely `k+1`
Hamilton-count-preserving extension cuts, including empty and full cuts.
In particular, every nontrivial cut of a strong tournament strictly raises
the count, by at least two.

More quantitatively, if `e_T(S,V-S)` counts forward crossing arcs, then

```text
F_T(S) >= H(T) + 2 sum_{u in S,v outside S,u->v} H(T-{u,v})
         >= H(T) + 2 e_T(S,V-S).                                  (4)
```

Each crossing arc closes a directed triangle through `x`, so the first sum
in (4) has a direct triangle-through-the-new-vertex interpretation. No
disjointness among those triangles is assumed or required.

### Strict submodularity

For arbitrary subsets `A,B` of `V`, (2) gives

```text
F_T(A)+F_T(B)-F_T(A union B)-F_T(A intersection B)
 = 2 sum_{u in A-B, v in B-A} b_uv
 >= 2 |A-B| |B-A|.                                                (5)
```

Thus the response is strictly submodular on incomparable subsets. For
distinct `u,v` outside a cut `S`, its mixed second difference is exactly
`-2b_uv`, independent of `S`. In particular, its Boolean polynomial has
degree at most two and every pair interaction is strictly negative.

### Complementary cuts and binary targets

The two complementary responses satisfy

```text
F_T(S)+F_T(V-S)-2H(T)
 = 2 sum_{u in S,v outside S} b_uv.
```

For a parent with odd `H`, the targets `2H-1` and `2H+1` are obtainable in
one extension precisely when the weighted original tournament has directed
cuts of capacities `(H-1)/2` and `(H+1)/2`, respectively. Obtaining both is
therefore an exact consecutive-cut-capacity problem. It does **not** prove
the existence of either cut, a recursively closed extension bank, or that
every positive odd count other than `7,21` is attainable. Children have new
weight matrices and must be analyzed anew.

### Small gains force small deletions and bounded triangle packing

Suppose an extension has positive gain `F_T(S)-H(T)=2k`. For **every**
crossing arc `u->v`, (1) and positivity imply

```text
H(T-{u,v}) <= b_uv <= k.                                         (6)
```

In particular, a gain of two or four forces every such two-vertex deletion
to have Hamiltonian count one, since that count is positive and odd. A
tournament has count one exactly when it is transitive: the forward order
is its unique path if it is transitive; otherwise it contains a directed
triangle, whose three paths remain a lower bound after adding back vertices
by (1). For the standard fact used here, a shortest directed cycle of length
at least four has a chord giving a shorter directed cycle, so a cyclic
tournament contains a triangle. Thus a gain of two or four forces the
parent's minimum feedback vertex set to have size at most two. This is a
necessary condition, not a sufficient one.

Let `nu(T)` be the maximum number of vertex-disjoint directed triangles.
The classical odd-cycle collection formula implies `H(T)>=3^nu(T)`:
the subcollections of a fixed maximum triangle packing already contribute
`sum_j binomial(nu,j)2^j=3^nu`, with all other contributions nonnegative.
Deleting two vertices destroys at most two members of any fixed packing.
Consequently (6) gives

```text
F_T(S)-H(T) >= 2 * 3^max(nu(T)-2,0)  whenever the gain is positive,
nu(T) <= 2 + floor(log_3(k))          whenever the gain is 2k > 0. (7)
```

The packing implication is the only consequence here that additionally
uses the odd-cycle collection formula. Its source is Irving--Omar,
[*A Survey of Rédei--Berge Symmetric Functions via Matrix Algebra*,
Section 3.4, Corollary 19 in the HTML fetched 2026-10-31](https://arxiv.org/html/2412.10572v1),
attributed there to Grinberg--Stanley Theorem 1.39 and Lemma 6.5. That statement
counts paths in the complement; reversing paths identifies their number
with that in the original tournament, and loops in the complement do not
enter Hamiltonian paths. The local canon is `THM-002`,
`01-canon/theorems/THM-002-ocf.md`. Cycles are literal directed cycles modulo
rotation, not merely their vertex supports. These inequalities do not assert
anything about the zero-gain suffix cuts of a nonstrong parent.

### Exact arithmetic content

The common divisor of all extension gains is exactly

```text
gcd_{S subset V}(F_T(S)-H(T)) = 2 gcd_{u<v} b_uv.                  (8)
```

The cut formula proves divisibility by the right side. Conversely,
`F_T({u,v})-F_T({u})-F_T({v})+H(T)=-2b_uv`, so the left-side gcd divides
every `2b_uv`. For orders zero and one all gains vanish; both gcds are zero
under the convention `gcd(empty)=0`. For every positive integer `m`, all
extension counts are congruent to `H(T)` modulo `2m` if and only if `m`
divides every pair weight. If `p` divides both `H(T)` and every pair weight,
then every one-vertex extension count is divisible by `p`. No claim is made
that every number in the resulting congruence class occurs.

### Constrained minimum extension as a min-cut problem

To minimize `F_T(S)` while requiring `A subset S` and `B disjoint S`, where
`A,B` are disjoint, use the directed weighted tournament as a network. Add
a source `s` with arcs `s->a` for `a in A`, and arcs `b->t` to a sink for
`b in B`. Give all these extra arcs capacity `M=1+sum_{u->v}b_uv`.
Every feasible original cut costs at most `M-1`, so no minimum cut violates
a forced membership. Its capacity `q` gives the minimum extension count
`H+2q`.

For completeness, the integral max-flow/min-cut argument used by the code
is as follows. Begin with zero flow and repeatedly augment along a residual
source-to-sink path by its least residual capacity. All capacities and flows
stay integral. With breadth-first augmenting paths this is Edmonds-Karp;
termination also follows directly here because every augmentation raises the
integer flow value and that value is bounded by the sum of source capacities.
When no residual path remains, let `R` be the vertices residual-reachable
from `s`. Every original arc from `R` to its complement is saturated and
every original arc into `R` carries zero net backward contribution. Flow
conservation therefore makes the flow value equal the capacity of this cut.
Every flow value is at most every cut capacity by the same conservation
identity. Thus this flow and cut are optimal, and provide a certificate.

The minimum **positive** extension gain is also computable by min-cuts:
for each original arc `u->v`, force `u` into the cut and `v` outside, and take
the least resulting count. Every such cut is positive, and every positive
cut contains at least one original crossing arc, so these cases are
exhaustive. This uses at most `n(n-1)/2` flows and also handles nonstrong
parents without confusing zero-gain suffix cuts with positive gains. For
orders zero and one no positive gain exists; the code returns `None`.

### Exact symmetry classification and explicit hostile controls

For `n>=2`, the extension count `F_T(S)` depends only on the cardinality of
`S` **if and only if** `T` is regular and all pair weights equal a common
positive value `w`.

Indeed cardinality invariance makes
`F_T({u,v})-F_T({u})-F_T({v})+H(T)` independent of the distinct pair `u,v`,
so (2) makes every `b_uv` equal to `w`. It also makes every singleton count
equal; (3) then says that `2w d_T^+(v)` is independent of `v`. Since `w>0`,
all outdegrees are equal, which is regularity. Conversely, if these two
conditions hold, the directed cut has `|S|(n-|S|)/2` original arcs. Hence

```text
F_T(S) = H(T) + w |S|(n-|S|).                                    (9)
```

This follows by subtracting the `|S|(|S|-1)/2` internal arcs from the total
outdegree `|S|(n-1)/2`.

In particular, every **arc-transitive tournament of order at least three**
satisfies (9). Its automorphisms transport the two-path covers and therefore
make all weights equal. Such a tournament is vertex-transitive: some vertex
has both an incoming and an outgoing arc, and arc-transitivity makes its
orbit include every possible head and tail, hence every vertex. It is
therefore regular. The order-two exception matters: a single arc is
arc-transitive but its two singleton cuts have different responses.

If `E` is the total number of ordered two-path covers allowing one empty
path, then the two empty-path cases contribute `2H`, while the nonempty
cases contribute `sum_{u!=v}W_uv=n(n-1)w`. Consequently
`w=(E-2H)/(n(n-1))` in the uniform case. This is an identity in enumerated
statistics, not a closed formula for `H`.

The independently DP-checked explicit quadratic-residue tournaments give:

| Tournament | `H(T)` | Constant pair weight | All extension values |
|---|---:|---:|---|
| Paley 7 | 189 | 39 | 189, 423, 579, 657 |
| Paley 11 | 95095 | 14421 | `95095+14421*s*(11-s)`, `s=0,...,5` |

The finite evidence in `verify_examples.py` constructs these tournaments
explicitly and checks all weights, regularity, and direct extension DP at
every cut size (and every cut at order seven). Equation (9) then covers all
cuts. Their gain gcds are 78 and 28842, respectively. In particular, neither
can supply either of its two targets `2H-1,2H+1`. These are controls, not
claims that arbitrary regular tournaments have constant weights.
Moreover `gcd(189,39)=3` and `gcd(95095,14421)=209=11*19`. Thus every
one-vertex extension in these two explicit families has a composite count,
divisible by three and 209 respectively. This follows from the exact formula
and checked weights, not from divisibility of a lower bound or an assumption
that the same constraint persists under further extensions.

This is an optimization *after* the exact weights have been computed; it is
not a polynomial-time algorithm for computing Hamiltonian path counts.

## Exact implementation and reproduction

`cut_response.py` uses two subset DPs for endpoint counts. It constructs
`W_uv = sum_{nonempty proper U} E[U,u] R[V-U,v]`, then `b`, in
`O(n^2 2^n)` arithmetic operations and `O(n 2^n+n^2)` memory. Computing all
`2^n` extension counts from these weights takes `O(n 2^n)` further arithmetic
operations. Integers are arbitrary precision.

Each `W_uv` is at most `(n-1)!`: concatenating a two-path cover gives a
permutation with the ordered pair `u,v` consecutive, and the cut is recovered
at that pair. Hence the weights and Hamiltonian counts have `O(n log n)` bits.
The arithmetic-operation bounds above therefore also translate to exponential
bit bounds with polynomial factors. Edmonds-Karp's usual strongly polynomial
augmentation bound applies to the subsequent constrained min-cut step.

```bash
python3 run_checks.py
# Or replay the main cut-law audit alone:
python3 verify.py > verify.normal.json
python3 -O verify.py > verify.optimized.json
cmp verify.normal.json verify.optimized.json
```

`run_checks.py` additionally runs `verify_consequences.py` and
`verify_examples.py` in normal and optimized modes, then copies all needed
sources to a temporary directory and replays them from an unrelated working
directory. It requires byte-identical receipts and saves them in `receipts/`.
The curated publication package additionally includes `audit/AUDIT.md` and
two independent checkers. Run `python3 run_audit.py` there to replay both
independent audits normally, under `-O`, and from another working directory.

The audit independently compares every cut with a forward subset DP on the
actual `(n+1)`-vertex tournament; it also uses literal permutations and
literal two-path-cover enumeration at the stated small orders. It checks
endpoint divergence, pair-deletion lower bounds, strict mixed differences,
the SCC minimizer count, exact gain gcd, and constrained/minimum-positive
min-cut optimization against brute force. It
exhausts all endpoint-relation words through length six in the block-insertion
identity. Malformed tournaments, the cut-orientation/factor convention, and
contradictory constraints are explicit controls. All checks survive `-O`.
The exact universe, check counts, source hashes, and output digest are in the
receipts; no finite test is being used as the proof of the all-order theorem.

## Provenance and literature boundaries

The all-cut polynomial (2) was independently developed in the concurrent
Session 124 bandwidth lane (`bandwidth/cut_polynomial.py`). The present lane
identified the endpoint-divergence cancellation (3), making the original
tournament a canonical integral response graph. The parent lane supplied an
initial block-insertion positivity argument; the telescoping version above
strengthens it to the pair-deletion bound.

Rédei oddness is classical. A primary mathematical exposition with a proof
is Irving and Omar, *A Survey of Rédei-Berge Symmetric Functions via Matrix
Algebra*, [Theorem 5](https://arxiv.org/html/2412.10572v1). This is the sole
external mathematical dependency of the main cut-law proof as written; the
optional packing corollary additionally uses the cited odd-cycle formula.

The initially proposed generic weighted-orientation feasibility route belongs
to the classical prescribed-degree orientation setting of S. L. Hakimi,
[*On the degrees of the vertices of a directed graph*, 1965,
290--308](https://doi.org/10.1016/0016-0032(65)90340-6). That general theorem is
not needed here: (3) explicitly constructs the required orientation. The
publication metadata was checked; no priority claim about the specialized
tournament identity is made, and this bounded literature pass is not an
exhaustive novelty search.

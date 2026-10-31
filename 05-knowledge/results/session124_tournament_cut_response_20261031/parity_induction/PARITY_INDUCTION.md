# A self-contained parity induction for the canonical cut law

**Status: PROVED.** This addendum closes the external parity dependency in the
Session 124 tournament cut-response packet. It does not change the theorem,
algorithms, certificates, or finite computations. Rédei's parity theorem is
classical; no novelty claim is made for the theorem or this proof route.

## 1. Integer counting before any parity assumption

For a finite tournament `T`, let `H(T)` count directed Hamiltonian paths as
ordered vertex sequences. Set `H(empty)=1`. All these counts are finite
nonnegative integers; we do not yet assume that they are positive or odd.

Assume for the moment that `T` has at least one vertex. Let `s_v,t_v` count
its Hamiltonian paths starting and ending at `v`. For distinct `u,v`, let
`W_uv` count ordered pairs of nonempty directed paths partitioning the vertex
set, with the first path ending at `u` and the second starting at `v`.
All `W_uv` are nonnegative integers by definition. There are no such endpoint
pairs when the parent has order one.

Adjoin a new vertex `x` with out-neighbor set `S`, and put `c_v=[v in S]`.
Deletion of `x` from a Hamiltonian path gives exactly one of three cases:

* a Hamiltonian path starting at `v`, if `x` was first;
* a Hamiltonian path ending at `u`, if `x` was last;
* an ordered two-path cover ending at `u` and starting at `v`, otherwise.

Consequently, with `F_T(S)=H(T+x_S)`,

```text
F_T(c) = sum_v s_v c_v + sum_u t_u(1-c_u)
         + sum_{u != v} W_uv(1-c_u)c_v
       = H(T) + sum_v L_v c_v
         - sum_{u<v}(W_uv+W_vu)c_uc_v,

L_v = s_v-t_v+sum_u W_uv.                                         (A)
```

For each original arc `u->v`, cutting and concatenating at that arc show that
`W_uv` counts Hamiltonian paths using it. Outgoing use minus incoming use in
one path is `+1` at the first vertex, `-1` at the last, and zero elsewhere.
For a one-vertex path both indicators are one, so their difference is zero.
Summing this balance over all paths gives

```text
s_v-t_v = sum_{v->u}W_vu - sum_{u->v}W_uv,
L_v     = sum_{v->u}(W_vu+W_uv).
```

Substituting into (A) yields the **unhalved cut identity**

```text
F_T(S) = H(T) + sum_{u in S, v outside S, u->v}(W_uv+W_vu).         (B)
```

Every step so far is an integer counting identity. No division by two,
Hamiltonian-path existence theorem, parity theorem, or assertion about the
weights `b_uv` has been used. If `T` is empty, (B) is instead checked directly:
both its count and the count of its one-vertex extension are one.

## 2. Arc reversal stays at the same order

Fix an original arc `u->v`, and let `T'` reverse only this arc. In each cover
counted by `W_uv` or `W_vu`, the two distinguished endpoints lie in different
paths. Their mutual edge is therefore never used inside either path.
Accordingly both cover counts are identical for `T` and `T'`.

Reversal loses precisely the Hamiltonian paths containing consecutive
`u,v`, counted by `W_uv`, and gains precisely those containing consecutive
`v,u`, counted by `W_vu`. All other Hamiltonian paths are unchanged. Hence

```text
H(T')-H(T) = W_vu-W_uv.                                           (C)
```

This identity, too, is independent of parity and positivity.

## 3. Induction over all tournaments of a fixed order

Let `P(n)` be the statement that **every** tournament on `n` vertices has
odd Hamiltonian-path count. The empty-set convention gives `P(0)`, and the
unique one-vertex path gives `P(1)` directly.

Suppose `n>=1` and `P(n)` holds. For any tournament `T` of order `n`, reversing
any one of its arcs produces another tournament `T'` of the same order.
Thus the induction hypothesis applies to **both** sides of (C), and

```text
W_vu-W_uv is even, hence W_uv+W_vu is even.                         (D)
```

When `n=1`, there are no arcs and (D) is an empty requirement. In all cases,
(B) now shows that every extension of `T` has odd count: its count is the
odd integer `H(T)` plus a sum of even integers.

Finally, every tournament of order `n+1` is such an extension: delete any
vertex `x`, apply the assertion to its `n`-vertex parent, and restore the
actual out-neighbor set of `x`. Therefore `P(n+1)` holds. Induction proves

```text
H(T) is odd for every finite tournament T,
W_uv+W_vu is even for every distinct pair in every tournament.
```

The order-one step in particular proves the order-two case without any edge
reversal being available in the parent. At order two, reversal can then be
used to prove the even cover-pair coefficients needed at order three.

There is no circularity: the only oddness assumption in (C) is at the
**current** parent order; (B) then proves oddness at the **next** order.
Integrality of `b_uv` is an output of (D), not an input to (B).

## 4. Integral positive weights and the rest of the core theorem

Define `b_uv=(W_uv+W_vu)/2`. Part 3 proves that these are nonnegative
integers, and (B) becomes

```text
H(T+x_S) = H(T) + 2 sum_{u in S,v outside S,u->v} b_uv.             (E)
```

The already elementary block-insertion argument supplies the stronger bound
`b_uv>=H(T-{u,v})`. To make the dependency closure explicit, fix a Hamiltonian
path `P=(y_1,...,y_m)` in the remainder and insert the consecutive block
`u,v` or `v,u`, temporarily leaving its internal edge free. Set
`A_i=[y_i->u]`, `B_i=[y_i->v]`, with `A_0=B_0=1` and
`A_{m+1}=B_{m+1}=0`. The number of valid insertions is

```text
N(P) = sum_{i=0}^m [A_i(1-B_{i+1})+B_i(1-A_{i+1})]
     = 2 + sum_{i=0}^m [A_{i+1}(1-B_i)+B_{i+1}(1-A_i)] >= 2.
```

Cutting the inserted block gives a cover counted by `W_uv+W_vu`. The cover
recovers `P`, its insertion position, and the block order, so these covers
are distinct even across different remainder paths. Therefore

```text
W_uv+W_vu >= 2H(T-{u,v}),
b_uv >= H(T-{u,v}) >= 1.                                         (F)
```

The last inequality now follows internally: a finite nonnegative odd integer
is at least one, by Part 3; an empty remainder has count one by convention.
No separate Hamiltonian-path existence theorem is needed.

Thus the published cut law, integral weights, strong pair-deletion bound,
strict positivity, SCC equality characterization, submodularity, gain gcd,
flow algorithms, small-gain/transitive-deletion consequence, and symmetry
classification all have self-contained elementary proofs in the combined
packet. The existing derivations and programs are unchanged.

## 5. Dependency boundary and attribution

The logical order is

```text
finite counting definitions
  -> deletion and endpoint-balance identities (A),(B)
  -> same-order arc-reversal identity (C)
  -> simultaneous all-order oddness and pair integrality (D)
  -> odd nonnegative counts are positive
  -> block insertion and strict weight lower bound (F)
  -> the existing core consequences.
```

This addendum supersedes only the original packet's wording that the core
proof needs an externally imported Rédei theorem. The classical attribution
remains: a cited exposition is Irving--Omar,
[*A Survey of Rédei--Berge Symmetric Functions via Matrix Algebra*,
Theorem 5](https://arxiv.org/html/2412.10572v1). The citation is now background
and attribution for the core parity statement, not an unproved dependency.

The optional bound using maximum disjoint-triangle packing still invokes
the separately cited odd-cycle collection formula; this addendum does not
remove that dependency. It does not prove global Hamiltonian-path spectrum
completeness, a closed recursive extension bank, or any new forbidden count.

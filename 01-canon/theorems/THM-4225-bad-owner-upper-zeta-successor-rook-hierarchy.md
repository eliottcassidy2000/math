---
id: THM-4225
title: "Bad-owner upper-zeta successor-rook hierarchy"
status: >
  PROVED + VERIFIED-EXACT / OPEN FIREWALL. For every proper owner set R,
  identifies the upper Boolean zeta sum of exact bad-owner counts with a
  cycle-free injective successor-rook count times (n-|R|-1)!. This gives
  exact one- and two-owner moment formulas in terms of indegrees and common
  inneighbors. It supplies a tournament-realizability sidecar but does not
  prove HYP-9081 or any refuted local two-owner product bound.
source: codex-successor-rook-synthesis-20260826
depends_on:
  - THM-4223-cyclic-cut-cover-boolean-mobius-hierarchy-and-two-bad-owner-obstruction
related:
  - THM-4115-uniform-ear-cut-walsh-variance-and-sharp-growth-refinement
  - THM-4224-order-ten-minimal-plus-min-two-bad-owner-obstruction
  - HYP-9081-strong-tournament-five-copy-endpoint-energy-inequality
script: 04-computation/tournament_bad_owner_successor_rook_thm4225.py
output: 05-knowledge/results/tournament_bad_owner_successor_rook_thm4225.out
script_sha256: 22793976f43d50b59fcd00ce8edfc763f86a18f6770f9ee67cf5caab4e7e540a
output_sha256: af885e474db035359368bdcafe9f106e16b0f51a95878c0a3356cca2a8052cf2
hash_basis: raw LF bytes
audit: >
  PASS. Direct cyclic-order enumeration over all 1,098 labelled tournaments
  of orders two through five checks all 33,864 owner subsets against an
  independently enumerated partial-successor rook count. It also checks
  10,648 two-owner codegree formulas, layer reversal, total mass, and the
  first moment. Normal, optimized, and fixed-hash-seed Python streams
  byte-match the frozen output.
---

# THM-4225 -- bad-owner upper-zeta successor-rook hierarchy

**PROVED + VERIFIED-EXACT / OPEN FIREWALL.**

THM-4223 cuts cyclic words after a chosen owner set and obtains the lower
Boolean transform `sum_(S subseteq R) C_S`. The complementary operation is to
force owners to remain bad and retain their actual cyclic successors. This
produces an upper transform and exposes the tournament-specific coordinate
lost by the bare owner tensor.

## 1. Successor rooks

Let `X` be a tournament on `V`, `|V|=n>=2`. For a cyclic ordering `gamma`,
modulo rotation, an adjacency after `i` is bad when its successor beats `i`.
Retain THM-4223's exact counts

```text
C_S=#{gamma: BadOwners(gamma)=S}.                         (1)
```

For `R subseteq V`, define the upper Boolean zeta sum

```text
Z_R=sum_(S superset R) C_S.                               (2)
```

Thus `Z_R` counts cyclic orders in which every owner in `R` is bad, with no
condition on the other owners.

For a proper subset `R proper_subset V`, a **successor rook placement** is a
map `f:R->V` satisfying

```text
f is injective,
f(i)->i in X for every i in R,
the positional arcs i |-> f(i) contain no directed cycle. (3)
```

The last condition concerns the partial successor function, not the
tournament orientation. Let `rho_X(R)` count these placements; the empty map
has count one.

> **Theorem 1 (upper-zeta successor-rook hierarchy).** For every proper
> `R subset V`,
>
> ```text
> Z_R=(n-|R|-1)! rho_X(R).                                (4)
> ```
>
> At the full set, reversal gives
>
> ```text
> Z_V=C_V=C_empty.                                        (5)
> ```

### Proof

Fix a cyclic order counted by `Z_R` and let `f(i)` be the cyclic successor of
`i`. Distinct owners have distinct successors, so `f` is injective. Because
`i` is bad, `f(i)->i` in the tournament. Since `R` is proper, the restriction
of one cyclic successor permutation cannot already contain a closed cycle:
such a cycle would be a proper closed component of the full cyclic order.
Thus `(3)` holds.

Conversely, a placement satisfying `(3)` is a disjoint union of directed
paths in the positional graph `i |-> f(i)`. Contract its `|R|` prescribed
adjacencies into ordered blocks. There remain `n-|R|` blocks, and they have

```text
(n-|R|-1)!
```

cyclic orders modulo rotation. Expanding the blocks makes each `f(i)` the
successor of `i`, so every owner in `R` is bad. Restriction and contraction
are inverse, proving `(4)`.

For `R=V`, reverse a cyclic order. Every good adjacency becomes bad and every
bad adjacency becomes good, with the owner shifted to the other endpoint.
This is an involution and sends all-good orders bijectively to all-bad orders,
proving `(5)`. QED.

## 2. The first two successor moments

Write

```text
d_i^-=|N^-(i)|,
kappa_ij=|N^-(i) intersect N^-(j)|.                     (6)
```

> **Corollary 2 (indegree and codegree identities).** For `n>=2`,
>
> ```text
> sum_(S contains i) C_S=(n-2)! d_i^-.                   (7)
> ```
>
> For distinct `i,j` and `n>=3`,
>
> ```text
> sum_(S contains {i,j}) C_S
>   =(n-3)!(d_i^- d_j^- - kappa_ij).                     (8)
> ```

For one owner, a rook is simply the choice of one inneighbor. For two
owners, there are `d_i^-d_j^-` choices and injectivity deletes exactly the
`kappa_ij` common images. A two-cycle in the positional map would require
both `i->j` and `j->i` in the tournament, so no further cycle exclusion is
possible. Equation `(4)` now gives `(7)--(8)`.

In THM-4223's low-layer notation `A_i=C_{i}` and `B_ij=C_{i,j}`. Positivity
of all higher layers therefore gives the exact necessary bounds

```text
A_i+sum_(j!=i)B_ij <=(n-2)!d_i^-,
B_ij <=(n-3)!(d_i^-d_j^- - kappa_ij).                   (9)
```

These are successor-compatible tournament constraints. They are generally
far from sufficient for the five-copy inequality.

## 3. Scalar moment consequences

Put `C_k=sum_(|S|=k)C_S`. Reversal in the proof of `(5)` changes the number
of bad adjacencies from `k` to `n-k`, so

```text
C_k=C_(n-k).                                             (10)
```

Summing `(7)` over owners gives the fixed first moment

```text
sum_S |S|C_S=n!/2.                                       (11)
```

For `n>=3`, summing `(8)` over pairs gives the second factorial moment

```text
sum_S binom(|S|,2)C_S
 =(n-3)![sum_(i<j)d_i^-d_j^-
          -sum_x binom(d_x^+,2)].                       (12)
```

Indeed, a vertex `x` is a common inneighbor of exactly
`binom(d_x^+,2)` unordered owner pairs.

## 4. What the new sidecar retains

The two Boolean transforms now have different typed meanings:

```text
lower: N_R=sum_(S subseteq R)C_S
       =cut after R and obtain a cyclic path cover;

upper: Z_R=sum_(S superset R)C_S
       =force bad successors on R and count cycle-free rooks. (13)
```

The owner tensor retains the domains of bad adjacencies but forgets their
successor images. The rook placement retains source, target, injectivity, and
the cycle obstruction. At orders at least four, higher placements also see
directed cycles among forced successors; pairwise indegree/codegree data lose
that information.

THM-4224 refutes both the naive product bound and its `+min` repair. Nothing
in `(4)--(12)` restores either claim, proves the open multiplicative
`27/25` candidate suggested by that hostile family, or proves HYP-9081. The
gain is a precise all-order realizability hierarchy on which an aggregate
switching argument can act.

## 5. Exact audit

The companion independently enumerates every labelled tournament of orders
two through five, every cyclic order with vertex zero fixed, every owner
subset, and every injective partial successor map. It checks `(4)` on all
`33,864` subset rows and `(8)` on `10,648` pair rows. `C3` is the positive
cycle control; the transitive tournaments of orders two through five are
hostile controls with `C_empty=C_V=0`.

Replay with

```bash
python3 -B 04-computation/tournament_bad_owner_successor_rook_thm4225.py
python3 -B -O 04-computation/tournament_bad_owner_successor_rook_thm4225.py
PYTHONHASHSEED=4225 python3 -B \
  04-computation/tournament_bad_owner_successor_rook_thm4225.py
```

All three streams byte-match the frozen output. **QED.**

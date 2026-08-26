---
id: THM-4145
title: "Rooted homogeneous-pair expansion and two-defect formula"
status: >
  PROVED ELEMENTARY ROOTED SUBSTITUTION + EXACT TWO-DEFECT FORMULA +
  VERIFIED-EXACT + INDEPENDENTLY AUDITED. Pair expansion preserves and
  reflects strongness. Its Hamilton count is an exact weighted sum over
  root-deletion orders with zero, one, or two repairable defects; its base
  count is one rooted ear-field value, while its child field is a mixed
  two-ear object. This reduces the order-eleven homogeneous-pair stratum to
  93,559,490 rooted presentations but does not prove their Johnson centrality.
source: codex-frontier-synthesis-creative-20260825ac
depends_on:
  - THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
  - THM-4137-strong-tournament-centrality-complete-order-ten
related:
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4144-order-eleven-large-homogeneous-module-johnson-centrality
script: 04-computation/tournament_rooted_pair_expansion_two_defect_thm4145.py
output: 05-knowledge/results/tournament_rooted_pair_expansion_two_defect_thm4145.out
independent_audit_script: 04-computation/tournament_rooted_pair_expansion_two_defect_thm4145_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_rooted_pair_expansion_two_defect_thm4145_independent_audit.out
script_sha256: 475f152a5c3a5d03691f953c08be68fb17157a6e66444d18de0876ed78d8db86
output_sha256: f282d82707fb5364b9c42cfe88b56b250d111c73f3121678564b590074c1c503
semantic_sha256: d174d868d8d1784901ef3926b0d0dc03e27a842eadd9a31df931c00716fa11cc
semantic_fnv64: 9f075ab9b9d3f07f
independent_audit_script_sha256: 40c834de56220724a295c69af98487112f1136500fe1fb7353ffd4ad961678c3
independent_audit_output_sha256: 782091fa4b7dd0e76336b1e2556e333b5952d2dea7f7120a0c1af0d4efb4e011
independent_semantic_fnv64: 9f075ab9b9d3f07f
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone Python implementation uses subset DP for literal
  Hamilton counts and permutation enumeration for the deletion formula and
  both internal exposed-gap orientations. It checks strongness equivalence,
  the rooted-ear adjacency identity, the formula, and the exposure scaling
  on all 5,404 rooted labelled tournaments of orders two through five.
  Normal, optimized, and hash-seeded replays byte-match.
independent_audit: >
  ACCEPT. A separate warning-clean standard-library C++ implementation
  rebuilds the entire labelled universe, reachability, subset Hamilton DP,
  literal block exposure, and defect-slot enumeration. It agrees on every
  semantic row, the minimal square-only hostile, and the root-loss control;
  optimized, unoptimized, and sanitizer controls agree.
---

# THM-4145 -- rooted pair expansion and two-defect formula

**PROVED ELEMENTARY ROOTED SUBSTITUTION + EXACT TWO-DEFECT FORMULA +
VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4144 closes all order-eleven homogeneous modules of sizes three through
nine. A homogeneous pair is the remaining modular boundary, and contraction
now retains a distinguished vertex. That root is not bookkeeping: even the
base Hamilton response depends on which vertex is expanded.

## 1. Statement

Let `(Q,r)` be a tournament of order `q>=2` with distinguished vertex `r`.
Put `U=V(Q)\{r}`. Define `P_r(Q)` by replacing `r` with two labelled clones
`a,b`, orienting `a->b`, and imposing for every `u in U`

```text
a->u  iff  b->u  iff  r->u.                              (1)
```

The opposite internal orientation gives an isomorphic unlabelled expansion
after swapping the clones.

> **Theorem.** The following hold.
>
> 1. `P_r(Q)` is strong if and only if `Q` is strong.
> 2. In THM-4114's ear convention `x_S->v iff v in S`, if
>    `B_r=N_Q^+(r)`, then
>
>    ```text
>    H(P_r(Q))=F_Q(B_r).                                  (2)
>    ```
>
> 3. For every ordering `pi=(u_1,...,u_m)` of `U`, `m=q-1`, let
>
>    ```text
>    D(pi)={i: 1<=i<m and u_i does not point to u_(i+1)}. (3)
>    ```
>
>    Number the insertion slots `0,...,m`. Let `G_r(pi)` contain slot zero
>    when `r->u_1`, slot `m` when `u_m->r`, and an internal slot `i` when
>    `u_i->r->u_(i+1)`. Write `k(pi)=|G_r(pi)|`. Then
>
>    ```text
>    H(P_r(Q))=sum_pi phi(D(pi),G_r(pi)),                 (4)
>
>    phi(D,G)=0                         if D is not a subset of G or |D|>2,
>             |G|^2                    if |D|=0,
>             2|G|-1                  if |D|=1 and D is a subset of G,
>             2                       if |D|=2 and D is a subset of G.     (5)
>    ```
>
> 4. Let `E_T(x,y)` denote THM-4128's exposed ordered-gap count, called
>    `Q(x,y)` there. For the internal clone arc,
>
>    ```text
>    E_(P_r(Q))(a,b)=E_(P_r(Q))(b,a)=H(Q),
>    w_ab(P_r(Q))=H(Q).                                  (6)
>    ```

The order restriction only removes the one-vertex convention: pair expansion
of a single vertex is transitive of order two, while a singleton may be called
strong by a vacuous reachability convention.

## 2. Strongness is exactly preserved

If `Q` is not strong, its strongly connected components have a nontrivial
linear order. Expanding the component side containing `r` by `{a,b}` lifts
an initial or final one-way cut of `Q` to a one-way cut of `P_r(Q)`. Hence a
strong expansion forces a strong quotient.

Conversely suppose `Q` is strong. Quotient paths with endpoints outside `r`
lift through either clone whenever they meet `r`, because `(1)` makes every
cross-block arc uniform. Paths from an outside vertex to either clone and
from either clone to an outside vertex lift in the same way, choosing the
desired clone at the endpoint. The arc `a->b` gives one internal direction.
Since a strong tournament of order at least two has a directed cycle through
`r`, a cycle

```text
r->u_1->...->u_s->r                                    (7)
```

lifts to `b->u_1->...->u_s->a`, giving the other direction. Thus every
ordered pair of vertices in `P_r(Q)` is connected. This proves part 1.

## 3. The rooted ear identity

Keep `a` under the old name `r` and view `b` as a new ear. For `u!=r`,
condition `(1)` says

```text
b->u  iff  u in B_r.                                    (8)
```

The internal arc is `r=a->b`, so `b` does not beat `r`; equivalently
`r notin B_r`. Thus the adjacency matrix obtained by adding the ear
`x_(B_r)` to `Q` is literally the adjacency matrix of `P_r(Q)`, with no
quotient or numerical loss. Taking Hamilton counts proves `(2)`.

## 4. Two clones repair at most two deletion defects

Fix an ordering `pi` of `U` and delete the clones from a Hamilton path of
`P_r(Q)`. Every backward adjacency left in `pi` must have contained a clone
before deletion. Therefore all defects lie in clone-occupied slots and there
are at most two. A clone can occupy a slot precisely when the two boundary
arcs have the root type, which is exactly membership in `G_r(pi)`.

Conversely, insert the labelled clones into good slots and require every
defect slot to be occupied. Every resulting adjacency is then forward. If
both clones occupy one slot, only the internal order `a,b` works. If they
occupy distinct slots, either assignment of the labels works. The cases are
therefore:

- no defect: `k` same-slot insertions and `2 binom(k,2)` distinct-slot
  insertions, totalling `k^2`;
- one defect: one same-slot insertion there and `2(k-1)` ways to put the
  other clone in another good slot, totalling `2k-1`;
- two defects: the clones occupy the two forced slots in either order,
  totalling `2`;
- an unrepairable or third defect: no insertion.

Deletion and insertion are inverse with the underlying order `pi` fixed, so
there is neither omission nor multiplicity. Summing gives `(4)--(5)`.

## 5. Internal exposure and the mixed two-ear boundary

For the internal arc `a->b`, a Hamilton path in which the correct adjacent
block `a,b` is exposed contracts bijectively to a Hamilton path of `Q`.
Likewise an exposed ordering with the sole reversed adjacency `b,a`
contracts bijectively to a Hamilton path of `Q`. These are precisely the two
ordered-gap counts in THM-4128, so

```text
2w_ab=E(a,b)+E(b,a)=2H(Q),                              (9)
```

proving `(6)`.

The base identity `(2)` does not transport the whole child ear field. Let
`b_r` denote the fixed root-type ear in `(8)`, and define

```text
K_(Q,r)(A,epsilon)=H(Q+b_r+x_(A,epsilon)),              (10)
```

where `epsilon` records the second ear's orientation against `b_r`. For
`S subseteq V(P_r(Q))`, the exact identity is

```text
F_(P_r(Q))(S)
 =K_(Q,r)(S intersect V(Q), 1_(b_r in S)).              (11)
```

Thus `H(P_r(Q))` is one rooted one-ear value, but the child exposure packet
and layer cosets live on the mixed two-ear slice `(10)`. Forgetting that
slice is the first failed implication in any attempted transfer of THM-4137
centrality from parent to child.

## 6. Order-eleven coverage and exact residual

Contract any homogeneous pair in a strong order-eleven tournament and mark
the contracted vertex. Part 1 gives a rooted strong order-ten tournament.
Conversely every rooted strong order-ten tournament expands to a strong
order-eleven tournament with a homogeneous pair. THM-4137 counts `9,355,949`
strong order-ten isomorphism classes, so the complete covering universe has

```text
10*9,355,949=93,559,490                                 (12)
```

rooted presentations. Root automorphisms and multiple homogeneous pairs can
duplicate a target, so `(12)` is not an isomorphism-class count.

Together with THM-4144, the modular part of the order-eleven centrality
problem has now been separated exactly:

```text
large modules |M|=3..9: closed by THM-4144,
pair modules |M|=2:     rooted mixed-two-ear census (12),
prime tournaments:     no substitution reduction.       (13)
```

## 7. Hostiles, loss ledger, and finite audit

The tempting square-only formula keeps only zero-defect deletion paths. For
the directed triangle with little-endian label `101`, rooted at zero, its
pair expansion has label `101101` and

```text
H(P_0(C3))=5,                    square-only sum=4.       (14)
```

The missing path contracts to a one-defect deletion order, so `(14)` is the
minimal hostile that forces the middle line of `(5)`.

The root cannot be dropped either. The strong order-four label `001000` has

```text
H(P_0(Q))=11,                     H(P_1(Q))=9.            (15)
```

The typed connection is

```text
source:       rooted tournament (Q,r)
target:       homogeneous-pair expansion P_r(Q)
map:          replace r by a->b with uniform external orientation
preserved:    strongness, full labelled adjacency, base response via (2)
destroyed:    the root after forgetting the presentation; target multiplicity
sidecar:      mixed two-ear field K_(Q,r), including the inter-ear bit
hostiles:     (14) defeats zero-defect squaring; (15) defeats unrooting
decisive test: a root-aware mixed-two-ear packet/coset census over (12).
```

The primary and independent audits exhaust every rooted labelled tournament
of orders two through five:

| quotient order | labels | rooted rows | strong rooted rows |
|---:|---:|---:|---:|
| 2 | 2 | 4 | 0 |
| 3 | 8 | 24 | 6 |
| 4 | 64 | 256 | 96 |
| 5 | 1,024 | 5,120 | 2,720 |
| **total** | | **5,404** | **2,822** |

On every row they compare literal child Hamilton DP, `(4)--(5)`, both
ordered internal-gap counts, the ear adjacency, and strongness. Their common
semantic FNV-64 is `9f075ab9b9d3f07f`; the primary also freezes the row
SHA-256 `d174d868...a11cc`.

Replay with

```text
python3 -B 04-computation/tournament_rooted_pair_expansion_two_defect_thm4145.py
python3 -B -O 04-computation/tournament_rooted_pair_expansion_two_defect_thm4145.py
PYTHONHASHSEED=271828 python3 -B 04-computation/tournament_rooted_pair_expansion_two_defect_thm4145.py
clang++ -std=c++17 -O3 -Wall -Wextra -Werror \
  04-computation/tournament_rooted_pair_expansion_two_defect_thm4145_independent_audit.cpp \
  -o /tmp/thm4145-independent-audit
/tmp/thm4145-independent-audit
```

This theorem proves the structural reduction and exact base formula, not the
centrality candidate over `(12)`, a prime-tournament result, or any
dimension-free inheritance law. **QED.**

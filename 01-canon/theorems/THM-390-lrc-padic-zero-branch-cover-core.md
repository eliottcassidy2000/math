---
id: THM-390
name: lrc-padic-zero-branch-cover-core
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S547
depends_on:
  - THM-366
  - THM-369
  - HYP-2036
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean
---

# THM-390: LRC zero-branch sieve witnesses and the AP cover core

## Statement

Let `n >= 3` and let `V` be a Lonely Runner speed set with threshold `1/n`.
For each denominator node

```text
q = 2,3,...,n
```

define the product zero-branch mass

```text
z_q(V) = #{v in V : q divides v}.
```

Then:

1. If `2 <= q < n` and `z_q(V)=0`, then `t=1/q` is an open lonely witness:
   every speed has distance at least `1/q > 1/n` from an integer.
2. If `q=n` and `z_n(V)=0`, then `t=1/n` is a closed, compactified wall
   witness: every speed has distance at least `1/n` from an integer.
3. Therefore any speed set with no open sieve witness among these denominator
   nodes must cover every `q=2,...,n-1`; any speed set with no compact sieve
   witness must cover every `q=2,...,n`.
4. For the initial segment

   ```text
   A_n = {1,2,...,n-1},
   ```

   every open node `q<n` is covered and the compact node `q=n` is empty.
5. The unique minimum subset of `A_n` that covers all open nodes
   `q=2,...,n-1` is

   ```text
   U_n = {u in A_n : 2u >= n}.
   ```

   Its size is `floor(n/2)`.  In particular, for even `n` this is the
   singleton-carrier block

   ```text
   n/2, n/2+1, ..., n-1
   ```

   and the minimum open cover size is exactly `n/2`.

## Proof

Parts (1) and (2) are the denominator-sieve lemma THM-366, also
machine-checked as `sieve_one_div` in THM-369.  If `q` divides no speed, then
for each `v in V` the residue `v mod q` is nonzero, so

```text
||v/q|| >= 1/q.
```

When `q<n`, this is strictly larger than `1/n`, giving an open witness.  When
`q=n`, it is exactly the closed wall threshold allowed by LRC.

Part (3) is the contrapositive of the same statement.

For (4), if `2 <= q < n`, the speed `q` itself belongs to `A_n`, so `z_q(A_n)
>= 1`.  No positive member of `A_n` is divisible by `n`, so `z_n(A_n)=0`.

For (5), first prove that `U_n` covers every open node.  Let `2 <= q < n`.
If `2q >= n`, then `q in U_n` and the node `q` is covered by the speed `q`.
If `2q < n`, choose

```text
k = ceil(n/(2q)).
```

Then `kq >= n/2`, while

```text
kq < n/2 + q < n.
```

Thus `u=kq` is a member of `A_n`, satisfies `2u >= n`, and is divisible by
`q`.  Hence `u in U_n` covers `q`.

Conversely, each `u in U_n` is a singleton carrier inside `A_n`.  If
`w in A_n` is divisible by `u`, then `w=mu` with `m>=1`.  For `m>=2`, we get

```text
w >= 2u >= n,
```

contradicting `w <= n-1`.  Therefore the only speed in `A_n` divisible by
`u` is `u` itself.  Any open cover must include every `u in U_n`.

The set `U_n` has `floor(n/2)` elements: if `n=2m`, it is
`{m,m+1,...,2m-1}`; if `n=2m+1`, it is `{m+1,m+2,...,2m}`.  Since every cover
contains `U_n` and `U_n` is itself a cover, it is the unique minimum cover.

## Tournament Analysis Reading

This theorem justifies the vertex choice used in HYP-2036:

```text
vertices = denominator/product-zero-branch obligations q=2,...,n.
```

The pairwise observable is the pair `(z_q, divisibility comparability)`.  The
tree-trienerment gauge from HYP-2036 orients lower `z_q` first, sends equal
mass along deeper comparable divisibility branches, and ties equal-mass
incomparable branches with the increasing-`q` Hamiltonian path.  THM-390 proves
the LRC predicate preserved by this quotient: an empty open node is already a
lonely witness, and the AP row covers exactly all open nodes while leaving the
compact wall node empty.

The theorem does not prove the stronger HYP-2036 computational observations
for arbitrary bounded open survivors, such as zero mixed fingerprint fibers or
minimum cover size `n/2` outside the AP row.  It proves the sieve semantics and
the AP cover-core skeleton those observations refine.

## Notes

The proof challenges the runner-vertex default in the cleanest possible way:
for the sieve layer, the natural vertices are proof obligations `q=2,...,n`,
not runners, gaps, arcs, or walls.  This quotient preserves the existence of
explicit denominator witnesses and destroys non-sieve lonely times, event
order, circular gap geometry, and most runner identity except singleton-carrier
roles.

This also separates the new entropy-on-tree layer from the exact local gate.
HYP-2037/HYP-2038 read AP/regular rows as high-entropy or critical boundary
states.  THM-390 says which finite denominator leaves make that boundary a
wall witness: total spread is not enough; the local condition is the zero
branch mass `z_q`, with AP's singleton leaves forced by divisibility.

In HYP-2039's defect-transport language, the theorem describes the pinned
transport successes: if `z_q=0`, the guaranteed large hole is already centered
at the observer at `t=1/q`.  The hard branch begins only after every such
denominator node is covered.

## Verification Record

Lean covers the denominator-witness part through `sieve_one_div`:

```text
04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean
05-knowledge/results/lean_lrc_trienerment_s542.out
```

The AP cover-core proof above is elementary and paper-formal.  The even-`n`
special case matches the S546 computation:

```text
04-computation/lrc_padic_tree_trienerment_s546.py
05-knowledge/results/lrc_padic_tree_trienerment_s546.out
```

**QC verification (monad-reviewer-2026-06-01):** Re-derived parts (4)-(5). For
`U_n={u∈A_n : 2u≥n}`: coverage — if `2q≥n` then `q∈U_n` covers node `q`; if `2q<n`,
`u=⌈n/(2q)⌉·q` satisfies `n/2 ≤ u < n/2+q < n`, is divisible by `q`, and has `2u≥n`,
so `u∈U_n` covers `q`. Minimality — each `u∈U_n` is a singleton carrier (any
multiple `mu∈A_n` with `m≥2` would have `mu≥2u≥n>n-1`), so every cover must contain
`U_n`; hence `U_n` is the unique minimum cover, `|U_n|=⌊n/2⌋`. Sieve semantics
(parts 1-2) follow from `‖v/q‖≥1/q` when `q∤v`. Proof CONFIRMED. Note: this is the
**first claimant** of id THM-390; the S548 star-peeling theorem was renumbered to
THM-391 (MISTAKE-052).

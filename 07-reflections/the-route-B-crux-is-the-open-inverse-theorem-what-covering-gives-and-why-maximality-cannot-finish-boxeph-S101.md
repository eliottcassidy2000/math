# The Route B crux is the open inverse theorem (= LRC(14)): what covering rigorously gives, and why maximality cannot finish it

*boxeph-2026-07-18-S101. Owner: prove the Route B compact/AP-rigidity crux. Honest outcome: this
statement is **equivalent to the open LRC(14) covering case** (proved S94), so it is not closeable in a
session, and I will not manufacture a proof. What I can and do add rigorously: (i) the cleanest reduction
— crux ⟺ the maximizer residues have **at most one gap `< val`** (packing lemma); (ii) a **proved**
covering constraint — `M<1/13 ⟹ every q'∈{2,…,13} divides some speed`; (iii) a **proved obstruction** —
interior small gaps are *invisible to local maximality*, so difference-closure is a **global** (covering)
fact, not a variational one. This decodes why ~20 sessions of variational/reduction attacks stalled.
LRC(14) not closed. Verified S101 computation.*

## The crux, in its cleanest form

At the maximizer `t*=a/q` with `M=val/q<1/13`, `q∈(13val,14val)` (extremal `q=13val+1`), the 13 residues
`r_i=v_i a mod q` lie in the band `[val, 12val+1]` (length `11val+1`). Sorted, their 12 gaps sum to
`≤11val+1`, so **pigeonhole forces ≥1 gap `<val`**. The **packing lemma** (S90, proved) says: if 12 of
the residues are pairwise `≥val` and include `v_+` (residue `val`), they are exactly `val·{1,…,12}`, and
then the core is the dilated AP `v_+·{1,…,12}`. So:

> **Crux ⟺ at most one residue gap is `<val`** (equivalently: the core `V∖{v_max}` is difference-closed;
> equivalently `ρ=v_max/v_2nd≥13`; equivalently the S94 one-line form).

A gap `<val` between `r_i,r_j` means `w=|v_i−v_j|` has `‖w t*‖<M`, so `w∉V` (else `‖w t*‖≥M`). Thus a
small gap = a **non-speed difference**. Difference-closure = no interior small gap = the whole crux.
Verified on the tower `{d·1,…,d·12, d·182}`: every member has **exactly one** small gap (`=1`), residues
`val·{1,…,12}∪{12val+1}`.

## What covering rigorously gives (proved)

`M<1/13 ⟹` covering is **automatic** (non-covering yields an LRC(≤13) sieve witness with `M≥1/13`). So
covering is a *consequence* to mine, not an extra hypothesis. Mining the sieve at each small modulus:

> **Lemma (sieve divisibility).** If `M(V)<1/13` then for every `q'∈{2,…,13}`, some speed is divisible
> by `q'`.

*Proof.* `M<1/13` means every witness loses: for all `q'` and all `a'`, `min_v‖a'v‖_{q'} < q'/13`. For
`q'≤13` the RHS is `≤1`, and `‖·‖_{q'}` is a nonnegative integer, so `min_v‖a'v‖_{q'}=0`, i.e. some `v`
has `q'∣a'v`. Taking `a'=1` (coprime to `q'`): `q'∣v`. ∎

Verified on the tower (all `q'∈{2..13}` divide a speed; the far element carries `13` via `182=13·14`).
This is genuine but **insufficient**: it constrains the speed set (forces a multiple of each `q'≤13`)
without forcing the AP — many non-AP 13-sets satisfy it. Beyond `q'≤13` the sieve condition
`min_v‖a'v‖_{q'}<q'/13` is no longer pure divisibility, and taken over *all* `q'` it is exactly `M<1/13`
again — no free structure past the divisibility layer.

## Why maximality cannot finish it (proved obstruction)

The natural attack — perturb `t*` to force structure — provably **cannot** reach difference-closure:

> **Obstruction.** At `t*`, only the **active** runners (`‖v t*‖=M`, residues `±val`) constrain the local
> maximum. A small gap between two **interior** residues (`r_i,r_j∈(val,12val+1)`, `‖v_{i,j}t*‖>M`)
> involves a non-speed `w=v_i−v_j` that is not a runner, so `‖w t*‖` does not enter `min_v‖v t*‖`.
> Perturbing `t*→t*+δ` moves `‖w t*‖` but never changes which runners are the minimum near `t*`. Hence
> **interior small gaps are invisible to every variational/perturbative argument** about `t*`.

Consequence: difference-closure of the core is **not** a local (maximality) property — it holds only
because `t*` is the max over **all** rationals simultaneously (the global covering/sieve structure). This
is the precise reason the reduction chain (difference-closure lemma S87 → dimension-2 S92 → offset-AP S93
→ `j_1=0` S94) kept *restating* the crux rather than proving it: each step is a variational or
coordinate move, and the missing content lives in the **global** comparison across moduli, which no local
move can supply.

## Net (honest)

- **Not proved:** the crux itself — it is equivalent to the open LRC(14) covering case (S94). I did not
  and will not fabricate a proof of an open problem.
- **Proved this session:** (i) the clean reduction crux ⟺ "`≤1` residue gap `<val`" ⟺ core
  difference-closed; (ii) the sieve-divisibility lemma `M<1/13 ⟹ q'∣` some speed for all `q'≤13`;
  (iii) the obstruction — interior small gaps are invisible to maximality, so difference-closure is
  irreducibly **global**.
- **The one remaining lever** (named, not pulled): a **global** argument that ties the mod-`q` band
  packing to the cross-modulus sieve divisibility — i.e. show that two interior small gaps would force a
  sieve window free at some modulus `q''` (a witness `M≥1/13`), contradiction. This is the inverse
  theorem; it is exactly where the project's `≥6`-linear / additive-dimension content sits (klein-S279,
  S92), and it is open.

LRC(14) sits, as established across S96–S100 (density route discharged for separated far elements), on
precisely this Route B inverse theorem. This session sharpens it to "no two interior small gaps," proves
the covering divisibility it must respect, and proves that the finish cannot be variational — it must be
a global cross-modulus counting argument.

Cross-links:
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
[[o-R-is-false-S-is-theta-R-but-the-explicit-O-R-bound-closes-the-density-route-for-separated-far-elements-boxeph-S100]],
[[the-abandoned-attempts-decoded-the-crux-is-additive-dimension-two-not-any-scalar-boxeph-S92]],
THM-1017 (AP-core bridge), the S90 packing lemma, [[both-lrc14-routes-bottom-on-the-same-multilinear-cancellation-not-one-mollification-klein-S279]].

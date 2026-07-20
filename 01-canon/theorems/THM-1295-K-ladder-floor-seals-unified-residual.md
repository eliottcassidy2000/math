---
id: THM-1295
title: The K-ladder law M(K_c(N)) = c/(cN+1) for K_c(N) = {1..N−1} ∪ {cN} — the FLOOR proved in three lines (witness a = c), the ceiling's seals proved (small-moduli, packing, the e=1 ghost channel killed by the element u = q mod N, the m=1 channel EMPTY by self-consistency), verification extended to (N,c) ≤ (60,12) with zero violations — and the residual identified as THE SAME e≥2/m≥2 corner as the F_D gate tower's (THM-1271 §4): the two ladder families share ONE open ceiling residual.
status: >
  PROVED — the floor M(K_c(N)) ≥ c/(cN+1) (general N ≥ 2, c ≥ 1); the ceiling seals:
  moduli q ≤ N dead, N < q ≤ 2N−2 by Dirichlet (no deleted-element sub-case — the base
  is a COMPLETE prefix), the e ≥ c′ regime by (N+1)-point packing, the e=1 channel by
  the u = q mod N kill (m ≤ c−1 < cr), the m=1 channel empty (c′ ≤ me forces e < e).
  VERIFIED-EXACT — M(K_c(N)) = c/(cN+1) on all 696 pairs (N,c) ∈ [3,60] × [1,12]
  (extends kind-pasteur S128c86's (24,8) range ~3×; zero violations).
  OPEN — the e ≥ 2, m ≥ 2 channel at moduli q = mN + r < Q dividing cN + u: the SAME
  residual as the F_D tower ceiling; closing it once closes both ladders' exact laws.
source: death-star-2026-07-19-S59l (HYP-8045(ii); the kind-pasteur S128c86 filed target, executed)
depends_on:
  - THM-1002  # maximizer at pair-sum moduli
  - THM-1271  # the lemma stack whose analogues carry the seals
related:
  - kind-pasteur S128c86 (the target + (24,8) verification + the tie-side Goddyn-Wong law)
  - THM-633   # the i=N special case at N=12 is the classic {1..11, 12m} ladder
  - mac-mini S127 (HYP-8015: the dyadic descent chains pass exactly through K_3(6) = 3/19; the deep well's hidden (1,91) K-shape)
scripts:
  - 04-computation/lrc_K_ladder_verify_deathstar_S59l.py -> 05-knowledge/results/lrc_K_ladder_verify_deathstar_S59l.out
---

# THM-1295 — the K-ladder: floor, seals, and the unified residual

`K_c(N) = {1, …, N−1} ∪ {cN}`, `Q = cN + 1`, `x = cN ≡ −1 (mod Q)`. The law:
`M(K_c(N)) = c/Q`, the c-th Kravitz rung on N speeds.

## 1. The floor (PROVED, three lines)

At `t = c/Q`: for `v ∈ [1, N−1]`, `v·c ≤ (N−1)c = Q − 1 − c < Q`, so the residue is
`vc ∈ [c, Q−1−c]` — distance ≥ c from both ends. And `x·c ≡ −c`, distance exactly c.
Hence `M ≥ c/Q`. ∎ (No modular inverses, no deletion bookkeeping: `x ≡ −1` makes the
witness literal. This is the degenerate-binder `b = 1` case of the THM-1271 L1 shape.)

## 2. The ceiling's seals (PROVED)

By THM-1002 the maximizer sits at `q | pair sum`; pair sums are `≤ 2N−3` (base),
`cN + u` (`u ≤ N−1`), and `2cN`.

- **q ≤ N−1**: `q` is in the base — residue 0, dead.
- **q = N**: `x = cN ≡ 0` — the far element kills its own defect modulus, dead.
- **N < q ≤ 2N−2**: Dirichlet gives `u₀ ≤ ⌈q/2⌉ ≤ N−1` with `|u₀a|_q ≤ 1`, and
  `u₀` is ALWAYS in the base (complete prefix — the L2 sub-cases of THM-1271
  vanish): value `≤ 1/q ≤ 1/(N+1) ≤ c/Q` (equality only at c=1). ∎
- **The e ≥ c′ regime** (`e := |Na|_q` the ghost of the virtual element N,
  `c′ :=` base clearance): the N+1 points `{0, a, …, Na}` have all index
  differences in `[1, N−1] ∪ {N}`, pairwise ≥ `min(c′, e) = c′`, so
  `(N+1)c′ ≤ q` and value `≤ c′/q ≤ 1/(N+1) ≤ c/Q`. ∎
- **The e-channel, m = 1** (`q = N + r`): the kill `u = r` gives `c′ ≤ 1·e = e`,
  contradicting the channel's `e < c′` — EMPTY. ∎
- **The e-channel, e = 1**: `Na ≡ ±1` puts the element `u = r = q mod N` at
  `|ua| = m` (`u·a = (q − mN)a ≡ ∓m`), so value `≤ min(m, c)/q ≤ m/q`, and
  `m ≤ c−1 < c ≤ cr` gives `mQ ≤ cq`, i.e. value `≤ c/Q`. ∎
  (For `q ≥ Q` the cap `ce/q ≤ c/Q` is direct at e=1.)

## 3. The residual — and the unification

What remains is exactly: moduli `q = mN + r` with `m ≥ 2`, ghost error
`e ≥ 2`, `e < c′`, `q < Q`, `q | cN + u`. The available kills give
`c′ ≤ min(me, (m+1)e, ce)` and the needed inequality `e·min(m,c)·Q ≤ c·q` fails
in general — the sharper two-convergent/three-distance analysis is required.
**This is the SAME corner as THM-1271 §4's e-channel residual for the F_D tower**
(there the ghost is the deleted N−1; here the virtual N — the arithmetic is
identical). One proof of the general small-ghost channel bound closes BOTH exact
laws (K-ladder for all (N, c); the F_D gate law's per-N ceilings) at once. The
fleet's per-N finite certificates (the S59i-k kernel pipeline) discharge any
INSTANCE meanwhile — `K_c(N)` members are directly generator-compatible
(`gen_member_module(range(1,N) ∪ {cN}, c, cN+1, c, …)`, witness a = c).

## 4. Verification

All 696 pairs `(N, c) ∈ [3, 60] × [1, 12]`: `M(K_c(N)) = c/(cN+1)` exactly, zero
violations (22 s, the exact pair-sum evaluator). Extends S128c86's `(24, 8)`.
Convergent evidence of the ladder's centrality: mac-mini S127's dyadic descent
chains pass exactly through `K_3(6) = 3/19`, and the deep well `{1..12,182}`
exposes a hidden `(1, 91)` K-shape one σ-level down.

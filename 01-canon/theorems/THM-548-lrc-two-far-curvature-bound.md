---
id: THM-548
title: The two-far curvature bound for LRC(14) — the mixed second difference I_B(u,v)=p0(B∪{u,v})−p0(B∪{u})−p0(B∪{v})+p0(B) has decorrelated limit Φ₂(B)=(2p₂(B)−p₁(B))/49 and a signed deviation controlled by the resonance frequency mu+nv (extends THM-546's one-far (6/49) bound to two far runners; the analytic complement to codex's HYP-2679 atlas)
status: STUB / namespace claim (mac-mini-2026-06-20-S3). Deriving the signed Fourier/Abel bound on I_B(u,v)−Φ₂(B); the limit Φ₂(B)=(2p₂−p₁)/49 is computed; the deviation Fourier form and the resonance-set reduction are in progress. Complements codex HYP-2679 (the exact two-far atlas).
source: mac-mini-2026-06-20-S3
depends_on:
  - THM-546   # one-far signed (6/49) bound
  - THM-547   # boundary-collar closure (bounded base)
  - HYP-2679  # codex's two-far curvature ledger (the exact atlas this bounds analytically)
  - THM-531   # scale-invariance (resonant far pairs reduce)
related:
  - HYP-2678  # true-wide Ruzsa/Plünnecke program
  - HYP-2657  # QR reality
  - HYP-2606  # covolume / relation-height
  - OPEN-Q-108
external: Lonely Runner Conjecture; boundary functions & curvilinear convergence (Kaczynski, Bagemihl, McMillan); Fatou nontangential limits.
---

# THM-548 — The two-far curvature bound (claim/derivation in progress)

## The object (codex HYP-2679)

`I_B(u,v) := p0(B∪{u,v}) − p0(B∪{u}) − p0(B∪{v}) + p0(B)` — the mixed second difference of
two far runners `u<v` over a base `B`. The far-element recursion's second order.

## The decorrelated limit (computed)

Let `A_j = {x : B misses exactly sector j}`, `A_{jj'} = {x : B misses exactly {j,j'}}`,
`p1(B)=Σ_j meas(A_j)`, `p2(B)=Σ_{j<j'} meas(A_{jj'})`. Then exactly
> `I_B(u,v) = −Σ_j ∫1_{A_j}(x)1_j(ux)1_j(vx)dx + Σ_{j<j'}∫1_{A_{jj'}}(x)[1_j(ux)1_{j'}(vx)+1_{j'}(ux)1_j(vx)]dx`.
As `u,v→∞` independently (2-D Weyl), each joint indicator decorrelates to `1/49`, giving the limit
> **`Φ₂(B) := lim I_B(u,v) = (2·p2(B) − p1(B)) / 49`.**

## The deviation = resonance (the crux)

Fourier-expanding (`1_j(y)=Σ_m ŝ_j(m)e(my)`):
> `I_B(u,v) − Φ₂(B) = Σ_{j,j'}(±) Σ_{(m,n)≠(0,0)} ŝ_j(m)ŝ_{j'}(n)·1̂_{A_{·}}(−(mu+nv))`.
The controlling frequency is **`mu+nv`**. For `(u,v)` with no small relation, `|mu+nv|` is large
for all small `(m,n)`, so by BV decay `|1̂_A(L)|≤#arcs/(π|L|)` the deviation `→0`. The deviation is
large only at **resonances** `mu+nv≈0` for small `(m,n)` — i.e. `u/v ≈ −n/m`, a small additive
relation between the two far runners. This is exactly:
- the **boundary-function** picture: `u,v` are two approach arcs; resonance = a shared asymptotic
  direction (Bagemihl ambiguous point); the resonant set is structured/small (the exceptional set);
- the **Freiman** picture: a small relation `mu+nv=O(1)` ⟹ `u,v` lie in a common low-dim GAP ⟹
  scale-invariance (THM-531) reduces the resonant pair to a bounded model.

## Plan to close true-wide (region III)

`p0(B∪{u,v}) = p0(B) + [one-far over bounded B, THM-547-style] + I_B(u,v)`, and
`I_B(u,v) = Φ₂(B) + deviation`. `Φ₂(B)` is a finite-base quantity (bounded `B`); the deviation is
`O(1/dist)` off resonance and Freiman-reducible on resonance. Target inequality: `p0(B∪{u,v})≤cap_k`
for all true-wide rows, splitting non-resonant (analytic bound) vs resonant (scale-invariance to a
finite check). Derivation + exact verification ongoing; merges with codex HYP-2679's atlas.

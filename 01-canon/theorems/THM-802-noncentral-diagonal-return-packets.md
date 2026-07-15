---
id: THM-802
title: Noncentral diagonal-return packets persist inside a five-core-safe interval
status: CLAIMED — exact all-height construction and proof identified; full proof and replay in progress
source: codex-2026-07-14-S10 reduced-holonomy continuation
depends_on:
  - THM-779   # normalized collision transition and prefix-legality equation
related: [THM-773, THM-778, THM-783, THM-786, THM-788, THM-794, HYP-6840]
---

# THM-802 — Noncentral diagonal-return packets

> **Namespace reservation.**  For `L=105H`, the eight-owner family
>
> ```text
> W=(L+1,L+22,L-20,L-6,L+15,L-13,L+8,2L+2)
> ```
>
> has `H` consecutive prefix-legal packets with owner word
> `(7,2,5,3,0,6,4,1,7)`.  Its count vector `(1,1,1,1,1,1,1,2)` is not the
> once-per-owner packet of THM-794, but every count times its inverse residue
> is `1 mod 7`, so the return is again diagonal.  The packet chain lies over
> `[2/15,1/7]`, which is safe for the five-core
> `P={1,2,11,12,13}`.  Exact wall-order inequalities, the collision-prefix
> certificate, tournament fingerprints, and an all-height replay are being
> written before this claim is promoted to `PROVED`.


---
id: HYP-3733
title: The covering-min as an ILP over the danger-circulant (set-cover on the breakpoint universe + primitivity constraints) PINS the primitive odd-n covering-min EXACTLY -- n=7:2/13, n=8:2/15, n=9:4/33 (all matching prior exact/exhaustive), n=11:3/31 (NEW), n=13:1/12 (<=, V=56); the margin M_prim(n)-1/n = 1/91, 1/120, 1/99, 2/341, 1/156 is IRREGULAR (= 1/(n(2n-1)) only at n=7,8, then deviates), confirming the covering-min is genuinely n-dependent. The set-cover (cover all breakpoints tau=k/d by radius-r danger arcs, min r) is correct ONLY over the full breakpoint universe d<=2V (the coarse single-modulus Z/m grid is WRONG -- it lets the solver pick Z/m-flat-but-spiky sets); its LP-dual is a PACKING of pairwise-incompatible lonely witnesses = an INDEPENDENT SET in the danger conflict graph, and the oriented danger circulant is a ROTATIONAL TOURNAMENT R_m whose Hamiltonian-path count is ODD (Redei: 15,175,3267 for m=5,7,9) and whose OCF H(R_m)=I(Omega,2) is the chromatic<->OCF bridge -- now COMPUTATIONAL. Four reframings of 'LRC parity = bipartiteness of C_n' (HYP-3729): (obs) the GEOMETRY (parity/C_n) is observer-invariant but the ARITHMETIC covering-min is observer-ANCHORED (changing observer breaks divisibility-covering, changes M); (translate) speed-shift = orbit ROTATION (breaks covering), t-shift = SHEAR; (all-n-points) = the view-obstruction / vertex-transitive C_n; (Ham paths) = the rotational tournament R_m + Redei + OCF
status: VERIFIED -- ILP (scipy.milp/HiGHS) pins primitive covering-min n=7,8,9,11 exactly (n=9 cross-checked vs exhaustive), n=13<=1/12 (V=56); margins exact; Redei odd Ham-paths verified m=5,7,9; observer/translate reframings verified. The OCF<->covering-min identity is a STRUCTURAL bridge (set-cover/packing duality + rotational tournament), not a proven equality. n=13 and larger are V-bounded (upper bounds if the true min needs speeds > V).
source: mac-mini-2026-06-30-S51
related:
  - HYP-3729  # LRC parity = bipartiteness of C_n (the result being reframed)
  - HYP-3727  # primitivity: full covering-min=1/n (non-primitive blocks), primitive >1/n (this pins the primitive values)
  - HYP-3725  # the construction is not the covering-min (the spread family is -- now pinned by ILP)
  - HYP-2566  # uniform looseness (primitive covering-min > 1/n) -- these exact values are evidence
references:
  - klein-2026-06-29-S34   # Ramanujan=Ihara-RH, metazeta (OCF tournament-side); Redei
  - opus-2026-06-30-S1     # even block 1/n (the non-primitive full-min the ILP finds WITHOUT primitivity)
results:
  - 04-computation/covering_min_ip_v2_macmini_20260630.py
  - 05-knowledge/results/covering_min_ip_primitive_macmini_20260630.out
  - 05-knowledge/results/covering_min_ip_v2_macmini_20260630.out
---

# HYP-3731 -- the covering-min ILP, the OCF bridge, and four reframings

Two asks: (B) implement the covering-min as an IP over danger-circulant independent sets (pin odd-n + make
the chromatic<->OCF bridge computational); (A) creatively reframe "LRC parity = bipartiteness of C_n"
(HYP-3729) -- change the observer, translate, copy to all n points, relate to Hamiltonian paths.

## (B) The ILP -- and the formulation that actually works
**Correct formulation.** `M(S) = max_t min_v ||vt||` is attained at a breakpoint `tau=k/d` with `d` a pairwise
sum/diff of speeds (`<= 2V` if speeds `<= V`). So: `M(S) <= r` iff for every `tau` in the **full breakpoint
universe** `{k/d : 1<=k<d<=2V}`, some `v in S` has `||v tau|| <= r` (the radius-`r` danger arcs cover the
universe). Binary-search the smallest feasible `r`; each feasibility is a **set-cover ILP** (scipy.milp/HiGHS):
`sum x_v = n-1`, divisibility-cover (`sum_{q|v} x_v >= 1`), **primitivity** (for each prime `p`,
`sum_{p does not divide v} x_v >= 1`, forcing `gcd=1`), and a cover row per `tau`.
*Pitfall (fixed):* gridding a single modulus `Z/m` is WRONG -- the solver returns sets flat on `Z/m` but
spiky at the true breakpoints (it reported `2/11` at n=9 vs the true `4/33`). The universe must be all `k/d`.

**Pinned primitive covering-min** (matches prior exact/exhaustive where known; n=9 cross-checked vs the
2M-set exhaustive):

| n | M_prim | witness t* | margin M-1/n | a worst-case set |
|---|--------|-----------|--------------|------------------|
| 7 | **2/13** | 4/13 | 1/91 | {1,6,7,8,9,15} |
| 8 | **2/15** | 11/15 | 1/120 | {2,6,7,8,10,13,14} |
| 9 | **4/33** | 29/33 | 1/99 | {1,4,5,6,7,11,32,36} |
| 11 | **3/31** | 29/31 | 2/341 | {2,6,8,9,10,11,13,14,17,19} |
| 13 | **<=1/12** | 17/36 | 1/156 | {3,6,..,33,26} (V=56) |

**Margin deviation (ask 1 of the prior turn, now exact):** `1/91, 1/120, 1/99, 2/341, 1/156` -- **irregular**;
it equals `1/(n(2n-1))` only at n=7,8 and deviates after (n=9 `1/99` not `1/153`; n=11 numerator 2). The
covering-min is genuinely `n`-dependent -- no clean closed form. (Drop the primitivity rows and the ILP
returns the *full* covering-min `1/n` via the non-primitive scaled blocks -- even block at n=8, `3*{1..8}` at
n=9 -- confirming HYP-3727's easy/hard split.)

**The chromatic<->OCF bridge, computational.** The set-cover's **LP-dual** is a packing of lonely witnesses
no single speed can simultaneously danger = an **independent set** in the danger conflict graph. The oriented
danger circulant `C_m(1,..,j)` is a **rotational tournament** `R_m`; its **Hamiltonian-path count is ODD**
(Redei -- verified `15, 175, 3267` for `m=5,7,9`), and its `H(R_m)=I(Omega,2)` is the OCF. So solving the
covering-min ILP and computing the OCF are independent-set computations on the *same* circulant/tournament --
the chromatic<->OCF bridge is now a concrete shared object (the danger circulant), not just an analogy.
(Structural bridge; the exact OCF=covering-min equality is not claimed.)

## (A) Four creative reframings of "LRC parity = bipartiteness of C_n"
- **Change the observer (same way).** The *geometry* (the equally-spaced orbit = `C_n`, and its
  bipartite/non-bipartite parity) is **observer-invariant** (vertex-transitive cycle). But the *arithmetic*
  covering-min is **observer-ANCHORED**: re-referencing to another runner subtracts its speed, which **breaks
  the divisibility-covering** (the new |speeds| omit some `q`) and changes `M` (even block from runner-2:
  `M=1/7`, and it no longer covers `7`). So the analogy persists at the level of the parity/geometry, but the
  worst-case *set and value* are anchored to the divisibility origin (observer 0). The worst case is the most
  symmetric (observer-blind) geometry meeting the least symmetric (origin-anchored) arithmetic.
- **Translate.** Translating the speeds (`v -> v+c`) **rotates** the orbit rigidly (circulant-invariant, `C_n`
  preserved) but breaks covering; translating the witness (`t -> t+c`) **shears** the orbit. The cycle
  geometry is translation-covariant; the covering arithmetic is not -- the same geometry/arithmetic tension.
- **Copy it to all n points.** Imposing the lonely condition at all `n` vertices of `C_n` is the
  **view-obstruction** reframe (Cusick): the `n` trajectories must miss the centre. For the extremal the `n`
  holes are vertex-transitive; the even/odd parity is whether they pair antipodally (bipartite) or not.
- **Hamiltonian paths.** The oriented danger circulant is the **rotational tournament** `R_m`; Redei's odd
  Hamiltonian-path count and the OCF `H(R_m)` tie the covering-min to the project's core tournament invariants
  (above). For odd `m`, `R_m` is a regular (Paley-like) tournament; for even `m`, none exists -- the same
  parity, now as tournament regularity.

## What it buys
A *correct, exact* computational handle on the primitive covering-min (the LRC hard core, THM-523/HYP-2566):
the ILP pins n=7,8,9,11 exactly and n=13 up to a speed bound, gives the exact (irregular) margin trajectory,
and -- via set-cover/packing duality and the rotational tournament -- makes the chromatic<->OCF bridge a
shared computational object. The four reframings locate exactly which structure is observer/translation
invariant (the `C_n` geometry/parity) and which is anchored (the divisibility-covering).

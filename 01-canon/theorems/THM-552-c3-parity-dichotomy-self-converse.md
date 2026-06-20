---
id: THM-552
title: The c3-parity dichotomy for phi-self-converse tournaments -- a grid-symmetric (half-tiling) tournament has c3 forced EVEN when n is even, and reaching both parities when n is odd, via the phi-invariant 3-set count (= the apex/diagonal fixed-vertex mechanism)
status: PROVED (anti-automorphism orbit argument); VERIFIED exhaustively over all grid-symmetric tilings n=3..8 (independent re-verification, 0 mismatches)
source: kind-pasteur-2026-06-20
depends_on:
  - THM-280   # grid reflection = converse+relabel (phi-self-converse <=> grid-symmetric tiling)
related:
  - THM-549   # mac-mini: half-tiling = complement-quotient fundamental domain
  - THM-550   # codex: half-tiling parity recurrence
  - THM-016   # Claim B alternating-sum even-odd split (DIFFERENT object; see Distinctions)
  - THM-017   # B(Li,Rj)=B(Lj,Ri) even-odd split
  - HYP-2687  # gnomon shells / even-odd dichotomy (this theorem is its proved core)
external: Moon, self-converse tournaments; anti-automorphisms.
---

# THM-552 — The c3-parity dichotomy for phi-self-converse tournaments

## Statement

Let `phi(i) = n+1-i`.  Let `T` be a `phi`-self-converse tournament on `{1,...,n}`
(equivalently, by THM-280, a GRID-SYMMETRIC tiling: `phi` is an anti-automorphism,
`T(phi u, phi v) = T(v,u)`).  Let `c3(T)` be the number of directed 3-cycles.  Then:

- **`n` even**  ⟹  `c3(T)` is **EVEN** (for every such `T`).
- **`n` odd**   ⟹  `c3(T) ≡ #{phi-invariant directed 3-cycles} (mod 2)`, and BOTH parities
  occur.  The `phi`-invariant 3-vertex-sets are exactly
  `{ (n+1)/2 , x , n+1-x }` for `x = 1, ..., (n-1)/2`, numbering `(n-1)/2 = d`
  (= the number of fixed diagonal cells of the half-tiling).

This is a genuine constraint: general (non-self-converse) tournaments on even `n` CAN have
odd `c3` (e.g. a 4-tournament with a single 3-cycle), but `phi`-self-converse ones cannot.

## Proof

Because `phi` is an anti-automorphism, it sends directed 3-cycles to directed 3-cycles:
a directed cycle `u->v->w->u` is reversed by arc-reversal to `u<-v<-w<-u`, i.e.
`u->w->v->u`, and relabeling by `phi` (which restores `T`) carries it to the directed
3-cycle on `{phi u, phi v, phi w}`.  So `phi` acts on the set of directed 3-cycles of `T`
as an involution.  Directed 3-cycles NOT fixed by `phi` come in 2-element orbits and
contribute an EVEN count; hence
```
        c3(T) ≡ #{ directed 3-cycles fixed by phi }   (mod 2).
```
A directed 3-cycle is `phi`-fixed **iff its vertex set is `phi`-invariant** (a `phi`-fixed
3-cycle must have `phi`-invariant support; conversely on a `phi`-invariant 3-set the unique
orientation respecting the anti-automorphism is fixed).  `phi` is an involution whose orbits
on vertices have size `<=2`; an ODD-size invariant set must contain a fixed point of `phi`.
`phi` has a fixed vertex iff `n` is **odd** (the median `(n+1)/2`).  Therefore:

- `n` even: `phi` is fixed-point-free, so there is **no** `phi`-invariant 3-set ⟹ the fixed
  3-cycle count is `0` ⟹ `c3(T)` is even.
- `n` odd: the `phi`-invariant 3-sets are exactly `{(n+1)/2, x, n+1-x}`, `x=1..(n-1)/2`,
  giving `(n-1)/2` of them; their cyclic-orientation count can be odd, so `c3(T)` can be odd.

∎

## Verification

`04-computation/half_tiling_gnomon_kpswf.py` and the independent
`04-computation/half_tiling_verify_contested_kps.py` (output
`05-knowledge/results/half_tiling_verify_contested_kps.out`): over ALL grid-symmetric tilings
(`2^half(n)` of them: 2,4,16,64,512,4096 for n=3..8) the c3 parities present are
`n=3:[0,1]`, `n=4:[0]`, `n=5:[0,1]`, `n=6:[0]`, `n=7:[0,1]`, `n=8:[0]` — even-only at even `n`,
both at odd `n` — and the per-tiling identity `c3 mod 2 = #phi-invariant-3-cycles mod 2`
holds with **0 mismatches** for n=3..8.

## Distinctions (avoid conflation)

- **NOT** the trivial complement-invariance `c3(T)=c3(T^op)` (mac-mini THM-549 / the 2x
  half-region computation): that says the *count* is fold-invariant; THM-552 is about its
  *parity* being pinned by the parity of `n`.
- **NOT** the THM-016/017 even-odd split, which concerns the alternating sum
  `sum_S (-1)^{|S|} ...` (Claim B).  Both ultimately stem from the converse/`phi` symmetry
  (THM-016/017 use the arc-reversal half), but the objects differ.

## Geometric source (half-tiling)

The `phi`-fixed vertex `(n+1)/2` exists only at odd `n`.  In the half-tiling it is the
apex/diagonal of the staircase in vertex-space; its presence is exactly why the odd half-region
is a centered square `k^2` (with a central fixed cell) while the even one is a pronic `k(k-1)`
with no center (THM-549/550, HYP-2687).  So the c3-parity dichotomy is the cycle-space face of
the same odd/even shape split, and a clean instance of the even-odd tournament distinction
(cf. THM-016/017) emerging from the converse fold.

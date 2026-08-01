---
source: kind-pasteur-2026-07-24-S157 (Opus 4.8)
status: RESULT (exact Z[omega] computation + S3 framework + reformulation). Attacks FC(3) via the triangle's S3
  symmetry with an EXACT Z[omega] engine. Finds three character-vanishing levers (C3-weight, sign-rep, invariant
  reduction) -- all Conway-Jones -- but shows S3 gives NO free reality/capacity boost (moment coeffs are genuinely
  complex). Exact D=2 capacity: the cyclic-weight family kills exactly 2 leaks and is FORCED to leave the 3rd
  nonzero (17 solutions of M3=M6=0, all with L(f^9)!=0). Tower TERMINATES at capacity=#params, like FC(2).
  Recasts FC(3) as a TRANSVERSALITY statement about the moment-leak equations. Higher-D open (needs the invariant
  reduction to be tractable).
tags: [factorial-conjecture, S3-symmetry, roots-of-unity, conway-jones, exact-arithmetic, capacity, transversality]
related: [kps-S154, kps-S156, THM-415, THM-2801]
---

# FC(3) via the triangle's S3 symmetry: exact Z[omega] capacity

Target: nonzero `f in C[X1,X2,X3]` with all `L(f^k)=0` (`L(X^a)=prod a_i!`). Cyclic-weight ansatz
(kps-S156): `f` in the `C3` eigenspace of weight `omega^2` -> `L(f^k)=0` for `3 !| k` automatically; solve the
`3|k` leaks `L(f^{3j})`.

## 1. The S3 framework -- three character-vanishing levers (all Conway-Jones)
`L` is **S3-invariant**, so `L(g)=L(Pi_triv g)` (trivial-isotypic projection). On the triangle `{X1,X2,X3}`:
1. **C3 rotation weight** (`rho=(123)`): `f` of weight `omega^2` => `L(f^k)=omega^{2k}L(f^k)=0` unless `3|k`.
   (Verified exact: `L(L1^3)=6`, `L(L1^6)=720`, `L1=X1+omega X2+omega^2 X3`.)
2. **Sign rep** (`tau=(12)`): the Vandermonde `Delta=(X1-X2)(X2-X3)(X3-X1)` is alternating, `L` symmetric, so
   **`L(Delta * g)=0` for every symmetric `g`** (verified: `L(Delta), L(Delta e1), L(Delta e1^2), L(Delta e1^3)=0`).
   But `Delta^2` is symmetric, `L(Delta^2)=24 != 0`: like C3, it kills a *subset* of powers, not the tower.
3. **Invariant reduction**: `L(f^{3j})` depends only on the S3-trivial part of `f^{3j}`, a polynomial in the two
   reflection invariants `theta_2=L1 barL1` (deg 2), `theta_3=L1^3+barL1^3` (deg 3). (The tractable route for D>=3.)

**All three are representation-theoretic (character-orthogonality / vanishing sums of roots of unity)** -- the same
Conway-Jones governor as the series thread (kps-S156) and LRC (THM-415). None closes the tower by itself.

## 2. No free reality/capacity boost from S3
Hoped: S3 forces `L(f^{3j})` real, halving equations. **False.** The exact moment polynomial (Z[omega] engine)
`M3=L(f^3)` for `f=L1+a P2+b Q2` has genuinely complex coefficients, e.g. terms `(36+36 omega) b^2`,
`(-72 omega) ab`, `(-432 omega) a^2 b`. So `L(f^{3j})` is complex for complex/general `(a,b)`; reality holds only
at the pure seed `a=b=0`. Each leak stays **one complex** equation (two real). No capacity doubling.

## 3. Exact D=2 capacity: the tower TERMINATES
`f=L1+a P2+b Q2` (2 complex params). Solving `M3=M6=0` (exact polynomials -> 4 real eqs, 4 real unknowns) gives
**17 nonzero solutions**, and at every one
```
|L(f^9)|/9!  ranges 0.63 .. 2.3e4   (nonzero on ALL solutions).
```
> **2 parameters kill 2 leaks (`j=1,2`) but are FORCED to leave the 3rd (`j=3, k=9`) nonzero: capacity = #params,
> tower terminates -- identical shape to the FC(2) control (kps-S156 sec 3).**
So the extra dimension gives **no** capacity advantage *at degree 2* (the deg-2 weight-`omega^2` basis in 3 vars
has the same size, 2, as the deg-2 antisymmetric basis in 2 vars). The dimensional difference is only asymptotic:
weight-`omega^2` corrections number `~d^3/6` in 3 vars vs `~d/2` in 2 vars -- more parameters at high degree, but
(so far) still `capacity = #params`.

## 4. Reformulation: FC(3) as a transversality statement
The pattern `capacity = #params` at every degree would say: the leak equations `L(f^{3j})=0` cut the parameter
space **transversally** -- each new parameter buys exactly one more killed leak, never two. Then no finite
cyclic-weight `f` kills the whole infinite tower, i.e. **FC(3) holds in the cyclic-weight family**.
> **Conjecture (transversality form of FC(3)):** for the cyclic-weight family of every degree `D`, the leak map
> `f |-> (L(f^3),L(f^6),...)` has capacity exactly `#params(D)`. D=2 verified.
A genuine FC(3) **counterexample** would be a *non-transversal coincidence* -- an `f` where the leaks collapse
below generic rank, killing more leaks than parameters. That is exactly a **period-rigidity failure** (KZ,
kps-S154): a cancellation beyond the character/Conway-Jones ones of Sec 1. So "FC(3) may be false" (friend) =
"the moment-leaks can be non-transversal"; "FC(3) true" = "they are always transversal."

## 5. Status (honest)
Delivered: exact Z[omega] engine (seed + M3 exact); the S3 three-lever framework; the refutation of the reality
shortcut; exact D=2 capacity (terminates, capacity=2). NOT resolved: D>=3 capacity -- the powers `f^{15}, f^{18}`
are too large for brute force; the **S3-invariant reduction (Sec 1.3)** is the route to push it, computing
`L(f^{3j})` in the 2-variable invariant ring `C[theta_2,theta_3]` instead of degree-`3jD` in 3 vars. That is the
concrete next computation. No FC(3) counterexample found; none excluded beyond D=2.

Files: `/tmp/{zw_engine,signrep,capacity}.py`. Builds on kps-S154/S156.

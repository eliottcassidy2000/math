---
id: HYP-3728
title: THE BASS FACTORIZATION IS THE GF(2) CUT(+)CYCLE SPLIT (the metazeta) -- Ihara zeta_G(u)^-1 = (1-u^2)^(r-1) . det(I-Au+Qu^2) splits EXACTLY into the CYCLE/even-graph half (1-u^2)^(r-1), r=cycle rank=dim cycle space, and the CUT/sandpile half det(I-Au+Qu^2) (at u=1 -> Laplacian L, det=spanning trees=sandpile group); verified K_4,K_5,C_5,K_3,3,Petersen. RAMANUJAN <=> Ihara-RH (nontrivial poles on |u|=1/sqrt(k-1)): K_4,K_5,C_5,Petersen Ramanujan, K_3,3 bipartite. THE IHARA/RAMANUJAN THREAD MEETS THE HEEGNER THREAD: the complete-graph metagraph K_n has Ihara cut-factor discriminant 9-4n, hitting the Heegner numbers -3,-7,-11,-19,-43 for n=3,4,5,7,13 (the n=4 Klein-four metagraph carries sqrt-7 = the apex); C_5 (pentagon) carries the GOLDEN RATIO (Fibonacci/H_2). The 'metazeta' = even-graph(cycle) (+) sandpile(cut) as one Euler product; the observer's 1 = the baseline, prime cycles accumulate above = the H-gradient
status: VERIFIED (Bass factorization splits into (1-u^2)^(r-1) cycle half + Laplacian-pencil cut half, small graphs; Ramanujan/Ihara-RH poles on |u|=1/sqrt(k-1); K_n disc=9-4n Heegner at n=3,4,5,7,13; C_5 golden ratio). Synthesis of the owner's Ihara/Bass/Ramanujan prompt with the project's GF(2) cut+cycle (tiling), the even-graph/sandpile avenues, and the Heegner thread.
source: klein-2026-06-29-S34
depends_on:
  - HYP-3730   # the Heegner / class-number-1 thread (this Ihara thread meets it)
related:
  - HYP-3726   # mac-mini: the floor margin = 1/hexagonal, Sum = ln4 (Borel-Cantelli budget)
  - THM-588    # metagraph algebraic connectivity = 4 (the spectral-gap / Ramanujan question)
  - HYP-3604   # the least-eigenvalue certificate (the spectral gap)
results:
  - 04-computation/bass_cut_cycle_ramanujan_klein.py
  - 05-knowledge/results/bass_cut_cycle_ramanujan_klein.out
---

# HYP-3728 — Bass = cut(+)cycle; Ihara/Ramanujan meets Heegner

## The Bass factorization IS the cut(+)cycle split (verified)
The Ihara zeta `zeta_G(u)^-1 = (1-u^2)^(r-1) . det(I - A u + Q u^2)` (`r=|E|-|V|+1` cycle rank, `Q=D-I`)
splits exactly into:
- **CYCLE / even-graph half** `(1-u^2)^(r-1)`: the exponent `r-1` is `dim(cycle space) - 1` -- the
  even-graph / `H^1` (the wiggly/cycle side of the project's GF(2) `E = Cut (+) Cycle`);
- **CUT / sandpile half** `det(I - A u + Q u^2)`: at `u=1`, `I - A + Q = D - A = L` (the Laplacian), whose
  reduced determinant is the number of spanning trees = the order of the sandpile group (the cut / score /
  matrix-tree side).
So **the owner's claim is exact: the Bass factorization of the Ihara zeta IS the GF(2) cut(+)cycle Hodge
split, made into one Euler product** -- the "metazeta." Verified `K_4, K_5, C_5, K_3,3, Petersen`.

## Ramanujan <=> Ihara-RH (verified)
A `k`-regular graph is Ramanujan (`|nontrivial A-eigenvalue| <= 2 sqrt(k-1)`) IFF its Ihara zeta's
nontrivial poles lie on `|u| = 1/sqrt(k-1)` (the RH analogue). Verified: `K_4` (poles `|u|=0.7071=1/sqrt2`),
`K_5` (`0.5774=1/sqrt3`), `C_5` (`1.0`), `Petersen` (`0.7071`) are Ramanujan; `K_3,3` is bipartite (the
`-k` eigenvalue, RH on `|u|=1/sqrt(k-1)` with the bipartite `+-`). The metagraph `G_n` (alg. connectivity 4,
THM-588) is the natural target -- its Ramanujan/Ihara-RH = its expansion = the least-eigenvalue certificate
(HYP-3604); `G_n` is irregular, so the irregular-Ramanujan (Lubotzky) form applies.

## The Ihara/Ramanujan thread MEETS the Heegner thread
The complete-graph metagraph `K_n` (`k=n-1`-regular, nontrivial eigenvalue `-1`) has Ihara cut-factor
`1 + u + (n-2)u^2`-type with **discriminant `9 - 4n`**:
```
n :  3   4   5   6    7   ...  13
9-4n: -3  -7  -11 -15  -19 ... -43
Heegner? Y   Y   Y   N    Y       Y
```
So `K_3, K_4, K_5, K_7, K_13` have **Heegner discriminants** `-3,-7,-11,-19,-43` (class number 1). The `n=4`
case (`K_4` = the Klein-four metagraph) carries `sqrt-7` -- the apex / Heegner `-7` / the conductor of
`X_0(14)`. And `C_5` (the pentagon) carries the **golden ratio** (`2cos(2pi/5) = 1/phi`, Fibonacci / `H_2`).
So the Ihara/Ramanujan structure of the small metagraphs hits exactly the project's number-theoretic
constants: the Heegner discriminants (`-3,-7,...`, HYP-3727) and the golden ratio (Fibonacci, HYP-3710).

## The H-gradient and the observer's 1
`zeta_G(u) = prod_{prime cycles}(1 - u^{length})^{-1}` and `log zeta(u) = sum_m N_m u^m/m` (`N_m` = closed
walks). The **observer's `1`** is the `u^0` baseline (the trivial / the transitive sink); the **prime
cycles** (closed geodesics = the metagraph's odd cycles / intransitivity) accumulate above it -- and that
accumulation IS the **H-gradient** (the metagraph's `H`-spectrum, which counts the cyclic structure). The
metazeta packages the even-graph (cycle), the sandpile (cut), and the `H`-gradient (prime cycles) as one
spectrum.

## The other threads (referenced)
- **ln 4 budget** (mac-mini HYP-3726): the floor margin `1/(n(2n-1))` sums to `2 ln 2 = ln 4` (log-det of the
  square Cartan), a finite Borel-Cantelli budget for the safe measure. The Ihara zeta connects: `ln 4` is
  also the central-binomial sum `sum C(2n,n)/(n 4^n)` (the `1/4` random-walk radius) -- the same `4` as the
  zeta's growth `N_m ~ k^m`.
- **odd/even-n covering-min** (post-MISTAKE-087): the covering-min is the SPREAD family (irregular: `n=7`
  `2/13`, `n=8` `2/15`, `n=9` `4/33`); even `n=2p` (`p<=7`) folds to the proven `LRC(p)` (the even-fold,
  `14->7`); the odd-`n` covering-min is the open frontier (mac-mini owns it). Not resolved here.

## Net
The Bass factorization IS the cut(+)cycle split (the metazeta = even-graph (+) sandpile Euler product);
Ramanujan = Ihara-RH; and the small-metagraph Ihara zetas carry the Heegner discriminants (`-3,-7,-11,...`)
and the golden ratio -- so the Ihara/Ramanujan thread meets the Heegner (HYP-3727) and Fibonacci (HYP-3710)
threads on the same metagraph. The H-gradient is the prime-cycle accumulation above the observer's baseline 1.

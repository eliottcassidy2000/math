---
id: THM-611
title: The tight-locus tower — local-max transfer + binding-pair divisibility. At a tight point of a family with hiding denominator q*=14m, the m-divisible sub-family E=m·U carries the whole local geometry (f_S = g_U(m·) near t*), the binding pair is divisible by m, and (via LRC≤13) U is STRICTLY loose for a primitive S — so a primitive tight family with q*>14 would be a loose U tightened by ≤12 non-E runners. Search (~200k q*=28 candidates, 0 hits) supports the CONFINEMENT primitive-tight ⟹ q*=14, reducing the rigidity to a finite mod-14 residue problem.
status: PROVED (Lemmas A,B,C + the LRC≤13 corollary — elementary; verified exactly). Lemma C proves the confinement for one tightener (m=2, |F|=1). CONJECTURE: the FULL confinement primitive-tight ⟹ q*=14 (Lemma C + search ~200k q*=28 candidates + AP perturbations, 0 primitive tight with q*>14; the multi-tightener case is not proved).
source: mac-mini-2026-07-03-S31
depends_on:
  - THM-610   # Lemma 2: tight ⟹ 14|q* (so q*=14m); this theorem is the fine structure of that
  - LRCUpTo13 # citation node: LRC proven for ≤13 runners — used to force U strictly loose
related:
  - HYP-4060  # kps-S36: tight-locus mechanism (exact-landing ⟹ dilated AP). Lemma B is its RIGOROUS local form; confinement sharpens the "principal branch" to a theorem-conjecture.
  - HYP-2913  # three-gap/Steinhaus + ±units-covering necessary condition; confinement ⟹ their mod-14 setting
  - HYP-2561  # tight-locus finiteness (inf L=1/1260, locus={AP,GW}); this is a geometric route to it
  - HYP-2909  # apex binding-pair (Lemma A is the sharpened, divisibility-carrying form)
results:
  - 04-computation/tight_locus_geometry_macmini_20260703.py
  - 04-computation/tight_q28_search_macmini_20260703.py
  - 05-knowledge/results/tight_locus_geometry_macmini_20260703.out
  - 05-knowledge/results/tight_q28_search_macmini_20260703.out
external: Lonely Runner Conjecture n=14; LRC proven ≤13 (Sungkawichai–Trakulthongchai arXiv:2604.23906); Goddyn–Wong Integers 6 (2006) #A38.
---

# THM-611 — the tight-locus tower (local-max transfer)

**Setup.** `n=14`; `S` a primitive 13-set, **tight**: `M(S)=max_t min_{v∈S}||vt|| = 1/14`, attained at
`t*=a/q*` (lowest terms). By THM-610 (Lemma 2), `14 | q*`; write `q* = 14m`, `m≥1`. Residues
`r_i = v_i·a mod q* ∈ [m, 13m]` (circular distance `≥ m`), with the min `=m`.
Let **`E = {v∈S : m|v}`** (the `m`-divisible runners) and **`F = S∖E`** (the rest). Write `E = m·U`.

## Lemma A (binding pair, divisibility-carrying). *sharpens HYP-2909*
There is a **binding pair** `(i,j)`: `v_i a ≡ m`, `v_j a ≡ −m (mod q*)`. Hence `q* | (v_i+v_j)`,
and — since `gcd(a,q*)=1` and `m ≡ 0 (mod m)` — **`m | v_i` and `m | v_j`**. More generally every
runner at exact binding distance (`r_k∈{m,13m}`) is divisible by `m`, so the **binding pair ⊆ E**.
*(Verified: AP `q*=14,m=1` pair `(1,13)`; even block `q*=28,m=2` pair `(2,26)`, both even.)*

**Proof.** `t*` is a local max of `f_S=min_v||vt||`, a min of tents; for a local max the active
(distance-`1/14`) set must contain a runner with `+` slope (phase `=m/q*`, i.e. `r=m`) and one with
`−` slope (phase `=1−m/q*`, i.e. `r=13m`). Those give `v_i a≡m`, `v_j a≡−m`; add/reduce mod `m`. ∎

## Lemma B (local-max transfer). *the tower step*
In a neighborhood of `t*`, the full view equals the `m`-divisible sub-family's view, rescaled:
> **`f_S(t) = g_U(m·t)`** for `t` near `t*`, where `g_U(τ)=min_{u∈U}||uτ||` is the view of `U=E/m`.

Consequently `a/14 = m·t*` is a **local max of `g_U` at height `1/14`**, and `M(U) ≥ 1/14`.

**Proof.** At `t*` the `F`-runners are *strictly* safe (`r_k∉{m,13m}` ⟹ distance `>1/14`), so by
continuity they stay `>1/14` on a neighborhood and never attain the min there. Thus locally
`f_S = min_{v∈E}||vt|| = min_{u∈U}||(mu)t|| = min_u||u(mt)|| = g_U(mt)`. ∎
*(Verified exactly on the even block `E=2·{1..13}`, `U={1..13}`: `f_E(t)=g_U(2t)` at all sampled `t`.)*

## Corollary (LRC≤13 forces the tightener picture).
If `S` is **primitive** and `m≥2`, then not all `v_i` are divisible by `m`, so `|E|=e≤12`, i.e.
`U` has `≤12` runners. By **LRC(≤13)** (citation), `M(U) ≥ 1/(e+1) ≥ 1/13 > 1/14`: **`U` is strictly
loose.** So `a/14` is a *non-global* local max of `g_U` (height `1/14 < M(U)`), and the `13−e ≥ 1`
runners of `F` must **suppress `U`'s entire super-`1/14` region** (in the `τ=mt` coordinate) down to
the global value `1/14`. A primitive tight family with `q*>14` is therefore precisely *a strictly-loose
`≤12`-runner set `U`, dilated by `m`, re-tightened to `1/14` by `≤12` non-`m`-divisible runners.*

## Lemma C (shift obstruction — the confinement, PROVED for one tightener). *`m=2`, `|F|=1` is impossible.*
Let `S` be primitive tight with `m=2` (`q*=28`) and exactly one non-`E` runner `F={w}` (`w` odd, `e=12`).
By the Corollary `U=E/2` is strictly loose, so `R={t: g_U(2t)>1/14}` has positive measure; `R` is
`(+1/2)`-invariant (`g_U(2(t+1/2))=g_U(2t+1)=g_U(2t)`). Tightness `M(S)=1/14` forces `w` to cover `R`
up to points: for a.e. `t∈R`, `||wt||<1/14`. But `w` is **odd**, so `||w(t+1/2)|| = ||wt+1/2|| =
1/2−||wt|| > 1/2−1/14 = 6/14` — `w` is *strictly safe* at `t+1/2`. Since `t∈R ⟹ t+1/2∈R`, the point
`t+1/2` is uncovered: `f_S(t+1/2)=min(g_U(2(t+1/2))>1/14,\ ||w(t+1/2)||>1/14)>1/14`, contradicting
`M(S)=1/14`. ∎ *(Verified: `E=2·{1..12}` + any single odd `w∈{1,…,299}` gives 0 tight families.)*

The same shift `t↦t+1/m` kills a lone tightener for any `m` with `||{w/m}||>1/7` (all `m=2`; and
`m≤7` unless `w≡±1 (mod m)`) — the obstruction is that one off-grid runner cannot cover a
`(1/m)`-periodic region it is forced to shift off. Two or more tighteners can in principle split the
`σ`-pairs, which is why the general confinement remains a conjecture (below).

## Confinement conjecture (Lemma C + search): *primitive tight ⟹ `q*=14`.*
The tightener picture above is never realized: **no primitive tight family has `q*>14`.**
Evidence: exact-`M` searches found **0** primitive tight families with `q*>14` over
- AP single/double perturbations (`~2.4k`), and
- **`q*=28`-targeted**: `2.6k` single- and `138k` double-primitivizations of the even block, plus
  `60k` random mixed-parity families — **0** primitive tight (the even block's tightness does not
  survive making it primitive). The only `q*=28` tight family is the **imprimitive** even block
  (`e=13`, `U=`AP, `M(U)=1/14` — the boundary case where the corollary's strict inequality lapses).

## Why this matters — the rigidity becomes finite mod-14 combinatorics
Granting confinement, every primitive tight family sits at `q*=14`: its phases lie on the **14th-root
grid** `{k/14}`, so it is a multiset of residues `R⊆{1,…,13}` (13 runners). The rigidity
`M=1/14 ⟹ {AP,GW}` then splits into:
1. **the mod-14 shell** (which residue multisets `R` admit a globally-tight speed lift) — a FINITE
   problem: `R` must cover the `±units {1,3,5,9,11,13}` (HYP-2913) and, by the three-gap conjecture
   `g(14)≤3` (HYP-2913), be a `≤3`-gap pattern ⟹ `R∈{AP-residues, GW-residues}`;
2. **the lift rigidity** (given `R`, which integer speeds `v_i≡r_i (mod 14)` stay globally tight) —
   controlled by Lemma B applied at `m=1`.

So the geometric skeleton of the open core is: **THM-610 (14|q*) → THM-611 confinement (q*=14) →
finite mod-14 shell + three-gap `g(14)≤3`**. The two remaining gaps are the confinement (this note's
conjecture) and `g(14)≤3` (HYP-2913) — both now sharply isolated on the 14th-root grid.

## Not claimed
Confinement is not proved (it is an existence/non-existence statement about the odd-tightener
construction, not reducible to LRC≤13 by the corollary alone). `g(14)≤3` remains open (HYP-2913).
This theorem contributes the proved tower lemmas (A,B + corollary) and the confinement reduction.

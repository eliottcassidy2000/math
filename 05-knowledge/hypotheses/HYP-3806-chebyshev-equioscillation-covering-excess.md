---
id: HYP-3806
title: CHEBYSHEV EQUIOSCILLATION + COVERING EXCESS + RIGIDITY -- creative abstract lenses synthesized and applied to the covering-min proof. (I) The covering-min M(S)=max_t min_v||vt|| is a CHEBYSHEV MINIMAX; its extremizer (the construction) is pinned by its ALTERNATION SET = the runners binding at t*. VERIFIED (exact, n=7..14): the alternation set is EXACTLY {1, n(n-1)} = {v : v ≡ ±1 mod Phi6}, multiplicity 2 for ALL n; and the killer n(n-1) ≡ -1 mod Phi6 IS the v≡-1 binding runner (the killer identity = the equioscillation). (II) COVERING EXCESS (a unifying invariant across the repo's two halves): M_C - 1/n = (n-1)/(n*Phi6) EXACTLY -- the price a COVERING constraint pays above the free LRC floor 1/n; the numerator (n-1) = the DROPPED speed = the CF partial quotient of t*=n/Phi6=[0;n-1,n]. Tournament analogue: flip-rank rho(n) - ceil(log2|G_n|) (HYP-3803/3804). (III) RIGIDITY: t* and 1-t* are the ONLY global maximizers (alternation length 2 = the MINIMAL possible), so the extremizer is an ISOLATED corner of the minimax = the deep-well isolation (mac-mini). SYNTHESIS: the alternation length (2) = the #atoms of the lonely measure (2) = the OPUC Verblunsky termination |alpha_1|=1 (2 atoms) = the #global maximizers (2) = the two solutions of v≡±1 mod Phi6 -- FIVE "2"s that are ONE fact, all carried by the killer identity n(n-1)≡-1 mod Phi6. APPLIED: the covering-min lower bound becomes a Chebyshev/LP-DUALITY problem with a 2-POINT-supported dual (on {1,killer} at {t*,1-t*}); a beater must break a length-2 alternation forced by v≡±1 mod Phi6 while covering all q<=n (THM-523) -- OPEN-Q-108 in Chebyshev/LP-dual form
status: MIXED (exact verification + clean synthesis + a proof TEMPLATE, not a proof). VERIFIED (exact Fraction, n=7..14): alternation set = {1,n(n-1)} = {v≡±1 mod Phi6} mult 2; M_C-1/n=(n-1)/(n*Phi6) exact; global maximizers = {t*,1-t*} only (grid N=2e6, n=7,10,14, rigid). The Chebyshev/equioscillation, covering-excess, and rigidity are established/elementary facts, here identified and unified. HONEST: this REFRAMES the covering-min lower bound as a 2-point Chebyshev/LP-dual construction and unifies the extremal structure; it does NOT close OPEN-Q-108 (no beater bound proved) -- a proof template + a synthesis, applied but not completed.
source: klein-2026-07-01-S73
depends_on:
  - HYP-3800   # S68: phase-residue p(w)=nw mod Phi6, killer identity, [0;n-1,n] (the arithmetic of the alternation)
  - HYP-3801   # S69: OPUC/Verblunsky (|alpha_1|=1 = 2 atoms = the alternation length)
related:
  - HYP-3803   # tournament flip-rank; the covering-excess analogue rho-ceil(log2)
  - HYP-3804   # packing/covering asymmetry (covering excess in the tournament cube)
  - HYP-3245   # equioscillation/autocorrelation atlas (prior equioscillation work)
  - HYP-3805   # opus-S15: flip-rank(7)=12 (tournament covering-excess 0,0,0,1,3 for n=3..7; obstruction = the Paley heptagon = LRC extremal object) -- reinforces "covering >> free"
  - HYP-3768   # the Dedekind/looseness M_C-1/n (= covering excess, now identified)
  - OPEN-Q-108 # the multi-far / no-beater bound (this is its Chebyshev/LP-dual form)
  - HYP-3715   # t*=n/Phi6 hexagonal binding
results:
  - 04-computation/chebyshev_covering_excess_klein.py
  - 05-knowledge/results/chebyshev_covering_excess_klein.out
---

# HYP-3806 — Chebyshev equioscillation, covering excess, rigidity: abstract lenses applied

## (I) The covering-min is a Chebyshev minimax; the extremizer equioscillates on 2 runners
`M(S) = max_t min_{v in S} ||vt||` is a **minimax** (a Chebyshev/best-approximation object). By Chebyshev
theory the extremizer is characterized by its **alternation set** — the constraints active at the optimum.
For the construction, the runners binding at `t*` (those with `||v t*|| = M_C`) are, **exactly and for all
n (verified n=7..14)**:
> **alternation set = `{1, n(n-1)}` = `{v : v ≡ ±1 (mod Phi6)}`, multiplicity 2.**
The condition `||v t*|| = M_C` is `n v ≡ ±n (mod Phi6)` i.e. `v ≡ ±1 (mod Phi6)`; the two solutions in range
are `v=1` and `v = Phi6-1 = n(n-1)` — **the killer, via the killer identity `n(n-1) ≡ -1 (mod Phi6)`.** So
the equioscillation set and the killer identity are the same fact.

## (II) Covering excess — a unifying invariant across the project's two halves
> **`M_C - 1/n = (n-1)/(n*Phi6)`** (exact, all n).
This is the **covering excess**: how far the covering-min sits *above* the free LRC floor `1/n` — the price
the covering constraint charges. The numerator `n-1` is the **dropped speed** (the CF partial quotient of
`t* = n/Phi6 = [0; n-1, n]`); the denominator is `n*Phi6`. (This equals the S56/HYP-3768 "looseness" /
Dedekind value `13/2562` at `n=14`.) The **tournament analogue** is the flip-rank excess `rho(n) -
ceil(log2|G_n|)` (HYP-3803/3804): in both halves a *covering* constraint forces the extremal value above
the *free/packing* floor, by an excess set by the symmetry/arithmetic.

## (III) Rigidity — alternation length 2 is minimal, so the extremizer is isolated
The global maximizers of `min_v ||vt||` are **only `t*` and `1-t*`** (verified, grid `N=2e6`, `n=7,10,14`).
Alternation length 2 is the *minimal* nondegenerate length, so the extremizer is an **isolated corner** of
the minimax — the potential-theoretic / deep-well isolation (mac-mini's 5001-set finding). A minimal
alternation is a rigid, corner solution: no continuous deformation preserves optimality.

## The synthesis: five "2"s are one fact
The number **2** appears as: the **alternation length** (2 binding runners); the **#atoms** of the lonely
measure (`t*, 1-t*`); the **OPUC Verblunsky termination** `|alpha_1| = 1` (a 2-atom measure, HYP-3801); the
**#global maximizers** (2); and the **two solutions of `v ≡ ±1 (mod Phi6)`**. These are not five
coincidences but one structure — the length-2 Chebyshev alternation at the `Phi6`-binding — all carried by
the **killer identity `n(n-1) ≡ -1 (mod Phi6)`**.

## Applied to the proof (the template)
The covering-min lower bound "every covering set has `M >= n/Phi6`" is, in this frame, a **Chebyshev /
LP-duality** statement with a **2-point-supported dual** (on the runners `{1, killer}` at the times
`{t*, 1-t*}`). Concretely, a **beater** (`M(S') < M_C`) must **break a length-2 alternation** that is
*forced* to `v ≡ ±1 (mod Phi6)` — while simultaneously **covering** every `q <= n` (THM-523). The two
obligations conflict on the `Phi6` phase-lattice (S68): covering pins the reachable binding to `Phi6`
(HYP-3802's "covering forces `q*: n -> Phi6`"), and the alternation there is rigidly length-2. Closing this
— constructing the 2-point dual certificate, or proving the forced alternation cannot be broken under the
covering constraint — is **OPEN-Q-108 in Chebyshev/LP-dual form**. This is the sharpest reframe yet of the
open crux: a *finite, 2-dimensional* dual object on a *finite* phase-lattice.

## Honest status
All facts verified exactly (alternation set, covering-excess identity, rigidity). The Chebyshev / covering-
excess / rigidity lenses are established mathematics, here identified, unified, and applied to reframe the
covering-min lower bound. This is a proof TEMPLATE + a structural synthesis — it does NOT close OPEN-Q-108.

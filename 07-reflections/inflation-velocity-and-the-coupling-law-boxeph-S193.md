# Inflation velocity and the coupling law: why WOWII witnesses work, and the strong-core reduction

*boxeph-2026-07-21-S193. Companion to klein-S395 (the-wowii-103-refutation…), opus-S437
(inflation-decoupling-counterexamples…), kind-pasteur THM-1845, klein THM-1850. Object: THM-1855.*

## The question behind the WOWII-103 refutation

WOWII-103 was `α(G) ≤ ⌊b(G) − ln(ecc_avg(G))⌋`, killed by a triangle with four leaves on two
of its vertices (α=9, b=10, ecc_avg=30/11 ⇒ ⌊8.997⌋=8 < 9). The **mechanism** (death-star-S77 and
I both isolated it): the leaf family pins `b − α = 1` *constant*, and pushes the weakly-coupled
correction `ecc_avg` **just past e** — `30/11 = 2.7273 > 2.71828 = e` — so `ln > 1` drops the floor
by one. The refutation is not luck: it is **velocity arbitrage**. Under "add a leaf," α and b move
in lockstep (both +1), so their gap is frozen; the *correction* moves independently. A conjecture
dies exactly when the correction's velocity is decoupled from the primary gap and can be driven
across an integer floor.

This reflection turns that observation into a **computable law for tournaments** (THM-1855): give
every iso-invariant a *velocity* under a fixed set of inflation operations; the velocity vectors
then (a) **predict** which conjectured inequalities are fragile and hand you the witness, and
(b) **prove** the survivors by reducing them to the strongly-connected core.

## The engine: order-join is an invariant algebra

Write `T = T₁ ▷ T₂` for the **order-join** (every vertex of T₁ beats every vertex of T₂). Because
all cross-arcs point one way, nothing cyclic or Hamiltonian crosses the cut, and (verified exactly
on all iso-class pairs with n₁+n₂ ≤ 7):

| invariant | behaviour under `▷` | one-line reason |
|---|---|---|
| `c3` (# 3-cycles) | **additive** `c3₁+c3₂` | a 3-cycle needs two arcs across the cut, impossible |
| `tr` (largest transitive subtournament) | **additive** `tr₁+tr₂` | concatenate the two chains |
| `scc` (# strong components) | **additive** | the cut never merges components |
| `n` | additive | — |
| **`H` (Rédei Hamiltonian-path count)** | **MULTIPLICATIVE** `H₁·H₂` | a Ham path traverses all of T₁, then all of T₂ |

The multiplicativity of the Hamiltonian-path count under order-join is the clean structural fact
that makes the whole method run.

**Atoms.** Every tournament factors uniquely as an ordered join of **strongly connected**
tournaments (its condensation is transitive — standard). So the strong tournaments are the
`▷`-atoms. Their iso-class counts: `n=3:1, 4:1, 5:6, 6:35, 7:353` (396 strong classes for n≤7).

## The coupling law (THM-1855)

Call an inequality `Φ` **join-monotone** if `Φ(T₁) ∧ Φ(T₂) ⟹ Φ(T₁ ▷ T₂)`.

> **Reduction theorem.** A join-monotone tournament-invariant inequality holds for *all*
> tournaments **iff** it holds for all *strongly connected* tournaments.

Proof: induct on the `▷`-factorization. The content of any join-monotone conjecture lives entirely
in the irreducible (strong) core — the same "reducible = block-triangular, SCCs are the atoms"
philosophy as THM-1830.

**Velocity corollary (single-vertex joins).** `D+` (add a source/dominator) `= v ▷ T` and `D−`
(add a sink/dominated) `= T ▷ v`. Their velocity table (constant over all 530 classes n≤7):

```
            Δc3   ΔH   Δtr   Δscc   Δsmin   Δsmax
  D+ (src)   0     0    +1    +1      0      var
  D- (sink)  0     0    +1    +1     var     +1
  complement fixes {c3, H, tr, scc, srange}
```

The decisive line: **under both D±, `tr` climbs by +1 while `c3` and `H` stay frozen.** `tr` is
*decoupled upward* from `{c3, H}`. This is the tournament version of "leaves pump α while ecc is
free."

**Fragility predictor.** `Φ: LHS ≤ RHS` is *inflation-fragile* if some operation `O` has
`Δ_O(LHS) > Δ_O(RHS)`; the witness is the `O`-orbit of a tight base tournament, iterated to the
first floor crossing. Because `tr` (and `srange`) outrun `c3, H` under `D±`, **any bound of the
shape `[tr- or srange-side] ≤ [function of frozen c3, H]` is fragile**, and its witness family is
the 3-cycle atom + a stack of transitive singletons — precisely THM-1830, the tournament
triangle+leaves.

## What the law does to the fleet's conjectures (all confirmed exhaustively, n≤7)

- **kind-pasteur's `srange ≤ tr`** — the engine's headline refutation. *Predicted fragile*
  (`srange` velocity up to +3 under `D−`, `tr` only +1 ⇒ decoupled); breaks at n=7
  (witness c3=4, tr=5, srange=6). My two **repairs are join-monotone survivors**, verified n≤7:
  `srange ≤ tr + c3` and `srange ≤ 2(tr−1)`. Both reduce to the strong core.
- **kind-pasteur's open `c3 ≤ H`** — *join-monotone* (airtight: if `H=1` then transitive so
  `c3=0`; else `H₁,H₂≥2 ⟹ H₁H₂ ≥ H₁+H₂ ≥ c3₁+c3₂`). ⇒ **reduces to the strong core**, and is
  verified on all 396 strong classes n≤7. kps's candidate is now a *strong-tournament* statement.
- **klein's `dom + tr ≤ n+1` (THM-1850)** — join-monotone, re-derived by velocity: under `D−` it
  is tight-preserving (`Δdom=0, Δtr=+1, Δ(n+1)=+1`), matching klein's "tight at the transitive
  tournament." Reduces to the strong core.
- **kind-pasteur sandwich `n − c3 ≤ tr ≤ smax+1` (THM-1845)** — both faces join-monotone; the
  lower face is `D±`-tight-preserving (`n−c3` and `tr` both +1). Reduces to the strong core.

## WOWII-103 proper: king-eccentricity straddling e

The faithful directed port uses `tr` (the α-analog) and **average king-eccentricity** (ecc(v) =
max out-distance; finite exactly on strong tournaments). The naïve `tr ≤ ⌊(n−1) − ln ecc_avg⌋`
breaks immediately (its RHS is simply too loose — the same triviality as klein's directed-103 (G)).
The *interesting* object is the **transcendental threshold itself**: among strong tournaments n≤7,
the achievable avg king-eccentricity straddling e is

- closest **below** e: **19/7 = 2.71429** (strong tournaments at n=7),
- closest **above** e: **17/6 = 2.83333** (a strong tournament at n=6).

So a WOWII-103-shaped bound with a `−ln(ecc_avg)` correction tips **between two explicit strong
tournaments exactly at e** — the tournament analog of `30/11` straddling e. This is HYP-8641's
concrete anchor; calibrating a *tight* directed-103 whose only slack is the `ln`-term (so the
straddle is decisive) is the open follow-up.

## The meta-point

The fleet already had *generate → exhaustive-check → Lean* (kps, klein, death-star). The missing
piece was **why**: a conjecture is fragile precisely where two invariants it couples are decoupled
by an inflation operation, and it is provable precisely when it is join-monotone (then its whole
content sits in the strong core). Velocity is the bridge from blind enumeration to prediction, and
from prediction to a reduction proof. Related threads: [[the-wowii-103-refutation-and-what-it-lends-the-repo-klein-S395]],
[[inflation-decoupling-counterexamples-the-wowii-motif]], THM-1830, THM-1845, THM-1850.

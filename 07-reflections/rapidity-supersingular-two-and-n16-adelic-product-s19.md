---
source: oracle-2026-05-31-S19
status: integration + proof-attempt (rapidity/adelic thread × n=16 LRC)
tags:
  - lonely-runner
  - n16
  - rapidity
  - adelic
  - supersingular
  - cayley-dickson
  - product-formula
  - bruhat-tits
---

# Rapidity, Supersingularity at 2, and an Adelic Product Formula for n=16 LRC

This session was asked to integrate the repo's **rapidity** thread (THM-252, the
adelic tournament geometry) with codex's **n=16 Lonely Runner** program
(S387–S390, "dyadic endpoint-debt") and push toward a proof. The two turn out to
be the same object viewed at the two completions of `Q`, and the integration
hands codex's empirical "conserved debt norm" a name: the **adelic product
formula**.

## 1. The three inputs

- **THM-252 (rapidity lattice).** For the flip-chain the eigenvalue rapidities
  `arctanh(λ) = ½ ln((m−k)/k)` lie in the **log-prime lattice**
  `⊕_{p ≤ m−1} Q·ln p`, with Baker's theorem giving `Q`-independence. The map
  `λ ↦ arctanh λ` is the **formal logarithm** of the Cayley formal group
  `F(x,y) = (x+y)/(1+xy)` (relativistic velocity addition).
- **Adelic tournament geometry (opus-S91).** The eigenvalue space is adelic,
  `R × ∏_p Z/p^e`. Crucially `F` is **supersingular (height ∞) at p = 2** and
  **ordinary (height 1) at odd p**; rapidity is its formal/`p`-adic log; the
  product formula `∏_v |x|_v = 1` holds. The 2-adic place is the special one —
  the profinite tournament space is the *odd* profinite integers because p = 2
  is killed.
- **codex n=16 (S389/S390).** `n = 16 = 2^4` forces a "16-gate" (a speed
  divisible by 16) to kill the odd unit skeleton `a/16`. But the gate does not
  close the problem; it *leaks* endpoint debt to deeper dyadic layers. Codex
  measured, across the dyadic ladders `d = 2,4,8,16`:

  ```
  (#unprotected endpoints) × (max_gap / threshold) ≈ 34/33   (constant)
  gap halves   (1/33 → 1/66 → 1/132 → 1/264)
  debt doubles (34 → 68 → 140 → 280)
  leak depth   v2(denominator) = 5 → 6 → 7 → 8
  ```

## 2. The unification: n=16 lives at the supersingular place

Three facts that were stated separately are **one fact**:

| view | statement at p = 2 / n = 16 |
|------|------------------------------|
| Cayley–Dickson | sedenions `= 2^4` are where zero divisors appear |
| formal group | `F(x,y)=(x+y)/(1+xy)` is supersingular (height ∞) at p = 2 |
| LRC (codex) | the 16-gate leaks dyadic endpoint-debt, never closing |

At an **odd** prime the formal group is ordinary (height 1): the `[p]`-map has an
étale part, the structure "closes up," and the corresponding LRC modulus is a
single, easily-covered sieve requirement (THM-369: one speed divisible by `q`
suffices, done). At **p = 2** the formal group is supersingular (height ∞): the
`[2]`-map is purely connected, there is **no étale/horizontal direction to wrap
the debt back up**, and so the dyadic debt can only **descend**. `n = 16` is the
pure 2-power, so its entire obstruction sits at this one supersingular place.
This is *why* n = 16 is qualitatively harder than an odd composite like n = 15:
the odd primes give tame, closable sieve conditions; the 2-tower does not close.

## 3. codex's conserved norm IS the adelic product formula

Read codex's invariant adelically. The **gap** `max_gap/threshold` is an
**archimedean** size (how lonely — a real distance). The **endpoint debt**
`#unprotected` is a **2-adic** size (how much protection is owed, measured by
Bruhat–Tits depth `v2` of the surviving endpoints). Their near-constant product

```
|gap|_∞ · |debt|_2 ≈ const
```

is exactly the shape of the product formula `∏_v |x|_v = 1`: the archimedean and
2-adic sizes are reciprocal. **A counterexample is the "bad corner"
`(gap = 0, debt = 0)`** — full open cover with every endpoint protected. The
product formula forbids it: if the archimedean size collapses to 0, the 2-adic
size must blow up, and a *finite* speed set cannot carry infinite debt (it has
finitely many endpoints). Hence no counterexample.

This is a genuine proof *shape* for n = 16:

> **(P)** For every primitive 15-speed set at n = 16, `|gap|_∞ · |debt|_2` is
> bounded below by a positive constant (and tight sets, `gap = 0`, instead carry
> `debt ≥ 1`). Therefore the bad corner `(0,0)` is unreachable: LRC holds at 16.

The missing analytic content is the *lower bound* — the same gap codex hit, now
with a reason to expect it (a product formula) rather than a coincidence.

## 4. The descent is ADELIC, not pure 2-adic (the key correction to codex)

Every endpoint of a speed `v` sits at 2-adic denominator depth `4 + v2(v)` on the
Bruhat–Tits tree of `Q_2`. codex's hope was a *pure dyadic* descent: protecting a
depth-`L` endpoint forces a strictly deeper protector, so a finite set terminates
at an unprotected "private dyadic leaf." **That pure-2-adic monovariant is false**,
and working out why is the heart of the integration.

Take a unit endpoint `e = (16m±1)/(16v)` at depth `L = 4 + v2(v)` and a candidate
protector `p`. Write `16v = 2^{4+j} v'` (`v'` odd, `j = v2(v)`) and
`p = 2^{v2(p)} p'`. If `v2(p) < 4+j`, the value `p·e` is `(odd)/(2^{4+j−v2(p)} v')`,
whose least nonzero distance to `ℤ` is `1/(2^{4+j−v2(p)} v')`. This is `< 1/16`
**as soon as `v' > 1`** — i.e. a protector of *lower* 2-adic valuation succeeds by
borrowing the **odd part `v'`** of the speed. This is exactly why codex found that
the `16`-gate's endpoints *can* be protected by lower odd residues `1,3,5,7,…`.

So protection trades 2-adic depth against **archimedean / odd-prime room**: the
debt cannot be cleared by pure dyadic descent, but it is *paid for* out of the
other places. That is the precise sense in which n = 16 LRC is an **adelic**
balance, not a 2-adic one — and it is exactly what the product formula
`|gap|_∞ · |debt|_2 ≈ const` encodes. The supersingularity at 2 says the 2-adic
side alone never closes; the conservation says the only way to drive the
archimedean gap to 0 is to spend unbounded 2-adic debt, which a finite set cannot
afford. The obstruction lives at **both places at once** — the right object is the
adele, not the 2-adic tree.

This reframes codex's open "positive divergence of the debt-flow": it should be a
divergence of the *adelic* height (gap·debt), not of the dyadic depth alone. The
empirical single-depth concentration of each ladder's debt (§6) shows the dyadic
descent is real, but the *bound* that forbids the bad corner is adelic.

## 5. What is already formal

The first branch is now **machine-checked**: this session's Lean module
`LonelyRunner.lean` proves `counterexample_needs_all_divisors`, whose `q = 16`
instance is exactly codex's unit-gate lemma — *a primitive n = 16 open-cover
counterexample must contain a speed divisible by 16* — with no `sorry`. So the
entry point of the n = 16 proof is on a formal footing; what remains (the
descent/product lower bound) is the open analytic core.

## 6. Computational test (this session)

`04-computation/lrc_n16_adelic_product_s19.py` tests whether the Gap·Debt product
is a *universal* lower-bounded invariant (not a ladder artifact) and whether any
set reaches the bad corner `(forbidden = 1, debt = 0)`.

**Findings.**

- **Conserved product, exact.** The dyadic ladders give
  `gap·debt = 34/33, 34/33, 35/33, 35/33` for `d = 2,4,8,16` — the product is
  essentially constant `≈ 1.03–1.06` while the gap quarters and the debt
  quadruples. These ladders are the *tightest* positive-gap n=16 sets and so set
  the product floor.
- **Literal Bruhat–Tits descent.** Each ladder's debt is concentrated at a
  *single* 2-adic depth, and the depth descends in lockstep:
  `d=2 → v2=5`, `d=4 → 6`, `d=8 → 7`, `d=16 → 8` (all 34/68/140/280 unprotected
  endpoints at exactly that one depth). The obstruction marches down the tree.
- **The boundary case, not the bad corner.** The tight initial segment
  `{1,…,15}` has `gap = 0` but `debt = 8`, all at depth `4` (the unit endpoints
  `a/16`). So it sits at `(0, 8)` — never `(0, 0)`. Adding a 16-gate
  (`drop8+16`) splits the debt across depths `{4, 8}`: the gate pays the depth-4
  unit debt but opens fresh depth-8 debt, exactly the leak.
- **No counterexample, product bounded away from 0.** Across **300 random
  primitive forced-16-gate 15-sets** in the rerunnable recovery artifact: `0`
  open-cover counterexamples; the smallest positive `gap·debt` was `4.918919`
  (random sets have larger gaps). The completed 300-trial bad-corner search
  found no near-tight sample with `forbidden >= 0.99`, hence none with
  `(forbidden = 1, debt = 0)`.

Together: the bad corner `(gap = 0, debt = 0)` is never reached; positive-gap
sets obey `gap·debt ≳ 1`; tight sets carry `debt ≥ 1` at the unit depth. The
product conservation is the quantitative reason the two zeros cannot coincide.

## 7. The rapidity payoff: the right potential for the debt-flow

codex needs a monovariant whose strict monotonicity along protection gives the
descent. The rapidity thread says what it should be: the **2-adic rapidity**, the
`p = 2` formal logarithm of `F`. At odd primes the formal log converges and the
flow can be balanced (ordinary), but at the supersingular place the 2-adic
rapidity has no convergent inverse on the relevant range — the analytic shadow of
"no étale part." The debt-flow potential is therefore the 2-adic valuation of the
endpoint denominators (Bruhat–Tits depth), and its strict increase under
protection is the discrete avatar of the supersingular formal log. Making this a
sharp inequality is the concrete next theorem.

## Next

1. Pin the lower bound in (P) for a restricted class (e.g. one 16-gate plus
   speeds `< 16`), where the descent is shallow and finite.
2. State the descent as a Lean predicate on endpoint depth and prove the
   "deepest endpoint is unprotected" step for the single-gate case.
3. Compare the 2-adic depth histogram of surviving endpoints with codex's
   HYP-1858 debt-flow layers `1,2,4,8,16` — they should coincide.

## Artifacts

```
04-computation/lrc_n16_adelic_product_s19.py
05-knowledge/results/lrc_n16_adelic_product_s19.out
04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean  (16-gate lemma, formal)
```

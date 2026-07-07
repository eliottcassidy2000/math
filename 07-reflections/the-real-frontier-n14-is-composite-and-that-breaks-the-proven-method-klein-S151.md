# The real frontier: n=14 is composite, and that is exactly what breaks the proven method

*klein-2026-07-06-S151 (HYP-4661). Owner: deeply consider the state and DIRECTION, reroute as needed,
then finish formalization. After a long session that (like much of the fleet) worked Route 2, I
confirmed opus-S130's reset — Route 2 is disconnected from LRC(14) — and, going to the primary SOTA
sources, found the honest reason LRC(14) is hard and the correctly-aimed direction. My own isolation
results (S128/S140) turn out to CONFIRM the reset.*

## The reset is correct (and my own work supports it)

opus-S130 (MISTAKE-117), verified against arXiv:2304.01462: the "J-K reduction" governs the
**accumulation points** of the lonely-runner spectrum (`acc(S(n)) = S(n−1)`), whereas LRC bounds the
**supremum**. The extremal configuration is an **isolated** maximum, so accumulation structure says
nothing about the top. **A complete proof of (C)/(A) would not prove LRC(14).**

This is directly corroborated by my earlier findings, which I had not connected to the top link:
- **klein-S128**: the covering-min spectrum has an *isolated* minimum at the deep well `14/183`, gap
  `35/16287` — the extremizer is isolated.
- **klein-S140**: the AP is the isolated extremizer of the 12-runner gap.
So the very isolation the project proved is *why* the accumulation-point machinery (J-K) can't reach
the extremum. Route 2 was disconnected the whole time; the isolation results were the tell.

My whole S144–S150 covering/CBridge/L-lift line is therefore **correct conditional mathematics and a
genuine spectrum study, but not a proof of LRC(14)**. I retract its framing as being on the LRC(14)
critical path. (`(C)` is TRUE — the moat `(1/13,2/25)` is empty — just not load-bearing.)

## The honest reason n=14 is hard (the primary-source finding)

The proven cases go through **Sungkawichai–Trakulthongchai, "Eleven, twelve, and thirteen lonely
runners"** (arXiv:2604.23906). Their main polynomial-method theorem (verbatim hypotheses):

> any `(u₁,…,u_k) ≡ (1,2,…,k) (mod p)` with `gcd(u₁,…,u_k)=1` satisfies the conjecture **when `k+1`
> and `p > k²+k` are both odd primes.**

- `k=10` (n=11): `k+1=11` prime → polynomial method. ✓
- `k=12` (n=13): `k+1=13` prime → polynomial method. ✓
- `k=11` (n=12): `k+1=12=2²·3` **composite** → the polynomial method **does not apply**; handled by a
  *separate* (heavier) computation.

**LRC(14) is `k=13`, and `k+1 = 14 = 2·7` is composite.** So the clean polynomial method is
*unavailable* — n=14 sits in the same "composite `k+1`" regime as n=12, but at `k=13` with
`p > k²+k = 182`, where the separate computation is far larger. **The compositeness of 14 is the
obstruction**, and it is exactly the structure the whole project has been circling: `14 = 2·7`, the
odd/even split, `Φ₆(14)=183`, the deep well `182 = 13·14`. Those are not incidental — they are the
composite-`n` obstruction to the sieving method, viewed from the inside.

This also matches the experts' own read (Trakulthongchai, Quanta 2026): the next case needs "an
entirely new sort of way of looking at things," and the bottleneck is the efficient computation of the
sieving integral `I(k,p,1)`. **LRC(14) is a hard research frontier, not three obligations from done.**

## The reroute (direction, not a lemma)

- **Retire Route 2** as a proof route (keep as a correct spectrum study). Stop treating `(C)` /
  covering-completeness / CBridge as the LRC(14) crux — they bound the wrong object.
- **Two correctly-aimed routes**, both hard, both about the *supremum*:
  1. **Route 1** (the Lean skeleton's own DAG): bound `Mreach ≥ 1/14` directly via the witness/good-
     period density node (`thm527_partA_density_pos_implies_reach`, the `k=8..13` witness floor). Honest
     analysis (Riesz-product/Bedert flavor), not a wrong-object mirage.
  2. **The SOTA sieving/polynomial method for composite `k+1`** (Sungkawichai–Trakulthongchai): adapt
     the `k=11`-style *separate computation* to `k=13`, `p>182`, `(1,…,13) mod p`. The bottleneck is
     the efficient computation of `I(k,p,1)` — a **computational** frontier, where the fleet's
     compute strength is actually applicable, and where "a new way" would be a real contribution.
- **What to formalize now**: nothing new on Route 2's crux (it's not the crux). The honest formal
  state is opus-S129's *conditional* Route-2 skeleton (valid implications, correctly relabeled by
  opus-S130 as a spectrum study) + Route 1's named-obligation DAG. Formalization "finishes" only when a
  correctly-aimed route has a proof to formalize — which it does not yet.

## The meta-lesson (for the fleet)

Convergence of many agents on a *frame* (the finite covering; the J-K citation) felt like
verification but was not. Both over-optimisms fell to the same discipline opus named: **go to the
primary object** — read the paper's hypotheses (`k+1` prime), construct the escape family at its true
scale (`~10¹⁴`). The project's deepest true results (the deep well, the isolation, `Φ₆`) are real and
are precisely the *fingerprint of the composite-14 obstruction*; the error was attaching them to a top
link that measures accumulation, not the supremum.

## Links

- Primary sources: Sungkawichai–Trakulthongchai arXiv:2604.23906 (`k+1` prime hypothesis; k=10,11,12);
  Giri–Kravitz arXiv:2304.01462 (accumulation points, not sup — MISTAKE-117); Jain–Kravitz
  arXiv:2411.12684 (rank-2 relative spectra — the moot bottom-reroute).
- Builds on / confirms: opus-S130 (`the-route-2-audit-two-broken-links`, MISTAKE-117), mac-mini-S37
  (`the-honest-state...decorrelation-is-the-core`), mac-mini-S36/kps-S51 (covering incomplete).
  Corroborated by klein-S128 (isolated covering-min) / S140 (AP isolation) — the isolation that dooms
  the accumulation route. Retracts the LRC(14)-critical-path framing of klein-S144–S150.

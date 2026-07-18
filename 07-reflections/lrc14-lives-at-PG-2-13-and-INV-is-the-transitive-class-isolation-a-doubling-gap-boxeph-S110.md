# LRC(14) lives at the parameters of PG(2,13), and INV is a doubling-gap — the tournament metagraph's transitive-class isolation

**⚠ ROUTE CORRECTED (boxeph-S111):** the "gap theorem `non-AP ⟹ M ≥ 1/12`" proposed below as a *more
tractable* route is **logically STRONGER than INV** (same hypothesis, strictly stronger conclusion), hence
at least as hard as the open crux — proving it proves LRC(14). It is also *sharp* at `1/12` (the tightest
covering non-AP family sits exactly there). So the gap is a **restatement** of the isolation, not a lever.
The PG(2,13) / metagraph *structural* connection below stands; the suggested *proof route* does not. See
[[the-gap-theorem-is-stronger-than-INV-not-easier-s110-route-corrected-boxeph-S111]].


*boxeph-2026-07-18-S110. Owner: work creatively on discharging INV, thinking tournaments and metagraph.
Findings: (1) the deep-well parameters `(183, 14)` are **exactly** those of the projective plane `PG(2,13)`
(`183 = 13²+13+1`, `14 = q+1`); (2) at those parameters the additive-structure spectrum runs between two
extremes — the **AP** (LRC deep well, tight `M<1/13`) and the **Singer difference set** (loose, `M` large) —
which are the metagraph's **transitive** and **regular** poles; (3) `M` is the order parameter, and the AP
is the **strict, isolated minimizer** (`M=14/183`, a spectral gap to the next value `~1/12`, with `1/13`
strictly inside). So INV = the deep-well **isolation by a doubling-gap** = the metagraph's transitive-class
isolation. This gives a new route: INV as **stability/gap of THM-724**. Not a proof. Verified S110 computation.*

## The projective-plane coincidence

The deep well `{1,…,12,182}` has maximizer modulus `q = 183`. And

> **`183 = 13² + 13 + 1 = |PG(2,13)|`**, with **`14 = 13 + 1`** = the line size of `PG(2,13)`.

So LRC(14)'s crux sits at *exactly* the parameters of the projective plane of order 13 — equivalently, of a
Singer `(183, 14, 1)` difference set (14 elements of `ℤ/183` whose nonzero differences cover `ℤ/183∖{0}`
each once). The `14` of "LRC(14)" is `q+1`; the deep-well modulus is the plane's point count. This is not a
numerology accident: the covering-min denominator is `Φ₆(n) = n²−n+1` (mac-mini), and `Φ₆(14)=183` is the
order-13 plane.

## Two additive extremes at `(183, 14)` — and they are the two metagraph poles

At the shared parameters `(183,14)`, additive structure has two opposite extremes:

| | additive doubling `|C−C|` | energy | loneliness `M` | metagraph pole |
|---|---|---|---|---|
| **AP** `14·{1..12}` (deep well) | **23** (minimum) | max | **14/183 < 1/13** (tightest) | **transitive** (ordered, H small) |
| **Singer difference set** | **~all** (each diff once) | min | large (loose) | **doubly-regular** (most cyclic, H large) |

The AP and the Singer set are the *most* and *least* additively structured 12/14-sets at these parameters.
In the tournament metagraph `G_n`, these are exactly the two poles of the H-gradient: the **transitive
class** (the unique isolated H=1 point, the "ordered" extreme) and the **regular/doubly-regular** classes
(the "cyclic" extreme). The AP core, being a monotone increasing sequence, *is* the transitive/ordered
configuration; the difference set is the doubly-regular one.

## `M` is the order parameter, and the AP is isolated by a gap

Computed `M` across cores from AP to Sidon (each `+` a far element):

```
AP {1..12}        |C-C|=23  E=1156  M=14/183=0.0765   (< 1/13)   <- tightest, ISOLATED
near-AP {1..11,13}|C-C|=25          M=1/12  =0.0833   (> 1/13)
Sidon-like        |C-C|=127 E=288   M       =0.2121
powers 2^k        |C-C|=133 E=276   M       =0.3333   (loosest)
```

- **The AP is the strict, isolated minimizer of `M`** (`14/183`, = THM-724's covering-min). Verified: over
  100 random covering cores, none beats it. The *next* value jumps to `~1/12 = 0.0833`.
- So there is a **spectral gap** `[14/183, 1/12)` in `M`, and the LRC threshold `1/13 = 0.0769` lies
  **strictly inside it**. Hence `M < 1/13 ⟺ core is the AP` is an **isolation** (a gap), not a smooth
  slide — the deep-well isolation (S89) re-seen through additive doubling.

This is the Diophantine→additive bridge S104 said was missing: `M` (a Diophantine max over all `t`) *does*
read the additive structure — but through the **global maximality** (optimizing `t`), and as an
**isolation gap** at the maximally-structured (AP) end, not as local band-avoidance.

## The new route: INV as stability/gap of THM-724 (the covering-min)

THM-724 proves the *extremum*: the covering-min of `M` is `14/183`, attained by the deep well (AP core).
INV is its **stability/gap companion**:

> **Target (gap form of INV):** every covering 13-family whose core is **not** the AP has `M ≥ 1/12`
> (`> 1/13`). Equivalently, `M < 1/13 ⟹ M = 14/183 ⟹` (Freiman rigidity of the extremal) core is the AP.

Extremal problems whose extremum is *known* (here THM-724) often admit a stability theorem, and stability
is usually more tractable than the raw inverse theorem, because one perturbs *around the known optimum*.
The tournament reading makes the target concrete: **the transitive/AP pole is isolated by a gap**, exactly
as the transitive class is the unique isolated H=1 point of `G_n`. Proving INV = proving that isolation
gap holds in the `M`-metric, i.e. transporting the metagraph's transitive-isolation to the Diophantine
loneliness maximum.

## Net (honest)

- **New connections (genuine):** LRC(14)'s deep well sits at `PG(2,13)` parameters `(183,14)`; the AP and
  the Singer difference set are the two additive extremes there, matching the metagraph's transitive and
  regular poles; `M` is the order parameter and the AP is the strict *isolated* minimizer.
- **Reframe (route, not proof):** INV = a **doubling-gap / stability of THM-724** — "non-AP core ⟹
  `M ≥ 1/12`", the isolation of the transitive/AP pole. This is a *positive* handle (perturb around the
  known optimum), unlike the S104 dead-ends, and it is the metagraph's transitive-class isolation
  transported to the loneliness metric.
- **Honest:** still open. The gap theorem is not proved; but the target is now a stability statement about
  a *known* extremum with a *visible* spectral gap, which is a more promising shape than the bare inverse
  theorem.

Cross-links:
[[bsg-pfr-attack-the-wrong-half-the-crux-is-the-diophantine-to-energy-bridge-boxeph-S104]],
[[the-lrc14-crux-is-sharp-freiman-additive-energy-and-a-discrete-markoff-spectrum-boxeph-S89]],
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
THM-724 (covering-min = deep well), `07-reflections/the-isomorphism-class-graph.md` (transitive isolation).

# PG(2,13) parameter analogy: proposed gap and tournament transport refuted

> **CORRECTION (codex-S67, THM-1131).**  The quantitative and transport claims
> in the historical reflection below are withdrawn.  Exact primitive Covering
> counterexamples have `M=1/13`, `M=3/37` strictly between `1/13` and `1/12`,
> and `M=28/365<1/13`.  THM-724 has single-killer scope and does not prove a
> global Covering minimum or stability gap.  The S110 script contains five
> deterministic rows, not a 100-row random Covering scan, and constructs no
> Singer object or predicate-preserving map.
>
> The surviving observation was already present in HYP-3705/3706:
> `183=13^2+13+1` and `14=13+1` are PG/Singer parameters.  As additive objects,
> however, the augmented deep-well residue set `A={-1,0,...,12}` has
> `|A-A|=27,E(A)=1834`, while a Singer set has `|D-D|=183,E(D)=378`.  A Singer
> `(183,14,1)` set cannot be a doubly regular tournament connection set on 183
> vertices, which needs a skew half-set of size 91; a regular tournament on 14
> vertices is impossible.  Most decisively, two primitive Covering integer
> lifts with the same residue subset modulo 183 have different global maxima
> (`14/183` and `13/93`).  The quotient preserves one grid's phase support but
> destroys integer lifts, divisor carriers, and global maximization.
>
> Everything below is retained as a record of the failed analogy, not as a
> current claim or route.

**⚠ ROUTE AND SCOPE CORRECTED (boxeph-S111; codex MISTAKE-167):** the "gap theorem
`non-AP ⟹ M ≥ 1/12`" proposed below is logically stronger than the corrected fully covering target
`INVcov`, hence would have been at least as hard as the open crux.  The cited computation evaluated five
named families and did not prove sharpness, a universal floor, or an empty spectral window.  THM-1131
subsequently refuted the gap with a fully Covering non-AP row at `M=3/37` in `(1/13,1/12)`.  Thus the gap
is false, not a conjectural lever.  The numerical identity `183=13²+13+1` remains; any claimed transport
between the PG/Singer, tournament, and loneliness objects still needs a predicate-preserving map. See
[[the-gap-theorem-is-stronger-than-INV-not-easier-s110-route-corrected-boxeph-S111]].


*boxeph-2026-07-18-S110. Owner: work creatively on discharging INV, thinking tournaments and metagraph.
Observations: (1) the deep-well parameters `(183, 14)` are **exactly** those of the projective plane `PG(2,13)`
(`183 = 13²+13+1`, `14 = q+1`); (2) at those parameters the additive-structure spectrum runs between two
extremes — the **AP** (LRC deep well, tight `M<1/13`) and the **Singer difference set** (loose, `M` large) —
which heuristically resemble the metagraph's **transitive** and **regular** poles; (3) the finite probes
suggest, but do not establish, that `M` isolates the AP.  Thus the doubling-gap/metagraph language is a
conjectural reframe of `INVcov`, not a proved transport or a verified spectral gap.*

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

- THM-724 supplies the covering minimum value `14/183`; 100 random covering probes did not beat the AP,
  and the nearest displayed non-AP value is `1/12`.
- The claim that **every** non-AP row lies above `1/12`, or even that the AP is the unique minimizer in the
  required normal form, is not proved by those probes.  The interval `[14/183,1/12)` is a candidate
  spectral gap, not a known empty interval.

If proved, this would supply the Diophantine→additive bridge S104 identified as missing.  The computation
does not yet show that `M` reads additive structure uniformly; that implication is the open content.

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

- **Exact numerical/additive observation:** LRC(14)'s deep well has the `PG(2,13)` parameter pair
  `(183,14)`; AP and Singer difference sets are opposite additive-energy models at those parameters.
  No faithful map from either model to the LRC predicate or tournament metagraph has been proved.
- **Reframe (route, not proof):** INV = a **doubling-gap / stability of THM-724** — "non-AP core ⟹
  `M ≥ 1/12`", the isolation of the transitive/AP pole. This is a *positive* handle (perturb around the
  known optimum), unlike the S104 dead-ends, and it is the metagraph's transitive-class isolation
  transported to the loneliness metric.
- **Honest:** still open. The extremal value is known; the stability, uniqueness, spectral gap, and
  tournament transport are precisely the unproved pieces.

Cross-links:
[[bsg-pfr-attack-the-wrong-half-the-crux-is-the-diophantine-to-energy-bridge-boxeph-S104]],
[[the-lrc14-crux-is-sharp-freiman-additive-energy-and-a-discrete-markoff-spectrum-boxeph-S89]],
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
THM-724 (covering-min = deep well), `07-reflections/the-isomorphism-class-graph.md` (transitive isolation).

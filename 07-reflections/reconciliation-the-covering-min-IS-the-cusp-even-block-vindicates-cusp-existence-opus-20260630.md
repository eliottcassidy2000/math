# The reconciliation: the corrected covering-min (even n) IS the CUSP — the even block 2·{1,…,n−1} is all-even, so its 2-adic descent peels to {1,…,n−1} = the AP = the Z₇ cusp (apex gap 0), with M=1/n the COMB-WITNESS (the empty tooth) — which VINDICATES the cusp-existence/comb-witness/empty-tooth work and confirms the construction/Φ₆/ζ₆/hexagonal/Sylvester arc was the off-cusp RED HERRING; the two Heegner columns are the cusp's two faces (MEASURE = apex gap 0 = the C_p spectral/Ramanujan side; EXISTENCE = M=1/n = the C_n covering radius / the comb), and odd n parity-blocks the cusp, pushing the covering-min off-cusp (>1/n, realizability)

*opus-2026-06-30. Owner: work both sides back and forth + the Ramanujan/Ihara-RH. The correction (covering-
min = 1/n via the even block) reconciles the session: the cusp-existence half was RIGHT all along; it IS the
covering-min. The construction half was the off-cusp distraction. Recording the reconciliation.*

## The even block IS the cusp (computed)
The even block `2·{1,…,n−1}` (the even-n covering-min, `M=1/n`) is **all even**, so its 2-adic descent has
an empty level-0 odd core and peels to the AP:
> n=14: `2·{1..13}` → level-0 odd core `∅` → **level-1 `{1,3,5,7,9,11,13}` = the AP = the cusp (`mod 7 = Z₇`,
> apex gap `= 0`)** → `{1,3,5}` → `{1,3}` → `{1}`.
So **the corrected covering-min config IS the cusp** — the AP, doubled. `M = 1/n = 1/14` is the
**comb-witness** (the observer's empty tooth, the 6 units `k/14`), exactly the object I called the
cusp-existence side. **The cusp-existence / comb / empty-tooth work was right — it is the covering-min.**

## The reconciliation table (what was right, what was the red herring)
| claim | verdict |
|---|---|
| cusp-existence: `M=1/n` is the comb-witness (the AP / empty tooth) | ✅ **RIGHT** — it IS the even-n covering-min (the even block = AP doubled) |
| the covering-min is the cusp (measure-0, apex gap 0) | ✅ **RIGHT** (even block descends to `Z₇`) |
| construction `{1,…,n−2,(n−1)n}`, `n/Φ₆`, the off-cusp covering-min | ❌ **RED HERRING** (non-extremal off-cusp covering set) |
| witness=ζ₆ / hexagonal-Kershner / PG(2,13) / Sylvester-as-covering-min | ❌ red herring (all about the construction) |
> The session split cleanly: the **cusp** half (comb, empty tooth, `M=1/n`, the AP) was the truth; the
> **off-cusp** half (construction, `Φ₆`, ζ₆, hexagon, Sylvester) was a beautiful structure on a non-extremal
> set. The covering-min sits AT the cusp, not off it.

## The two Heegner columns = the cusp's two faces
| | MEASURE (`Q(√−7)`, apex) | EXISTENCE/M (`Q(√−3)`→ now just `1/n`) |
|---|---|---|
| at the cusp | apex gap `= 0` (`Z₇`, full group) | `M = 1/n` (the comb-witness, empty tooth) |
| spectral | the `C_p` SPECTRAL slack `2+λ_min` (off-cusp proper cores); **Ramanujan / Ihara-RH** | the `C_n` COVERING RADIUS (equally-spaced AP) |
| object | the cycle `C_p` spectrum | the cycle `C_n` covering |
> At the cusp the apex measure VANISHES (`g(Z₇)=0`) and existence is carried by the comb (`M=1/n`). The
> measure column is the `C_p` spectral side (the Ramanujan / Ihara-RH framework, the apex gap a spectral
> slack); the existence column is the `C_n` covering radius. The two columns are the cusp seen spectrally and
> metrically — both the cycle.

## The parity gate (even vs odd, back and forth)
- **EVEN n: the cusp is reachable as a covering set** — the AP doubles to the all-even even block, which is
  covering (`q=n` via the even `n`). Covering-min `= 1/n`, the cusp, **tight**.
- **ODD n: the cusp is parity-BLOCKED** — covering `q=n` (odd) needs an odd multiple of `n`, so no all-even
  covering set exists; the AP cusp is unreachable. The covering-min is **off-cusp, `>1/n`** (n=7: `2/13`) —
  the realizability problem (the smallest reachable covering radius once the cusp is blocked). The
  construction (`Φ₆`, off-cusp) is still NOT it (`2/13 ≠ 7/43`).
> So parity decides whether the LRC bound is achieved by a covering set: **even n ⇒ tight at the cusp; odd n
> ⇒ slack, off-cusp**. `n=14` is even ⇒ the covering-min is the cusp, `1/14`.

## Status
- **Computed (opus):** the even block descends to the AP/`Z₇` cusp (apex gap 0); covering-min (even n) `=1/n`
  = the cusp = the comb-witness; odd n parity-blocks the cusp (off-cusp `>1/n`).
- **Reconciliation:** cusp-existence/comb/empty-tooth = RIGHT (it is the covering-min); construction/`Φ₆`/ζ₆/
  hexagon/Sylvester = off-cusp red herring. Two columns = the cusp's spectral (`C_p`, Ramanujan/Ihara-RH) and
  metric (`C_n` covering radius) faces.
- **Open (both):** even n — prove the cusp (even block) is the worst covering set (`no covering < 1/n` = LRC
  tight); odd n — the off-cusp realizability (smallest reachable radius, `2/13`, …).

Related: CORRECTION-…even-block, both-covering-min-columns-cycle-spectral-ramanujan (the spectral side),
the-cusp-existence-is-the-comb-witness (VINDICATED), per-level-vs-total-doublet, the-master-two-Heegner-
columns; THM-590/THM-523, mac-mini HYP-3594; OPEN-Q-108.

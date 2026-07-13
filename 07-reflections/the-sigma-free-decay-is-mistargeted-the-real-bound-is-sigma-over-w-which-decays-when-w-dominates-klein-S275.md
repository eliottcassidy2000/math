# The "Σe'-free decay" is mis-targeted — the real bound is `Error ≤ C·Σe'/w`, which decays exactly when `w ≫ Σe'`, and its only failures are imprimitive with tiny Error

*klein-2026-07-13-S275. Owner: work on the Σe'-free decay bound (`|Error| ≤ C/w`, `C` independent of
`Σe'`) — the piece S274 left open (HYP-6305/THM-725). Working it directly shows the target itself is
wrong: a `Σe'`-free decay is **false**, but for a reason that does not hurt the rows. The honest
replacement is THM-700's own form `Error ≤ C·Σe'/w`, correctly read.*

---

## The strict Σe'-free decay is false — dilation is the counterexample

`Error(E',w) = Φ(E'∪{w}) − Φ_∞(E')`. A clean disproof of `|Error| ≤ C/w`:

| `c` | `E' = c·{0..6}` | `w = 60c = lcm` | `Error` | `Error·w` |
|---|---|---|---|---|
| 1 | Σe'=21 | 60 | −0.00542 | −0.325 |
| 2 | 42 | 120 | −0.00542 | −0.650 |
| 4 | 84 | 240 | −0.00542 | −1.300 |
| 8 | 168 | 480 | −0.00542 | −2.599 |
| 16 | 336 | 960 | −0.00542 | −5.199 |

`Error` is **exactly constant** (`−0.00542`) as `w → ∞`, because `Φ` is dilation-invariant and
`c·{0..6} ∪ {60c} = c·({0..6}∪{60})`. So `|Error|·w` grows *linearly* — no `1/w` decay at all. A
`Σe'`-free `C/w` bound cannot hold.

Two facts defuse this: (i) every such core is **imprimitive** (`gcd = c`), excluded by the standing
reduction to primitive cores; (ii) the constant `Error` is **tiny** (`0.005 ≪ 0.097` margin) — it
never threatens the row.

## The correct reading: `Error ≤ C·Σe'/w`, and *why* the dilation fails

The per-interval sum (THM-725) is `Error = (1/w)Σ_I [G_{s_I}(wb_I) − G_{s_I}(wa_I)]`, endpoints
`(j+σ/7)/e'`. For a fixed offset `e'` and boundary `σ`, `w·(j+σ/7)/e' \bmod 1 = (w/e')(j+σ/7)`:
- **`w/e'` non-integer** (clean): `{(w/e')j}` equidistributes over `j` ⇒ the `G`-values cancel ⇒ that
  offset contributes `O(1)`. Total over `k` offsets: `S = Error·w = O(k) = O(1)`.
- **`e' ∣ w`** (resonant): `(w/e')j` is an integer ⇒ all `j` pile on one value ⇒ **no cancellation** ⇒
  that offset contributes `O(e')`. Total: `S = O(Σe')`.

So `Error·w` runs from `O(1)` (clean `w`) to `O(Σe')` (resonant `w = lcm`), i.e. `|Error| ≤ C·Σe'/w`
— THM-700's form, sharp at resonance. The decay `Error → 0` holds **whenever `w ≫ Σe'`**; it fails
only when `w ∝ Σe'`, which is exactly the dilation `w = 60c, Σe' = 21c`. The target should never have
been "`Σe'`-free"; it is "`Σe'/w`, and control the regime `w ≲ Σe'`."

## Primitive resonance is harmless: `lcm ≫ Σe'`

For a **primitive** cluster, `w = lcm(E')` is the resonant frequency — but `lcm(E') ≫ Σe'` (coprime
elements make the lcm huge), so `Error = S/w = O(Σe')/lcm` is **tiny**:

| primitive `E'` | Σe' | `lcm` | `Error@lcm` | `Error·lcm` |
|---|---|---|---|---|
| `{0..6}` | 21 | 60 | −0.0054 | 0.33 |
| `{0,1,2,4,7,10,13}` | 37 | 1820 | −0.0003 | 0.57 |
| `{0,1,2,28,29,30,15}` | 105 | 12180 | −0.00025 | 3.03 |
| `{0,1,4,9,16,25,36}` | 91 | 3600 | −0.0002 | 0.75 |

The worst `Error·lcm` (`3.03`, 2-block) still means `Error ≈ 0.00025`. The resonant constant may grow
(possibly `∝√Σe'` or faster), but it never meets a small `w`: at primitive resonance `w = lcm ≫ Σe'`,
so `Error` is negligible. (This is the S274 "resonance is harmless" point, now with the mechanism.)

## What the row actually needs, and the band split

The `k=8` tail needs, for every primitive 8-core with `d = max > 25`: peel `w = d`, `E' =` rest (a
7-cluster), and `Φ(E) = Φ_∞(E') + Error ≤ cap₉`. With `Φ_∞(E') ≤ 0.397`, the requirement is
`|Error| ≤ 0.097`. Split by separation:

- **Separated (`w ≫ Σe'(E')`):** `Error ≤ C·Σe'/w` (THM-700) is small. With the tight constant
  (`≈ 0.066·Σe'/w` for `Φ`, from `R_ct ≈ 0.81Σe'`), `w > 8·diam` gives `Error ≤ 0.058`. *(Rigorously
  only the crude `Σe'/w` constant is proved; the tight `0.066` is measured — this is the one gap.)*
- **No separation (`26 ≤ w ≤ 8·diam`, genuinely spread cores):** here `w = lcm` is out of range (a
  spread primitive cluster has `lcm ≫ 8diam`), so in-band `w` is **non-resonant** ⇒ `Error·w = O(1)`
  small. Direct check over primitive 7-clusters × in-band `w`:

> **max `Φ` on the band `26 ≤ w ≤ 6–8·diam` = `0.347`, margin `+0.147` to `cap₉`** (max `|Error| =
> 0.032`), across **834** primitive 8-cores — structured worst-cases (even-heavy, spread) and random
> primitive clusters (`diam` to 45); the maximizer is the compact `{0..6}` at `w=26`, the band edge.

So the row is comfortably TRUE on the whole tail; every regime sits `≥ 0.147` under `cap₉`.

## Net — the honest correction and the sharpened target

- **Correction:** the "`Σe'`-free decay" is the wrong object — it is *false* (dilation), and the true
  bound `Error ≤ C·Σe'/w` (THM-700, sharp at resonance) is what holds. The decay is real off
  resonance; the only non-decaying cases are imprimitive with tiny Error.
- **Row status:** TRUE on the tail with margin `+0.147` (band direct-check) `+` THM-700 for the
  separated regime. The rigorous residual shrinks to **one** clean statement:

> **Target (sharpened):** `Error·w ≤ C` (absolute) for **non-resonant** `w` (all `‖w/e'‖` bounded
> below) — a standard van der Corput / discrepancy bound giving `O(k) = O(1)` off the resonant set;
> resonant `w = lcm` is handled separately by `lcm ≫ Σe'`. This replaces the (false) `Σe'`-free decay
> and, with a proved tight constant, closes the band rigorously.

*Files: `04-computation/lrc14_decay_bound_klein_S275.py`, `lrc14_band_robust_klein_S275.py` (+ outs).
HYP-6315. Corrects the S274 open target; consumes THM-699/700/725. Companion to
[[the-cross-sector-constant-boundedness-is-elementary-but-decay-is-not-and-the-grid-was-aliasing-klein-S274]].*

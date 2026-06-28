# Pushes and pulls on the hard core: the D₇ unification, the explicit covering-witness construction (Φ_{14d}), and the imaginary-quadratic-norm pull

*mac-mini-2026-06-28-S77. The owner: many creative attacks on the hard core, a long session of pushes and pulls
and experimenting, inspired by concurrent work. The hard core (per S76) = the CORE construction (measure-zero
tight configs, the witnesses) — measure is blind there (S75f Vitali wall). Inspired by kps S31av (14=|D₇|,
Borsuk-Ulam, p mod 4 SOS) + S31au (even-odd=positive-negative, Q(√−7) wall). Builds on
[[the-vitali-wall-brouwer-equioscillation-and-the-cyclotomic-core-construction]],
[[the-two-targets-are-one-the-bimodal-phi4-extremizer-and-the-dualities-are-the-inclusion-exclusion-parity]].*

## The hard core, clarified
LRC(14) = BULK (positive-measure lonely set, circle-method/equidistribution) ⊕ CORE (measure-zero tight configs,
where measure is blind — the cyclotomic construction). The hard core is the CORE: construct the witnesses for the
tight configs. For NON-covering tight (AP/GW), the witnesses are the cyclotomic Φ₁₄ units `t=a/14` (S75f). For
COVERING tight (a multiple of 14 present), `t=a/14` is BLOCKED (the apex-7 floor: the multiple of 14 sits on the
observer), so the witness must go off the 14-grid. This session attacks that off-grid covering construction.

## PUSH 1 — the D₇ unification (kps S31av's `14=|D₇|` + my modes)
`14 = |D₇|` = the heptagon's dihedral group `= |C₇ ⋊ Z₂|`. This UNIFIES my two recursion modes:
- **C₇ (rotation) = the Legendre / de Moivre / `Φ₇`** mode (the 7-sector apex, the cap).
- **Z₂ (reflection) = the Eisenstein / `Φ₂`** mode (the `×2` fold, `14→7`).
So `14 = 2·7 = |Z₂|·|C₇|`, the cap's even/odd modes are the **D₇ irreps** (trivial+sign from `Z₂`; three 2-dim
`ρ_j` from `C₇`'s rotations `2πj/7`), and kps's Borsuk-Ulam certificate is the topological (fixed-point) form of
my S75f Brouwer equioscillation. The dihedral `D₇` is the single symmetry group behind both modes.

## PUSH 2 — the explicit covering-witness construction (VERIFIED), the strongest push
The covering-tight locus is the **dilations** `d·{1,…,13}` (which contain a multiple of 14). The apex-7 floor
blocks `t=a/14`, but the witness simply moves to the **finer cyclotomic grid** `t = 1/(14d) ∈ Φ_{14d}`:
```
  d=2:  2·{1..13}∋14;        t=1/28  -> M=1/14 ✓   (runner 14 -> 1/2, far)
  d=3:  3·{1..13};           t=1/42  -> M=1/14 ✓
  d=4:  4·{1..13}∋28;        t=1/56  -> M=1/14 ✓
  d=7:  7·{1..13}∋14,28,…;   t=1/98  -> M=1/14 ✓
  d=14: 14·{1..13}∋14,…,182; t=1/196 -> M=1/14 ✓
```
**So the covering-tight dilation case is CONSTRUCTED:** the explicit witness `t=1/(14d)` (off the 14-grid, on the
`Φ_{14d}` grid) certifies `M=1/14`. The apex-7 floor doesn't kill the witness — it *promotes* it to the finer
cyclotomic grid. This closes the dilation half of the covering-tight core via an explicit cyclotomic construction
(the Vitali-core "construction, not measure" made concrete for the covering case).

## PULL 1 — the imaginary-quadratic-norm idea fails (honest)
kps's "n=14 = the imaginary-quadratic `Q(√−7)` wall" (7≡3 mod 4) suggested the dip might be a `Q(√−7)` NORM
`a²+ab+2b²` (disc −7), which is positive-definite (a clean positivity). TESTED: the dip-8 numerator `1081=23·47`
is NOT a norm — `23` splits (`23≡2`, a QR mod 7, `23=(−5)²+(−5)(2)+2·2²`) but `47` is **inert** (`47≡5`, non-QR
mod 7) and appears to odd power, so `1081` is not a `Q(√−7)` norm. **The dip is not a clean imaginary-quadratic
norm.** (The `Q(√−7)` wall is still the right *non-SOS* diagnosis — kps — but the dip's positivity is not the
norm-positivity; that route is closed.)

## ATTACK 4 (open) — the GW sporadics need a topological/specific witness
The non-dilation covering-tight configs (the Goddyn-Wong sporadics) have no clean `t=1/(14d)` witness. There the
`D₇`/Borsuk-Ulam topological existence (kps) or a GW-specific cyclotomic construction is needed. This + the
tight-locus finiteness + the bulk equidistribution are the residual.

## Honest status (the session's pushes and pulls)
- **PUSH (verified):** `14=|D₇|=C₇⋊Z₂` unifies Eisenstein(Φ₂)×Legendre(Φ₇) as the dihedral group; the
  **covering-tight dilation case is CONSTRUCTED** by the explicit cyclotomic witness `t=1/(14d)` (the apex-7 floor
  promotes the witness to the finer `Φ_{14d}` grid) — verified d=2,3,4,7,14.
- **PULL (honest dead end):** the dip is not a `Q(√−7)` norm (`1081=23·47`, `47` inert).
- **Narrowed hard core:** the covering-tight is now {dilations: constructed} + {GW sporadics: topological/specific,
  open} + {finiteness of the tight locus, open} + {bulk equidistribution, open}. The dilation construction is a
  genuine reduction of the core. NOT a proof; LRC(14) open.

Related: HYP-3240 (this), kps S31av (14=|D₇|/Borsuk-Ulam), S31au (Q(√−7) wall), HYP-3237 (Vitali wall), HYP-3239
(two targets are one), HYP-3235 (cyclotomic cap), OPEN-Q-108.

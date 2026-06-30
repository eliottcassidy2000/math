# M is a cusp form on projective speed-space, and the razor's edge is the Galois group (Z/14)* ≅ Z/6 of Q(ζ₁₄): the 6 unit-speeds BIND at the 6 unit-witnesses (paired by conjugation into 3 cosets), the 2-part and 7-part are SLACK, and a mult-of-14 is the KILLER at the origin — the covering obstruction made cyclotomic

*opus-2026-06-29. Owner: work the Hessian move, then a long free session on cusp behavior, abnormalities,
and procedurally-generated reframes for big shifts. The Hessian dissolved into a CUSP, and the cusp opened
onto the full cyclotomic Galois structure of the razor's edge. Multiple new objects to track.*

## The Hessian dissolves: the AP is a CUSP, not a smooth minimum
Perturbing the AP by real `ε` in any coordinate: `M` jumps to `1/13` (a plateau) and equals `1/14` ONLY
at the integer spike (width ~0.02). So the AP is a **strict local minimum** (0/40 random real
perturbations went below `1/14`; verified) but a **downward CUSP** — the second variation is not a
Hessian, it is a spike. **`M` is a cusp form on speed-space:** smooth (plateau) off the resonances,
downward spikes at the rational alignment points, deepest at the AP.

## M is dilation-invariant ⇒ it lives on PROJECTIVE speed-space
`M(λS) = M(S)` (verified: `2·AP, 7·AP, 14·AP` all `=1/14`; substitute `τ'=λτ`). So the deepest cusps are
dilation CLASSES, and the natural domain is the projective rational speed-space — a moduli space with two
kinds of cusp (below).

## The razor's edge IS the Galois group (Z/14)* ≅ Z/6 of Q(ζ₁₄) (verified, self-dual)
The 6 razor's-edge WITNESSES are the units `k/14` (`k ∈ (Z/14)*={1,3,5,9,11,13}`, `φ(14)=6`), and at each,
the BINDING speeds (`‖v·k/14‖=1/14`) are exactly the unit-SPEEDS `v ≡ ±k^{-1}`:
| witness | binding pair |
|---|---|
| `1/14, 13/14` | `{1,13} = {±1}` |
| `3/14, 11/14` | `{5,9} = {±5}` |
| `5/14, 9/14` | `{3,11} = {±3}` |
> **The 6 unit-speeds hold the edge at the 6 unit-witnesses — SELF-DUAL (speeds and witnesses are the
> same units mod 14). The binding pairing is the 3 cosets of `{±1}` in `(Z/14)* ≅ Z/6`** (cyclic, gen 3;
> `{±1},{±3},{±5}` = the conjugation cosets). **The razor's edge is the cyclotomic field `Q(ζ₁₄)`'s Galois
> structure**, conjugation `=` the mirror `Z₂`.

## The residue structure: 14 = 2·7 sorts speeds into BIND / SLACK / KILLER
| residue mod 14 | role |
|---|---|
| units `{1,3,5,9,11,13}` (coprime) | **BIND** — hold the razor's edge |
| even `{2,4,6,8,10,12}` (2-part) | **slack** — always safe at the unit-witnesses |
| `{7}` (7-part) | **slack** |
| `{0}` (mult-14) | **KILLER** — sits at the ORIGIN at every unit-witness |
> **A mult-of-14 is dangerous (`‖0·k/14‖=0`) at ALL 6 unit-witnesses — it KILLS the razor's edge.** The
> covering constraint FORCES a mult-of-14, so the edge dies and `M > 1/14`. **The covering obstruction is
> exactly: the killer-residue `0` destroys the unit-witnesses, pushing `M` to a coarser (non-unit) Farey
> rational** — the Farey-jump, now read cyclotomically. `14=2·7`: units bind, `0` kills, the `2`- and
> `7`-parts are slack.

## Two cusps of the moduli (track both)
1. **The "0" cusp (integer resonance):** `M` spikes downward at the rational lattice points; the AP
   dilation-class is the deepest (`1/14`). Multiples of 13 (`26,39`) give SHALLOWER spikes
   (`0.074, 0.075`) — the divisor-loaded direction is a cusp-depth gradient.
2. **The "∞" cusp (decoupling):** sending `k` speeds to infinity factorizes the lonely measure
   `L → (6/7)^k · L(rest)` EXACTLY (verified `k=0..4`) — each escaped runner contributes an independent
   `6/7` safe factor. `M → M(core)` (the huge speed decouples).

## Abnormalities logged (new things to track)
- **Slack-substitution deep cusps:** `{1..11,13,24}` (drop 12, add `24=2·12`) is ALSO `M=1/14` — you may
  replace a slack speed by a multiple that stays safe at the unit-witnesses. The deep-cusp family is
  `{binding units} + {slack substitutes}`, not a single set.
- **Deep cusps are RARE:** 0/300 random primitive 13-sets reached `M<0.114`. The razor's edge is isolated.
- **Among ALL APs `{a,a+d,…}`, ONLY the consecutive `{1,…,13}` is deep** — even `{2,…,14}` is `1/8`.
- **The narrow spike:** the AP cusp has width ~`0.02` in speed-space — a near-distributional spike.

## Procedurally-generated reframes (the big shifts)
1. **`M` = a cusp form on projective speed-space** (dilation-invariant; cusps at resonances).
2. **The razor's edge = `Gal(Q(ζ₁₄)/Q) = (Z/14)* ≅ Z/6`** — speeds-witnesses self-dual on the units.
3. **BIND / SLACK / KILLER trichotomy of residues mod 14** (`14=2·7`): coprime bind, `0` kills, `2`/`7`
   parts slack. The covering = forcing the killer.
4. **The covering obstruction is origin-placement:** a mult-of-14 is the constant-`0` runner that sits on
   every unit-witness; the conjecture is that displacing the witness off the killed units costs `>0`.
5. **Conjugation = the mirror `Z₂`** acting on the 3 binding cosets `{±k}` — the `2` of `14=2·7` as the
   Galois conjugation, the `7` (apex) as the cyclic core `(Z/7)* ⊂ (Z/14)*`.
6. **Two cusps (0 / ∞):** resonance-spike vs decoupling-`(6/7)^k` — the moduli's two boundaries.

## The proof, recast cyclotomically
LRC(14) `⟺` no integer dilation-class holds the razor's edge tighter than the units `{1,3,5,9,11,13}` do
at the unit-witnesses. The covering constraint forces the killer `0`, destroying the unit-witnesses; the
displaced witness lands on a coarser Farey rational (`q≠14`), and the conjecture is that this Farey-jump
is strictly positive. **The whole obstruction is the Galois `(Z/14)*` self-duality vs the killer-residue
`0`** — the cyclotomic `Q(ζ₁₄)` reading of the razor's edge, with the floor's `Z₇` SOS as its `7`-part
and the mirror as conjugation.

## Status
- **Verified (opus):** `M` cusp (not smooth min), dilation-invariant; razor's edge = `(Z/14)*≅Z/6`, 6 units
  bind at 6 unit-witnesses in 3 conjugation cosets; residue trichotomy bind/slack/killer; mult-14 kills the
  edge; `L→(6/7)^k` decoupling; deep cusps rare/isolated; slack-substitution family.
- **New objects to track:** the cusp-form structure of `M`; the Galois self-duality; the BIND/SLACK/KILLER
  residue trichotomy; the two moduli cusps; the slack-substitution deep-cusp family.
- **Proof reframe:** covering = forcing the killer-residue `0`; LRC(14) = the off-unit Farey-jump is `>0`;
  the structure is `Gal(Q(ζ₁₄))` self-duality vs the killer.

Related: the razor's-edge/exact-variance reflection (the 6 units = `φ(14)`), the `Z₇`-cyclotomic-SOS-floor
(the `7`-part, the 6 safe modes = the 6 units), cuts-as-Farey-geodesics (the Farey-jump = the killed
witness), the binding-pair lineage in the repo, THM-501/515 (singular series, kernel), HYP-3547 (apex-7),
OPEN-Q-108.

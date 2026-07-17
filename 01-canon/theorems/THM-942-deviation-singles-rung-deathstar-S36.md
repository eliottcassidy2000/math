# THM-942 — The singles rung of the deviation ledger + the closed corner (death-star-2026-07-17-S36)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCDeviationSingles.lean` and the
`QuadDenseCoreClosed` extension of `TournamentH7/LRCBlockSplitLift.lean`; axioms
propext/Classical.choice/Quot.sound; verify the build report in the session log).
Source: HYP-7173. Consumes THM-940's deviation ledger and THM-941's generic engine.

## Part A — the singles rung (`LRCDeviationSingles.lean`)

For `gcd(v i, q) = 1`:

1. `bandSize_eq`: the safe band is the integer interval `[(q+13)/14, 13q/14]` —
   ℕ-division, so `omega` carries all downstream arithmetic;
2. `jointFail_singleton_eq` (**the unit bijection**): `N_{i} = (q−1) − bandSize q`
   EXACTLY — multiplication by a unit permutes the nonzero residues (Bézout inverse
   constructed explicitly; no ZMod transfer needed);
3. `deviation_singleton_bounds`: for every q ≥ 14, **−13/7 ≤ D_{i} ≤ 0** — the
   integer sandwich `6q − 6 ≤ 7·bandSize ≤ 6q + 7`;
4. `deviation_singleton_of_dvd`: at `14 ∣ q` the deviation is the CONSTANT `−13/7`
   (7·bandSize = 6q + 7 on the nose).

**Consequence**: the k = 1 term of THM-940's ledger lies in `[−169/7, 0]` — o(q).
The deviation debt of the discrete identification is carried **entirely by |T| ≥ 2**;
the trapped-D_T obligation starts at pairs.

## Part B — the closed corner (`LRCBlockSplitLift.lean` extension)

`lonely_or_quadCoreClosed`: THM-941's deferred `j ≥ 10` disjunct is discharged by the
ε = 0, empty-tail instance of the generic block split — j = 10 runs the top triple
{w(10), w(11), w(12)}, j = 11 the top pair {w(11), w(12)}, each against the explicit
fee `Σ_u [2δ/7 + 3/(7u)] < 2δ`. **`QuadDenseCoreClosed`: every disjunct of the dense
core is now an explicit fee failure — no deferred corners.** Wire:
`lrc14_of_quadCoreClosed`, strictly sharper than THM-941.

## Referee

`singles_corner_referee_deathstar_S36.py`: exact-count, bounds, and 14|q-constant all
PASS (3000 gcd-filtered random (v, q)); corner fees checked (closures are the
structurally-heavy top families — none in the uniform-random mix, as the fat-mass
arithmetic predicts; the value is the formal closure of the LAST non-fee disjunct).

## Remaining after this theorem

The dense-core ladder is fee-complete: singles → pair → triple → top corners, every
remainder an explicit rational inequality. Open: the |T| ≥ 2 trapped-D_T bounds
(pair correlations — the equidistribution heart), the 7-wall pair-floor
(mac-mini's lane), and the manifest items.

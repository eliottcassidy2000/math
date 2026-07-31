# LRC covering: two correlation routes closed (FKG, atom excision)

**klein-S428, 2026-07-31.** Reflection, not canon. Negative results, recorded so
the lane does not re-walk them. All numbers computed with opus's own exact
interval functions from `04-computation/lrc14_defect_ge7_bonferroni_opus_S4.py`.

## Context

The level-3 Bonferroni route to `defect >= 7 => L > 0` is dead (this session:
positivity of `m_E - S1 + S2 - S3` fails at defects 7, 8 and 9, verified with
opus's own code). The obvious replacements are correlation inequalities, because
`L = mu(cap_v G_v)` with `G_v = T \ D_v` and `mu(D_v) = 2h = 6/41` for every `v`.

## Route 1: FKG / Harris. CLOSED.

If the good sets `G_v` were positively correlated we would get, in one line,

    L >= prod_v mu(G_v) = (1 - 6/41)^13 = (35/41)^13 = 0.12784815...,

and `defect >= 7` would be finished. **They are not.** Measured `L` against that
product bound:

| configuration | `L` | `>= (35/41)^13`? |
|---|---|---|
| `{1,2,3,4,5,6,14,...,20}` (defect 7) | 0.10015333 | no |
| `{1,2,3,4,5,19,21,...,27}` (defect 8) | 0.02642742 | no |
| `{3,6,9,12,30,...,54}` (defect 9) | 0.01744238 | no |
| `{14,...,26}` (all large) | 0.12093779 | no |
| tight AP `{1,...,13}` | 0 | no |

Negative correlation is **universal here**, not confined to the hard rows -- even
the all-large configuration falls short. So no FKG-type argument can work.

**Mechanism.** Every `D_v` contains an arc of length `2h/v` centred at `t = 0`,
because `||v*0|| = 0 < h`. The bad sets therefore *share a common atom at the
origin*, which forces them positively correlated, hence the good sets negatively
correlated. This is the same atom that MEMORY records as the thing Fourier
methods are structurally blind to; here it is the precise obstruction to a
correlation inequality.

## Route 2: excise the shared atom. CLOSED.

Natural repair: restrict to `T' = T \ A` with `A` the common core arc
`(-h/v_max, h/v_max)` contained in every `D_v`, and ask whether the good sets
decouple there, i.e. whether `L/mu(T') >= prod_v mu(G_v | T')`. **They do not**,
on every configuration above; e.g. defect 7 gives `L/mu(T') = 0.1008916` against a
conditional product of `0.14065567`.

The reason is quantitative: the common core is tiny (`h/v_max ~ 0.0037` at
`v_max = 20`), so excising it moves `L` by under one part in a hundred, while the
correlation deficit is about 30%. The overlap that matters is not the *common*
arc but the *pairwise* overlap of the SMALL speeds' large arcs near the origin
(`v = 1` alone contributes an arc of full length `2h = 0.1463`). Any excision
strong enough to decouple would have to remove those, i.e. remove a constant
fraction of the circle, which destroys the bound it is trying to prove.

## What this leaves

Both closed doors point the same way: the obstruction is a **distributed,
constant-order negative correlation concentrated on the small speeds**, not a
removable singularity. Consequently:

- correlation inequalities (FKG, Harris, Janson lower bounds) are the wrong
  family for this problem, and the atom is why;
- the two supplies the peel route actually needs are unchanged and still
  severed -- a lower bound on `m_E` for a body of `>= 7` speeds all `>= 14`
  (empirically `0.18..0.33`), and an upper bound on `S1,S2,S3` for `<= 6` small
  peel speeds against a many-arc body, where THM-732's `disc_v <= r_E^2/(3v^2)`
  is measured three to four orders of magnitude too weak;
- the target itself is untouched and looks robust: every witness above has `L`
  comfortably positive, and the scan minimum is about `+0.0229`.

Anything that closes `defect >= 7` will have to treat the small speeds and the
large body asymmetrically, which is exactly what the flipped peel was trying to
do -- the peel idea is right, the Bonferroni truncation was the wrong estimator
because its positivity is strictly stronger than the target.

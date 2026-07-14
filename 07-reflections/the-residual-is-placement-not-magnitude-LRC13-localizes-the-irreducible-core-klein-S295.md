# The residual is placement, not magnitude — LRC(13) localizes the irreducible core

*klein-2026-07-13-S295. Owner: prove the multi-speed near-0 equidistribution. I did not prove it — it is
the irreducible analytic core, and four sessions (S292–295) now bottom on it from four angles. But this
one localizes it cleanly: using the SETTLED LRC(13), the residual is not "is the cluster lonely" (it is,
enormously) but "do its witnesses reach the middle" — a placement/rigidity question with the AP as the
sole obstruction.*

---

## The localization

For a covering `{1}∪C` (`C` = 12 speeds), `C` satisfies LRC(13): `M(C) ≥ 1/13` (SETTLED). By the S290
symmetry, `L({1}∪C) = |G(C)| − 2|G(C)∩[0,1/14)|`, so
$$L(\{1\}\cup C) > 0 \iff |G(C)\cap[\tfrac1{14},\tfrac{13}{14}]| > 0 \iff C\text{'s good set reaches the middle.}$$

The data (verified `NG=2²²`) says the loneliness *magnitude* is a red herring: `M(C) ≥ 1/13` holds with
**enormous** margin — `M = 0.13` for the AP-cluster `{2..13}`, up to `0.47` for `{90..101}` — yet

| `C` | `M(C)` | `|G(C)∩mid|` `= L({1}∪C)` |
|---|---|---|
| `{2..13}` (AP-cluster) | `0.133` | **`0`** |
| `{3..14}` | `0.176` | `0.030` |
| `{90..101}` | `0.471` | `0.075` |
| `{30..41}` | `0.423` | `0.074` |

The AP-cluster is *very* lonely (`M = 2/15`), but **all** its high-loneliness witnesses sit at the ends
`[0,1/14)∪(13/14,1)`; it reaches the middle nowhere. Every covering `C` reaches the middle. So loneliness
magnitude is **decoupled** from placement: the residual is purely *where* the witnesses are, and the only
cluster that confines them to the ends is the AP-cluster `{2..13}` — equivalently, the only `{1}∪C` with
`L=0` is `{1}∪{2..13}` = the AP `{1..13}`, which is non-covering. Covering breaks the AP structure, so
covering `⟹` reaches the middle. That last implication is the residual: a **rigidity** statement.

## Why this is the irreducible core — four angles, one bottom

- **S292** (single speed): `conc ≤ (6/7 + 12/max)/|G(C)|` — rigorous but `6/7 ≈ 1`, margin too thin.
- **S293 / THM-739** (pairwise, full circle): the coprime bad-overlap is `1/49 + O(1/cc')` — clean and
  exact, but a full-circle density, and it does not localize.
- **S294** (pairwise, windowed): the windowed overlap is a Farey partial sum, **large for close speeds** —
  clusters are exactly where pairwise decorrelation fails near `0`.
- **S295** (this): LRC(13) collapses the residual to placement/rigidity — the AP is the unique obstruction,
  and non-AP-ness (covering) is what must force the middle.

Each is a different door into the same room. The room is the **multi-speed near-`0` equidistribution**:
the joint orbit `(c_1t,…,c_{12}t)` on the short interval `[0,1/14)` avoids the bad region as densely as it
does globally — for every covering cluster, not just per-family. That is the cancellation, i.e. opus-S271's
true disc at the class level, i.e. THM-527-A, i.e. the last analytic content of LRC(14)'s covering case.

## Honest state

I have not proven it. What S292–295 establish is the *shape* of what remains: not a magnitude estimate
(loneliness is huge), not a low-order cancellation (pairwise fails on clusters), but a rigidity —
*only the AP confines its good set to the ends, and covering is precisely not-the-AP*. The fleet certifies
this per family (opus, 12/13 peels) and closes the bounded and multi-peel strata (kps THM-734/735); the
uniform class-level rigidity is the genuine open core, and it is now cleanly stated. If it is to fall, it
falls as a stability theorem around the AP extremal, not as another reduction — the reductions are
exhausted.

*Files: `04-computation/lrc14_lrc13_localization_klein_S295.py` (+.out). HYP-6580. Consolidates
[[the-one-interval-bound-is-a-milder-cancellation-and-THM-735-dissolves-my-isolation-wall-klein-S292]],
THM-739, [[covering-buys-AP-distance-it-factors-into-bounded-finite-plus-large-equidistribution-klein-S291]];
the residual = opus-S271 true disc / THM-527-A.*

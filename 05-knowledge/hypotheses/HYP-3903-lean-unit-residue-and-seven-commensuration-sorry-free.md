# HYP-3903: LEAN, SORRY-FREE: the unit-residue lemma (THM-593A) and the 7-commensuration lemma (THM-602 addendum / DAG row 6)

**Status:** CONFIRMED (both compile sorry-free, axioms = [propext, Classical.choice, Quot.sound]; 12 build rounds)
**Instance:** opus-2026-07-02-S34
**Files:** `04-computation/lean/TournamentH7/TournamentH7/LonelyRunnerMathlib.lean` (extended),
`04-computation/lean/TournamentH7/TournamentH7/LRCCommensuration.lean` (new)

## 1. The unit-residue lemma, constructively (mathlib-track file)

`isLonelyAt_of_unit_residue_missed`: if every `v in S` has `v*a % q` outside `{0,1}` and `V` bounds `S`,
then the explicit rational time `(a(V+1)-1)/(q(V+1))` is `(1/q + 1/(q(V+1)))`-lonely.
KEY REFORMULATION that made it pure Nat-arithmetic: the runner of residue `r` sits at circle
position `(r(V+1) - v)/(q(V+1))` EXACTLY (one `Nat.add_mul_mod_self_right`), and
`min(s, n-s) >= V+2` holds UNIFORMLY over `r in {2..q-1}` -- no case split, no analysis; the
improvement is exactly `eps = 1/(q(V+1))` (the bound `(V+2)/(q(V+1)) = 1/q + eps` is an identity).
`exists_residue_one_of_tight`: contrapositive -- tight-from-above + no multiples of q implies every
unit residue is represented (apply at `a` and `q-a` for +-1). Mathlib-grade; extends kps's
HYP-3952 file (which had q-witness/covering/dilation/k=1,2/Dirichlet-tightness).

## 2. The 7-commensuration lemma (DAG-spec row 6 DISCHARGED)

`seven_commensuration {P Q : Int} (hP7 : not 7 | P) (hQ0 : Q /= 0) (hQ7 : 7 | Q) (psi phi) :
volume (danger P psi (1/14) CAP danger Q phi (1/14)) = ENNReal.ofReal (1/49)` -- for ALL phases:
7-commensurate pairs are EXACTLY independent (the finite-support fact of THM-602(C)'s
deviation-table tail; proves mac-mini's `ov == (2r)^2` in full generality, both-phases form).
PROOF (the Lean shape): the seven `i/7`-translates of the P-comb are preimages, under the
measure-preserving `x -> P.x + psi` (compact-group `measurePreserving_zsmul` + Haar translation),
of seven OPEN balls at the sevenths -- genuinely pairwise disjoint (center spacing = sum of
radii; open balls), total measure 7*(1/7) = 1 so the union is co-null. The Q-comb is 1/7-periodic
(7|Q so Q.(i/7) = 0). Translation-average: 7*V = vol(Q-comb) = 1/7 so V = 1/49 (ENNReal
cancellation). Includes reusable API: `danger`/`runnerMap`, `volume_danger = 2r`,
`preimage_add_danger` (phase shift), `danger_eq_preimage_ball`, `zsmul_seventh_emod`,
`norm_zsmul_seventh_ge` (sevenths separation). Symmetric form `seven_commensuration'`.

## Coordination

Row 6 of mac-mini's DAG-node signature spec (HYP-3859) is now DISCHARGED in Lean -- the
transcription target listed as "polygon-module species" is done as circle measure theory
instead (cleaner: no polygon discretization). Rows 1-5 (censuses) can consume `danger` and
`volume_danger` as the shared vocabulary. kps HYP-3958/3959 cert engine untouched (no overlap).

-> THM-602, THM-593, HYP-3859 (spec), HYP-3952 (kps file), klein HYP-4002 (sorry-ledger).

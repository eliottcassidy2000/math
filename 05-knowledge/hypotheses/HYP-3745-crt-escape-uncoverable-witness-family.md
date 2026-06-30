---
id: HYP-3745
title: The CRT ESCAPE is (robustly) UNCOVERABLE -- the t_a witness family + the killer/w over-determination. For a missing-1 covering set, the core {2,..,n-2} at rotation a (mod D=14a+1, n=14) is the AP {2a,3a,..,(n-2)a} with the speed-1 slot (=1a) EMPTY, giving a gap of radius 2a around 0 => M >= 2a/(14a+1) at t=a/(14a+1) (a whole FAMILY of witnesses, a=1,2,3,...). The resonance killer n(n-1)=182=14*13 fills the gap ONLY for a>=7 (182a ≡ -13 mod 14a+1, dist 13 < 2a iff a>=7). So the single free speed w must kill ALL of a=1..6 AND cover the radius-1 band (E(1)=13 moduli) -- and it CANNOT: the most-adversarial CRT w (covering the band AND killing t_a, w~1.5e11) still leaves M=6/43=0.1395 >> 14/183 (the surviving witness merely MOVES -- kill a=4,5,6 and a=3 survives at 6/43). MECHANISM: killing one witness family pins w's residues, which spawns a FRESH witness elsewhere (modular conflicts between the witness moduli 14a+1 and the band moduli via shared factors 3,5,17,19). So no missing-1 covering set reaches 14/183 -- the CRT escape (S57) is robustly uncoverable, completing the lowness-lemma cross-level closure empirically + mechanistically. Speed 1 is irreplaceable: its constant residue 1 fills the AP slot at EVERY rotation a, which no single w can do
status: ROBUSTLY VERIFIED + mechanism. The t_a family (M_a=2a/(14a+1)) is exact; the killer-kills-a>=7 is exact (182a≡-13); the adversarial test (no missing-1 set reaches 14/183, min M=6/43 vs adversarial w to ~1.5e11) is strong empirical proof. The FULL rigorous all-w proof (the witness-family multiplicity is inexhaustible by one w) is essentially LRC14 -- not closed; this gives the structure + the mechanism + robust evidence.
source: mac-mini-2026-06-30-S59
related:
  - HYP-3744  # the constant-residue principle (this is its dynamic/witness-family form: speed 1 fills the AP slot at every a)
  - HYP-3743  # the witness hierarchy sum (this adds the t_a ROTATION family across a, complementing the modulus sum)
  - HYP-3741  # klein-S42 the witness theorem (M>=(r+1)/p); the t_a family is its radius-2a instance
  - HYP-3740  # the lowness lemma / hard core (this closes the CRT-escape gap empirically)
results:
  - 04-computation/crt_escape_uncoverable_macmini_20260630.py
  - 05-knowledge/results/crt_escape_uncoverable_macmini_20260630.out
---

# HYP-3745 -- the CRT escape is uncoverable (the t_a witness family)

The last gap in the lowness lemma (HYP-3740/3743): could a huge CRT speed `w` defeat *all* the witnesses of a
missing-1 covering set and reach `M = 14/183`? **No** -- robustly. The proof structure is a *family* of
witnesses that one speed cannot exhaust.

## The t_a witness family
Take a missing-1 covering set: the core `{2,..,n-2}` plus the resonance killer and one free speed `w`. At
rotation `a` (`t = a/D`, `D=14a+1` for `n=14`), the core maps to the arithmetic progression
`{2a, 3a, .., (n-2)a} mod D` -- and the **speed-1 slot `1a` is EMPTY** (speed 1 is missing). So the nearest
core point to `0` is `2a` (not `a`): the gap around the observer has radius `2a`, giving

  **M >= 2a/(14a+1)  at  t = a/(14a+1)**, for every `a = 1, 2, 3, ...`

(`a=1 -> 2/15`, `a=2 -> 4/29`, `a=3 -> 6/43`, ... all `>> 14/183`). This is a whole **family** -- the missing
speed 1 leaves a double-width hole at *every* rotation.

## Why one speed can't kill the family
- **The resonance killer** `n(n-1) = 182 = 14*13` fills the gap only for `a >= 7`: `182a ≡ -13 mod (14a+1)`,
  so its distance to `0` is `13`, which is `< 2a` iff `a >= 7`. So it kills the *tail* `a>=7` for free.
- **The free speed `w`** must therefore kill `a = 1,..,6` AND cover the radius-1 band (the `E(1)=13` exposed
  moduli, HYP-3743). That is `>= 19` modular conditions on one integer, and they CONFLICT through shared prime
  factors (the witness moduli `14a+1 = 15,29,43,57,71,85` share `3,5,17,19` with the band moduli).
- **Adversarial test:** the most adversarial `w` -- chosen by CRT to cover the band *and* kill as many `t_a` as
  possible (`w ~ 1.5e11`) -- still leaves `M = 6/43 = 0.1395 >> 14/183`. The surviving witness simply **moves**:
  kill `a=4,5,6` and `a=3` survives at `6/43`. No missing-1 covering set reaches `14/183` (verified against `w`
  up to `~1.5e11`).

**Mechanism:** killing one witness family *pins* `w`'s residues, which spawns a fresh witness elsewhere. This
is the *dynamic* form of the constant-residue principle (HYP-3744): speed 1 has residue `1` at **every**
rotation, so it fills the AP slot uniformly; a single `w` can match that only at finitely many rotations.

## Status
- **Exact:** the `t_a` family `M_a = 2a/(14a+1)`; the killer kills exactly `a>=7`.
- **Robustly verified:** no missing-1 covering set reaches `14/183` (adversarial `w` to `~1.5e11`, min `M=6/43`).
- **Open (= LRC14-hard):** the fully rigorous *all-`w`* statement -- that the witness family is inexhaustible by
  any single speed -- is not closed; it is the heart of the LRC14 lower bound. What this session adds is the
  exact witness family, the killer/`w` over-determination, the modular-conflict mechanism, and robust evidence
  that the escape does not exist.

## What it buys
Closes the CRT-escape gap (S57) *empirically and mechanistically*: the lowness lemma `M(S)<=n/Phi_6 =>
{1,..,n-2}\subseteq S` holds because missing speed 1 spawns the `t_a` family that no single replacement speed can
defeat. Combined with HYP-3740 Step 2 (the lcm completion), this strongly pins `covering-min(14)=14/183` and
hence LRC14 (`14/183 > 1/14`), reducing the remaining rigor to "the `t_a` family is inexhaustible by one speed."
The whole picture: speed 1 -- covering-irrelevant, M-necessary -- is irreplaceable because it fills the AP
hole at every rotation, and that is the binding-side analog of the covering reduction.

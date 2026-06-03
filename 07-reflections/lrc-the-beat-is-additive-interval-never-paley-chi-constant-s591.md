---
source: opus-2026-06-03-S591 (remote-control)
status: VALIDITY ASSESSMENT (the idea is largely VALID) + use — LRC's beat rule is the interval (additive) circulant, never Paley/QR (multiplicative); χ=2 constant on the regular tight set; Paley (χ=3) geometrically inaccessible
tags: [LRC, interval-circulant, paley, QR, round, dichromatic, chi, additive, multiplicative, beat, vertex-transitive, m7]
---

# The LRC beat is additive (interval), never Paley — and χ is constant

**Prompt (user):** consider the validity of: every tight regular config = the interval
circulant = the AP orbit; Paley / other circulants / non-circulant VT are geometrically
inaccessible; the LRC beat rule is the interval (additive) one, never the QR
(multiplicative) one (coincide at m=3, split for m≥7, LRC takes interval); hence χ is
constant on the LRC-tight set.

**Verdict: largely VALID, and it sharpens THM-401.** Verified below; the one honest
qualifier is the word "regular."

## 1. The core mechanism — why only round (interval) is accessible

LRC's comparator is the **half-turn / circular-position** rule: at time `t` runner `i`
beats `j` iff `j` is within the forward half-circle of `i` (positions `v_i t`,
`v_j t`). An out-neighbourhood is therefore always a **contiguous arc** — the tournament
is **round** (locally transitive). So the LRC-accessible tournaments are exactly the
**round** ones (`A000016`), a vanishing fraction of all tournaments (`A000568`).

## 2. Interval vs Paley — verified

| m | interval `S` | Paley `S=QR` | interval round | Paley round | iso |
|---|---|---|---|---|---|
| 3 | `{1}` | `{1}` | ✓ | ✓ | **same** |
| 7 | `{1,2,3}` | `{1,2,4}` | ✓ | **✗** | different |
| 11 | `{1,2,3,4,5}` | `{1,3,4,5,9}` | ✓ | **✗** | different |

> **The interval circulant is round; Paley is NOT round for m≥7.** They coincide only at
> `m=3` (`QR(3)={1}` is an interval), and split for every Paley prime `m≥7`. So the LRC
> beat (round) can realise the **interval** circulant but **never** Paley — exactly the
> user's claim. Other circulants (non-interval connection sets) and non-circulant VT are
> likewise non-round → inaccessible.

## 3. χ is constant — verified, with the value

Dichromatic number (min transitive parts): **`χ(R_m) = 2` for m=3,5,7 (constant)**, while
**`χ(Paley_7) = 3`**. So:

> **χ ≡ 2 on the interval circulants** (the regular LRC-tight configs = the AP orbit);
> Paley sits at `χ=3`, confirming it is a *different, more tangled* object that LRC
> cannot reach. Since the regular tight set is the single iso-class `R_m`, **every
> iso-invariant — χ, `H`, score sequence — is constant on it.**

**Honest qualifier.** The full LRC-tight set is `R_m` (the AP/regular) **plus
non-regular sporadics** (e.g. the `n=8` rows), which are also round but not
vertex-transitive. The claim "= the interval circulant" is exact for the **regular**
tight configs; the sporadics are *other* round tournaments. χ is constant on the
regular part for sure (one class, χ=2); whether the sporadics are also χ=2 (i.e. χ
constant on the *whole* round/tight set) is the residual — round tournaments are
generally low-χ, and `R_m` is χ=2, so it is plausible but not proven here.

## 4. The use — this IS THM-401, at the tournament level

The idea cleanly separates the two algebraic faces (S586–S590) onto two tournament
roles:

```
ADDITIVE  (interval circulant, connection set {1,…,(m-1)/2})  =  the BEAT RULE (the tournament itself)
MULTIPLICATIVE  ((ℤ/n)^* / QR, inverse clocks)                =  the SYMMETRY (the unit witness orbit, HYP-2124)
```
- The **beat is additive** — the interval circulant — i.e. THM-401's additive modulus
  `2n−1` *is* the beat (the antipodal-shell / round structure).
- The **multiplicative units do NOT change the beat**; they act as the **symmetry**
  permuting the witness orbit (the primitive roots, S588). Paley would be a
  *multiplicative beat* — and it never occurs.

So "LRC is additive, not multiplicative" (THM-401) gets a tournament-theoretic proof:
the comparator is round = interval = additive; the multiplicative QR comparator is
non-round = inaccessible. **The additive face is the beat; the multiplicative face is
the symmetry.** This is the cleanest statement yet of the add/mult division of labour.

## 5. Further uses

- **Drastic restriction of the search space.** LRC explores only round tournaments
  (`A000016`), and the tight/regular core is the *single* class `R_m` (χ=2). The proof
  may restrict to round tournaments and, for the worry-set, to `R_m` + the round
  sporadics — never the `χ≥3` combinatorial tournaments.
- **Rigidity, restated.** "χ (and every invariant) constant on the regular tight set" =
  the worry-set is a *single rigid iso-class* — the cyclotomic/unit-orbit rigidity
  (HYP-2124) seen as iso-class invariance.
- **Paley is the foil.** Paley `m≥7` (`χ=3`, non-round, the maximally-multiplicative
  tournament) is precisely what LRC *cannot* build — a clean witness that LRC is an
  additive, not a quadratic-residue, phenomenon.

## 6. Honest status

- **Verified:** interval round / Paley not-round (m=7,11); interval≠Paley (m≥7),
  coincide (m=3); `χ(R_m)=2` (m=3,5,7), `χ(Paley_7)=3`.
- **Valid (rigorous):** LRC comparator → round only; the regular tight config = `R_m`;
  χ constant on the regular tight set; THM-401's additive face = the beat.
- **Residual (the qualifier):** whether χ is constant on the *whole* tight set
  (including non-VT round sporadics) — plausible (round ⟹ low χ), not proven here; and
  "regular round = uniquely `R_m`" assumed standard (VT round = interval circulant).

**Artifacts:** `04-computation/lrc_interval_vs_paley_s591.py`,
`lrc_interval_vs_paley_light_s591.out`. Builds on THM-401 (additive modulus), HYP-2124
(unit-orbit symmetry), S566/S570 (worry-set = round/regular), S588 (cyclotomic). New:
**HYP-2141**.

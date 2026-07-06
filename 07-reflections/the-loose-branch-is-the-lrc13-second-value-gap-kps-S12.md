# The loose branch is the LRC(13) second-value gap — and the "big pair ≥ 38" filter is pure Farey

**kind-pasteur-2026-07-05-S12 (HYP-4147).** A reframe of the LRC(14) loose branch,
after HYP-4137 killed every finite-template route. What remains is a clean spectral
statement about the *level below*, and its main "profile filter" turns out to be a
one-line Farey fact.

## The object

Peel the argmax killer from a primitive compressed covering Fin-13 family; the base
is 12 speeds. Write `M(base) = max_t min_k ‖v_k t‖` (the covering-min / max-min).
The dichotomy `TightLooseDichotomy` splits on `M(base)`:

- **tight**: base `= c·{1,…,12}` (a dilated AP), or
- **loose**: `M(base) ≥ 2/25` (a real 2/25-witness).

## Fact 1 — the floor is FREE (cited LRC(13))

Every 12 distinct speeds satisfy `M ≥ 1/13` by the Lonely Runner Conjecture at 13
runners — **which is cited as settled** (Sungkawichai–Trakulthongchai, ≤13). No
covering/compression needed. So `M(base) ≥ 1/13` costs nothing.

The AP `{1,…,12}` attains it: at `t = 1/13`, `{k/13 : k=1..12}` are all nonzero
residues mod 13, min distance `1/13`. The AP is the *equality case* of LRC(13).

## Fact 2 — the loose branch is exactly a SECOND-VALUE GAP

Given Fact 1, the dichotomy's real content is the jump:

> **`M(base) > 1/13  ⟹  M(base) ≥ 2/25`** for primitive compressed covering bases,
> and `M = 1/13` only for the AP.

Equivalently: the second-smallest value in the compressed-covering max-min spectrum
is `σ₂ = 2/25`, with an empty gap `(1/13, 2/25)` above the cited floor. `σ₂ ≤ 2/25`
is trivial (`{1,…,11,24}` attains `2/25` at `t = 2/25`). `σ₂ ≥ 2/25` — gap emptiness
— is the whole open content. Exhaustively confirmed for all **1,531,052** primitive
compressed covering 12-bases in `[1,26]` (`lrc_spectral_gap_kps_S12.c`): the AP is the
**unique** floor (count 1 at `1/13`), the gap is **empty**, and `σ₂ = 2/25` (count 2).

## Fact 3 — the "big pair ≥ 38" filter is PURE FAREY

`1/13` and `2/25` are **Farey neighbors**: `1·25 − 2·13 = −1`. So NO fraction with
denominator `≤ 37` lies strictly between them (their mediant is `3/38`). By grid
attainment (THM-592 / HYP-4108), `M(base)` is a value `p/q` attained at a merge point
`m/(v_i + v_j)`, so its reduced denominator divides `v_i + v_j`. Hence:

> a gap value forces reduced denominator `≥ 38`, hence a pair `v_i + v_j ≥ 38`.

This is exactly the fleet's merge-exclusion / `gap_forces_big_pair` filter (`38`,
THM-592, formalized in `LRCMergeExclusion`) — but it is not an ad-hoc computational
observation: it is the Farey-neighbor relation of the two spectrum values. The `d ≥ 3`
and `w_max ≥ 19` companions fall out the same way (`3/38` is the mediant; `38 = 2·19`).
Gap-emptiness for bases with `max ≤ 18` is therefore immediate (no pair reaches 38).

## Fact 4 — nothing FINITE closes the gap (why it is genuinely analytic)

Three natural finite/effective schemes all provably fail on the remaining (big-pair)
case:

- **Census / fixed template** (HYP-4137, S11): the witness modulus grows like
  `log(height)`; free-modulus witnesses are CRT-killable and a runner `≡ 0 mod L`
  blocks all pinned-only witnesses. No bounded modulus set suffices.
- **Measure / peeling** (this session): peeling one runner leaves an 11-sub whose
  `2/25`-lonely measure is `≈ 0.085 ≪ 4/25 = 0.16` for *every* base tested, loose or
  not — so "a loose 11-sub survives the 12th runner" never fires. The witness is an
  arithmetic point (a merge point), not a positive-measure event (cf. MISTAKE-072).
- **Finite harmonic gates** (S2 + S11): the `k=1,2` reciprocal/second-harmonic gates
  cover narrow bases and the `{1..11,24}` value, but the uncovered set is exactly the
  unbounded profile-survivor family — no finite set of gates reaches it.

So the loose branch is irreducibly a statement about the *arithmetic* location of the
best merge point, uniform over unbounded scale.

## The precise open target

> **CONJECTURE (compressed-covering LRC(13) second-value gap).** For every primitive
> compressed covering 12-speed base other than the dilated AP, `M ≥ 2/25`; equivalently
> `σ₂ = 2/25` and `(1/13, 2/25)` is spectrum-empty.

It is a **rigidity + gap** pair: (a) `M = 1/13 ⟹` AP (uniqueness of the LRC(13)
extremizer among compressed covering families); (b) `M > 1/13 ⟹ M ≥ 2/25` (the jump).
Both live one level down from LRC(14), in the extremal spectrum of the *proven* LRC(13)
— which is where the real analytic work now sits. Literature lead: the Lonely
Runner / view-obstruction "spectrum" and its gaps (Cusick, Schoenberg, Perarnau–Serra
"tight locus widely open") — the second-value gap at a specific level may be addressable
by their rigidity methods, which the ad-hoc profile never tapped.

## Why this reframe helps

The fleet's "gap-violator profile" (covering, spread, compression, big-pair, pinning,
ladder) is a pile of necessary conditions assembled to make a census tractable. Naming
the target as *the LRC(13) second-value gap* (i) shows the floor is already paid for by
citation, (ii) derives the load-bearing `≥ 38` filter from one Farey identity rather than
from THM-592's machinery, and (iii) points the remaining work at a known genre of
problem (extremal-spectrum rigidity) instead of an ever-larger census that HYP-4137
proved can never close.

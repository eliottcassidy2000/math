# The AP is self-protecting: only a multiple of 12 breaks it

**opus-2026-07-06-S103** (HYP-4356). A creative, clean piece of the (C)/hdich crux,
formalized (LRCAPProtection.lean, GREEN).

## The mechanism

The 12-runner Farey gap (C) reduces (S100) to the AP-completion step: adding a 12th runner
to the dilated AP `{1,…,11}` must give the tight completion or a loose family, never an
in-window value. The reason is almost embarrassingly simple once seen:

> At `t = 1/12`, every speed `s` with `12 ∤ s` sits at circle-distance `≥ 1/12` from `0`,
> because `dist(s/12, ℤ) = |s − 12m|/12 ≥ 1/12` (the nearest integer multiple of 12 is at
> least 1 away). The AP `{1,…,11}` is entirely made of such speeds — so it is LONELY at
> `1/12` with margin exactly `1/12`, and it stays lonely when you add ANY runner that is
> not a multiple of 12.

`1/12 > 2/25`, so `{1,…,11, v}` is LOOSE for every `12 ∤ v`. The ONLY way a 12th runner can
pull `M` down into the window is to be `≡ 0 (mod 12)` — to land exactly on `0` at `t = 1/12`
and destroy the protection. Those `v = 12w` are precisely the l = 1 lift stratum
(`LRCLiftRigidityRows`: `w = 1` is the tight completion `M = 1/13`, every other `w` gives
`M ≥ 2/25`). So the AP-completion step is two clean cases, both settled:

* `12 ∤ v`: loose at `1/12` (this file, `ap11_loose_of_not_dvd`, GREEN);
* `12 ∣ v`: the l = 1 lift rigidity (kernel rows, GREEN).

## Why this is the inductive core of (C)

(C)/hdich is "a 12-family with `M < 2/25` is the dilated AP." The induction (kps-S13's
recursive AP tower): if some 11-subfamily is tight (a dilated AP `c·{1..11}`, forced by
LRC(12) rigidity), then AP-completion forces the 12th runner to complete the AP (tight) or
the family is loose. The self-protection lemma is the "add-a-runner" heart of that step,
and it exposes the arithmetic cleanly: the AP's loneliness lives at the single point
`1/12`, and divisibility by 12 is the exact knife that cuts it.

## The wider echo

This is the "free action / divisibility" meta-mechanism (S100b) in miniature: the AP's
protection is a divisibility condition (`12 ∤ s`), and breaking it requires hitting the
`12ℤ` sublattice — the same "regularity forces divisibility forces a gap" shape that runs
through the whole project. And `12` here is `k` (the AP length is `k−1 = 11`, the protecting
denominator is `k = 12`); the lemma generalizes verbatim to every `k`, giving the universal
AP-completion behind the Farey ladder `j/(kj+1)`.

# The dominant/compressed split is a value split — the killer is an offset-forcer (M↓14/183) or a free-rider (M=1/13)

*mac-mini-2026-07-04-S45, **corrected and deepened S46**. The clean statement is: the `13×` line
(dominant vs compressed) separates killers that **kill the base optimum** (→ offset → `14/183`) from killers
that **ride safely through it** (→ `M = 1/13`). The compressed floor is `1/13`, not `7/89`.*

> **⚠ S46 CORRECTION.** S45 (below) claimed the compressed floor is `7/89 = {1..11,13,84}` (the "lcm(12,14)
> shadow"). That is only the **minimal-tightener / non-dilated** floor. The TRUE compressed floor is **`1/13`**:
> the **dilated** deep-well `c·{1,…,12} ∪ {182}` (`c≥3`) is compressed and covering with `M = 1/13` exactly
> (`< 7/89`). `7/89` is the lowest *non-dilated* compressed rung; `1/13` is the real floor. Both are `> 14/183`,
> so the covering-min is untouched — but the clean tight statement is `compressed ⟹ M ≥ 1/13`, not `≥ 7/89`.

## The real image: offset-forcer vs free-rider
A covering family's killer `a` (with `13∣a`) has, at the base's max-min `t*`, a phase `a·t*`. Everything turns
on whether that phase is an **integer**:
- **DOMINANT** killer (deep well `{1..12,182}`, `t*=1/13`): `a t* = 182/13 = 14 ∈ ℤ` → killer at `0` → it
  **kills the base optimum**, hiding must offset by `1/2379` → `M = 14/183`. *(Offset-forcer. Discharged: the
  giant runner is peelable, kps HYP-4087.)*
- **COMPRESSED** killer (dilated `{3,…,36,182}`, `t*=1/39`): `a t* = 182/39 = 14/3 ∉ ℤ`, `‖·‖=1/3 ≥ 1/13` →
  the killer **rides through** the base optimum untouched → `M = M(base) = 1/13`. *(Free-rider.)*

The dilation is exactly what turns the offset-forcer into a free-rider: spreading the base by `c` moves its
optimum from `1/13` to `1/(13c)`, and `182/(13c)` stops being an integer. So the compressed floor `1/13` is
the base's own value, with the killer contributing nothing — the family hides like its 12-runner base.

## Why `1/13` can't be beaten (the CRT peel)
Peel `v* = max`; `W = V∖{v*}` has `M(W) ≥ 1/13` (LRC(13)). The floor `1/13` is attained only when `W` is
LRC(13)-**tight**, i.e. a dilated AP `c·{1,…,12}` with optima `{k/(13c) : 13∤k}`. Can the killer be unsafe at
*all* of them (forcing `M<1/13`)? **No — CRT forbids it.** Primitivity `gcd(V)=1` with `gcd(W)=c` forces
`gcd(c,v*)=1`, so `v*/(13c)` carries the factor `c` **coprime to `13`**: `v*`'s safety is indexed mod `c`, the
base optima mod `13`, and `13 ⊥ c` makes them independent — some `k` is simultaneously a base optimum (`13∤k`)
and `v*`-safe. And `c>1` is guaranteed because `c=1` *is* the deep well (`182>13·12`, dominant, not compressed).
I built the sharpest adversary — `157·{1,…,12} ∪ {18382}`, the killer a near-multiple of `13·157` engineered
unsafe at `j=1,…,12` — and it lands at `M=1/13` exactly, saved by the optimum at `j=14`.

## Why it matters (the operational sentence)
kps's sole open leaf `hcomp` (compressed ⟹ lonely) is implied by the **tight, structural** bound
`compressed ⟹ M ≥ 1/13` — a `1/13 vs 1/14` margin, from LRC(13), with no razor-thin `14/183` anywhere. If the
peel route closes (loose-base width bound + the `n=13` tight locus), it discharges `hcomp` **uniformly** — the
compressed twin of opus's dominant peel. The two sides of the `13×` line, both peeled to LRC(13): that is the
whole covering side, without a census.

## The clean open question
Complete the compressed peel: (i) loose base — `W`'s good interval is wider than `v*`'s danger-arc spacing
(needs the compressed bound `v* ≤ 13·(2nd)`); (ii) the `n=13` tight locus — is it only dilated APs, or is there
a GW-analog that also needs the CRT argument? Close both and `compressed ⟹ M ≥ 1/13` is a theorem.

---
*(S45 original text below — the `7/89`/lcm-shadow framing is the non-dilated special case; superseded as the
global floor by the `1/13` correction above.)*

## [S45] The two extremizers are the same construction, one lcm apart
The Lean covering dispatch (opus `covering_lonely_of_dominant_or_compressed`) cuts every covering family at
`largest = 13 × second-largest`:

| branch | how it covers the hard pair via ONE killer | killer | tight base | `M` | CF |
|---|---|---|---|---|---|
| **DOMINANT** (discharged, kps HYP-4087) | `{13,14}` together | `182 = lcm(13,14)` | `{1,…,12}` | `14/183` | `[0;13,14]` |
| **COMPRESSED** (open leaf `hcomp`) | `{12,14}` together | `84 = lcm(12,14)` | `{1,…,11,13}` | `7/89` | `[0;12,7]` |

Both are *one big-lcm killer + the tightest base that leaves*. The deep well spends its one killer on the
**hardest** pair `{13,14}`; the compressed extremizer can't — and that "can't" is the whole point.

## Why compressed is *forced* one lcm down
To cover `{13,14}` with a single runner you need `lcm(13,14)=182 ∣ killer`. With a base whose second-largest
is `≤ 12`, that killer is `≥ 182 > 156 = 13·12` — **necessarily dominant**. So *the moment a covering family
tries the deep-well move (one killer for `{13,14}`, base `{1..12}`), it lands in the dominant branch.* A
compressed family is barred from `182`; its best single-killer move is the **next** lcm down that still
covers `14`: `lcm(12,14)=84` (`84 ≤ 13·13`, compressed), covering `{12,14}` and freeing `13` as its own
small killer over `{1,…,11}`. That is `{1,…,11,13,84}`, `M = 7/89`.

So the `14/183`↔`7/89` gap is not numerology: it is `lcm(13,14)` vs `lcm(12,14)`, i.e. the **cost of being
denied the `182` move.** The razor-thin rung belongs to the branch that *can* play `182` — which is exactly
the branch already closed (a giant killer is peelable, kps HYP-4087). The open branch is the one that
*can't*, and it sits a definite `35/16287` higher.

## Why this matters for the endgame (the one operational sentence)
**All remaining razor-thinness is on the closed side of the `13×` line.** kps's sole open obligation `hcomp`
(compressed covering ⟹ lonely) never has to clear the sharp `14/183` — its own families are `≥ 7/89`, so any
loss-of-constants argument (folding, pigeonhole, census) that clears a *margin* suffices. This is why opus-S72
can close the `m=2,f=2` residual at `M ≥ 1/12` and klein-S129 can call the target "non-sharp": the sharp
point was quarantined into the dominant branch by the very construction that makes it sharp.

## The clean open question it leaves
Is there a **compressed peel** mirroring opus's dominant peel: *compressed covering ⟹ `M ≥ 1/13`* (the
12-runner LRC floor), on the grounds that being denied the single `{13,14}` killer costs one runner-slot,
dropping the family to the `(n-1)`-runner problem? That would discharge `hcomp` from `LRC(13)` with a clean
margin — the compressed twin of the dominant peel — and finish the covering side without ever touching
`14/183`. (`7/89 > 1/13 > 14/183`, so even the weaker `≥ 1/13` clears the leaf.)

See HYP-4089; and the results it confirms: HYP-4090 (klein), HYP-4087 (kps dominant peel), HYP-4091
(kps/opus compressed leaf). Ladder framing: [[the-covering-min-is-an-ostrowski-ladder-and-the-ap-and-deep-well-are-its-ends]].

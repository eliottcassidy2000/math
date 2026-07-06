# One mechanism, many surfaces: the LRC gap is a forbidden band, and the whole project has been circling it

**opus-2026-07-06-S100b** (HYP-4306 continued). A breadth-first fan-out across the repo's
Farey-adjacent corpus (6 parallel explorers over ~275 mediant + ~192 three-gap/renormal
reflections) returns one conclusion: the project has been discovering ONE mechanism under
a dozen surface names, and the LRC(14) spectral gap is its cleanest instance.

## The one mechanism

> **A regularity (equal spacing / a free group action / continued-fraction extremality)
> forces an arithmetic DIVISIBILITY, which QUANTIZES an achievable quantity onto a LADDER
> whose inter-rung windows are FORBIDDEN.**

It appears, provably, at least nine times:

| surface | regularity | divisibility forced | ladder / forbidden band |
|---|---|---|---|
| density quantization (THM-412, S703) | free rotation `μ` on norm-shells | `w ∣ count` | density in `(w/2)ℤ`; 4,5 forbidden |
| totient ladder (THM-416, S706) | free `μ(K)` in any dim | `w ∣ r(D)` | quantum `w/2`; nontotient dims drop |
| **LRC Farey minimax (HYP-4306)** | grid attainment `t=a/(v_i+v_j)` | denom ∣ a pair-sum | rungs `j/(kj+1)`; `(1/(k+1),2/(2k+1))` empty |
| lonely-measure (THM-522) | danger-arc endpoints | `L ∈ (1/(n·lcm))ℤ` | `L>0 ⟹ L ≥ 1/(n·lcm)` |
| 12-runner AP rigidity (klein-S140) | optimality ⟹ equal spacing | `Q = 13r` | only `1/13` realized; gap void |
| covering = congruence subgroup (HYP-3553) | level-`N` structure | `Γ₀(N)` index | floor depends only on `N` |
| Ostrowski covering-min (HYP-4078) | continued-fraction extremality | convergents only | rungs `k/(k(n-1)+1)` |
| modular units (HYP-3793) | unit-orbit extremizers | `φ(N)` atoms | maximizers = `(Z/N)*` |
| Farey-cell-void (THM-622) | Farey-neighbour `1/13, 2/25` | reduced denom ≥ 38 | cell interior void; mediant `3/38` |

The project already PROVED the last six as separate theorems. What the S100 Farey-ladder
reframe adds is that they are one statement seen at different scales.

## The three skeletons that keep reappearing (and how they fit)

Everything above decomposes into three interlocking structures, each of which the repo
found independently:

1. **Three-distance / Steinhaus — the RIGIDITY ENGINE.**
   `the-LRC-for-the-AP-IS-the-three-distance-theorem` (opus-S630): the LRC for the AP is
   literally Sós–Steinhaus on `{i·t}`. Extremal configs are `{kα mod 1}` progressions,
   which Steinhaus forces to have ≤ 3 gaps. VERIFIED: the tight `{1..k}` at `1/(k+1)` has
   2 gaps; the second-best `{1..k-1,2k}` at `2/(2k+1)` has 3 gaps. The open core (HYP-2913,
   the `±units` covering condition `g(n) ≤ 3`) is EXACTLY what the crux hinges on — proving
   the extremizer is `{kα}`-structured.

2. **Stern-Brocot / continued fractions — the ADDRESS SYSTEM.**
   `the-loneliness-integral-limit...regimes-are-the-Stern-Brocot-tree` (opus-S630): the
   Mode-A recursion `O(n)-O(n-1) = φ(n-1)` births exactly the new Farey mediants; the cap
   kernel `K(a,b)=g(a,b)/(7ab)` walks the continued fraction of `a/b` with breakpoints at
   convergents (HYP-3230). The "hard" covering cap and the "easy" additive census are the
   SAME recursion at different scales.

3. **Free action / divisibility — the ALGEBRAIC HEART.**
   At `k=12`, `13` is PRIME, so the tight speeds `{1..12}` ARE the unit group `(Z/13)*` —
   one free multiplicative orbit (verified). Perturbing the orbit breaks the free action
   and jumps to the next Farey rung `2/25`. This is `density-quant`/`totient`/`modular-units`
   realized as the LRC gap's rigidity. (For composite `k+1` the tight set is all nonzero
   residues, not the unit group, and the rigidity routes through CRT instead — the reason
   the mechanism is CLEANEST at the primes, and `n=14 → 13`-prime is a lucky level.)

## The doubly-prime peel (the fan-out's sharpest new insight)

The fan-out surfaced THM-539 / HYP-2623 (kps-S9, `lrc-spectral-gap-dips-along-primorials`):
the collar just above the floor `1/(k+1)` carries Stern-Brocot dips `a/(a(k+1)-1)` that
descend TOWARD the floor, and how deep they go is governed by the PRIME FACTORIZATION of
`k-1` — `a=3` dips appear only when `6 ∣ (k-1)`, `a=4` only when `30 ∣ (k-1)`, etc. (a
primorial gate; the large speed `a(k-1)` lands on integers at every clock dividing `k-1`,
so more prime factors ⟹ deeper dip, saturating at `a=4`).

This EXPLAINS the peel that the whole reduction performs. The raw `n=14` problem is the
13-runner level (`k=13`), where `k-1 = 12 = 2²·3` is COMPOSITE — so the window
`(1/14, 2/27)` carries a real dip at `3/41` (the Farey neighbour flagged in
`anatomy-of-a-tight-runner-set` condition 7). The project peels to the 12-runner rigidity
level (`k=12`), where `k-1 = 11` is PRIME — no `a≥3` dips, so the window `(1/13, 2/25)` is
CLEAN (only tight + mediant). And simultaneously `k+1 = 13` is PRIME, so the tight speeds
`{1..12}` are exactly the unit group `(Z/13)*` (verified). **The peel is not bookkeeping:
it moves the problem onto the DOUBLY-PRIME level where both `k-1 = 11` (no primorial dips,
window clean) and `k+1 = 13` (tight config = free unit orbit) are prime.** `13` being a
lucky doubly-adjacent-prime-flanked value is why the `n=14` LRC rigidity is provable at all.
(Verified: my 12-runner candidate families all land exactly at `1/13`, `2/25`, or above —
none in the window; the `3/41` dip at the 13-runner level is THM-539's, cited not
re-derived.)

## The tower of ladders (the connective synthesis)

There is ONE ladder per runner-count: `ladder(k) = { j/(kj+1) }`. The LRC(14) reduction
DESCENDS this tower — `n=14` (13 runners, ladder(13): `1/14, 2/27, …, 14/183` deep well)
peels to the `n=12` rigidity (12 runners, ladder(12): `1/13, 2/25, …, 14/169` deep well).
Every "magic number" is a rung of some ladder(k):
- `1/13` = ladder(12) rung 1 (tight) = my S99 projection floor;
- `2/25` = ladder(12) rung 2 (mediant) = the dichotomy threshold;
- `3/38` = the first Stern-Brocot mediant-CHILD inside the window = mac-mini's HYP-4252 cell;
- `14/169` = ladder(12) rung 14 = the deep well (my S77);
- `14/183` = ladder(13) rung 14 = the n=14 construction deep well;
- `1/6` = ladder(5) rung 1 = mac-mini's 7-spread census infimum (the ≤5-class ceiling).
The census/torus machinery is verifying that ladder(12)'s FIRST inter-rung window is empty,
one level of a self-similar tower.

## The proof architecture the fan-out suggests

The renormalization cluster surfaced a coherent inductive skeleton (not yet a proof):
- **Base:** the fixed point of the scale flow is the AP (`union-bound-dies-at-7`,
  HYP-3900); the renormalization sends a deep `n/2`-cluster to its difference core (an AP),
  contracting ~`n/2` runners per level, terminating in ≤ 2 levels — so the infimum is a
  FINITE census at height 1.
- **Rigidity step:** three-gap forces the extremizer to be `{kα}`; the `±units` condition
  (HYP-2913) forces its residues to be the unit orbit; the mediant is the smallest
  Farey step off it. This is `second-best = mediant`, hdich's twin at every `k`.
- **Descent:** the cyclotomic moment-ladder × 2-adic fold (kps-S31aq) bottoms both
  recursions at the Fejér/de Moivre kernel; the `p`-adic tree (n=14 is rank-1, one 7-adic
  tower) gives a single-tower descent; the cluster-gcd ladder + k-stratification bound the
  `k≥2` heights so only `k=1` needs the census.

Induction on the runner-count `k`, with three-gap rigidity as the step and renormalization
as the contraction, is the natural shape — and it is the same shape as the density-ladder
proof (free action → divisibility → gap), only with the free action replaced by
Steinhaus rigidity.

## Niche threads worth reviving (flagged by the fan-out)

- **Beatty–Pell "one layer below"** (`beatty-pell-crossover-word`): the proof-bearing
  coordinate is one layer under the visible token (address + carry, not just the word).
  The LRC analog: a "blocked denominator" row must record shell class + carry residue +
  endpoint atom, not just "blocked" — a discipline that keeps recurring.
- **Eisenstein cusp = three-distance** (`the-eisenstein-cusp-dichotomy-is-the-three-distance
  -theorem`): ties `Q(√-3)` cusps to the gap; never connected to the ladder.
- **7/89 binding modulus** (`seven-over-89-...anti-Littlewood-floor`): a Stern-Brocot ray
  treated as a curiosity — but `89` is a prime core of the window's deep Stern-Brocot
  descendants (`8/101, 7/89` appeared in my window enumeration).
- **The 21-resonance at dimension 12** (S706): quantum 21 = `3·7` = the forbidden
  Hamiltonian-path count, via `Q(ζ_42)`; flagged as "a resonance, not a theorem."

## The takeaway

The gap is not a wall the project keeps failing to climb; it is a FORBIDDEN BAND that the
project keeps rediscovering the physics of. Density quantization, the Ostrowski covering
ladder, the cap-kernel three-gap recursion, the congruence subgroup, the modular units,
and the Farey-cell void are the same forbidden-band mechanism read on different instruments.
Naming it once — regularity ⟹ divisibility ⟹ quantized ladder ⟹ empty window — turns the
crux from "prove no covering family lands in (1/13, 2/25]" into "prove the first rung of the
LRC Farey ladder is empty," a statement whose engine (Steinhaus three-gap), coordinate
system (Stern-Brocot), and heart (`(Z/13)*` free orbit) are all already in the repo.

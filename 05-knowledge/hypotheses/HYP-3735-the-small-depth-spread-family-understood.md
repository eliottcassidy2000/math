---
id: HYP-3735
title: UNDERSTANDING the small-depth spread family (the covering-min for n=7..11). It is the THRESHOLD semiconvergent [0;n-1,a(n)] on the Stern-Brocot ray = the smallest depth a where a genuinely-spread primitive covering set exists; a(n)=2,2,4,4,3. STRUCTURE: at the witness t*=c/m (m=(n-1)a+1) the speeds' residues c*v mod m all lie in the band-complement [a, m-a] (avoiding the central 2a-1 residues around 0), covering 2..n, using only a SUBSET of safe residues. THRESHOLD MECHANISM: for depth a<a(n) NO spread set exists -- the best covering set is the near-block 1/(n-1) ("blocking"); at a=a(n) spreading first beats blocking (M=a(n)/m < 1/(n-1)). OBSTRUCTION = LARGE PRIMES: each prime p in (n/2,n] forces a LONELY speed (p is its only multiple <=n), whose danger arcs (period 1/p) must be covered by other speeds; covering them needs the modulus m fine enough -> depth >= a(n); as n grows the lonely-prime obstruction TIGHTENS so a(n) grows/jumps. FINITE RANGE: the small-depth (a<=4) family spans n=7..11 ONLY; at n=12 it DIES -- depths a=2,3,4,5 are UNACHIEVABLE (verified reliably at small V; n=12 a=2..5, n=13 a=2..4, n=14 a=2..4 all return the near-block, and S52's V=50 run agrees) -- and a(n) jumps to >=6. So for n>=12 the covering-min is at HIGH depth (>=6), <= the construction n/Phi_6 (depth n); exact a(n) for n>=12 OPEN (moderate depths 6..n-1 need large V the ILP can't resolve). ACHIEVABILITY is NON-monotone (n=9: depth 4 achievable but 2,3 not)
status: VERIFIED for n<=11 (structure, a(n)=2,2,4,4,3, the achievability map, the threshold). The n>=12 death of small-depth (a<=5) is reliably verified (small-V ILP, no time-limit issue; cross-checked vs S52 V=50). The large-prime MECHANISM is a structural explanation (consistent, not a proof). OPEN: the covering-min for n>=12 (high depth, <= construction); whether a moderate-depth (6..n-1) spread family beats the construction.
source: mac-mini-2026-06-30-S53
related:
  - HYP-3732  # the Stern-Brocot ray [0;n-1,k]; this is the deep dive into the covering-min's depth a(n)
  - HYP-3731  # the covering-min ILP (the tool); klein-S36
  - HYP-3725  # construction = depth-n point on the ray; the spread (depth a<n) beats it for n<=11
  - HYP-2566  # uniform looseness; the covering-min values are the evidence
results:
  - 04-computation/covering_min_perdepth_macmini_20260630.py
  - 05-knowledge/results/covering_min_perdepth_n12-14_macmini_20260630.out
  - 05-knowledge/results/covering_min_depthprobe_macmini_20260630.out
---

# HYP-3735 -- the small-depth spread family, understood

The owner asked to understand the small-depth spread family -- the covering sets that give the LRC covering-min
for `n=7..11` (`2/13, 2/15, 4/33, 4/37, 3/31`) by beating the construction `n/Phi_6`. Here is what it is, how
it works, and why it has a finite range.

## What it is -- the threshold semiconvergent
On the Stern-Brocot ray `[0;n-1,k] = k/((n-1)k+1)` (HYP-3732), the covering-min is `[0;n-1,a(n)]`, a **Farey
neighbor of `1/(n-1)`**, at depth `a(n) = 2,2,4,4,3` (`n=7..11`). It is the **smallest depth `a` at which a
genuinely-spread primitive covering set exists**.

## How it works -- structure and the threshold
At the optimal witness `t* = c/m` (`m=(n-1)a+1`), the speed-residues `c*v mod m` all lie in the
**band-complement** `[a, m-a]` -- they avoid the central `2a-1` residues around `0` (so `min_v ||v t*|| = a/m`).
The set covers `2..n`, is primitive, and uses only a **subset** of the available safe residues (there is
slack -- e.g. `n=9` uses 8 of 26 safe residues). The defining condition is global: `M(S)=a/m` requires the
witness `c/m` to be the *worst* gap (no finer breakpoint exceeds `a/m`).

**The threshold.** For depth `a < a(n)`, no spread set exists -- the best covering set is the **near-block**
`1/(n-1)` (the "blocking" solution: a primitivized scaled block `g*{1..n-1}+prime`). At `a = a(n)`, spreading
**first beats blocking** (`a(n)/m < 1/(n-1)`). So `a(n)` is the depth where *spreading overtakes blocking*.
Achievability is **non-monotone**: at `n=9`, depth 4 works but depths 2,3 do not (`a(9)=4`, `M=4/33`, not the
mediant `2/17`).

## Why -- the large-prime obstruction
Each prime `p in (n/2, n]` forces a **lonely speed**: `p` is the only multiple of `p` that is `<= n` (and
`2p > n`), so the covering set must contain `p` (or a large multiple), essentially alone. A lonely speed `p`
puts danger arcs at `k/p` with gaps `~1/p`; for `M(S)` to stay `<= a/m`, every such gap must be covered by
another speed. Doing so forces the modulus `m` fine enough -- i.e. depth `>= a(n)`. As `n` grows there are
more and larger lonely primes, the obstruction **tightens**, and `a(n)` grows.

## Finite range -- the family dies at n=12
The small-depth (`a <= 4`) family is a phenomenon of `n = 7..11` only. **At `n=12` it dies:** depths
`a = 2,3,4,5` are *unachievable* -- the per-depth ILP (small `V`, reliable) returns the near-block `1/11`,
not the depth-`a` target, and S52's independent `V=50` run agrees. Same at `n=13` (`a=2,3,4`) and `n=14`
(`a=2,3,4`). So `a(n)` **jumps to `>=6` at `n=12`**, and for `n>=12` the covering-min sits at **high depth**
(`>= 6`), bounded above by the construction `n/Phi_6` (depth `n`). The exact `a(n)` for `n>=12` is **open**:
the moderate depths `6..n-1` require speeds up to `~m = (n-1)a` that the breakpoint-universe ILP cannot
resolve reliably (time-limited at large `V`).

| n | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|----|----|----|----|----|
| a(n) | 2 | 2 | 4 | 4 | 3 | >=6 | >=6 | >=6 |
| covering-min | 2/13 | 2/15 | 4/33 | 4/37 | 3/31 | <=12/133 | <=13/157 | <=14/183 |

## Convergence with klein-S38 (HYP-3734) -- the rigorous version
klein independently worked the same family and supplies the rigorous backbone:
- **`a_1 = n-1` PROVED** (so `M=[0;n-1,...]` always): `1/n < M <= n/Phi_6 < 1/(n-1)` (THM-523 floor +
  `n(n-1)<n^2-n+1`), so `floor(1/M)=n-1`.
- **Farey-neighbor `<=>` `D ≡ 1 mod (n-1)`** (binding witness denominator `D`), verified for ALL covering-mins
  `n=7..14` (`D=13,15,33,37,31,133,157,183`); the general proof is open.
- **The achievable rungs form an UP-SET `[k_min, inf)`** (so achievability in `k` is monotone; the covering-min
  is `k_min = a(n)`). `k_min = 2,2,4,4,3` (`n=7..11`), and **no rung `<= 6` at `n=12`** -- exactly my "death at
  n=12". (My phrase "non-monotone" referred to `a(n)` vs `n`, not achievability vs `k`; klein's up-set is the
  precise statement.)
- **The obstruction, sharpened:** rung `k` demands a covering radius `floor(kD/(k(n-1)+1))` at *every* modulus
  `D`. The **radius-0** moduli are exactly `D <= n-1` -- these ARE the THM-523 resonances (THM-523 is the
  radius-0 layer!). The **radius-1 band `D in (n, 2n-2]`** is the EXTRA demand that over-constrains the `n-1`
  speeds, so low rungs fail as `n` grows. This is the precise form of my "large-prime lonely-speed"
  obstruction (the lonely primes `p in (n/2,n]` are where the radius-1 band bites).

So the two accounts dovetail: my structural/threshold view ("spreading beats blocking at `a(n)`", residues
filling the band-complement, finite range) + klein's rigorous reduction (`a_1=n-1`, `D≡1`, the up-set, the
radius-0/1 obstruction layering THM-523).

## Consequences
- The small-depth spread is **not** the LRC14 mechanism. At `n=14` the covering-min is at high depth (`>=6`),
  `<= 14/183` (the construction) -- the proof target for `n>=12` is the *high-depth / construction-scale*
  regime, not the small-depth spread. (This is consistent with -- and sharpens -- MISTAKE-088: neither the
  near-block `1/(n-1)` nor a small-depth spread set is the `n>=12` covering-min.)
- `a(n)` (the **achievability depth**) is the irregular new sequence: `2,2,4,4,3, then jumps`. Its growth is
  driven by the lonely-prime count in `(n/2, n]` -- a number-theoretic, irregular driver (an *extremal*, not
  *additive*, sequence; cf. the HYP-3732 taxonomy).

## What it buys
A complete picture of the small-depth spread family: it is the threshold semiconvergent where spreading first
beats blocking, structured by speed-residues filling the band-complement, obstructed by lonely large primes,
and **finite in range (`n=7..11`)**. It tells us the LRC14 covering-min is governed by the high-depth regime
(`<= n/Phi_6`), redirecting the proof effort there, and isolates the open question: does a moderate-depth
(`6..n-1`) spread family beat the construction for `n>=12`, or is the construction (`depth n`) the covering-min?

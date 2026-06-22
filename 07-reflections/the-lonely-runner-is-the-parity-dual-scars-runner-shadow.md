# The Lonely Runner is the parity-dual scar's runner-shadow

*kind-pasteur-2026-06-22-S31i. Applying the parity-dual clique-scar structure to LRC(14).
Third of a trilogy with `seven-is-the-unique-forbidden-clique-value.md` and
`the-parity-dual-clique-scars.md`.*

## The structure we built

The clique spine `I(K_r,2) = 2r+1` carries three selection rules:
- **tournaments** realize every clique *except* the odd `K_3` (forbidden conflict graph, forces
  `C_5`) -> the UNIQUE permanent gap at `7`;
- **even graphs** realize *only the odd* cliques (parity: `K_r` even iff `r` odd) -> gaps at the
  even cliques `5,9,13`, and they HEAL the tournament's `7`-scar;
- `7` is the unique value where orientation-*realizability* and degree-*parity* disagree.

## The Lonely Runner factors the scar

`LRC(14)` lives at `14 = 2 * 7`, and these are exactly the two scar factors:

> `14 = (parity factor 2) x (the unique forbidden odd-clique value 7 = I(K_3,2))`.

The repo's whole LRC(14) reduction runs through **7 inner sectors**, obtained by the
**complement-halving** `x -> -x` (the `Z/2` that takes `14` to `7`, THM-280/549). That `Z/2` is
the SAME involution as the even graph's "forget orientation": the apex prime `7` is the
odd-clique-scar side, the halving `2` is the parity side. The threshold `1/14 = 1/(2*7)` is the
product of the two scars.

## The obstruction IS the cyclic / K_3 content (grounded)

By HYP-2605 the Lonely Runner is the winding tournament `T(x)`: loneliness at phase `x` <=> `T(x)`
has a scale-`1/7` LOCAL SINK (a `>1/7` empty arc). A sink is HIERARCHICAL (few 3-cycles);
no sink is CYCLIC (many 3-cycles), the direction of the purely-cyclic `K_3` apex. Computed on AP
clusters (`lrc_winding_cyclic_scar_kps.py`):

| k | lonely fraction | mean #3-cycles (lonely) | mean #3-cycles (non-lonely) | ratio |
|---|---|---|---|---|
| 8 | 0.94 | 15.9 | 20.0 | **1.25x** |
| 10 | 0.78 | 32.3 | 39.3 | **1.22x** |
| 12 | 0.57 | 54.7 | 68.4 | **1.25x** |

**Non-lonely phases carry ~1.25x the cyclic content of lonely phases**, robustly. So the Lonely
Runner OBSTRUCTION (the non-lonely set whose measure must be killed) is precisely the winding
tournament pushed toward the cyclic / odd-clique `K_3` direction -- the same `K_3` that scars the
H-spectrum at `7`. This makes the owner's old lead literal: *non-loneliness is the odd-clique
direction of the conflict structure that forbids `H in {7,21}`.*

## Why 14, and why uniquely hard

`7` is the UNIQUE permanent prime gap (the one forbidden odd clique). Therefore `14 = 2*7` is the
UNIQUE parity-doubling of the scar -- the one Lonely Runner case sitting on the orientation-forbidden
odd clique. The proven cases `<=13` and the other composite cases `>=15` are hard for ordinary
reasons (compositeness, size); `14` is hard because it is the runner-shadow of the single defective
prime. The Lonely Runner, the H-spectrum, and the even-graph metagraph are three projections of one
object -- the `K_3 <-> C_5` odd-clique/odd-hole imperfection at the apex prime `7` -- and `LRC(14)`
is its shadow on the circle, read at scale `1/(2*7)`.

-> HYP-2605, HYP-2878, HYP-2879, HYP-2880, THM-200, THM-280,
`seven-is-the-unique-forbidden-clique-value.md`, `the-parity-dual-clique-scars.md`.

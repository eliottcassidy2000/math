# The Kakeya-needle obstruction to M = 3/38: medium-modulus needles cover the mod-38 band

*boxeph-2026-07-19-S124. Owner: work the 3/38 residue system for achievability; think Kakeya needles.
Result: a concrete obstruction map for the depth-minimal gap value 3/38. Filling the mod-38 safe band (a
necessary condition for M=3/38) makes a family LOOSE — the q=38 hole is real but never the deepest; a
deeper hole always sits at a MEDIUM modulus q'∈(13,38), giving M≥2/25. The obstruction is a Kakeya-style
NEEDLE-COVERING: no single medium modulus universally beats the q=38 hole (a family can block the mod-19
needle by containing 19), but the UNION of medium needles covers every band-filled covering family — you can
evade any one, not all. This builds the (previously unbuilt) Kakeya bridge to the q=38 maximizer. Not a
proof; the escape tail is untouched. Verified S124.*

## The target and the residue system

3/38 is the **unique** numerator-3 value in the open gap (1/13, 2/25) (boxeph-S123), the depth-minimal
mediant. A family at M=3/38 must (opus-S117 + S123):
- be **covering** (since 3/38 < 1/12, the S121 cascade forces a multiple of every q=2..12);
- have a **determinant-3 maximizing pair at s=38**, active at residues ±3;
- place all 12 residues mod 38 in the **safe band [3,35]** (avoid the hole {0,±1,±2}), with 3/38 the
  **global** deepest hole.

The Kakeya lens: each modulus q is a **needle direction**, and the witness t=p/q measures the family's
deepest hole along it. M=3/38 demands blocking every needle to depth ≤ 3/38, with q=38 the deepest at
exactly 3/38. The obstruction is that the needles cannot all be blocked at once.

## (A) Band-filling alone makes the family loose

Constructed covering families with residues mod 38 in [3,35] and the det-3 pair (3,35): the q=38 hole is
genuinely at 3/38, but the **true M is much larger**, realized at a medium modulus:

| family (covering, band-filled, pair (3,35)) | hole @1/38 | true M | at q' |
|---|---|---|---|
| `{3,5,7,8,9,10,11,12,13,15,21,35}` | 3/38 | 1/8 | 24 |
| `{3,5,7,8,9,10,11,12,17,21,24,35}` | 3/38 | 5/29 | 29 |
| `{3,7,8,9,10,11,12,15,20,22,33,35}` | 3/38 | 4/23 | 23 |
| `{3,5,7,9,10,11,12,14,16,24,33,35}` | 3/38 | 5/21 | 21 |

The q=38 hole is never the deepest — a **medium modulus q'∈(13,38)** always supplies a deeper hole, pushing
M above 2/25 (the loose side). To force M = 3/38, one would have to *also* close every medium-modulus hole,
which collapses the family to the AP {1,…,12} (M=1/13). The gap value 3/38 is squeezed out between "loose
(medium hole, M≥2/25)" and "the AP (M=1/13)."

## (B) The mod-19 parity split (38 = 2·19)

At t=m/38 an even speed `2u` has `‖2u·m/38‖ = ‖u·m/19‖`, a **multiple of 1/19 = 2/38**. Since the band
demands distance ≥ 3/38 and the multiples of 1/19 skip from 2/38 to 4/38, a band-satisfying even speed sits
at `≥ 2/19 = 4/38 > 3/38`. So:

> **The 3/38 hole is carried entirely by ODD speeds** (the pair (3,35) is odd, sum 38); even speeds are
> pushed to `≥ 2/19` and never touch it.

Consequently the **mod-19 needle** `t=n/19` (min-distance is a multiple of 1/19) must be blocked to `≤ 1/19`
for every rotation n — otherwise a `2/19 = 4/38` hole beats 3/38 and M ≥ 2/19. Blocking it forces the speed
residues mod 19 to cover all `±`-unit classes (or the family to contain a multiple of 19), with even speeds
barred from `±1`. This is the finite mod-19 feasibility system macmini-S27 posed, now with the parity and
determinant structure attached.

## (C) The needle-covering: adaptive, not universal

No **single** medium needle universally beats the q=38 hole. In particular a family can **block the mod-19
needle by containing 19** (then its mod-19 witness is 0). But such families are caught by a **different**
medium needle. Over **1066** band-filled covering families tested (subsets of [3,35] with the pair (3,35)):
- the mod-19 needle alone beats 3/38 for **772/1066** (72%, those not blocking it);
- **271/1066** contain 19 and so evade the mod-19 needle — but are deeper at another modulus;
- **the union of medium needles q'∈[14,37] beats the q=38 hole for ALL 1066/1066** — every band-filled
  covering family has M > 3/38.

This is the Kakeya signature: not one needle in a fixed direction, but a family-**adaptive** needle — you
can defeat any chosen direction, yet the union of directions leaves no gap. It matches codex-S16's reading
of the n=12 object as an **adaptive graphic rank**, not a fixed dimension: the "direction bush" that beats a
family varies with the family. Here the bush is the medium moduli (13,38), and blocking one (e.g. adjoining
19) merely shifts the deepest hole to another.

## Why the needles are unavoidable — the covering tension

The needles are forced by **covering + tightness**. Covering (kps-S45: M<2/25 ⟹ a multiple of every
q∈{7,…,12}) plants many small/medium multiples; each multiple `q·w` interacts with the mod-q' structure at
medium q', and keeping the deepest hole at a *large* modulus (38) while covering all the medium ones
over-constrains 12 speeds. This is boxeph-S102's finding that the crux forces at **medium moduli
13<q'<q**, made concrete at q=38: the medium needles are exactly where the forcing lands, and they are
adaptive.

## Honest status

- **New:** the Kakeya-needle lens is now pointed at the q=38 maximizer (the bridge was unbuilt per the
  survey). The concrete obstruction — band-filling ⟹ loose via a deeper medium-modulus needle; the mod-19
  parity split isolating odd speeds; the adaptive needle-covering — is a genuine mechanism for why 3/38 is
  squeezed out.
- **Not a proof.** 3/38 is verified unachievable for all primitive covering bases in [1,26] (kps-S12); the
  needle-covering is observed on band-filled families with small elements. The **unbounded-modulus escape
  tail** (macmini-S36 / HYP-4667: `{i+L·k_i}`, L=lcm(2..Q₀)) uses large elements and is untouched — those
  families approach 2/25⁺ from above, never entering the gap. The obstruction is mapped, not closed.
- **Consistent with the analytic-core status:** the reason 3/38 resists a finite proof is exactly that the
  beating needle is adaptive and its modulus is unbounded across the escape tail.

Cross-links:
[[the-determinant-stratified-gap-numerator-two-is-excluded-and-3-38-is-the-depth-minimal-target-boxeph-S123]],
[[kakeya-needles-are-an-adaptive-graphic-rank-not-a-dimension-analogy-codex-S16]],
[[the-free-sieve-window-lever-is-refuted-the-crux-forces-at-medium-moduli-not-the-sieve-boxeph-S102]],
[[the-finite-covering-is-incomplete-escape-at-the-next-prime-macmini-S36]],
opus-S117 (q=38 residue system), macmini-S27 (mod-19 parity descent), kps-S45 (covering ⟹ mult of 7..12),
`lrc14_kakeya_needles_3over38_boxeph_S124.py`.

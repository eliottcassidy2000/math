# The AP-stability theorem IS the covering rigidity — shadow gaps, and the THM-523 boundary

*klein-2026-07-13-S296. Owner: attempt the AP-stability theorem directly. I attempted it, found the
mechanism, and the attempt resolved into a clean characterization rather than a proof: the AP-stability's
difficulty is **exactly** the covering condition. The clean modular witness `t=1/q` works precisely for
non-covering sets (it is THM-523); covering removes every such witness, and what remains — the good set
reaching the middle through subtle resonance-shadow gaps — is the LRC extremal rigidity itself.*

---

## The formulation

`L({1}∪C) > 0 ⟺ G(C)` reaches the middle `[1/14,13/14] ⟺` the bad sets `∪_{c∈C} D_c` fail to cover the
middle. The AP-cluster `{2..13}` exactly **tiles** the middle (`L=0`, the rigid extremal); a covering `C`
must leave a **gap** (a good point). So the stability theorem is: *covering forces a gap in the
middle-tiling.*

## The mechanism: resonance-shadow gaps

Where do the middle good-arcs of near-AP covering clusters sit? At the **resonances `j/k` of a released
speed** (verified `NG=2²³`): `{2..12,14}` (no 13) is good at `j/13`; `{2..10,12,13,14}` (no 11) at `j/11`;
`{2..14}\{6}` near `j/6`. And the exact arithmetic is clean and provable:

> `t = j/k` (`gcd(j,k)=1`, `k≤13`) is good for `C` `⟺` **no `c∈C` is a multiple of `k`**
> (`‖c·j/k‖ ≥ 1/k > 1/14` unless `k|cj`, and `gcd(j,k)=1 ⟹ k|cj ⟺ k|c`).

So the clean middle-witness `t=1/k` exists **iff `C` misses a multiple of `k`** — i.e. iff `{1}∪C` is
**non-covering** at `q=k`. That is exactly **THM-523** (the non-covering witness `t=1/q`).

## The characterization: the difficulty *is* covering

This pins the AP-stability precisely. Split by the reduction:
- **Non-covering** `{1}∪C` (misses a multiple of some `q≤14`): the clean witness `t=1/q` lands in the
  middle — `G(C)` reaches it, `L>0`. **Done, and it is THM-523.**
- **Covering** `{1}∪C` (a multiple of *every* `q≤14`): **no** `t=1/q` is good — every clean modular
  witness is destroyed. `G(C)` reaches the middle only through the **resonance-shadow gaps**: at `j/k` a
  covering set replaces the *wide* bad arc of the small speed `k` (width `1/(7k)`) by the *narrower* arcs
  of its larger multiples, and a sliver of good opens on either side *if* the remaining speeds clear it.
  Whether such a gap always opens, uniformly, **is** the LRC extremal rigidity.

So the AP-stability theorem, attempted head-on, **equals the covering rigidity** — not easier, but now
exactly located: *covering is precisely the condition that removes all clean modular witnesses.* Every
non-covering set is handled by a single rational point; the covering sets are the ones with no such point,
where the middle-reaching becomes the delicate shadow-gap question. This is the same wall as S285–295,
seen once more — but this view says *why* it is a wall: the covering hypothesis is defined by the very
property (hitting every modulus) that annihilates the elementary witnesses.

## Honest state

Not a proof — the shadow-gap existence for all covering `C` is the open rigidity. What the attempt gives:
the mechanism (shadow gaps), the exact arithmetic (`t=1/k` ⟺ misses `k`), and the clean characterization
that the AP-stability difficulty **is** the covering condition, dovetailing with THM-523. (Note: LRC(14)
needs only `M≥1/14`, a point; but covering `⟹ M≥14/183>1/14` (THM-724/726) `⟹ L>0`, so `L>0 ⟺` the needed
`M>1/14` there.) If the residual falls it is as a shadow-gap / stability theorem for covering sets; the
elementary and reduction routes are exhausted.

*Files: `04-computation/lrc14_ap_stability_klein_S296.py` (+.out). HYP-6590. Consumes THM-523/366, opus-S271,
THM-527-A/663; continues
[[the-residual-is-placement-not-magnitude-LRC13-localizes-the-irreducible-core-klein-S295]].*

# "0+0=1" is the vacuum carrying the unit — Rédei is its sharpest instance, and the forbidden values are NOT in this world

**Source:** kind-pasteur-2026-06-14-S5. Dispatch: spend a long session pursuing
the "number system where 0+0=1" idea (T800), specifically whether the forbidden
values H∈{7,21} live in it.

## What a "0+0=1" number system actually is

Stripped of mystique, `0+0=1` says: **the additive-zero object (the empty / zero
configuration) is assigned the multiplicative unit `1`.** That is the universal
convention of a *generating function / partition function*: the empty
configuration has weight `1` (it's the empty product), while configurations
*combine* additively (you stack them). A "0+0=1" system is a semiring evaluation
read multiplicatively on the vacuum and additively on the bulk — the log/exp
seam. The repo lives in three such systems at once:

- **Tiling cube.** The XOR-identity (the all-base-orientation tiling, `z=0`) is the
  **transitive tournament**, whose `H = 1`. So the additive zero of the tiling
  hypercube *is* the OCF unit. Moving away from the origin (flipping tiles) grows
  `H` by creating cycles. `0 ↦ 1` literally.
- **OCF.** `H = I(Ω,2) = Σ_{independent S} 2^{|S|}`; the empty independent set
  `α_0` contributes `2^0 = 1` — the vacuum unit.
- **LRC.** The `c=0` term of the resonance-lattice Poisson sum is the independence
  density `(1−2δ)^n` — the vacuum, with all other lattice terms corrections.

So "0+0=1" is not a strange new arithmetic; it is the partition-function vacuum,
and the repo's three master functionals (`H`, the LRC density, the tiling cube) are
all evaluated in it.

## Rédei is the sharpest instance — the unit floor

For `H` the vacuum is not merely the unit, it is the *only odd contributor*:
`H = 1 + Σ_{|S|≥1} 2^{|S|} ≡ 1 (mod 2)`. That is **Rédei's theorem**, and it is the
ground floor of THM-466's 2-adic tower (`H ≡ Σ_{k<m} α_k 2^k mod 2^m`, `α_0=1`). In
2-adic language `H` is always a *unit* (`v_2(H)=0`). So:

> "H is a unit in the 0+0=1 world" = Rédei's theorem = the vacuum carries the unit.
> This is TRUE and it is exactly the unit floor — nothing more, nothing less.

## But the forbidden values are NOT in this world (tested, falsified pattern)

The dispatch's sharp question: is the unit/0+0=1 framing the *home of the forbidden
values* `H∈{7,21}`? A tempting 2-adic pattern says yes:
`7 = 1+2+4 = (111)_2`, `21 = 1+4+16 = (111)_4` — both "111" in a doubling base, so
the next forbidden value should be `73 = 1+8+64 = (111)_8`.

**It is false.** Computed (`04-computation/h_spectrum_0plus0_kps.py`): `73` is
ACHIEVABLE (at `n=7`), and `273=(111)_16`, `1057=(111)_32` are achievable at
`n=8,9`. So `7,21` being `(111)_{2^k}` is a **two-point coincidence**; the forbidden
set is not a 2-adic/unit pattern. (The exhaustive small-`n` spectrum confirms `7,21`
are the only *permanent* gaps; the larger gaps at finite `n` — `35,39` at `n=6`,
the `≥107` cluster at `n=7` — are "not-yet-achievable-at-this-`n`" and fill in as
`n` grows.)

Why `7,21` are forbidden is a **conflict-graph realizability** fact, one layer above
the unit floor. `H=7` needs the OCF-coefficient vector `α_1 + 2α_2 + ⋯ = 3`, whose
only candidate `(α_1,α_2)=(3,0)` is killed by THM-029 (three pairwise-intersecting
odd cycles force a fourth — a Helly/intersection obstruction in `Ω`), and `(1,1)` is
impossible (`α_2 ≤ C(α_1,2)=0`). `H=21` needs `α_1+2α_2=10`, blocked six ways
(HYP-1081). These are **inequalities and Helly facts about the intersection
structure of `Ω`**, not parity or 2-adic-digit facts.

## The verdict, and why it was predictable

> The "0+0=1" world is the right home for **Rédei** (the unit floor: `H` odd, the
> vacuum carries the unit) but **not** for the forbidden values. The forbidden set
> is a conflict-graph realizability obstruction (the "baby-Hodge" cone of achievable
> OCF-coefficient vectors), which the 2-adic/unit structure cannot see.

This is the *same boundary* THM-499 drew last week: `H` is strictly finer than the
spectrum, and its deep facts (the gaps, the `H∈{7,21}` exclusions) require the
conflict graph `Ω` (the `α_2` disjointness / Helly layer), not the spectral or
parity data. The 2-adic "0+0=1" world is the parity/unit side of that boundary; the
forbidden values live on the conflict-graph side. So the dispatch's question has a
clean answer consistent with the established boundary: **the unit floor is 0+0=1
(Rédei); the gaps are realizability, and the two are genuinely different layers** —
the discipline of testing the falsifiable `(111)_{2^k}` pattern (and watching it die
at 73) is what makes the separation sharp rather than rhetorical.

## What lifts

The genuinely useful crystallization: the three master functionals are
partition-function evaluations with the vacuum = unit (`0↦1`), and Rédei is the
unique place where the vacuum is the *sole* unit-residue (everything else even). The
forbidden-value program should therefore be pursued in the conflict-graph
realizability cone (THM-499's non-spectral side, the OCF-coefficient `(α_k)` vector),
**not** by chasing 2-adic patterns in `H` — a tempting blind alley this session
closed off with a single counterexample (`73`).

Cross-links: T800 (the seed), THM-466 (the 2-adic tower / unit floor), THM-499 (the
spectral/conflict-graph boundary — the same divide), THM-002 (OCF), THM-029/075
(the 7/21 mechanisms), THM-477 (the complement/blue-code side of 0+0=1),
[[the-triangular-number-is-the-n4-metagraph-kps]] (the prior session's `a=+1`, `b=/2`).

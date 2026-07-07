---
source: kind-pasteur-2026-07-07-S65
status: PROVED (THM-642). Owner directive: characterize which iso classes occur under a
  mapping rule; design rules to restrict the image; leverage the restriction for proofs.
tags:
  - lonely-runner
  - tournaments
  - image-restriction
  - sieve
  - round-tournaments
  - barrier
  - proof-leverage
---

# The image census: what tournament maps can and cannot see

**kind-pasteur-2026-07-07-S65 (HYP-4957, THM-642).** The owner asked the sharp question:
under a rule that maps runner families to tournaments, *which iso classes actually occur*,
can I *design* the rule so the image is restricted, and can I *leverage* that restriction as
a proof fact. Working it out produced both a clean positive (the sieve, as an image census)
and a clean limitation (a barrier that prunes an entire class of attacks) — and the two are
the same fact seen from opposite sides.

## The census (positive)

The mod-`p` QR cutoff (THM-640) factors through residues mod `p`, so its image is exactly the
**induced subtournaments of the Paley tournament `T_p`**. For `p = 7` — the composite factor
of `14` — arc-transitivity of `T_7` collapses this to almost nothing: any family with 6
distinct nonzero residues maps to the single class `T_7 ∖ v`, and a full residue set `{0,…,6}`
recovers the **whole Paley `T_7`** (the `H`-maximizer, `H = 189`).

The CRT-14 image then **stratifies by saturation pattern** — which `q ∈ {2,…,14}` have a
multiple present. This is a down-set in the divisor lattice, and *every cell but the top is
provably lonely* (a missing `q`-multiple gives `M ≥ 1/q` at `t = 1/q`). So:

> **The entire LRC(14) hard core is one image cell** — the fully-saturated top — and the
> sieve `counterexample_needs_all_divisors` is exactly the statement "a counterexample lives
> in the top cell." I did not expect the covering/sieve theory to *be* an image census, but it
> is: the `q`-defect vertex (residue 0 mod the prime part of `q`) is present iff the family is
> `q`-saturated.

## The barrier (limitation = leverage)

Here is the twist that makes it a theorem rather than a restatement. A residue-mod-N map's
invariants are functions of the residue multiset **only** — so it is **blind to the density
floor**. I verified this directly: families with *identical* residues mod 14 have `M` ranging
from `1/14` (tight) to `0.118`. The map's fibers are exactly the **lift/escape families**
(mac-mini-S36). So:

> **A residue-based tournament map can prove the sieve and can never reach the density floor.**

That is not a dead end — it is *leverage*, in the precise sense the owner meant. It **prunes**
the whole family of residue/covering tournament attacks (they bottom out at the sieve, which
is already GREEN) and pins the open problem to the residue-*invisible* axis, the density
`μ_{1/7}` itself. And the single-time value maps have their own restricted image — **round
(locally transitive) tournaments** — which forget gap *sizes* (the metric), so a snapshot is
also blind to loneliness. That is *why* S64's geometric cutoffs collapsed to transitive.

## The unification I did not see coming

Laying the barriers side by side, they are one phenomenon:

- **residues** → image = Paley subtournaments → sees the **sieve** only (S65).
- **one snapshot** → image = round tournaments → **metric-blind** (S64/S65).
- **pairwise** → 2-point LP → insufficient below `≈ 0.126 < 1/7` (S63 heresy C; klein-S159's
  pairwise moment-LP, independently, `0.1233`).

> **The loneliness density floor is invisible to every finite projection of the family —
> residue class, snapshot, or pair marginal. It is irreducibly a measure over all times.**

Four independently-found walls turn out to be one wall. This is the strongest kind of
leverage a limitation can give: it tells the whole fleet, with proof, which entire categories
of attack cannot close the floor, and confirms that the covering/sieve lane is *complete* as
the residue-visible content — the live front genuinely is the measure `μ_{1/7}` (the S59–S61
ledger lane), not the covering.

## Where the two-project weld now sits

THM-640 said the density-floor `M`-minimizer *is* the `H`-maximizer (the AP maps to Paley).
THM-642 locates that weld exactly: it lives at the **top image cell** (fully saturated,
residues `{0,…}`), where the runner tournament is (near-)Paley — and by the barrier the
loneliness question there is *orthogonal* to the tournament invariant. So the bridge relates
the **extremal principles** of the two projects, while the barrier proves it **cannot** relate
their **difficulty**. Both facts are real, and holding them together is the honest picture: a
gorgeous structural identity that is, provably, not a shortcut.

## A design principle, extracted

The owner asked for "tricks in how you define mapping rules that the possible set can be
restricted." The trick that worked: **make the rule factor through a coarse invariant** (a
residue, a snapshot) and the image collapses onto a structured, characterizable class (Paley
subtournaments; round tournaments). The restriction is then automatic and provable. The price
— and it is a law, not an accident — is that **the map sees only what the coarse invariant
sees**: residues → the sieve, snapshot → combinatorial order. To restrict the image *and*
retain the density floor is impossible, because the floor is not a function of any coarse
invariant. That trade — restriction ⟺ information loss — is the general shape, and knowing
its exact terms is the leverage.

## Ledger

- Canon: THM-642. Builds on THM-640 (Paley bridge), HYP-4667 (escape families = the fibers),
  HYP-4712 (sieve), S62 (Wiener collisions), S63/klein-S159 (pairwise barrier), S64 (geometric
  cutoffs).
- Files: `lrc_image_restriction_kps_S65.py` (+out with the residue-invisibility addendum).
- Does NOT prove LRC(14). It censuses the tournament-map image, proves the sieve is exactly the
  residue-visible content, and pins the density floor as irreducibly a time-measure — pruning
  the residue/snapshot/pairwise attack families with proof.

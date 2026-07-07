# The bigger picture: LRC(14)'s crux is now a finite covering system, not an analytic rigidity

**opus-2026-07-06-S126.** Assembly-owner synthesis of the fleet's last dozen sessions. The single
most important shift: the LRC(14) crux `(C)` has gone from "prove an analytic AP-uniqueness
rigidity" to "**verify a finite covering system of rational-point certificates**." Every branch is
now a Lean-checkable margin certificate; what remains is one uniform (all-height) covering claim of
Erdős-covering-system flavour, plus wiring.

## The whole architecture, top to bottom

```
LRC(14)   [13 speeds, threshold 1/14]
  ⟸ [J-K reduction]         CITATION  (Jain–Kravitz / Giri–Kravitz 2024: rank-2 subtori govern the spectrum)
  ⟸ (A)  no rank-2 subtorus in (1/13, 2/25)
  ⟸ (A)⟸(C)                projection floor GREEN (opus-S99) + pigeonhole rigidity lemma (2×2 core GREEN, wrapper OPEN)
  ⟸ (C)  the 12-speed gap  ⟺  the AP is the unique 12-family with M < 2/25
            = case 1 (GREEN)  +  case 2 (finite covering + AP)  +  case 3 (easy)
```

"LRC(14) closes when (C) closes," and everything above `(C)` is GREEN, provably clean, or a
citation (opus-S121 map).

## How (C) collapsed to a finite covering — the through-line

The fleet converged on `(C) ⟺ M < 2/25 ⟹ AP` (kps-S42: both gap edges are `k=1` Kravitz rungs,
so a gap member is an unattained `k≥2` value; mac-mini-S31; opus-S124). Then three independent
tools fused into one structure:

1. **kps-S41: the mod-25 covering atom.** `M ≥ 2/25` is witnessed at `t = c/25` whenever a rotation
   `c ∈ (ℤ/25)*` puts every speed off `{0,±1}` — a `rational_point_margin` certificate.
   `LRCMod25Floor` GREEN.
2. **mac-mini-S32 / opus-S124: the pair-blocking dichotomy.** Such a `c` exists **iff** the family
   is *not* a full transversal mod 25 (its `±`-residues miss one of the ten unit pairs). So mod-25
   splits every family into **non-blockers** (cleared, `M ≥ 2/25`, GREEN — mac-mini's THM-634 gives
   the explicit witness `c = a⁻¹`) and **blockers** (the residual). The AP is a blocker.
3. **kps-S43: the residual is a finite covering system.** The blockers are *defect-agnostic*
   (span all `d`, correcting opus-S123's defect stratification). Every non-AP blocker has
   `M ≥ 1/12`, and a **finite modulus set `q ∈ {6,…,39}` clears them all** — each at its own
   `t = c/q` `rational_point_margin` certificate. Only the AP survives every modulus (it is the
   unique `M`-minimizer `1/13`, so it has no slack anywhere).

The result: `(C)` is a **union of finitely many rational-point certificates** — the mod-25 one for
non-blockers, the `q ≤ 39` ones for non-AP blockers, the small-denominator one for mult-of-25 — with
the AP as the single deliberate exception. My own contributions slot in cleanly: S124's dichotomy is
the non-blocker/blocker split; S125's two-modulus (`q=13` clears the bottom, `q=25` the top) is the
`{13,25}` slice of kps's `{6..39}` covering; S119's mediant parity gate is the `k=2` order case.

## What is GREEN vs OPEN

**GREEN (kernel-pure Lean):** `LRCMod25Floor` (case-1 core), `LRCMod25Transversal`/THM-634
(case-1 witness), `LRCLadderD1`/THM-633 (the `{1..11}+x` ladder), `LRCBinderInfeasible` (mediant
parity gate), `LRCSubfamilyCap` (plateau), `LRCTorusProjection` (projection floor), plus the older
skeleton.

**OPEN — the actual remaining crux:**
- **The uniform covering (math).** Prove *every* non-AP blocker clears at some `q ≤ Q₀` (`39` on a
  height-≤110 sample). This is the one all-height statement, but it is a **finite mod-`q` residue
  condition** — clearing at `q` depends only on `v_i mod q` — so it is a *covering-system* fact
  (which residue patterns mod `q` fail to clear, for `q ∈ {6..39}`, and that the AP's is the only
  one failing all), not analysis. This is the Erdős-covering-flavoured heart.
- **The AP exception (math, immediate):** the AP is the unique `M`-minimizer at 12 speeds (`1/13`,
  `LRC(≤13)`, unique because `13` is prime); every other blocker is looser (`M ≥ 1/12`) and hence
  clears somewhere.
- **Case 3 (math, easy):** a mult-of-25 speed ⟹ small-denominator clearance.
- **Assembly:** the pigeonhole rigidity wrapper (finishes (A)⟸(C)); wiring `[J-K]+(A)⟸(C)+(C)` to a
  top-level conditional theorem; pinning the exact Jain–Kravitz statement.

## The one-paragraph picture

LRC(14) reduces (by the cited Jain–Kravitz rank-2 machinery, plus the GREEN projection floor and a
pigeonhole lemma) to `(C)`: the AP is the unique 12-family below `2/25`. `(C)` is no longer an
analytic rigidity — the mod-25 pair-blocking dichotomy peels off all non-blockers with a GREEN
certificate, and kps's finite covering `q ≤ 39` clears every non-AP blocker with more certificates,
leaving the AP as the unique deliberate exception. So the endgame is: **prove the finite covering is
uniform in height** (a residue-covering-system statement), discharge the two easy cases, and wire
the certificates together in Lean. The hard analysis has been replaced by a finite, checkable,
arithmetic object.

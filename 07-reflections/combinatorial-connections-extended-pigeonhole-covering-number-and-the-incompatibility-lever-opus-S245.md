---
source: opus-2026-07-11-S245
status: A combinatorial search session. NEW: the EXTENDED pigeonhole (wider-band composites, φ(q)/(2m)
  threshold) raises the provably-forced coverage 44%→48% of divisor-complete (rigorous, 0 violations); a
  clean DICHOTOMY (odd-prime-heavy → pigeonhole; even-heavy/~6-odd → clustering). CONNECTIONS (mined): the
  core is an exact COVERING-NUMBER problem (THM-718) but over-cover/unbounded/pinned ⟹ Hough/BBMST do NOT
  apply; the promising untried lever is the 26-necessary-conditions INCOMPATIBILITY (close without the
  inverse theorem).
tags:
  - lrc14
  - pigeonhole
  - covering-number
  - combinatorial
  - dichotomy
  - incompatibility-lever
  - search
---

# Combinatorial connections: the extended pigeonhole, the covering-number framing, and the incompatibility lever

**opus-2026-07-11-S245.** Owner: another session searching for combinatorial connections and simplifications.
One new proved extension, one clean dichotomy, and — via a repo-mining pass — the exact status of the
covering-system connection and the untried lever worth prioritizing.

## New: the extended pigeonhole (44% → 48%, rigorous)

The S242 pigeonhole clears a divisor-complete family at composite `q ∈ {15,21,25,27}` when
`#coprime-to-q < φ(q)/2`. It generalizes to the **wider band**: at a composite `q` (prime factors ≤ 13, danger
band `{0,±1,…,±m}`, `m = ⌈q/14⌉−1`, danger residues `{1..m}` coprime to `q`), the family clears when
`#coprime-to-q < φ(q)/(2m)` — each occupied unit fold-class blocks at most `m` unit multipliers. Adding
`q ∈ {33,35,39,49,55,65,77,91}` raises the provably-forced (no anti-concentration) coverage from **43.9% to
47.9%** of divisor-complete (verified 0 violations end-to-end). Lifts the total elementary-provable coverage
of LRC(14) a little above the ~95% of S242.

## New: the odd-prime-weight dichotomy

The extended-pigeonhole **remainder** (~52% of DC) is exactly the **odd-prime-light / even-heavy** families:
`#div-3 < 5` and `#div-5 < 4` (else forced at 27/25), so `a₃ ≈ 4`, `a₅ ≈ 2.5`, but `#even ≈ 7.4`. By klein's
*even-heavy ⟹ ~6 odd runners*, these **are** the ~6-odd-runner crux; they still clear, but by residue
**clustering** (a fold-class empty despite enough coprime speeds) at moderate `q ∈ [15,24]` — the
anti-concentration, not the pigeonhole. So the residual splits cleanly:

- **odd-prime-heavy** (`≥5 div by 3`, or `≥4 div by 5`, …) → **pigeonhole-forced** (48%, proved, no a.c.);
- **even-heavy / ~6-odd** → clears by **clustering** = the ~6-odd anti-concentration (klein's crux, 52%).

The pigeonhole handle scales exactly with odd-prime weight; the residual crux is the even-heavy tail = the
~6-odd-runner problem (S244: favorable, spread core misses `G'`).

## Connection: the core is a covering-number problem — but Hough/BBMST do NOT apply

Mining the repo pins the covering-system relationship precisely (THM-718, macmini-S24, HYP-3954):

- **It is genuinely a covering-number statement** — THM-718 (PROVED): `q` clears ⟺ the dilated-negated set
  `{±j·vᵢ : j=1..m}` **misses** a residue mod `q`; `clearing_count = (q−1) − |{±j·vᵢ mod q}|`. My auto-safe
  (S241) is the composite/coprime-fold-class refinement.
- **But it is the WRONG sub-theory for Hough / Mirsky–Newman / BBMST.** Three mismatches: (1) it is an
  **over-cover** (danger arcs cover with multiplicity ≥ 1), not an exact cover (Mirsky–Newman) nor a
  minimum-modulus question (Hough/BBMST); (2) the covering modulus is **unbounded and adaptive** — an
  explicit spread DC family blocks every `q ∈ [15,43]` and first clears at `q=44`, so no bounded window works
  (this also kills the Q50 finite census, MISTAKE-110); (3) covering-system **depth** (Guo–Sun: odd covering
  ⟹ ≥ 22 primes; Erdős–Selfridge) lives on the **shifted** side, and shifted-LRC is **false from n=5** — so
  importing covering-system depth risks proving a false statement. Net: the same *regime* (covering), the
  wrong *sub-theory*.

So the covering-system literature is a **documented dead end as a bound** — useful only as language (THM-718
makes it exact), not as a lever.

## The promising untried lever: 26-conditions incompatibility

The mining surfaced the strongest *combinatorial* (non-analytic) path, from `s554o`'s "twenty-six necessary
conditions for a counterexample": **seek two conditions that are provably incompatible** — which closes the
core *without* the inverse theorem. My extended pigeonhole adds a **new necessary condition**:

> **A counterexample is odd-prime-light** (`#div-3 < 5` ∧ `#div-5 < 4` ∧ the wider-band analogues) — else it
> is pigeonhole-forced to clear.

Combined with the standing necessary conditions — **A3**: a multiple of each maximal prime power `8,9,5,7,11,13`
(my "6 anchors", verbatim in s554o); mult of 14; multiplicand-maximal; etc. — this sharpens the counterexample
locus. The lever: show the odd-prime-light constraint (few div-by-3/5) is **incompatible** with divisor-
completeness + the energy/coverage conditions on the ~6-odd core. This is a *finite combinatorial*
incompatibility, not an inverse theorem, and it is the natural home for the extended pigeonhole.

## Net

- **Banked (proved):** the extended pigeonhole (44→48%) and the odd-prime-weight dichotomy.
- **Settled (connection):** the anti-concentration core is an exact covering-number problem (THM-718) but
  over-cover/unbounded/pinned — Hough/Mirsky–Newman/BBMST and the Q50 census are dead ends as bounds.
- **The three live levers** (mined, priority): (1) the three-gap/Steinhaus AP-inverse theorem (shared crux,
  S239); (2) push the pigeonhole into the clustering regime — *this session's extension is the first step*;
  (3) additive-energy/Sidon-dissociation on the ≤6 core (loose end is dissociated), respecting S181's
  necessary-not-sufficient caveat. **New candidate:** the 26-conditions incompatibility, with the extended
  pigeonhole's odd-prime-light condition as a fresh incompatibility partner.

→ opus-S242 (pigeonhole, the theorem this extends), opus-S241 (auto-safe), opus-S243/S244 (≤6 core, even-fold),
opus-S239 (shared crux), THM-718 (covering-number, exact), macmini-S24 / HYP-3954 (why Hough doesn't apply),
s554o (26 conditions + the incompatibility program), MISTAKE-110 (no finite census). Files:
`lrc14_extended_pigeonhole_dichotomy_opus_S245.py` (+`.out`).

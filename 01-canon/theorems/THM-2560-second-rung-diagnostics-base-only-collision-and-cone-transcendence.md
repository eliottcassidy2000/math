---
id: THM-2560
title: "Second-rung diagnostics on the canonical typed row: base-only first collision (r = 3) and cone-transcendence of the positive drift"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED (two exact
  companions; the realization branch is checked by an exact-value engine, an independent
  support-fold decision, and a second exact decomposition; the sidecar is
  gated by full 13^3-table reconstruction; ordinary/-O timing-normalized
  transcripts agree with the stored outputs).  Both results are mechanism
  NEGATIVES with
  exact redirects; neither contradicts THM-2550's positive drift.
  SCOPE: canonical TYPED row THM-2309 (25) only (not an asserted scalar
  cover); no physical current; no row removed; LRC(14) OPEN.
source: opus-2026-07-27 (executing the two finite tests specified by the
  post-MISTAKE-281 realization/transport maps)
depends_on:
  - THM-2471 (collision index (24), fibre product (31), sidecar (44)-(47))
  - THM-2368 (rooted words (9), bare tensor (29), sidecar (34)-(37), Sec 8)
  - THM-2365 (drift dichotomy)
  - THM-2550 (the positive-drift/non-replica inputs, with MISTAKE-283 scope)
related:
  - THM-2512 (Section 5's generic common-ancestry toothpick bridge test;
      only the specific u-slaved candidate tested here needs r = 5)
  - MISTAKE-281 (common-base discipline, obeyed by both companions)
  - MISTAKE-285 (cyclotomic sign-certificate representative repair)
  - MISTAKE-286 (generic THM-2512 bridge versus the specific r=5 ansatz)
scripts:
  - 04-computation/lrc14_ancestry_realization_stage_opus_20260727.py
  - 04-computation/lrc14_phase_cone_sidecar_opus_20260727.py
outputs:
  - 05-knowledge/results/lrc14_ancestry_realization_stage_opus_20260727.out
  - 05-knowledge/results/lrc14_phase_cone_sidecar_opus_20260727.out
script_sha256_ancestry: b060ada76076320330444386ef934f0ff02644b09e178708ee3f50ba601e63b1
output_sha256_ancestry: cc6634e617ea5288e32d68abebde596b608b30416cd832b73a90df57e571362f
script_sha256_phase_cone: a54a9bbc48fb2a1bb95c0ce1524461fc0d4f83a28fcef2fde42e336b32661f4f
output_sha256_phase_cone: d19d572f6c9d9ed8326d66ffc688ccedc083e47b18c875b30f077f3ba72dde0f
hash_basis: working-tree bytes (LF); replay comparison ignores elapsed-time fields only
---

# THM-2560 -- what the second rung looks like, exactly

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Row: `THM-2309 (25)`, `w = (1,14,27,40,53,66,13,2197,742586)`, owner
`j = 1`, word `sigma = {a}`.

## (A) The first collision is at depth 3 and the deep root is BASE-ONLY

The THM-2471 (24) collision index of the source-refined packet
(`e = 1_{E_1}`, `f = 1_Q P^2 1_{E_1}`, `c = c_1 = 13`) is

```text
r = 3:   I_1 = I_2 = 0 exactly,   I_3 = 9926558757352/109707098520974955,
```

the zero/positive collision branch confirmed by an independent support-fold
route and each displayed value cross-checked by a second exact decomposition;
service mass
`169 I_3 = 9926558757352/649154429118195`. The collision depth
`d = 13^3 = 2197` COINCIDES with the middle blocker scale `c_2` (a new
typed-row fact: `r = nu_13(c_2)`). The deep-root sidecar (44)-(47) at
`C = c_3 = 2*13^5` gives `delta = 2197`, `h = C/d = 338 = 2*13^2` with
`13 | h`: the deep phase is **base-only** -- independent of the
collision root `u` AND the sheet `a`. Consequently the `r = 5` affine
invariant `theta = t - 2u` does not exist here, and the companion's
specific proposal to realize one THM-2512 Sec 5 toothpick by slaving its
deep index through that invariant is NOT instantiable on this row's packet
at this clock.  THM-2512 itself neither assumes `r=5` nor rules out another
common-ancestry realization.  Redirects the files license: (i) a base-only
bridge (pair
the collision colours against a y-side cut probe); (ii) other lawful
clocks/word strata where the collision index may hit the `r = 5`
window. (Neither the `r >= 6` structural negative nor `r = 5` occurred.)

## (B) The positive drift is CONE-TRANSCENDENT

The THM-2368 (37) phase/event sidecar of the bare tensor `H_E` was
built exactly: `1,489,966` chambers (`212,845` contributing; `19`
distinct rooted snapshots; `10,998` aggregated states), GATED by exact
reconstruction of the full `13^3` table (all `2197` cells) and by the
row-sum match to the stored successor numerators. All `168` nonzero
target-mode walk endpoints are nonzero (consistent with
`D_{H_E} > 0`, THM-2550). But the Sec 8 phase-cone certificate **fails
for every one of the 168 characters**: after aggregation by positive rational
proportionality in `Z[zeta_13]`, each walk's direction classes
(`5,740--10,322` per character) span MORE than a half turn -- no
complex phase rotates any walk into a closed half-plane.

Here **cone-transcendent** is a deliberately narrow local term: for every
nonzero `a=0` target character, the set of nonzero aggregated chamber
contribution directions is not contained in any closed semicircle.  It does
**not** mean arithmetic transcendence, a nonzero winding number, chronological
full rotation, quotient-invariant cone failure, or failure of the untested
`a!=0` lines.  The canonical drift is therefore positive but not
half-plane-explainable at this exact level; its nonzero endpoints remain after
partial cancellation of non-semicircular phase support.

First cone-breaking event of the accumulated distinct-direction set in the
companion's increasing-`y` sweep starting at the owner edge `13/14` (the exact
diagnostic anchor in this fixed gauge, not a winding certificate): character
`(lam,mu) = (0,1)`, chamber left endpoint `y = 9653643/10396204`
(exact gap `25/10396204`, approximately `2.404724e-6`, from the owner-window
edge `13/14`), first direction making the accumulated set fail the phase-cone
certificate
`(-6,0,0,0,-1,-2,-2,-2,-2,-2,-1,0,0) in Z[zeta_13]`.

## Consequences for the uniform program

The THM-2368 (38) covering-row argument cannot be anchored on a plain
phase cone of the (37) sidecar: any uniform positivity mechanism must
(i) exploit the deep-coloured `a != 0` mode lines (untested here),
(ii) quotient by simultaneous-translation/factor-label symmetries
before the cone test, or (iii) break the rotating circulant word
BEFORE y-integration -- e.g. at the recorded first cancelling event on
the owner-window edge. Similarly, the ancestry/Boolean realization on
live rows must either go base-only or target row/clock/stratum
configurations with `r = 5`. Both companions obey MISTAKE-281's
common-base discipline (all pairings inside one integral / one chamber
partition; gates before interpretation).

## Independent hostile audit

The audit checked the two charts separately.

- In Part (A), `P^sA=P_(13^(s+1))f`; the adjoint/preimage reduction has the
  same normalization as the original integral.  The support-fold route fixes
  the first-positive index independently, while the correlation and `y`-side
  decompositions cross-check the exact values on complementary rungs.  With
  `r=3`, THM-2471 gives `d=c_1 13^(r-1)=13^3`, `d_0=1`, and
  `h=c_3/d=338`; both the collision-root and sheet terms are integers modulo
  one.  MISTAKE-286 records the repaired distinction between THM-2512's
  generic bridge question and the particular `r=5` slaving ansatz tested here.
- In Part (B), the root substitution, target-shift signs, common `c_3`
  safe/probe word, Jacobian, and `13*T_DEN` scaling reproduce the independently
  recomputed `2197`-cell tensor and its stored `169` row sums.  The integer
  chamber arrays are nonnegative and bounded by
  `13*T_DEN/7=553125667414320<2^63`.  MISTAKE-285 repairs the one
  genuine verifier gap: angular signs now evaluate the same conjugate
  sum/difference representative whose `l1` norm controls the error, after an
  exact outward-rounded trigonometric interval certificate.  The repaired
  run leaves every endpoint, span class, and witness unchanged.
- Noncontainment in a closed semicircle is invariant under a common complex
  phase and under the Galois relabelling of characters.  The *first* breaking
  boundary is not invariant: it belongs to the explicitly fixed standard
  embedding and increasing-`y` sweep.  No quotient-cone or winding statement
  is inferred from it.
- The audit also checked the negative boundary: `a=0` modes do not control the
  untested deep-coloured lines; a typed non-cover is not a scalar survivor;
  and neither a nonzero endpoint nor failure of this sufficient cone
  certificate supplies a physical ancestry current, future semantic root,
  Hall atom, row exclusion, or LRC(14) conclusion.

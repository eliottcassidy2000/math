---
source: claudebox-2026-06-02-S579
status: REFLECTION — geometric/sheaf reorganization of S559; one new verified
  statement (apex = unique whole-line section) + a verified lift test; residual open.
tags: [LRC, n14, 2q, apex, certificate, sheaf, line-arrangement, polynomial-method,
  P1, lift, S559, HYP-2063]
---

# The apex is the only certificate that is a whole line, not a point

**Prompt (user):** spend this math session thinking about an apex-lift certificate sheaf.

The phrase names three threads that turn out to be one object. The **apex** is
HYP-2063/S559's zero-divisor `q = n/2` where `ℤ_{2q}` stops being a field. A
**certificate** is the polynomial-method corrector (Sungkawichai–Trakulthongchai
Prop 4.1): a pair of units `(s,r)` that lands every runner in the strict interior,
witnessing loneliness. A **sheaf** is what you get when you stop asking "does a
global corrector exist?" and start asking "what are the *local* certificates, and
do they glue?"

## The reorganization

S559 already did the hard part: by CRT it reduced the `2q` corrector to a mod-`q`
condition, `∃ s',r' ∈ 𝔽_q^× : s'w_i + r'c_i ≠ f_i` for all runners `i`. Read that
line geometrically. Runner `i` is not a *constraint*; it is a **line** `L_i` in the
`(s,r)` plane over `𝔽_q` — the locus where runner `i`'s certificate *fails*. A
loneliness certificate is a single point sitting off every `L_i`. The correctors
are the complement of a **line arrangement**, and the universal-corrector question
is whether that arrangement covers the plane.

Make it a sheaf: `𝒞_v` = the constant sheaf on the certificate locus, extended by
zero. `H^0(𝒞_v)` is the set of correctors. "A corrector exists" is literally
"`𝒞_v` has a nonzero global section." For parity-matched tuples every line passes
through the origin, so each runner collapses to a single **point** of `ℙ¹(𝔽_q)` —
its forbidden slope `ρ_i` — and S559's Residual Theorem ("uncorrectable iff the
`ρ_i` cover `𝔽_q^×`") becomes "the deleted points fill the projective line."

## What the sheaf sees that the inequality didn't

Every non-apex runner deletes **one point** of `ℙ¹`. The apex `i=q` has `c_q ≡ 0`,
so its line is `s w_q = f_q` — horizontal, `r`-independent. And when the apex speed
is divisible by `q` (the tight tuple, `w_q = f_q = 0`) that "line" is `0 = 0`: it is
not a line at all, it is the **whole plane**. The apex is the unique runner whose
certificate is *non-transverse* — support of codimension 0, not 1. One section of
`𝒞_v` is the zero sheaf, and a single zero section annihilates `H^0`. That is the
whole obstruction, and it is visibly *local*: it lives at one runner, the apex.

This is the same fact HYP-2063 states as "the apex is the zero-divisor where `ℤ_{2q}`
stops being a field" — but the sheaf picture says *why that matters for certificates*:
a zero-divisor coefficient `c_q` is exactly a linear form that drops rank, i.e. a
hyperplane that degenerates to all-or-nothing. Field ⇒ every nonzero form cuts a
genuine hyperplane ⇒ every runner deletes only codimension 1 ⇒ finitely many
deletions can't fill the plane. The apex is the one place that argument has no field.

## The lift

S559's open lead (a) was to add the `r/p` freedom of the witness time
`t = s/(2q) + r/p`. In the sheaf language this is base change: work over
`𝔸²(𝔽_q) × 𝔽_p` with `p` the least prime not dividing `2q`. The apex speed mod `p`
is `q mod p ≠ 0` — a **unit** — so the degenerate horizontal line acquires a `p`-slot
and becomes a genuine hyperplane again. Verified: under the tight corrector the apex
lands at `2q·(s mod p)`, a nonzero even multiple of `2q`, strictly interior, so its
forbidden-fraction drops from `1.000` (whole space — unliftable in the original
plane) to **exactly `0.000`** at q=3,5,7,11. The lift doesn't *avoid* the apex; it
restores its transversality. The whole-line section becomes the zero section.

## The honest boundary

The lift clears the **apex**, for the **tight tuple**. It does not yet clear the
**ratio-spread residual** — the parity-matched tuples whose `ρ_i` cover `𝔽_q^×`,
which S559 identified as the real remaining hardness for `n=14`. Whether the *lifted*
arrangement ever fills `𝔸²(𝔽_q) × 𝔽_p` for those tuples is the next computation, and
I did not do it. Nor did I prove the tempting `H^1` statement — that the obstruction
is a cohomology class supported on the apex flat of the intersection poset. The apex
being the unique non-transverse flat is verified small-`q` evidence for it, not a
proof. The gain this session is organizational and one new verified local fact: the
field-failure of `ℤ_{2q}` is, certificate-by-certificate, the single runner whose
forbidden set is a whole line — and that is exactly the runner the lift repairs.

## The transcending pattern

Twice now the project has found that an *obstruction localizes to one object and one
degeneracy*: S556's wall = an LP degeneracy in the first spreading window; S559's
wall = the apex zero-divisor; here = the unique non-transverse section of a line
arrangement. Each time the move that helps is not a global re-derivation but adding
**one coordinate** that restores transversality (a second modulus, a lifted time, an
extra residue channel). The recurring lesson: LRC's hard cases are thin because they
are *non-generic*, and genericity is bought one prime at a time. Follow that — the
next obstruction will again be one degenerate stratum, and the next fix one more
coordinate.

**Artifacts:** `04-computation/lrc_apex_lift_certificate_sheaf_s579.py`,
`05-knowledge/results/lrc_apex_lift_certificate_sheaf_s579.out`; new **HYP-2101**.
Builds on HYP-2063/S559 (apex = field failure), HYP-2024 (certificates as
vertices/sections), THM-380 (pressure certificate), THM-396 (pinch/shield, the
apex `(q,q)` partner); arXiv:2604.23906 Prop 4.1/4.4.

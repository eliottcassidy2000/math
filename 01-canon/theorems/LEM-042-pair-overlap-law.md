---
id: LEM-042
title: THE PAIR-OVERLAP LAW (the full-period arithmetic core of the 7-wall pair-crumb). (A) THE EXACT FORMULA: for danger sets D_v = {t : ‖vt‖ < 1/14}, speeds a ≤ b, g = gcd(a,b), δ = g/(ab), W = (a+b)/(14ab), w_min = 1/(7b), μ(D_a ∩ D_b) = g·Σ_{j∈Z} min((W − |j|δ)₊, w_min). (B) THE EQUIDISTRIBUTION VALUE: the trapezoid's normalized integral is exactly 1/49. (C) THE EXACT CONSECUTIVE LAW: for a consecutive reduced pair (n,n+1), r = n mod 7, μ = 1/49 + r(6−r)/(49n(n+1)); hence μ ≥ 1/49, with equality iff r ∈ {0,6}. An integer AP u,u+d,… has consecutive-reduced adjacent pairs iff d divides u; equal gaps alone do not imply this. (D) FAILURES REACH THE DIAGONAL: for every q ≥ 1, the coprime pair (14q−1,14q+1) has μ = 1/49 − 6/(49(14q−1)(14q+1)) < 1/49, while its ratio tends to 1. (E) FULL-PERIOD SCOPE: the exact credits cross the seven-set Hunter wall on the full circle, but they do not furnish inherited-window credits; for (50,51), some half-circle windows have zero intersection mass.
status: PROVED ((A) arc counting; (B) integral identity; (C) exact residue-class collapse plus gcd identity; (D) specialization of the exact sum; (E) Hunter on the full circle) + REFEREED EXACT (204 formula-vs-brute pairs; consecutive values through n = 200; four full-period wall blocks; window collapse measured)
source: boxeph-2026-07-17-S69 (owner directive: finish the LRC 14 proof; integrate and extend incoming work — opus THM-956's one-item list, mac-mini HYP-3874, klein HYP-4021 hunter ledger)
depends_on: [klein's path_hunter_add_le (the ledger consuming these credits), independent of all my LEM-030..041 frame line]
related: [codex pair law (S33 usage: the D-matrix deviation — the box-hit face of the same object), THM-884(E)]
script: 04-computation/lrc14_pair_overlap_law_boxeph_S69.py -> 05-knowledge/results/lrc14_pair_overlap_law_boxeph_S69.out
---

# LEM-042 — the pair-overlap law

The 7-wall's pair-crumb, full-period face is an exact finite trapezoid sum.
Its continuous mean is `1/49`, but its discrete samples can lie on either
side of that mean.  The exact consecutive formula identifies one clean
nonnegative family; the near-diagonal counterfamily below shows why neither
"equal gaps" nor "ratios close to one" may be substituted for that arithmetic
hypothesis.

## Exact consecutive law and its scope

Let `(a,b) = g(n,n+1)` and write `n = 7h+r`, `0 <= r <= 6`.  Collapsing the
finite sum in (A) residue class by residue class gives

```text
mu(D_a cap D_b)
  = 1/49 + r(6-r)/(49 n(n+1)).
```

Thus the overlap is at least `1/49` for every consecutive-reduced pair.  It
equals `1/49` exactly in the two residue classes `n = 0,6 (mod 7)` and is
strictly larger in the other five classes.  In particular, the old wording
"equality only at n = 6,7" confused the first two examples with the two
infinite equality classes.

Equal gaps do not force consecutive reduction.  For an integer arithmetic
progression

```text
a_j = u + j d,
```

one has

```text
gcd(a_j,a_{j+1}) = gcd(u+jd,d) = gcd(u,d),
```

so the reduced adjacent difference is `d/gcd(u,d)`.  Therefore every (or,
equivalently, any) adjacent pair is consecutive after reduction if and only
if `d | u`.

## A counterfamily arbitrarily close to the diagonal

For every integer `q >= 1`, the pair

```text
(a,b) = (14q-1,14q+1)
```

is coprime, and the same exact trapezoid sum gives

```text
mu(D_a cap D_b)
  = 1/49 - 6/(49(14q-1)(14q+1))
  < 1/49.
```

Since `(14q+1)/(14q-1) -> 1`, failures are not confined to a shallow finite
region or bounded away from the diagonal.  Any usable pair-credit condition
must retain residue/arithmetic data, not merely a ratio cone.

## Inherited-window guardrail

The formula and Hunter credits above use Haar measure on the whole circle.
They cannot be inherited after restricting to an arbitrary interval `I`:
even for `(a,b)=(50,51)` there are half-circle windows with
`mu(D_a cap D_b cap I)=0`.  Consequently a block decomposition of a larger
LRC instance may use these credits only after proving a positioned-window or
aggregation lemma; it may not silently renormalize the full-period value on
each child window.

## Evidence log
- [x] formula exact ×204; exact consecutive law; AP divisibility criterion
- [x] infinite near-diagonal counterfamily `(14q-1,14q+1)`
- [x] full-period Hunter wall checked on four blocks
- [x] window collapse measured (the constraint on per-cell obligations)
- [x] named item ADVANCED (LEM-043(E), S70): the difference-comb law + beat census (mass in d beats, uniform in the gcd-7 case; single beat at d = 1; off-beat windows carry zero — 3/3 demo); the coverage synthesis assigns regimes to tools

---
id: THM-1237
title: POSITIONED SEVEN-WALL HUNTER TRANSFER — global strict-spectrum credit survives on an interval after an explicit harmonic and gcd discrepancy debt
status: PROVED (exact interval forest-Hunter theorem; strict-spectrum tree and protected-needle consumers; periodic BAD transfer; dependency-free referee; sorry-free Lean arithmetic core)
source: codex-2026-07-19-S78 continuation with path-inequality agent
depends_on: [THM-1166, THM-1221]
related: [THM-1203, THM-1218, THM-1219, THM-1226, THM-1242, HYP-7870]
script: 04-computation/lrc14_positioned_seven_wall_hunter_transfer_thm1237.py
output: 05-knowledge/results/lrc14_positioned_seven_wall_hunter_transfer_thm1237.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCPositionedSevenWallHunter.lean
script_sha256: 2995ea267cbd2a7cd222057d21297bc6651c9b59717246e1a4cfa75579bffff3
output_sha256: 6653250e834a8f27e69c8c67ac6a516033ff640a2f8bd52b42b2749641eb1dd2
formalization_sha256: 6d5eab12e2496ebe390ac5b0600684a414b1ac6b1b6131683b50f6eedcd92e93
---

# THM-1237 — positioned seven-wall Hunter transfer

## 1. Exact max-forest interval theorem

For seven distinct positive integer speeds `s1,...,s7`, put

```text
Di={t:||si t||<1/14},
Safe7=(R/Z) minus union_i Di,
H=sum_i 1/si.                                           (1)
```

Let `I` be any real interval of length `L`.  For every pair set

```text
rho_ij=mu_T(Di intersect Dj),
g_ij=gcd(si,sj),
e_ij=rho_ij(1-rho_ij)/g_ij,
w_ij(L)=L rho_ij-e_ij.                                  (2)
```

Then every forest `F` on the seven labels satisfies

```text
mu(I intersect Safe7) >= sum_(ij in F) w_ij(L)-(6/49)H. (3)
```

Consequently the strongest certificate of this form is the maximum-weight
forest:

```text
mu(I intersect Safe7)
 >= max_F sum_(ij in F) w_ij(L)-(6/49)H.                (4)
```

Negative positioned edges are omitted automatically.  In particular, if the
seven danger combs cover `I`, then

```text
max_F sum_(ij in F) w_ij(L) <=(6/49)H.                 (5)
```

Equation (5) strengthens the adaptive covered-interval certificate already
isolated in THM-1166 by identifying its exact safe-mass contrapositive and
max-forest form.

## 2. Proof

At a phase `t`, let `A(t)` be the active danger labels.  The graph induced by
`A(t)` inside a forest has at most `|A(t)|-1` edges when `A(t)` is nonempty.
Therefore, pointwise,

```text
1_(union_i Di)
 <=sum_i 1_Di-sum_(ij in F)1_(Di intersect Dj).         (6)
```

Taking complements and integrating over `I` gives

```text
mu(I intersect Safe7)
 >=L-sum_i mu(I intersect Di)
      +sum_(ij in F)mu(I intersect Di intersect Dj).    (7)
```

THM-1166's sharp density-weighted periodic discrepancy lemma supplies the
correctly signed bounds

```text
mu(I intersect Di)<=L/7+6/(49si),
mu(I intersect Di intersect Dj)>=L rho_ij-e_ij.         (8)
```

Substitution proves (3).  Since there are finitely many forests, choose a
global maximizer after the interval weights are computed; no invalid
pointwise tree switching occurs.  This proves (4)--(5).

## 3. Strict-spectrum tree and protected-needle debt

THM-1221 supplies a spanning tree `T0` with

```text
sum_(ij in T0) rho_ij >=15/154.                         (9)
```

Thus (3) gives

```text
mu(I intersect Safe7)
 >=(15/154)L-(6/49)H-sum_(ij in T0)e_ij.               (10)
```

There is a sharper universal upper endpoint for the pair mass than the
`1/7` estimate previously used in discrepancy projections:

```text
1/91<=rho_ij<=1/14,                                    (11)
```

and upper equality occurs only at reduced ratio `1:2`.  To prove the upper
bound, reduce a pair to coprime `a<b` and use THM-1166's folded formula

```text
rho(a,b)=1/49+[F(a+b)-F(b-a)]/(196ab),
F(r)=r(14-r), 0<=r<14.                                 (12)
```

The target is `F(a+b)-F(b-a)<=10ab`.  It is automatic for `ab>=5` because
the folded difference is at most 49.  The only coprime rows with `ab<5` are
`(1,2),(1,3),(1,4)`, whose differences are `20,16,12`; equality occurs only
in the first row.  Since `x(1-x)` is increasing on `[0,1/14]`,

```text
e_ij<=13/(196g_ij).                                    (13)
```

Suppose `I` is a protected needle with

```text
L>=1/(7m).                                             (14)
```

If the seven walls cover it, (10), (13), and (14) force the exact scalar debt

```text
24m H+13m sum_(ij in T0)1/g_ij >=30/11.                (15)
```

Strict reversal of (15) forces a seven-wall-safe point in `I`.  The weighted
coefficients matter: replacing `13` by the old symmetric `24` throws away the
sharp pair upper endpoint.

## 4. Localizing the deletion-incidence margin

The same proof transports not only `Safe7`, but the positive incidence
margin from THM-1221.  Let `B` be a measurable forbidden event on the same
phase circle, periodic with period `1/h`, of density `beta<=beta0<1/2`.
THM-1166's discrepancy lemma gives

```text
mu(I intersect B)<=beta0 L+beta0(1-beta0)/h.            (16)
```

Subtracting this from (4),

```text
mu(I intersect (Safe7 minus B))
 >=max_F sum_F w(L)-(6/49)H
      -beta0 L-beta0(1-beta0)/h.                       (17)
```

For a THM-1203 four-speed BAD event, `beta0=2/21`.  Using `T0`, (17) becomes

```text
mu(I intersect (Safe7 minus BAD4))
 >=L/462-(6/49)H-sum_(T0)e_ij-38/(441h).               (18)
```

Off THM-1218's arithmetic-progression branch, `beta0=60/637`, and

```text
mu(I intersect (Safe7 minus BAD4))
 >=45L/14014-(6/49)H-sum_(T0)e_ij
      -34620/(405769h).                                (19)
```

Thus HYP-7870's formerly unnamed transport loss is now fully explicit:
singleton fragmentation, pair-gcd positioning, and BAD-period erosion.  The
remaining all-packet problem is to prove that one of (18)--(19) is positive,
not to invent another global incidence margin.

## 5. Active-edge graphic rank

For `lambda>0`, let `G_lambda` contain precisely the edges with positioned
profit

```text
w_ij(L)>=lambda,
```

and let

```text
r_lambda=7-number_of_components(G_lambda)              (20)
```

be its graphic rank.  A spanning forest of `G_lambda` has `r_lambda` edges,
so (3) yields

```text
mu(I intersect Safe7)>=r_lambda lambda-(6/49)H.         (21)
```

Coverage therefore requires

```text
r_lambda lambda<=(6/49)H              for every lambda>0. (22)
```

If `G_lambda` is connected, this simplifies to `lambda<=H/49`.  This is a
rigorous active-edge blocking theorem: sufficiently many *positioned-profit*
edges cannot all be erased.  The adjective positioned is essential; global
strict-spectrum color alone is insufficient.

## 6. Guardrail and carrier/tournament audit

THM-1219 supplies the sharp guardrail.  For the consecutive seven-wall packet
`{a,...,a+6}` with `a=7m+1`, one selected `a`-gap has only `O(a^-2)` safe mass,
although the global safe mass is at least `15/154`.  For large `a` its global
strict-high graph is complete, while the six adjacent difference clocks are
all the blind `q=1` clock.  Hence no constant `theta>0` can guarantee

```text
local safe mass >=theta*(global safe mass)/a.           (23)
```

The pairwise observable for Tournament Analysis is the positioned profit
`w_ij(L)`.  Switching at `lambda` produces `G_lambda`; orienting its edges by
speed is only a gauge.  The meaningful fingerprint is its component count
and graphic rank, not directed cycles or a Redei path.  A runner tournament
retains neither `rho`, `gcd`, nor interval address.

We challenged runners, gaps, fixed sections, boundary walls, pair channels,
gcd periods, local overlap components, BAD quartets, and proof obligations as
vertices.  The faithful object is the positioned weighted forest together
with the deletion-circuit/period sidecar.  THM-1242 develops the complementary
per-edge beat-clock carrier; neither projection subsumes the other.

## 7. Verification and scope

The dependency-free referee checks the exact pair formula on `19,023` coprime
rows through 250, including the unique `1:2` upper equality, all `16,807`
labelled seven-vertex trees on all `127` nonempty active subsets (`2,134,489`
forest-Hunter checks), the threshold-rank consumer, and every rational
constant in (13)--(19).  Normal and optimized runs are byte-identical.

The Lean module kernel-checks the forest count, pair-error maximum,
protected-needle debt, periodic-forbidden transfer, BAD constants, and
graphic-rank consumer.  THM-1166's continuum discrepancy and THM-1221's
global tree selection remain explicit providers; there are no proof
placeholders.

Frozen hashes are

```text
source         2995ea267cbd2a7cd222057d21297bc6651c9b59717246e1a4cfa75579bffff3
output         6653250e834a8f27e69c8c67ac6a516033ff640a2f8bd52b42b2749641eb1dd2
formalization  6d5eab12e2496ebe390ac5b0600684a414b1ac6b1b6131683b50f6eedcd92e93
```

THM-1237 proves the requested transport inequality and exposes its exact debt.
It does not prove that every LRC(14) packet makes that debt small enough, and
it does not prove LRC(14).

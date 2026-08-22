---
id: THM-2594
title: "The realized theta-slaved contraction: S != 0 at the r = 5 window on a common ancestry base"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED (one common-base table, one
  integral per entry; exact Phi_91 reduction with 5/5 split-prime
  cross-checks; python/-O identical modulo elapsed fields; beta-zero
  and constant-column controls vanish exactly, while the genuine
  fixed-absolute-root control is nonzero).  MISTAKE-286 TYPING: this is a positive
  instance of THM-2471's deep-root-sidecar CANDIDATE realization on
  the canonical typed row -- NOT a generic THM-2512 bridge theorem,
  NOT a physical current; no scalar row is excluded; ledger 165;
  LRC(14) OPEN.  MISTAKE-281/293 discipline: all factors lie on one
  finite Boolean ancestry fibre and one base integral before every
  marginalization and DFT, but at linked nodes w_u, X_(u,a), and
  Y_(q,e'), not at one circle point.  The k=6 table is verified
  directly; no claim that k=6 exceeds THM-2449's eventual k_0.
source: opus-2026-07-28 (executing the Stage-2 step opened by
  THM-2575/THM-2577's r = 5 window)
depends_on:
  - THM-2471 (stalk, colours, deep-root sidecar (44)-(47))
  - THM-2449 ((14)-(16) response factors; (15) source-safe deletion)
  - THM-2512 ((12)-(15) contraction and factorization; (17)/(19) ignition)
  - THM-2575 / THM-2577 (the r = 5 window at sigma = {b})
related:
  - THM-2560 (the r = 3 base-only diagnostic this supersedes at {b})
  - MISTAKE-281 / MISTAKE-286 / MISTAKE-293 / MISTAKE-295 (typing repaired)
script: 04-computation/lrc14_stage2_theta_contraction_opus_20260728.py
output: 05-knowledge/results/lrc14_stage2_theta_contraction_opus_20260728.out
script_sha256: 09c43af0a0a56c7a0833bbfd13ed6a96bc5a7a3718aa1bc6b77a144bde101a06
output_sha256: bef4ee9a18ff3e2f455bad66a95252dd9989b2f60953e26e8ea0c2dc6ae7f5df
recovery_audit_script: 04-computation/lrc14_stage2_theta_contraction_r5_independent_audit_20260816.py
recovery_audit_output: 05-knowledge/results/lrc14_stage2_theta_contraction_r5_independent_audit_20260816.out
recovery_audit_script_sha256: 8be9c1b69b33ab51ac16ce2c2a7f836aae4b811e2817b90e25921c234578c568
recovery_audit_output_sha256: 087bb018408d3438494fb4ade461a4371e26bec60e31421540fea929316dd679
replay_semantic_sha256: 184c159dd4aedb02b38b2029846b9e88f301fcc63dd163682bbabddc438b825e
hash_basis: working-tree bytes (LF)
---

# THM-2594 -- the theta-slaved contraction survives

Row `THM-2309 (25)`, owner `j = 1`, word `sigma = {b}`, packet clock
`K = 2`, exact auxiliary response exponent `k = 6` (class
`k = 6 mod 16380`, `eps = +1`); collision at `d = 13^5` (`r = 5`, the
unique-invariant window).  The exponent is evaluated directly below; it is
not asserted to lie beyond THM-2449's unspecified eventual threshold `k_0`.

## Statement

Build the ONE common-base table
`N(u,q,ell,theta) = int A(y,u) F(y,q) cell_ell(y) DEEP_theta(y) QW_ell(y) dy`
(`U = P_{13^5} f`, `V = P_{13^5} e`, `f = 1_{Q_{1,{b}}} P^2 1_{E_1}`),
where the THM-2449-shaped factors are evaluated on their typed linked stalk
nodes (owner clock at `w_u`; deep probe and word at `X_(u,a)`).  They descend
to the common base because the former is root-blind and the latter two are
sheet-independent, realizing precisely `theta = t - 2u` per THM-2471 (45).
Doubly centre and form the THM-2512 (12)-(14) toothpick contraction with the
output index slaved through `theta`. Then at the fixed primitive quadruple
`(tau, a0, alpha, beta) = (1,1,1,1)`:

```text
S = Psi_{1,1}(1,1) != 0   exactly in Q(zeta_91)
```

(72-coefficient integer numerator over
`DEN_S = 2^4 3^3 5 7^3 11 13^19 53`; split primes 547/911/1093/2003/
2549 give 295, 213, 380, 1405, 883, all nonzero). By the exact
factorization `Psi = K_{alpha tau, beta} dtilde(alpha, -beta a0)`
(verified at 4 quadruples) and rational-table Galois transitivity,
**all 5,184 primitive cut coefficients of the slaved bundle fire.**

## Exact Boolean-fibre meaning

Use THM-2471's stalk coordinates

```text
w_u=(y+u)/13,
X_(u,a)=(w_u+a)/d,
Z_(u,a,b)=(X_(u,a)+b)/13^2,
Y_(q,e')=(w_q+e')/d.
```

Expanding `A(y,u)F(y,q)` as in THM-2471 (31) gives the exact formula

```text
N(u,q,ell,theta)
 =1/(13^2 d^2) integral_T
   sum_(a,e' mod d) sum_(b mod 13^2)
    1_Q(X_(u,a)) 1_E(Z_(u,a,b)) 1_E(Y_(q,e'))
    d_(1,ell)(w_u) Delta_(2u+theta)(X_(u,a))
    Q^(+1)_ell(13^6 X_(u,a)) dy.
```

At exponent `k=6`, direct identities give

```text
d_(1,ell)(w_u)=cell_ell(y),
Delta_(2u+theta)(X_(u,a))=DEEP_theta(y),
Q^(+1)_ell(13^6 X_(u,a))=QW_ell(y).
```

All three right sides are independent of the marginalized sheets.  Hence
`N(u,q,ell,theta)` is exactly the normalized sum of the corresponding
Boolean products with these three factors inserted *before* integration and
marginalization.  This proves one common ancestry fibre.  It does not put the
three factors at one circle point and does not identify `N` with THM-2449's
one-variable table `H^R`.

## What the controls do and do not localize

Two controls vanish identically: `beta = 0` by the row-zero law, and the
constant-column table obtained by copying the deep marginal into every
column.  The latter is only a deep-label-erasure control.  The genuine
fixed-absolute-root table

```text
A_abs(ell,t)=sum_(u,q)N(u,q,ell,t-2u)
```

has nonzero centred interaction and nonzero `Psi_(1,1)(1,1)`.  Therefore the
computation does **not** show that the affine slaving is the unique cause of
nonvanishing (MISTAKE-295).  What it does show is the first exact instance
where the collision root, deep phase, owner clock, word, and cut probe form a
nonzero `theta`-indexed contraction inside one integral on a live-branch row.

## Structural facts

- The slaving transports deep support with the root: per root `u`,
  only deep labels `t in {2u, 2u+1, 2u+2}` can fire
  (`DEEP_theta` empty for `theta in {3..12}`).
- `ell = 0` row vanishes (word guard-safe factor disjoint from
  `cell_0`); reflection symmetry `(ell,theta) <-> (7-ell, 2-theta)`.
- `34/169` ancestry fibres `(u,s)` carry nonzero interaction (none at
  `s = 0`).
- The genuine fixed-absolute-root control has `18/91` nonzero raw cells,
  `91/91` nonzero centred cells, and `66/72` nonzero reduced coordinates in
  `Psi_(1,1)(1,1)`; it is the sharp hostile to any causal-slaving reading.
- Gates: `J(0) = I_5 = 48602521488933856/337437093630814766589`
  reproduced from scratch; all colour/zero-sum/Hermitian laws;
  independent route-2 aggregate; toy-model pipeline validation;
  160-point slaving spot-check.

## What follows and what does not

The ancestry/Boolean realization program now has, on one row: the
owner-clock host array (THM-2575 B), the all-clock r = 5 window
(THM-2577), and a REALIZED nonzero theta-slaved contraction (this
theorem). Remaining, per MISTAKE-286/293's honest typing: relating this
candidate realization to THM-2512's generic bridge, to a target-active
`B(q)`-type object, to THM-2449's one-point lawful response table, and to
anything uniform over the 165 rows. NOT implied: a physical current, a
transplant theorem, any row exclusion, or LRC(14).

## Independent hostile audit

An independent audit rederived the `1/(13^2 d^2)` Boolean expansion from
THM-2471 (31), checked the three linked-node identities, and verified that
`k=6` is used only by direct evaluation rather than an inference about
THM-2449's `k_0`.  It separately rederived the THM-2512 factorization signs
and rational Galois argument.  Normal, optimized, and stored runs agree after
normalizing elapsed fields, and the stored hashes match the frontmatter.

As a hostile exhaustive check, the auditor directly enumerated all `5,184`
primitive quadruples of the theta table: all were nonzero (ordered reduced
vector digest
`0a1ba90d58076fea11f45690cd950092ee1da915167b99d3f473670de6d47ba0`).
The audit also found MISTAKE-295 and then independently verified the repaired
fixed-absolute-root table, including its nonzero primitive bundle.  Thus the
promoted result is the narrow auxiliary linked-node contraction stated here,
not a one-point lawful response, physical current, or row exclusion. **QED.**

## 2026-08-16 recovery audit

The recovery audit named in the frontmatter found that the parent-worktree
copies and these canonical artifacts are identical after CRLF-to-LF
normalization.  It reproduced the primary transcript in normal and optimized
modes, masking only elapsed fields.  Starting from the materialized exact
common-base table, a separately written ANOVA/cyclotomic implementation then
checked the THM-2512 factorization and directly exhausted all `5,184`
primitive quadruples: all are nonzero.  It also exhausted the genuine
fixed-absolute-root hostile, where all `5,184` primitive coefficients are
again nonzero.

The same audit rederived the linked-node identities

```text
13 w_u-y=u,
c_3 X_(u,a)-2(y+u)/13=2a,
13^6 X_(u,a)-y=u+13a,
```

and checked the deep-window law on every exact breakpoint cell.  These
integer differences prove descent to one outer base while preserving the
distinction between `w_u`, `X_(u,a)`, and `Y_(q,e')`.  The audit therefore
confirms the theorem's present narrow status and independently rejects any
claim that affine slaving is the unique cause of nonvanishing.  It does not
add a second independent constructor for the underlying interval table and
does not strengthen the theorem to a generic bridge or current.

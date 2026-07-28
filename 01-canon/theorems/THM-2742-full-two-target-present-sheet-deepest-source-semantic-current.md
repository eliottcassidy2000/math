---
id: THM-2742
title: "Full two-target present sheet and deepest-source semantic current"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the canonical
  typed LRC(14) row, THM-2720's deepest-source wall is exactly the frozen
  second-target slice t=0.  Restoring the lawful target coordinate t makes
  E3 intersect F_(ell,s,t) nonempty exactly when ell and t are nonzero,
  giving 936 positive sections.  Every one also reaches the prescribed D^6
  fork Q_(3,{1,2}), and every nontrivial deepest-target character of the
  resulting rational mass table is nonzero.  This is a target-active semantic
  current, not a private-root physical lift, endpoint current, row exclusion,
  or LRC(14) proof.
source: root/full-two-target-present-sheet-2026-07-28
audit: full-target-sheet-hostile-audit-2026-07-28 (independent periodic-antiderivative reconstruction, witness audit, character and replay checks)
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2407-owner-or-source-deletion-target-current-dichotomy
  - THM-2720-unshifted-deepest-source-present-packet-global-disjointness
related:
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
  - THM-2712-semantic-following-congruence-lock-and-address-coboundary-descent
script: 04-computation/lrc14_two_target_present_semantic_attachment_probe_20260728.py
output: 05-knowledge/results/lrc14_two_target_present_semantic_attachment_probe_20260728.out
script_sha256: d64a9e52db49c4ce4e391119d45e593453cbb3081c36ca2ae24d81bd7ab56b8f
output_sha256: d7159343e91b593d4be670cd7a53b89b5423ea077f78fc038eb7766c43939c03
hash_basis: LF-normalized bytes
---

# THM-2742 -- the full target sheet repairs the deepest-source present wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2720 proves that the old present packet cannot host an exclusive deepest
source because it is `c3`-safe while `E3` is `c3`-dangerous.  The old packet
is only one target-coordinate slice.  Restoring the discarded lawful target
coordinate moves the `q2/c3` graft together and turns the wall into a positive
two-target semantic current.

## 1. The lawful two-target sheet

Work on the canonical typed row

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5).                 (1)
```

Write `d_v(theta)` for the width-`1/14` danger comb of speed `v` centred at
`theta`, `g_v(theta)=1-d_v(theta)`, and `C_H` for the width-`1/7`
guard-safe set.  THM-2365's two target directions are

```text
eta=e_c2-e_q1,                    lambda=e_c3-e_q2.    (2)
```

For `ell in Z/7` and `s,t in F_13`, the full source-one present section is

```text
F_(ell,s,t)
 =d_c1(ell/7) * 1_(C_H) * product_(i=3)^5 g_qi(0)
  *g_q1(-s/13)g_c2(s/13)
  *g_q2(-t/13)g_c3(t/13).                             (3)
```

The last two paired shifts have exactly the signs in `(2)`.  The canonical
packet of THM-2720 is the slice

```text
F_(ell,s)=F_(ell,s,0).                                (4)
```

Thus restoring `t` changes no unrelated factor and does not relabel the
source; it restores the target coordinate suppressed in `(4)`.

Let

```text
A0=C_H intersect all five ordinary-unit safe sets,

E_j=A0 intersect d_cj(0) intersect product_(i!=j)g_ci(0).
                                                               (5)
```

Then the exact compatibility laws are

```text
E1 intersect F_(ell,s,0) nonempty  iff ell=0,
E2 intersect F_(ell,s,0) nonempty  iff ell!=0 and s!=0,
E3 intersect F_(ell,s,0) empty     for every ell,s,       (6)

E3 intersect F_(ell,s,t) nonempty iff ell!=0 and t!=0.   (7)
```

There is no restriction on `s` in `(7)`.  Hence the complete sheet has

```text
6*13*12=936                                             (8)
```

positive deepest-source sections among its `7*13^2=1183` labels.  Their
unnormalized exact grid masses range from

```text
242609986080 at (ell,s,t)=(1,9,1)

to

1362887826300 at (ell,s,t)=(3,0,2).                    (9)
```

This identifies the old wall as a one-target projection loss rather than a
uniform incompatibility of the lawful present object with `E3`.

## 2. Every repaired section reaches the prescribed fork

The deepest terminal word of THM-2305 is

```text
Q=Q_(3,{1,2})
 =A0 intersect d_c1(0) intersect d_c2(0) intersect g_c3(0).
                                                               (10)
```

For each label in `(7)`, define the exact semantic mass

```text
H_(ell,s,t)
 =measure(E3 intersect F_(ell,s,t) intersect D^(-6)Q),
D(x)={13x}.                                             (11)
```

All `936` values are positive.  Their exact range is

```text
114819491/12545122758259
 <=H_(ell,s,t)<=
672887730/12545122758259,                              (12)
```

with the extrema at the labels in `(9)`.  A strict rational witness is

```text
(ell,s,t)=(1,0,1),
x=15991693680925/100360982066072,
D^6x=3179229/20792408 in Q.                            (13)
```

Thus the repaired sheet attaches to the literal ordinary
`E3 -> D^6 -> Q_(3,{1,2})` cospan.  It is not merely a positive source
intersection or a changed affine chronology.

## 3. Every moving deepest-target character survives

Fix `ell!=0` and `s`.  The rational vector

```text
(H_(ell,s,t))_(t in F_13)                              (14)
```

is zero at `t=0` and strictly positive at all twelve nonzero `t`.  It is
therefore nonconstant.  Let `zeta_13` be primitive.  If one nontrivial
Fourier coefficient vanished, the rational polynomial

```text
P(X)=sum_(t=0)^12 H_(ell,s,t)X^t                       (15)
```

would be divisible by `Phi_13(X)=1+X+...+X^12`.  Both have degree at most
twelve, so every coefficient of `P` would be equal, contradicting `(14)`.
Galois conjugacy consequently gives

```text
sum_t H_(ell,s,t) zeta_13^(tau*t) !=0
                 for every tau in F_13^*.              (16)
```

All `6*13*12=936` marked coefficients survive.  Summing first over every
`(ell,s)` preserves the zero-at-zero/positive-off-zero profile, so all twelve
nontrivial deepest-target characters survive in the aggregate as well.

The aggregate profile satisfies `A_t=A_(-t)`, hence these coefficients are
real and occur in conjugate pairs.  Equation `(16)` proves target activity,
not an orientation of the two directions and not an identification with a
physical deck character.

## 4. Consequence and exact boundary

The proved chain is

```text
E3 semantic source and prescribed word
 + full lawful F_(ell,s,t) present sheet
 -> 936 positive semantic masses
 -> every nonzero deepest-target character.            (17)
```

The theorem does not intersect `(17)` with an affine half/`C_221` cycle, a
rail/carry/private-root/primitive-unit lift, the THM-2712 inner physical
triangle, or a separately covariant endpoint coefficient.  In particular it
does not prove a private-root edge, target/physical diagonal identification,
row exclusion, or LRC(14).  The residual ledger remains `165`.

The next object must retain `(ell,s,t)` while rebuilding those sidecars.
Projecting to `t=0` first is now proved to erase exactly the desired deepest-
source direction.

## 5. Exact reproduction and audit

Run

```bash
python3 04-computation/lrc14_two_target_present_semantic_attachment_probe_20260728.py
python3 -O 04-computation/lrc14_two_target_present_semantic_attachment_probe_20260728.py
```

Both modes byte-match the declared output.  The companion pins the audited
canonical interval implementation, reconstructs all three sources and the
terminal fork, exhausts the `1183` labels, uses the exact `D^6` prefix
identity, checks `(13)` factor by factor, and reduces every primitive target
coefficient in the power basis of `Q(zeta_13)`.  It contains no Python
`assert` nodes.

An independent hostile audit rederived the two shift signs, rebuilt every
cell through a separate periodic-antiderivative path, reproduced the
`936/1183` laws and both exact extrema, and checked the witness with strict
positive margins.  It matched all carrier and artifact hashes and replayed
the `936/936` marked and `12/12` aggregate character certificates in normal,
optimized, and stored modes.

QED.

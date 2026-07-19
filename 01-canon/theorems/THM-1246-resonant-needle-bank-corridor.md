---
id: THM-1246
title: THE RESONANT NEEDLE BANK HAS AN EXPLICIT TRANSVERSE CORRIDOR — up to twelve 14a-harmonics leave a sharp positive reciprocal-width lonely interval
status: PROVED (all-real projective-band corridor; sharp maximal needle-bank component; exact thirteen-speed family; reciprocal endpoint ladder; dependency-free exact referee; sorry-free Lean arithmetic core)
source: codex-2026-07-19-S78 continuation with address-cycle-multicell agent
depends_on: [THM-1239, THM-1243]
related: [THM-1197, THM-1240]
script: 04-computation/lrc14_resonant_needle_bank_corridor_thm1246.py
output: 05-knowledge/results/lrc14_resonant_needle_bank_corridor_thm1246.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCResonantNeedleBankCorridor.lean
script_sha256: 89412e52879369b75f694376371d06f1ba65fffd00e2446c17f2ae93a436aa27
output_sha256: a7e3a13abb8d9efeb782e5224b90e5906b9aa7e22b4ba7b661247c2d194779db
formalization_sha256: aa71f4ebd95836f7a2b022d29ddf3d7c61809c0127eab0f77e22bea797fd4d5c
---

# THM-1246 — resonant needle-bank corridor

## 1. The all-real corridor theorem

Let `a>0`, let `1<=H<=12` be an integer, and put

```text
L=15/196,
U_H=(14H+13)/(196H)=1/14+13/(196H),
R_H=182H/(14H+13).                                    (1)
```

In the scaled coordinate `x=at`, define

```text
C_H=[L,U_H],                  I_(a,H)=C_H/a.           (2)
```

Then every `t in I_(a,H)` is safe at radius `1/14` for

```text
14a,28a,...,14Ha                                                (3)
```

and for every additional real speed in the whole projective band

```text
a<=v<=R_H a.                                           (4)
```

The corridor has the exact positive length

```text
|I_(a,H)|=(13-H)/(196Ha).                              (5)
```

Every interior point is strictly safe for the needle bank.  Moreover (2) is
an exact maximal connected safe component for (3): the speed `14a` is
dangerous immediately to its left and `14Ha` is dangerous immediately to its
right.

## 2. Proof for the needle bank

For `1<=h<=H` and `x in C_H`,

```text
14hx-h in [h/14,13h/(14H)]
          subset [1/14,13/14].                         (6)
```

Thus the fractional part of `14hx` stays in the closed safe band, proving
the claim for every speed `14ha`.

The endpoints are sharp.  At the left endpoint,

```text
14L=1+1/14,                                            (7)
```

the right boundary of the address-one tooth of `14a`.  Moving left enters
that open tooth.  At the right endpoint,

```text
14H U_H=H+13/14=(H+1)-1/14,                            (8)
```

the left boundary of the address-`H+1` tooth of `14Ha`; moving right enters
that tooth.  Equations (6)--(8) prove maximality.  Subtracting the endpoints
gives (5).

## 3. The projective-band consumer

If `lambda=v/a` obeys `1<=lambda<=R_H`, then for every `x in C_H`,

```text
lambda x >=L>1/14,
lambda x <=R_H U_H=13/14.                              (9)
```

There is no wrap, so `||lambda x||>=1/14`.  This proves (4).  The theorem is
therefore not restricted to consecutive packets: any finite selection of
speeds from the real interval `[a,R_Ha]` can be added to the entire harmonic
needle bank without shrinking its safe corridor.

## 4. Exact thirteen-speed families

For an integer `a>=1`, suppose

```text
(12-H)(14H+13)<=(168H-13)a.                            (10)
```

Then the thirteen distinct integer speeds

```text
V_(a,H)={a,a+1,...,a+12-H}
         union {14a,28a,...,14Ha}                      (11)
```

all satisfy the projective bound (4).  Indeed (10) is exactly the cleared
inequality

```text
(a+12-H)/a<=R_H.                                      (12)
```

Thus the complete closed interval (2), not merely one rational phase, is
lonely for (11).  The integral condition is automatic for every `a>=2` and
every `1<=H<=12`.  For `a=1`, it holds exactly when `H>=4`.

At the balanced row `H=6`,

```text
V_(a,6)={a,...,a+6} union {14a,28a,42a,56a,70a,84a},
|I_(a,6)|=1/(168a).                                   (13)
```

This global corridor survives six simultaneous resonant needles even though
THM-1239 shows that a single `14a` needle can erase all local curvature
cracks in a selected gap.

## 5. Reciprocal toothpick self-similarity

The corridor left endpoint never moves as a new harmonic is added; only its
right endpoint descends.  For `H>=2`,

```text
U_(H-1)-U_H=13/[196H(H-1)].                            (14)
```

Thus the widths form a literal reciprocal toothpick ladder.  This is a
genuine self-similarity law, not an early-count coincidence: adding harmonic
`H` removes exactly the terminal interval in (14), and the new boundary is
the next tooth of that harmonic.  The corridor remains positive through
`H=12`; at `H=13` the two sharp boundary mechanisms meet.

## 6. Kakeya, Fano, and tournament audit

The local crack picture samples one carrier gap as a Kakeya needle and can
see a high-degree blocker star.  The corridor theorem instead retains a
transverse interval in the global `x=at` chart.  Its invariant is the nested
component with two labelled boundary owners `(14a,14Ha)`, not the number of
local cracks.  This explains recursively why resonant tooth replication does
not trap the phase circle.

On runner vertices, speed order gives the usual transitive tournament with
score histogram `(0,1,...,12)`, no directed triangles, singleton SCCs, and
one Hamiltonian path.  It loses the harmonic address, the common left
boundary, and the nested endpoint decrement (14).  We challenged runners,
gaps, needles, harmonic labels, boundary events, projective ratios, corridor
components, and proof obligations as vertices.  The faithful carrier is

```text
(C_H; left owner 14a, right owner 14Ha;
 projective band [1,R_H]; harmonic offset h).          (15)
```

## 7. Verification and scope

The dependency-free exact referee checks all `78` harmonic endpoint rows and
`119,988` thirteen-speed families for `H=1,...,12`, `a=2,...,10000`.  It
replays (5)--(14), the `a=1` threshold, strict midpoint witnesses, and the
`H=6` width.  Normal and optimized outputs are byte-identical.

The Lean module kernel-checks the rational corridor width and positivity,
both sharp bank boundaries, harmonic offset safety, the projective-band
consumer, the integer-family condition, reciprocal endpoint step, and `H=6`
width.  Interpreting the resulting fractional-part bands as circle norms is
the explicit paper layer; there are no proof placeholders or
`native_decide` calls.

Frozen hashes are

```text
source         89412e52879369b75f694376371d06f1ba65fffd00e2446c17f2ae93a436aa27
output         a7e3a13abb8d9efeb782e5224b90e5906b9aa7e22b4ba7b661247c2d194779db
formalization  aa71f4ebd95836f7a2b022d29ddf3d7c61809c0127eab0f77e22bea797fd4d5c
```

THM-1246 closes a broad resonant harmonic bank, not arbitrary blocker
address words or LRC(14).  Together with THM-1243 it shows that local
high-degree self-similarity must be charged by global address mass or an
external-blocker alphabet; blocker degree alone cannot be the missing
invariant.

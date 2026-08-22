---
id: THM-3660
title: "LRC exceptional leakage functional and fourteen-edge boundary"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The quotient of THM-3657's correction plane by its generic projective line
  has a normalized functional taking value 0 off the eight exceptional
  addresses and values +1/-1 on two reversal-paired four-sets.  Its vertical
  derivative is the signed boundary of that exceptional set: exactly 14
  edges, only 11 of them visible to ordinary base-13 carry.  The induced
  carry ledger has 1,092 nonzero ordered pairs, zero signed first moment, and
  exact squared ledger 1,560.  This is a static finite-field detector, not a
  physical positivity or LRC(14) theorem.
source: kps-s191 / THM-3659 exceptional-boundary continuation, 2026-08-21
audit: >
  PASS -- agent Poincare independently reconstructed the detector from the
  THM-3593/3657 tensor, reproduced the exact signed support, all 14/11 edges,
  seven reversal pairs, all 28,561 carry-pair multiplicities and the
  1092/0/1248/1560 ledger, both replays, and every hash.  It specifically
  audited the integer-lift and positivity wording; no characteristic-zero,
  PSD, physical-current, or scope inference was present.
depends_on:
  - THM-3657-lrc-two-current-quotient-address-atlas-and-reversal-gate
related:
  - THM-3659-lrc-mod169-carry-response-cochain-obstruction
  - THM-3534-lrc-two-current-endpoint-relative-cospan
script: 04-computation/lrc_exceptional_leakage_boundary_thm3660.py
output: 05-knowledge/results/lrc_exceptional_leakage_boundary_thm3660.out
script_sha256: 37a7c775096d80c4d8acd9d39fc0812e32545a963ee572ea773157d0ef6804c6
output_sha256: d24cda381808b5ec6009a47fc967f1293486e3ffb95d9153f4e4cd0bd3ed0730
semantic_sha256: 7811b56d46c9b6e0abdfcc86d6539b5b0169f5aaebc45b6e4431eed0b8c0b65e
hash_basis: raw LF bytes
---

# THM-3660 -- the eight exceptional addresses are an exact signed boundary

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  THM-3657's projective gate can be sharpened to a scalar detector.
The detector is exact, normalized, reversal-compatible, and supported on
precisely the eight addresses that the many-row cancellation gate requires.

## 1. The normalized quotient functional

Work in THM-3657's pinned RREF coordinate copy

```text
E=F_p^2,             p=755373809845391722745761,       (1)
```

and write

```text
b=371578917865089240854253                            (2)
```

for the slope of the generic 124-address line.  The functional

```text
ell_0(x,y)=y-bx                                       (3)
```

annihilates that line.  Its value at the exceptional address `(12,0)` is

```text
A=290125781062402394993132 !=0.                       (4)
```

Normalize

```text
ell=A^(-1) ell_0.                                     (5)
```

If `e(r)` is the correction row at address `r`, exact reconstruction gives

```text
g(r)=ell(e(r)) =
  +1, r in X_+,
  -1, r in X_-,
   0, otherwise,                                      (6)
```

where

```text
X_+={(12,0),(0,11),(6,5),(9,3)},
X_-={(0,12),(12,1),(6,7),(3,9)}.                     (7)
```

Thus `supp(g)=X=X_+ union X_-`, exactly THM-3657's exceptional locus.
This is stronger than merely saying that the eight rows occupy distinct
projective lines: after quotienting by the generic line they all have the
same absolute amplitude.

For simultaneous branch reversal

```text
j(r0,r1)=(12-r0,12-r1),                               (8)
```

one has

```text
g(jr)=-g(r).                                          (9)
```

Indeed THM-3657's row action sends `ell_0` to `-ell_0`, and reversal swaps
`X_+` with `X_-`.

## 2. The fourteen-edge vertical boundary

Let

```text
v(r0,r1)=(r0,r1+1 mod 13),
Dg(r)=g(vr)-g(r).                                     (10)
```

Then

```text
supp(Dg)={r:r in X or vr in X}.                       (11)
```

There is no accidental cancellation on an edge joining two exceptional
addresses.  The support consists of exactly the 14 labels

```text
(0,10),(0,11),(0,12),
(3,8),(3,9),
(6,4),(6,5),(6,6),(6,7),
(9,2),(9,3),
(12,0),(12,1),(12,12).                               (12)
```

With standard integer representatives, the value counts are

```text
Dg=-2 on 2 edges,
Dg=-1 on 4 edges,
Dg=+1 on 8 edges,
Dg= 0 on 155 edges.                                  (13)
```

Hence the exact combinatorial boundary ledgers are

```text
sum |Dg|=16,             sum (Dg)^2=20.               (14)
```

The affine edge reversal

```text
rho(r)=j(vr)=(12-r0,11-r1)                            (15)
```

preserves orientation after reversing both the edge and the odd potential,
so

```text
Dg(rho(r))=Dg(r).                                     (16)
```

The 14 edges form seven reversal pairs.  Also, on each of the 13 vertical
cycles,

```text
sum_(r1 in F13) Dg(r0,r1)=0.                          (17)
```

## 3. What ordinary carry sees

In THM-3659, a carry-one pair has split sum `t=(t0,t1)` with `0<=t0<=11`,
and its response is

```text
Gamma(x,y)=e(vt)-e(t).                                (18)
```

Applying `ell` gives

```text
ell(Gamma(x,y))=Dg(t).                                (19)
```

The ordinary carry chart therefore omits the three boundary edges in column
`t0=12` and sees exactly 11:

```text
(0,10),(0,11),(0,12),
(3,8),(3,9),
(6,4),(6,5),(6,6),(6,7),
(9,2),(9,3).                                         (20)
```

Their value counts are

```text
(-2:1), (-1:4), (+1:6),                              (21)
```

so their unweighted `L1` and square ledgers are `12` and `14`.

Each reachable split sum has exactly `13(12-t0)` carrying preimages.
Weighting (20) by those multiplicities gives the exact ordered-pair ledger

```text
nonzero quotient-leak pairs =1092,
signed first moment          =0,
absolute first moment        =1248,
squared moment               =1560.                  (22)
```

The zero signed moment is a hostile boundary: uniform linear averaging loses
the exceptional signal completely.  The nonzero squared ledger shows why a
lawful positive-semidefinite or variance lift would detect it, but (22)
does not itself supply such a lift.

## 4. Consequence and remaining type gap

The map `ell` is an exact exceptional-address detector:

```text
ell(e(r)) != 0  iff  r in X.                          (23)
```

Consequently, any typed current identity whose correction response can be
projected by `ell` reduces the eight-address gate to a signed scalar boundary
problem.  Moreover, any positive energy that controls `(Dg)^2` must charge an
edge incident to `X`.

What is still absent is exactly the predicate needed to use this statement
for LRC(14): no theorem identifies physical current coefficients with the
uniform carry multiplicities, makes their squared leakage positive, or
forces them to traverse the vertical response edges.  Arbitrary finite-field
or signed coefficients can cancel the first moment by (17) and (22).

## 5. Exact verification and scope

Reproduce with

```bash
python3 -B 04-computation/lrc_exceptional_leakage_boundary_thm3660.py
python3 -B -O 04-computation/lrc_exceptional_leakage_boundary_thm3660.py
```

The assertion-free companion source-pins THM-3657, reconstructs all 169
correction rows, verifies (6)--(23), both boundary supports, reversal,
all 13 telescoping identities, and the complete weighted pair ledger.
Normal and optimized streams are byte-identical.  The principal digests are

```text
detector: 974c4c8a41f96f2db1de947f21cd9ebad2161112a0e8a1e67832592ae2d6a004,
boundary: c2dac7a05d8b76c9569334848602ad728b865ad5e89de2f460a2cfd06f8c1b4a,
semantic: 7811b56d46c9b6e0abdfcc86d6539b5b0169f5aaebc45b6e4431eed0b8c0b65e. (24)
```

This is a static finite-field quotient-detector and signed-boundary theorem.
The integer ledgers in (13)--(22) use the canonical lifts of
`0,+1,-1,-2`; they do not constitute a characteristic-zero lift or a
physical positive form.  No chronology, current chain map, admissible
coefficient law, exceptional entry, row exclusion, or LRC(14) follows.
**QED.**

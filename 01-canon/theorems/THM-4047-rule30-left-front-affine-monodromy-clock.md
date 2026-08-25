---
id: THM-4047
title: "Rule 30 left-front affine monodromy and physical zero-column clock"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For the single-seed Rule 30
  orbit, every fixed left-front column has an exact affine one-period
  monodromy. A nonzero preceding tail resets the phase; an all-zero preceding
  tail makes the monodromy a translation whose bit is the parity of the prior
  column's period word. Odd parity doubles the least period, while even parity
  retains one physical phase bit. A physical all-future bank through r=100000
  has eventual-zero columns (2,7,28,399,53207,58286,87866), doubling successors
  (3,8,29,400,87867), and identity successors (53208,58287). A separate direct
  packed orbit reproduces the census and an exact global period 32 after time
  2000000. This is a fixed-column/source-depth theorem: the center is the
  moving diagonal c_t=ell_t(t), and the proved onset bounds do not uniformly
  reach that diagonal, so none of the three prizes follows.
source: codex/session-high-value-rule30-20260824, 2026-08-24
audit: >
  PASS. The primary path composes affine monodromies, queries the physical
  phase only at the seven multiplier-one columns, checks every wrapped
  recurrence, and verifies one exact packed-state repeat. A no-import audit
  evolves the closed 100001-bit left-front strip directly, with a scalar
  orientation control, and independently recovers the same zero, doubling,
  identity, and period histograms. Normal and optimized streams agree.
depends_on:
  - THM-3468-rule30-radial-green-fold-innovation-discrepancy-and-fixed-seed-carrier-boundaries
related:
  - THM-4050-rule30-half-arc-marked-cylinder-and-radius-nine-hostile
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
  - THM-3516-rule30-marked-van-der-put-carry-and-power-section-bridge
external:
  - Wolfram Rule 30 Prizes, https://rule30prize.org/ (CURRENT OFFICIAL LISTING checked 2026-08-24; problem statements and submission status only)
script: 04-computation/rule30_left_front_affine_monodromy_thm4047.py
output: 05-knowledge/results/rule30_left_front_affine_monodromy_thm4047.out
script_sha256: d454cc5b40315c02ebf486f29e227ebc9a79c78ede18fc449a7dd32f8dc21148
output_sha256: bae500127999c350ebeff77c3145fbb4abf7b5b7292d2882eeecc6492fba75e3
independent_audit_script: 04-computation/rule30_left_front_affine_monodromy_thm4047_independent_audit.py
independent_audit_output: 05-knowledge/results/rule30_left_front_affine_monodromy_thm4047_independent_audit.out
independent_audit_script_sha256: aa0cafcdf194b7073c50e02004af71a3fb91bcdb87b5cf042682e917286c303c
independent_audit_output_sha256: 4dd45b0073b15d18dfe6bc68c976fecf09fe7d86055c3bf673449b84bb521c1f
hash_basis: raw LF bytes
---

# THM-4047 -- the exact left-front monodromy clock

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** The universal part of
this theorem replaces a loose period bound by an exact reset/doubling/phase
criterion. The finite-exact part discovers the first late exceptions. The
exceptions also show why the scalar sequence of zero-column indices needs a
phase decoration.

## 1. Coordinates and recurrence

Let `a_t(j)` be the Rule 30 orbit from one seed:

```text
a_0(0)=1,  a_0(j)=0 for j!=0,
a_(t+1)(j)=a_t(j-1)+a_t(j)+a_t(j+1)+a_t(j)a_t(j+1)   (1)
```

over `F_2`. For `r>=0`, follow the left-moving front at fixed inward depth:

```text
ell_r(t)=a_t(-t+r),             ell_r(t)=0 for r<0.   (2)
```

Then `ell_r(0)=1` only for `r=0`. Substituting `j=-t+r-1` in (1) gives

```text
ell_r(t+1)=ell_(r-2)(t)+ell_(r-1)(t)
              +(1+ell_(r-1)(t))ell_r(t).             (3)
```

This is THM-3468's left-offset recurrence, now treated as a one-bit affine
system driven by the two lower columns.

## 2. Exact affine monodromy

Suppose `a_t=ell_(r-2)(t)` and `b_t=ell_(r-1)(t)` are `P`-periodic for all
`t>=S`. One step of (3) is

```text
x |-> a_t+b_t+(1+b_t)x.                              (4)
```

The composition of one driver block is an affine map

```text
Phi_(S,P)(x)=Mx+C,       M=product_(t=S)^(S+P-1)(1+b_t).  (5)
```

There are exactly three cases.

1. **Reset.** If the period word of `b` is nonzero, then `M=0`.
   Both incoming states have the same image `C`; hence `ell_r` is
   `P`-periodic from time `S+P` onward.
2. **Doubling.** If `b` is the zero word, then
   `C=sum_(t=S)^(S+P-1)a_t`. If `C=1`, every solution satisfies
   `x(t+P)=x(t)+1`. When `P` is the least period of `a`, the least period of
   `x` is exactly `2P`: a smaller period of `x` would also be a smaller period
   of its difference word `a_t=x(t+1)+x(t)`.
3. **Identity.** If `b` is zero and the weight of `a` is even, then
   `Phi` is the identity. The tail is `P`-periodic, but the lower tails do not
   determine which of the two complementary physical phases occurs. The
   missing sidecar is the single bit `x(S)`.

This proves the criterion for every `r`, with no computation. Induction from
the eventual base tails `ell_0(t)=1` for all `t` and `ell_1(t)=1` for `t>=1`
also reproves that every fixed left-front column is eventually
`2`-power-periodic. More sharply:

```text
period growth occurs only immediately behind an eventual-zero column,
and it doubles iff the preceding least-period word has odd weight.         (6)
```

An even-weight zero event is not a doubling event. It is precisely where
quotienting to period length or to the zero-column index loses the physical
phase.

## 3. Exact physical bank through depth 100000

The primary certificate recursively stores, for each `0<=r<=100000`, a
certified start, its least period word indexed by absolute time, its affine
monodromy `(M,C)`, and whether a physical phase query was needed. Reset
columns need no transient history. Multiplier-one columns query bit `r` of
the exact packed physical state

```text
L_(t+1)=L_t+(L_t<<1)+(L_t<<2)+((L_t<<1)&L_t),        (7)
```

masked to the closed strip. Every stored word is then checked around two
period wraps. The resulting eventual-zero depths are

```text
Z cap [0,100000]
  = (2,7,28,399,53207,58286,87866).                   (8)
```

Their ordinal encoding as ordinary natural numbers is

```text
n       1   2   3    4      5      6      7
z_n     2   7  28  399  53207  58286  87866.          (9)
```

The useful sequence object is decorated, not scalar:

```text
(successor, monodromy, physical phase)
(3,     doubling, 0)
(8,     doubling, 0)
(29,    doubling, 0)
(400,   doubling, 1)
(53208, identity, 0)
(58287, identity, 1)
(87867, doubling, 0).                                (10)
```

Thus the first four values do not license a triangular-number or other
low-order extrapolation. The two late identity events have opposite physical
phases even though their scalar monodromy type agrees. The decorated event is
the information-preserving discrete-to-natural map.

The bank has period histogram

```text
period       1    2    4    8     16     32
columns     16   10   56  668  87118  12133.          (11)
```

The recursive bank enters one global period-`32` state by time `1790730`.
The independent audit does not import that certificate or its starts. It
directly evolves the physical `100001`-bit closed strip and proves

```text
L_2000000=L_2000032,       L_2000000!=L_2000016.      (12)
```

Determinism makes (12) an all-future, least-global-period certificate; it is
not an observed tail fit. Extracting the 32 exact states independently gives
the same census (8)--(11).

## 4. The moving-diagonal firewall

The center is not any fixed column. Equation (2) gives the exact identity

```text
c_t=a_t(0)=ell_t(t).                                  (13)
```

So the center walks through the diagonal `r=t`. In the certified bank, only
`r=(0,1,2,3,4)` lies at or after the **stored recursive certificate start**.
Every `5<=r<=100000` lies before that sufficient start. These starts are not
claimed least: a column can enter its actual periodic tail earlier. The valid
firewall is that the fixed-column induction supplies no uniform bound
`S_r<=r`, and therefore does not control all diagonal values. This blocks
three tempting implications:

1. fixed-column eventual periodicity does not imply center periodicity;
2. the fixed-column zero/doubling clock does not control center discrepancy;
3. a precompiled fixed strip gives no lower bound for evaluating the moving
   diagonal bit.

THM-4050 gives the moving diagonal an adaptive nearest-zero address and an
exact marked-cylinder representation. Its finite radius-nine witness is a
hostile to a radius-eight shortcut, while its Haar law still does not bound
the deterministic temporal cylinder discrepancy needed here.

For THM-3471, the positive content is different: a fixed source depth uses
fixed left offsets, so (5)--(6) provide an exact procedural compiler for its
period tariff. The remaining Rule 30 tail begins at unbounded source depth
and Green slack. The clock locates where that compiler grows and where it
must retain a phase sidecar.

## 5. Interfaces with the three prizes

- **Prize 1 (non-eventual-periodicity): OPEN.** The new clock is proof-grade,
  but (13) changes column with `t`, outside the scope of the fixed-column
  limit. A proof must control the diagonal uniformly or transfer a marked
  invariant across increasing `r`.
- **Prize 2 (asymptotic balance): OPEN.** Fixed source-depth tails are now
  exactly compilable, but THM-3468's discrepancy invoice and THM-3471's
  `u>=3` completion still require an unbounded, marked tail estimate. Period
  counts alone discard the signed mass.
- **Prize 3 (fixed-seed query complexity): OPEN.** The finite left-front tail
  subsystem has a short periodic lookup after preprocessing; this is a
  hostile control, not hardness. The actual query is the moving diagonal
  (13), and no lower bound in a fixed formal machine/cost model is proved.

The official prize page remains an active listing accepting submissions; per
MISTAKE-403, this supports a dated inference that the problems remain open,
not a claim that the page literally contains the word "open."

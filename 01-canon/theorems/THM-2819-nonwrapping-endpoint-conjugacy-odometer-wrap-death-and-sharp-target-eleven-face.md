---
id: THM-2819
title: "Nonwrapping endpoint conjugacy, odometer wrap death, and sharp target eleven-face"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  Translation by tau=7/13^6
  conjugates source labels 1,...,12 at carry 12 to target labels
  0,...,11 at carry 6, including the rail, present, private-root,
  delayed-prefix, marked carrier, and carry-cell data.  It transfers the
  positive source eleven-face {2,...,12} to the positive target face
  {1,...,11}.  Canonical target label 12 is not the image of source label
  0: it pulls back to integer lift 13 and retains the odometer displacement
  91/13^6.  Since C3*(91/13^6)=14, its present C3-safe factor is invariant
  and is disjoint from the strict marked C3-danger band.  Thus the unique
  target twelve-face omitting label zero also dies.  All other target
  twelve-faces contain target label zero and die at the translated
  source-label-one future digit.  The target maximum is exactly eleven
  labels, or simplex dimension ten, on every first-fourteen marked rail.
  No outside-rail, row, or LRC(14) claim is made.
source: lrc-a12-carry-bridge/endpoint-translation-conjugacy-2026-07-28
depends_on:
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
  - THM-2809-sharp-labelwise-eleven-label-maximum-on-fourteen-marked-rails
related:
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2687-slope-seven-global-configuration-switching-positive-thirteenfold-no-go
script: 04-computation/lrc14_source_target_nonwrap_odometer_sharp_eleven_thm2819.py
output: 05-knowledge/results/lrc14_source_target_nonwrap_odometer_sharp_eleven_thm2819.out
script_sha256: 06cc68b2b962302513dc55f1612267459ba9a8a95b41ca30be4c1d7ff7592553
output_sha256: dfa30a53cfc17d00a5d6edc056952638d284614525e2c8fb1d983ba0ca528f7a
hash_basis: LF-normalized bytes
---

# THM-2819 -- endpoint translation is exact until the odometer wraps

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

THM-2809 proves that the full marked **source** endpoint on each of the first
fourteen THM-2749 rails has sharp label maximum eleven.  It deliberately
does not transfer that statement to the target endpoint.  The apparent
transfer is nearly an ordinary cyclic relabelling, but not quite: the
canonical integer-linear lift of one full `C13` turn is the nonzero odometer
translation

```text
13*(7/13^6)=91/13^6=7/13^5.                                    (1)
```

This theorem resolves both sides exactly.  The twelve nonwrapping charts
are conjugate.  The thirteenth chart retains `(1)`, but the displacement
does not repair the missing twelve-face: it stabilizes the deep `C3` phase,
so the wrapped present-safe packet is disjoint from the marked danger band.

## 1. Source and target chart conventions

Use the THM-2672 scales

```text
p=13,                  R=13^6,                  tau=7/R.          (2)
```

For an integer label `d`, source carry `12`, and a fixed THM-2640
configuration `e`, write

```text
c_s(d)=12+7d                                      in F_13,
U_d^(12)={x:x+d tau belongs to P_(e,c_s(d))}.                    (3)
```

For a canonical target label `epsilon in {0,...,12}`, put

```text
c_t(epsilon)=6+7epsilon                           in F_13,
U_epsilon^(6)={y:y+epsilon tau belongs to
                         P_(e,c_t(epsilon))}.                    (4)
```

Here `P_(e,c)` is the complete unintegrated THM-2640 packet: rail, present
factors, private deep half, predecessor carry, delayed word, future
half-digit, and exact unit flag.  The definitions use integer-linear lifts
before reduction modulo `13`.

Let `T_tau` denote literal physical translation from the marked source
coordinate to THM-2749's marked target coordinate.  For
`d=1,...,12`, set `epsilon=d-1`.  Then

```text
y=x+tau
implies
x+d tau=y+(d-1)tau,                                            (5)

c_s(d)=12+7d=6+7(d-1)=c_t(d-1)             in F_13.             (6)
```

Therefore

```text
T_tau U_d^(12)=U_(d-1)^(6),             d=1,...,12.              (7)
```

This is a physical packet identity, not only a carry or coefficient
identity.

## 2. Everything in the packet is genuinely conjugated

The exact companion checks `(7)` separately on:

1. all `14*12` translated rail supports;
2. all `7*12` present-factor banks;
3. the carry equality `(6)`;
4. the THM-2640 unit and private-root rows;
5. the future coordinate, which is fixed because `R tau=7` is integral;
6. all `14*7` complete marked source/target carriers; and
7. all `14*7` marked carry-cell restrictions.

For the last two items, THM-2749 gives literally

```text
T_tau C_(j,ell)^source=C_(j,ell)^target,
T_tau C_(j,ell,carry12)^source
       =C_(j,ell,carry6)^target.                                (8)
```

No independent target mask is inserted after translation.  This is exactly
the common-carrier convention under which THM-2749 proved equality of the
source and target delayed vectors.

## 3. The positive eleven-face transfers with its exact mass

THM-2809's positive source face is

```text
L_s={2,3,...,12}.                                                (9)
```

Applying `(7)` label by label gives

```text
L_t={1,2,...,11}.                                               (10)
```

The exact companion reconstructs `(10)` directly in the **target**
coordinate rather than inferring positivity from a coefficient equality.
On every first-fourteen rail it intersects the translated marked target,
all eleven rail and present pullbacks, and the target carry-six delayed
prefix.  Its positive clock support and exact numerator are:

```text
rail  0: (1)       399580256360672050023360
rail  1: (6)        74205644260590152069760
rail  2: (2,3)     724908063903933297548160
rail  3: (5)       565104521676801927300480
rail  4: (2,3)    1130171627188809393027840
rail  5: (4,6)     682117240653421081629120
rail  6: (2,3)    1267162127454119622485760
rail  7: (5,6)     941819893732588135224960
rail  8: (1,2,3)  1449825103908006680574720
rail  9: (5)       596479469905204957431360
rail 10: (1,3)     676409208657101856256320
rail 11: (5)       562231844838877400066880
rail 12: (2)       582228901121553500855040
rail 13: (5)       399555625773821502585600.                     (11)
```

These are exactly the source rows, as literal translation must preserve.
Both THM-2809 adjacent-edge choices remain honest for each of the eleven
labels, so all `2^11` assignments realize the same translated atom.

## 4. Twelve target faces containing zero die nonwrapping

Every target twelve-face except the one omitting target label zero contains
`epsilon=0`.  Under `(7)`, target label zero is the image of source label
`d=1`.

THM-2809's full unit/half scan gives the unique compatible source-label-one
row

```text
(sector,edge,kappa,h)=(1,0,1,12).                               (12)
```

Its future half-digit lies in

```text
(25,26)/26,                                                     (13)
```

whereas the complete marked source prefix, and hence its translated target
prefix by `(8)`, lies in

```text
(13,14)/26.                                                     (14)
```

Because `R tau=7`, endpoint translation fixes the future coordinate
`{Rx}`.  Thus `(13)` and `(14)` remain literal disjoint open intervals at
the target endpoint.  Every target face containing zero is empty after the
full delayed prefix is imposed.

## 5. The canonical label-twelve wrap is not label zero

The sole remaining twelve-face omits target label zero and contains target
labels `1,...,12`.  Labels `1,...,11` are the positive face `(10)`.
Target label `12` is where the tempting cyclic argument fails.

Pull it back to the source coordinate:

```text
y=x+tau,
y+12tau=x+13tau.                                                (15)
```

Hence

```text
T_(-tau) U_12^(6)=U_13^(12),                                   (16)
```

not `U_0^(12)`.  Although the carry has wrapped back to `12`, the physical
translation is

```text
13tau=91/R !=0                         mod 1.                    (17)
```

This is precisely the THM-2657/2672 odometer cocycle.  It cannot be silently
removed by reducing the label modulo `13`.

## 6. The odometer is killed by its C3 stabilizer

The wrapped chart does not survive.  The relevant deep speed is

```text
C3=742586=2*13^5.                                               (18)
```

Therefore

```text
C3*(91/R)=14                         in Z.                       (19)
```

The residual translation `(17)` fixes the `C3` phase exactly.  Every
`F_(ell,7)` present packet retains the `C3`-safe factor, while the marked
source deep overlap is the strict danger band

```text
H_mark=(169,181)/182.                                           (20)
```

Consequently

```text
T_(-91/R) F_(ell,7) intersect H_mark=empty,
ell=0,...,6.                                                    (21)
```

The sign in `(21)` follows the preimage convention in `(3)`; invariance in
`(19)` makes the conclusion sign-independent.

This is the first universal failed factor, not an inference from the final
zero.  Before applying `(20)`, the wrapped present packet still meets both
present-complement factors with exact positive masses

```text
ell : after source-safe, after translated-target-safe
 0  : 784513810080, 728825300280
 1  : 356159643840, 330902579700
 2  : 473353040160, 439916713920
 3  : 625685860800, 581278296060
 4  : 620815595400, 578179396872
 5  : 438008061491, 407334158691
 6  : 394663389120, 366653853720.                              (22)
```

Intersecting the next `C3` factor gives zero in all seven rows.  No rail,
delayed prefix, semantic section, or other label is needed for the
obstruction.  Direct target-coordinate reconstruction of labels
`1,...,12` is zero on all fourteen rails and all seven clocks.

## 7. Sharp target classification

Let `L` be any twelve-element target-label subset, with arbitrary compatible
THM-2640 configurations selected label by label.

- If `0 in L`, Section 4 kills the face at target label zero's future digit.
- If `0 notin L`, then `L={1,...,12}`, and Sections 5--6 kill target label
  twelve at the invariant `C3` safe/danger gate.

Thus all thirteen target twelve-faces are empty.  Section 3 exhibits the
positive eleven-face `{1,...,11}`.  Therefore, on every one of the first
fourteen fully marked THM-2749 target rails,

```text
maximum target-label cardinality=11,
maximum target simplex dimension=10.                            (23)
```

The exact failure census is

```text
12 faces: translated label-zero future digit;
 1 face:  odometer-wrapped label-twelve C3 gate.                 (24)
```

## 8. Information ledger and stopping boundary

| item | exact content |
|---|---|
| source | THM-2809 carry-12 marked source label family |
| target | THM-2749 carry-6 marked target label family |
| map | literal endpoint translation `x -> x+7/R` |
| preserved before wrap | rail, present, carry, private root, future digit, marked section, exact mass |
| destroyed by canonical reduction | integer lift of the target label |
| missing sidecar | odometer displacement `91/R` |
| wrap action on `C3` | trivial, because `C3*91/R=14` |
| wrap consequence | present-safe factor remains opposite the marked danger band |
| positive survivor | target labels `{1,...,11}` with all `2^11` edge assignments |

The result is a useful correction to the naive statement that target labels
are merely source labels cyclically shifted.  They are so only on the
nonwrapping interval of integer lifts.  The wrap retains a physical cocycle;
here its stabilizer is strong enough to kill the last twelve-face.
There is also no internal edge-cube substitute for that cocycle: THM-2809's
`2^11` edge assignments are all the same physical atom, so every internal
Boolean edge derivative vanishes.  The external translation jet `(17)` and
its retained `C3` phase are the information that the edge cube forgets.

The theorem does **not** prove a statement outside the first fourteen
THM-2749 rails, at a different marked band or carry, for a one-sided target
bank not transported by `(8)`, for a complete row, or for LRC(14).

## 9. Exact companion

Run

```bash
python 04-computation/lrc14_source_target_nonwrap_odometer_sharp_eleven_thm2819.py
python -O 04-computation/lrc14_source_target_nonwrap_odometer_sharp_eleven_thm2819.py
```

Both modes must byte-match the stored transcript.  The companion pins the
THM-2672, THM-2749, and THM-2809 scripts; checks every nonwrapping rail and
present chart identity; verifies all marked source/target carrier and carry
translations; reconstructs the target-label-zero delayed-prefix death;
tests the integer-lift-13 wrap and the exact `C3` stabilizer; records the
positive pre-deep masses `(22)` and zero deep masses; and directly rebuilds
both the positive target eleven-face and the empty target twelve-face
omitting zero.  It uses exact integer arithmetic, explicit exception gates,
and no truth-bearing Python assertions or floating point.

Promotion requires an immutable independent hostile audit of the
source/target coordinate convention, shift signs, canonical
integer lifts, carry relabelling, target-zero prefix gate, `C3` factor order,
strict endpoints, all fourteen direct target rows, all `2^11` edge choices,
normal/optimized replay, hashes, and documentation gates.

QED, pending independent hostile audit.

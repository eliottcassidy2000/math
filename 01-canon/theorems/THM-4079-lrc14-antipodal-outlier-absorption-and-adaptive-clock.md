---
id: THM-4079
title: "LRC(14) antipodal outlier absorption and adaptive clock"
status: >
  PROVED RELATIVE TO THM-2061/2066/2072 + VERIFIED-EXACT + INDEPENDENTLY
  VERIFIED-EXACT. A strict antipodal safe margin for a finite core absorbs
  every sufficiently large one-speed outlier. The same margin gives an
  adaptive owner clock 4B, compressed to 2B for even B, with empty eligible
  set. In particular every C_B={1,...,10,B}, B a positive multiple of 12012,
  is a primitive divisor-complete eleven-core whose dyadic two-tail seam is
  closed at clock 2B, although infinitely many such cores make any prescribed
  finite clock bank completely blind. A symbolic inclusion-minimal nine-speed
  obstruction proves that the antipodal premise is not universal. LRC(14)
  remains open.
source: codex-frontier-synthesis-creative-20260825d / LRC structural beyond-30 lane
audit: >
  PASS. The primary Fraction-exact construction tests 138 strict rational
  antipodal bodies, 2,898 cases in each of the continuous, 4B, and even-2B
  threshold regimes, 568,512 safe-phase inequalities, and 5,050 exact band-cover
  mesh points. The independent path checks 22,000 canonical-family phases,
  5,045,040 literal odd tail classes, 5,050 fixed-bank divisibility gates,
  the five symbolic obstruction intervals, 292 equality walls plus two
  endpoints and all 293 open cells, and 144 deletion-witness inequalities.
  Normal and optimized outputs byte-match; both scripts have zero assert
  nodes and zero floating literals.
depends_on:
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-2066-dyadic-seam-owner-word-crt-atlas
  - THM-2072-fixed-owner-clock-bank-no-go-and-half-shift-certificate
related:
  - THM-4075-lrc14-divisor-complete-dyadic-owner-word-closure-through-30
script: 04-computation/lrc14_antipodal_outlier_adaptive_clock_thm4079.py
output: 05-knowledge/results/lrc14_antipodal_outlier_adaptive_clock_thm4079.out
script_sha256: dff69f67bb0534cb3b60265aefab60282935a8f1477783668f350616116b5acd
output_sha256: 6dac3c17b85f428ef79b4569c627efeeef75a8e51341859bd887e0171abb90a7
independent_audit_script: 04-computation/lrc14_antipodal_outlier_adaptive_clock_thm4079_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_antipodal_outlier_adaptive_clock_thm4079_independent_audit.out
independent_audit_script_sha256: f86a5430d313262b46ce48181f7c476ffaf16154f40f2452d38a71522a60914b
independent_audit_output_sha256: 7099e4f5c31a9b6489ce9493ad3c44ae0b8c0e75b26e3a14a946a499ab24f0ce
hash_basis: raw LF bytes
---

# THM-4079 -- strict antipodal margin absorbs a large outlier

**PROVED RELATIVE TO THM-2061/2066/2072 + VERIFIED-EXACT + INDEPENDENTLY
VERIFIED-EXACT.** THM-4075 closes the divisor-complete seam only through
maximum 30, while THM-2072 proves that no core-independent finite clock bank
can be uniform. This theorem takes the complementary route: a strict
continuous margin absorbs an unbounded outlier and manufactures a clock from
that outlier itself.

## 1. Inherited seam and owner coordinates

Put `delta=1/14`. For a finite nonempty set `C` of positive integers, write

```text
G_C={theta in R/Z: ||c theta||>=delta for every c in C}.       (1)
```

For a positive clock `N`, retain THM-2066's labelled packet and tail relation

```text
A_N(C)={0<=r<N:14|cr|_N>=N for every c in C},                 (2)

E_N(C)={odd z mod 2N:7|zr|_N<N for every r in A_N(C)},       (3)

R_N(C)={(u,v) in E_N(C)^2:
          omega_v(r)=1-omega_u(r) for every r in A_N(C)}.     (4)
```

Here `omega_z(r)=nint(zr/N) mod 2`. THM-2061/2066 prove that

```text
R_N(C)=empty                                                 (5)
```

excludes a strict LRC(14) counterexample of the dyadic-seam form

```text
2C union {x,y},        x,y distinct positive odd.             (6)
```

THM-2072 gives the parallel continuous certificate: `(6)` is impossible if
`G_C` contains both `theta` and `theta+1/2`.

## 2. Adaptive antipodal outlier absorption

Let `D` be finite and nonempty, put `R=max D`, and suppose some `t_0 in R/Z`
and `eta>0` satisfy the strict antipodal margin

```text
||d(t_0+j/2)||>=delta+eta       for d in D and j=0,1.          (7)
```

Let `B` be a positive integer not in `D`. If

```text
boxed: B>=R/(14eta),                                          (8)
```

then there is `t` at circular distance at most `1/(14B)` from `t_0` such
that

```text
t,t+1/2 in G_(D union {B}).                                   (9)
```

To prove this, put `x_0=Bt_0 mod 1`. The phases simultaneously safe for the
new speed at both half-separated times are

```text
K_even=[delta,1-delta],
K_odd =[delta,3/7] union [4/7,1-delta],                       (10)
```

according to the parity of `B`. Every complementary component is an open
circular arc of length `2delta=1/7`. Choose `y` in the appropriate band and
`Delta in [-delta,delta]` with `y=x_0+Delta mod 1`, then put

```text
t=t_0+Delta/B.                                               (11)
```

The new speed is safe by `(10)`. Since distance to the nearest integer is
one-Lipschitz, for every old speed and both half-shifts,

```text
||d(t+j/2)||
 >=delta+eta-d/(14B)
 >=delta+eta-R/(14B)
 >=delta.                                                     (12)
```

This proves `(9)`, including equality in `(8)` because the target safe sets
are closed. THM-2072 then closes every seam `(6)` over the enlarged core.

## 3. The margin manufactures its own owner clock

The same hypothesis has an exact rational strengthening. If

```text
boxed: B>=R/(4eta),                                           (13)
```

choose an odd integer `r` nearest `4Bt_0`, so `|r-4Bt_0|<=1`. At clock
`N=4B`, the two residues

```text
r, r+2B mod 4B                                               (14)
```

belong to `A_(4B)(D union {B})`, and

```text
boxed: E_(4B)=R_(4B)=empty.                                  (15)
```

Indeed `t=r/(4B)` moves by at most `1/(4B)`, so `(13)` preserves every old
inequality. The new phase is `Bt=r/4`, a quarter or three-quarter; a
half-shift either fixes or interchanges those values. Hence both labels in
`(14)` are safe.

For any odd tail class `z`, the normalized distances at those labels are

```text
||zt||,       ||z(t+1/2)||=1/2-||zt||.                       (16)
```

They cannot both be strictly below `1/7`. Eligibility in `(3)` would require
exactly those two strict inequalities, so no odd class is eligible. This
proves `(15)` before owner words are even needed.

There is a sharper even-outlier grid. If `B` is even and

```text
boxed: B>=R/(2eta),                                           (17)
```

choose odd `r` nearest `2Bt_0`. At clock `N=2B`, the residues `r,r+B` are in
the safe packet and again

```text
E_(2B)=R_(2B)=empty.                                          (18)
```

Now the time drift is at most `1/(2B)`, while `Bt=r/2` is a half-integer and
the half-shift adds the integer `B/2`. Evenness is load-bearing: for odd `B`,
the same half-shift would turn that half-integer into an integer and make the
new speed unsafe.

## 4. An explicit all-height family beyond every fixed bank

For

```text
D={1,2,...,10},             t_0=1/12,                         (19)
```

the minimum of the twenty distances in `(7)` is `1/12`. Thus

```text
R=10,       eta=1/12-1/14=1/84,                              (20)
```

and the three sufficient thresholds `(8)`, `(13)`, `(17)` are exactly

```text
60, 210, 420.                                                (21)
```

In particular every `B>=60` is continuously absorbed, every `B>=210` has a
`4B` owner certificate, and every even `B>=420` has a `2B` certificate.

Now let `B` be any positive multiple of

```text
12012=lcm(11,12,13,14),
C_B={1,2,...,10,B}.                                          (22)
```

This eleven-core is primitive and divisor-complete through 14. At clock
`N=2B`, take

```text
r=B/6+1,
t=r/(2B)=1/12+1/(2B).                                       (23)
```

The integer `B/6=2002(B/12012)` is even, so `r` is odd. Each old speed loses
at most `5/B` from the margin `1/84`, while

```text
Bt=B/12+1/2,                                                 (24)
```

and both `B/12` and `B/2` are integers. Therefore

```text
r,r+B in A_(2B)(C_B),       E_(2B)(C_B)=R_(2B)(C_B)=empty.   (25)
```

Every dyadic two-tail seam over every core `(22)` is closed.

This adaptive success coexists with arbitrarily severe sensor blindness.
Given any finite nonempty clock bank `F`, put

```text
L=lcm(12012,{N:N in F}),       B=mL,       m>=1.              (26)
```

For each `N in F`, the speed `B` is zero modulo `N`, so

```text
A_N(C_B)=empty,
E_N(C_B)=O_N,                R_N(C_B)=O_N^2.                  (27)
```

The actual odd pair `(1,3)` is compatible across the whole generalized-CRT
bank, so the bank relation remains nonempty. Nevertheless the adaptive clock
`2B` closes the core by `(25)`. There are infinitely many such `B` for every
fixed bank. This is certificate blindness, not geometric hardness and not an
LRC counterexample.

## 5. Exact failure boundary of the antipodal method

The antipodal premise is not universal. Define

```text
D_0={1,3,5,7,8,9,13,21,24}.                                (28)
```

This set has no antipodal safe pair, and is inclusion-minimal with that
property.

For the symbolic proof put `t=2theta` and

```text
O={1,3,5,7,9,13,21},
S_o=union_(k=0)^(o-1)[(7k+1)/(7o),(7k+6)/(7o)].              (29)
```

For odd `o`, simultaneous safety at `theta,theta+1/2` is exactly
`||ot||>=1/7`, or `t in S_o`. Exact interval intersection gives

```text
intersection_(o in O) S_o
 =[15/91,6/35] union [12/49,13/49]
  union [71/147,76/147]
  union [36/49,37/49] union [29/35,76/91].                   (30)
```

On the first interval, `12t-2` lies in `[-2/91,2/35]`, so speed 24 fails
strictly; symmetry handles the fifth. On the second, `4t-1` lies in
`[-1/49,3/49]`, so speed 8 fails strictly; symmetry handles the fourth. On
the central interval, `4t-2` lies in `[-10/147,10/147]`, again below `1/14`.
Thus speeds 24 and 8 cover every odd-core possibility.

Every deletion is feasible, with the following exact witness `theta`:

| omitted speed | `1` | `3` | `5` | `7` | `8` | `9` | `13` | `21` | `24` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| witness | `4/147` | `11/63` | `29/336` | `15/112` | `43/336` | `11/49` | `29/126` | `43/182` | `15/182` |

Direct substitution gives all 144 remaining speed/phase inequalities. This
proves inclusion-minimality, not minimum cardinality or uniqueness.

Positive margin is also necessary for the perturbative theorem. The weak-safe
pair `D={1,13}`, `t_0=1/14` has no locally safe direction: for
`0<epsilon<1/182`, moving left fails speed 1 and moving right fails speed 13.

## 6. Audits, connection ledger, and scope

The primary audit generates 138 reduced rational anchors through denominator
30, takes every speed at most 40 that is strictly antipodal-safe there, and
tests 21 consecutive outliers at or above each threshold, with `B>40`. It
checks a 5,050-point exact mesh, including extremizers for both parity bands
in `(10)` that realize the symbolically proved covering radius `1/14`.

The second audit starts from the canonical family `(22)`, checks every odd
tail class at the two labels in `(25)` for the first 20 multiples, verifies a
thousand core phases and one hundred growing blind banks, reconstructs the
five intervals `(30)`, and separately sweeps the full `theta` line. That sweep
has 292 genuine equality walls plus endpoints `0,1`, hence 294 checkpoints
and 293 open cells; none is feasible.

Replay from the repository root:

```bash
python3 -B 04-computation/lrc14_antipodal_outlier_adaptive_clock_thm4079.py
python3 -B -O 04-computation/lrc14_antipodal_outlier_adaptive_clock_thm4079.py
python3 -B 04-computation/lrc14_antipodal_outlier_adaptive_clock_thm4079_independent_audit.py
python3 -B -O 04-computation/lrc14_antipodal_outlier_adaptive_clock_thm4079_independent_audit.py
```

The source is THM-2072's antipodal safe pair and fixed-bank hostile. The map
is a parity-appropriate phase projection followed by a core-dependent odd
grid. It preserves the strict old-core margin, the two antipodal labels, and
the owner eligibility predicate; the fixed bank loses the outlier because it
sees only zero residues. The result closes one unbounded structural family of
THM-2061 dyadic seams. It does not handle arbitrary simultaneous outliers,
prove universal antipodal feasibility, or prove LRC(14); nor does it claim
that `2B` or `4B` is the least certificate clock.

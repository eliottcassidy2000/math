---
id: THM-2831
title: "Cyclic resonance simple-zero and full-lower-bank exclusion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In every
  normalized R=3k+2 Faber bank with d=s=0, a common simple zero of V and
  M is incompatible with A_src K_Q=lambda M, even when all retained
  lower Faber rows are present.  Consequently the cyclic THM-2796
  carrier cannot realize any THM-2827 resonance in this chart.  Nonzero
  d or s and other chart-entry mechanisms remain open.
source: root/cyclic-resonance-simple-zero-2026-07-28
audit: >
  uniform-all-pole-audit independently rederived the zero-source Faber
  specialization, unique top/least active terms, both simple-zero
  valuation regimes, cyclic and UFD typing, and sharp c=0/v=0
  boundaries; normal, optimized, stored, LF hashes, AST, and docs pass.
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2827-uniform-pole-order-faber-nonresonance-atlas-and-double-pole-exclusion
related:
  - MISTAKE-317
script: 04-computation/jc_cyclic_resonance_simple_zero_thm2831.py
output: 05-knowledge/results/jc_cyclic_resonance_simple_zero_thm2831.out
script_sha256: 3826afcd985354a459b6cab2dac341dd60450461143823ae53357904510cc683
output_sha256: c0d88fed97f9fae85c752f27d024977129779b348f1a196997c29c9bf361cf72
hash_basis: LF-normalized bytes
---

# THM-2831 -- cyclic resonance simple-zero exclusion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2827 leaves a formal local resonance whenever `R=2 mod 3`.  Its
first witness has `(R,ord(V),ord(A_src))=(8,13,6)`.  The witness satisfies
every local valuation and flux equation, so another one-point valuation
comparison cannot remove it.  The missing global datum is that the
cyclic response carrier has a polynomial which is simultaneously a
simple factor of `V` and `M`.

This theorem uses that common factor.  In fact the proof does not need
the cyclic formula or the resonance divisibility: a common simple zero
already excludes the entire `d=s=0` Faber bank.

## 1. Statement

Fix

```text
R=3k+2>=8,                         k>=2,               (1)
```

and the normalized reduced Faber bank of THM-2827,

```text
Q=E_(4R-2)+sum_(j=1)^(R-2) a_j E_(4j-2).             (2)
```

Work in its nonsplit polynomial exact-square-prefix chart, with the
notation

```text
T=A_src^2/V,                       A_src K_Q=lambda M,
lambda!=0.                                               (3)
```

Assume

```text
d=s=0.                                                   (4)
```

If a nonconstant squarefree polynomial `S` has a root `alpha` satisfying

```text
ord_alpha(V)=ord_alpha(M)=1,                             (5)
```

then equations `(2)--(5)` are impossible.

The THM-2796 cyclic chamber has

```text
S=y^N+c,      c!=0,
V=v S y^(N+2),                    M=S y,                (6)
```

so `(5)` holds at every root of `S`.  Hence no cyclic carrier enters the
chart `(4)`.  In particular, none of the residual THM-2827 resonance
degrees

```text
R=3k+2,       D=4k+3,             N=D delta             (7)
```

can be realized there, even after arbitrary lower rows in `(2)` are
retained.

## 2. The exact odd Faber bank at `d=s=0`

THM-2827 gives

```text
K_j=-(d/2)phi_j+T H_j,                                  (8)

H_j =
 4 sum_(h,u,ell>=0; ell+2u=j-2-3h)
   binom(j-1/2,j+1+h) binom(j+1+h,ell)(-1)^ell
   binom(-2-2h,u) T^(2h)(Td)^u s^ell.                  (9)
```

Under `(4)`, the only surviving summands have `u=ell=0`.  Thus

```text
K_j=0                                      if j!=2 mod 3,

K_(3h+2)=c_h T^(2h+1),                    c_h!=0.       (10)
```

The coefficient is explicitly

```text
c_h=4 binom(3h+3/2,4h+3),                               (11)
```

which is nonzero.  Because the top row `R=3k+2` is present and every
retained lower row has index at most `R-2`,

```text
K_Q=T P(T^2),                                           (12)
```

where `P` is nonzero, its highest exponent is uniquely `2k+1`, and every
other active exponent is a smaller positive odd integer.

## 3. A simple zero cannot have response order one

Let

```text
a=ord_alpha(A_src)>=0.                                  (13)
```

The order is finite: `A_src=0` would make the right side of `(3)`
nonzero and the left side zero.  From `(3)` and `(5)`,

```text
ord_alpha(T)=2a-1.                                      (14)
```

There are two exhaustive cases.

If `a=0`, then `T` has order `-1`.  The unique largest exponent in
`(12)` is therefore the unique least-valuation term.  Put

```text
e=2k+1.
```

Then

```text
ord_alpha(A_src K_Q)=-e<0.                              (15)
```

If `a>=1`, let `ell>=1` be the least active odd exponent of `(12)`.
It is now the unique least-valuation term, and

```text
ord_alpha(A_src K_Q)
 =a+ell(2a-1)>=2.                                       (16)
```

But `(3)` and `(5)` require

```text
ord_alpha(lambda M)=1.                                  (17)
```

Neither `(15)` nor `(16)` equals `(17)`.  This proves the theorem.

The mechanism is not top-term domination in every valuation direction:
the controlling term switches from the highest odd exponent at `a=0`
to the lowest active odd exponent at `a>=1`.  What matters is the
missing order `1` between the two regimes.

## 4. Independent top-row UFD model

The top-only specialization exposes the resonance arithmetic directly.
Put

```text
D=4k+3,             e=(D-1)/2=2k+1.                    (18)
```

Then `(10)` gives `K_Q=c_R T^e`, so `(3)` becomes

```text
c_R A_src^D=lambda M V^e,             c_R!=0.          (19)
```

For the cyclic carrier `(6)`,

```text
A_src^D
 =unit * S^(e+1) y^[1+e(N+2)].                         (20)
```

Since `c!=0`, the polynomial `S=y^N+c` is squarefree and coprime to
`y`.  Every irreducible factor of `S` has valuation

```text
e+1=(D+1)/2,                                            (21)
```

strictly between `0` and `D`.  Equation `(20)` is therefore impossible
in the UFD `C[y]`.  This is an independent proof for the top row and
explains why the primitive Farey resonance of THM-2827 still cannot
globalize in the cyclic zero-source chart.

The assumption `c!=0` is load-bearing and sharp for this UFD argument.
If one formally puts `c=0`, then the two displayed `y` factors merge and
their total exponent is

```text
(e+1)N+1+e(N+2)=D(N+1),                                (22)
```

an exact multiple of `D`.  But `c=0` is not a THM-2796 cyclic carrier:
the response identity has `C=-Nc!=0`, and `S` would also cease to be
squarefree and disjoint from the pole factor.  Likewise `v=0` destroys
the nonsplit potential and makes `T=A_src^2/V` undefined.

## 5. Consequence and remaining chart boundary

THM-2827 reduced every balanced passport in residue class `R=3k+2` to
parts divisible by `D=4k+3`.  THM-2831 now removes the complete cyclic
`e=0` carrier from the normalized `d=s=0` source chart, including the
first `R=8,N=11` case.

This does **not** prove that the cyclic response map fails to enter every
Keller chart.  A nonzero `d` or `s` makes `K_Q` a polynomial in
`T,s,Td`; distinct terms can then have order `1` or collide.  Arbitrary
nonpolynomial prefixes, other source charts, and the full planar
Jacobian conjecture remain open.

## 6. Exact companion

The exact companion:

1. reconstructs every active coefficient `(11)` for `R=8,...,32` and
   checks the unique odd-exponent ladder;
2. rechecks the first hostile coefficient
   `c_R=-195/131072` at `R=8`;
3. exhausts the two simple-zero valuation regimes over a finite audit
   grid without using the grid as proof; and
4. verifies the UFD exponent `(21)` and formal boundary `(22)` for
   several resonance scales.

It uses only exact integer and rational arithmetic, contains no Python
`assert`, and has identical normal, optimized, and stored output.  Run

```text
python 04-computation/jc_cyclic_resonance_simple_zero_thm2831.py
python -O 04-computation/jc_cyclic_resonance_simple_zero_thm2831.py
```

The universal proof is Sections 2--4, not the finite replay.

## 7. Independent hostile audit

An independent audit rederived `(10)` directly from THM-2827, checked
that arbitrary lower coefficients cannot cancel either controlling
term, and separately verified the two simple-zero regimes,
`A_src=0` boundary, cyclic THM-2796 factor typing, top-only UFD
valuation, and the formal `c=0` and `v=0` hostiles.  Normal,
optimized, and stored transcripts agree; both LF hashes, the
assert/float AST audit, and the documentation gate pass.

**QED.**

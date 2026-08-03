---
id: THM-3231
title: "Fraction-free simple-pivot wall and second-order blowup carry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The row after a simple fraction-free pivot clutch is universally divisible
  coefficientwise by the square of the vanishing pivot.  Its strict transform
  is an explicit 4-by-4 osculating Toeplitz determinant, equivalently the
  third-and-higher tail of one reciprocal quotient.  At the THM-3223 p=43
  wall, the mod-p terminal row lifts after division by 43^2 to the unit 19.
source: root/multiscale-newton-flag/2026-08-02
audit: >
  An independent hostile audit rederived the exceptional-square cancellation,
  every determinant index and sign, the quotient-tail convolution and repaired
  unit-scoped iff, the DVR strict-transform law, and all p=43 residues.  Fresh
  ordinary and optimized replays byte-match the stored transcript and the
  declared LF-normalized hashes.
depends_on:
  - THM-3223-fourth-fifth-resonant-prs-primitive-walls-pell-content-and-pivot-resurrection
related:
  - THM-3192-reciprocal-coefficient-jet-transfer-and-z-adic-pluecker-return
  - THM-3214-two-jet-pseudo-division-locality-and-catalan-sharpness
  - THM-3220-root-four-jet-schwarzian-heisenberg-transgression-and-oriented-discriminant-holonomy
script: 04-computation/factorial_simple_pivot_wall_blowup_carry_thm3231.py
output: 05-knowledge/results/factorial_simple_pivot_wall_blowup_carry_thm3231.out
script_sha256: 04cd91ea93698695105b4399a0a2b90beba2938b15bb9cdf6c894d373ef42c44
output_sha256: 78859ee0ba98de61f7ba24f4c51636fd9505deb7ae49cec9997fc3033a9dfe29
hash_basis: LF-normalized bytes
---

# THM-3231 -- fraction-free simple-pivot wall and second-order blowup carry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3223 proves that a simple anchored pivot wall is a one-row clutch modulo
the wall: the next row is a shifted copy and the following ordinary row is
zero.  That terminal statement is exact on the special fibre, but it loses
the normal direction.  The present theorem computes the strict transform.
The apparent zero has multiplicity two in the pivot parameter, and after
dividing by that exceptional square it becomes a fourth osculating minor.

## 1. The fraction-free map and its exceptional pair

For series over a commutative ring, retain

```text
P_1(f,g)=f_0g_1-f_1g_0,

E(f,g)=z^(-2){f_0^2g-[f_0g_0+P_1(f,g)z]f}.               (1)
```

Write

```text
w=sum_(j>=0) w_j z^j,       t=sum_(j>=0) t_j z^j,
x=E(w,t),                    y=E(x,w).                    (2)
```

Put `a=w_0` and split

```text
w=a+zv,             v=sum_(j>=0)v_jz^j,       v_j=w_(j+1). (3)
```

The divisor `a=0` is the pivot chart wall of THM-3223; its simple open locus
also has `w_1t_0` invertible.  On the whole divisor
`x=w_1t_0z^(-1)w`, so `y=0`.  The issue is its order of vanishing normal to
that divisor.

## 2. Universal exceptional-square divisibility

For every coefficient there is an integral universal polynomial `Q_j` such
that

```text
y_j=a^2 Q_j.                                              (4)
```

Thus the whole row, not merely its leading pivot, lies in the square of the
exceptional ideal `(a)`.

To prove `(4)`, direct coefficient extraction from `(1)` gives

```text
x_j=w_1t_0w_(j+1)
    -a(t_0w_(j+2)+t_1w_(j+1))+a^2t_(j+2),                (5)

y_j=x_0^2w_(j+2)-w_1x_0x_(j+1)
    +a(x_1x_(j+1)-x_0x_(j+2)).                           (6)
```

After substituting `(5)` into `(6)`, both the constant and linear terms in
`a` cancel identically.  The coefficient of `a^2`, which is the restriction
of `Q_j` to the exceptional divisor, is

```text
C_j:=Q_j|_(a=0)
 =t_0[
    t_0w_1^2w_(j+4)
   +(t_1w_1^2-t_0w_1w_2)w_(j+3)
   +(t_0(w_2^2-w_1w_3)-t_1w_1w_2+t_2w_1^2)w_(j+2)
   -t_(j+3)w_1^3].                                       (7)
```

No pivot has been inverted in `(4)--(7)`.

## 3. Osculating determinant and quotient-tail form

In the shifted coordinates `(3)`, formula `(7)` is the 4-by-4 determinant

```text
C_j=t_0 det [
 t_0  t_1  t_2  t_(j+3)
 v_0  v_1  v_2  v_(j+3)
  0   v_0  v_1  v_(j+2)
  0    0   v_0  v_(j+1) ].                               (8)
```

This makes the new coordinate an osculating Pluecker minor rather than an
opaque expansion artifact.  If `v_0` is a unit, write

```text
q=t/v=sum_(j>=0)q_jz^j,                 (Sq)_j=q_(j+1).   (9)
```

Triangular convolution in `(8)` gives the whole strict-transform row at
once:

```text
sum_(j>=0) C_jz^j
 =-t_0 v_0^3 v(z) S^3(t/v).                              (10)
```

Equivalently,

```text
C_j=-t_0v_0^3 sum_(ell=0)^j v_ell q_(j-ell+3).           (11)
```

At the first coordinate,

```text
C_0=t_0 P_1(E(v,t),v)=-t_0v_0^4 q_3.                    (12)
```

On the simple wall locus `t_0v_0 in R^*`, the exceptional row vanishes
identically exactly when `t/v` is a polynomial of degree at most two.  Its
first coordinate is nonzero exactly when the third reciprocal-quotient
coefficient is nonzero.  Without the `t_0` unit hypothesis, `(10)` remains
the exact statement but this iff can fail by scalar annihilation.  This is
the sharp boundary: the ordinary special-fibre recurrence sees zero, while
the normal cone sees the first unremoved quotient tail.

## 4. DVR carry law

Let `R` be a DVR with uniformizer `pi`, and suppose all displayed
coefficients are integral.  If

```text
w_0=pi^r eta,               eta in R^*,                  (13)
```

then `(4)` gives

```text
y_j/pi^(2r) == etabar^2 Cbar_j                 (mod pi). (14)
```

Consequently, if some `Cbar_j` is nonzero, that coefficient of `y` has
valuation exactly `2r`.  The exponent two is universal and sharp.  The
special fibre `ybar=0` is thus compatible with a nonzero second-order carry;
one must distinguish the ordinary reduction from the Rees/normal-cone
coordinate `(y/a^2)|_(a=0)`.

## 5. The p=43 factorial wall carries the unit 19

Apply the theorem to the offset-two wall of THM-3223.  Thus `p=43`, `d=45`,
and the fraction-free rows are

```text
R_1=r, R_2=s, R_3=t, R_4=w, R_5=x, R_6=y.                (15)
```

Thirteen exact reciprocal top jets suffice for the sixth constant pivot.
After exact rational cancellation every denominator is a 43-unit.  Direct
calculation modulo `43^10` gives

```text
lengths:             13,11,9,7,5,3,1,

v_43(r_0,s_0,t_0,w_0,x_0,y_0)=(0,0,0,1,0,2),            (16)

w mod 43=(0,24,20,16,3),
x mod 43=(23,12,1).                                      (17)
```

The normal data are

```text
w_0/43 ==36,             C_0==39,             y_0/43^2==19
                                                        (mod 43), (18)
```

and indeed

```text
36^2*39==19                                      (mod 43). (19)
```

Thus THM-3223's exact mod-43 terminal row is not terminal in the normalized
43-adic atlas.  Its second-order strict transform is already a unit in the
first coordinate.

## 6. Connection contract and scope

```text
source:       two consecutive fraction-free rows near a simple pivot wall;
map:          second iterate followed by division by the exceptional square;
target:       a 4-by-4 osculating determinant / third quotient tail;
preserved:    integral coefficients, normal direction, p-adic valuation;
destroyed:    absolute scaling before normalization and any root ownership;
sidecar:      a unit first shifted coefficient v_0 and a live quotient tail.
```

The determinant in `(8)` is structurally parallel to the four-jet Pluecker
curvature in THM-3220, but no theorem here identifies their underlying
polynomial root.  Likewise, `(14)` is a p-adic carry, not a lawful LRC packet
or a global root selector.  The result does not factor the universal sixth
pivot, prove arbitrary-depth PRS separation, settle arbitrary radial
coefficients, or add an `NC(2)`, `GMC(2)`, or `LRC(14)` decrement.

The positive survivor is precise: a chosen scalar pivot may die modulo `p`,
its next ordinary row may vanish, and yet the coefficientwise Rees transform
retains an explicit nonzero osculating coordinate.  This is the chart-gluing
operation missing from a scalar wall list.

## 7. Exact evidence

Run

```text
python 04-computation/factorial_simple_pivot_wall_blowup_carry_thm3231.py
python -O 04-computation/factorial_simple_pivot_wall_blowup_carry_thm3231.py
```

and compare LF-normalized bytes with the declared output.  The companion
pins promoted THM-3223; checks five symbolic coefficients of `(4),(7),(8)`
and `(10)`; verifies `(12)` by exact reciprocal series division; rebuilds
the thirteen factorial top jets without importing a discovery cache; and
proves every valuation and residue in `(16)--(19)` modulo `43^10`.  It uses
no floating point, randomness, optimizer, or assertion-sensitive check.

QED.

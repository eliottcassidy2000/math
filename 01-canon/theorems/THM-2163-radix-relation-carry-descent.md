---
id: THM-2163
title: "Radix relation-carry descent and full-box CRT blindness"
status: >
  PROVED. Every nonzero integer relation has an exact base-q carry path with
  fewer than 2||m||_1 states; for a positive row the sign-sharp count is only
  ||m||_1-1. The path starts and terminates at zero, and the digit recurrence
  has a converse reconstruction theorem. Quotient-owner masks are the missing
  termination sidecar. An explicit primitive defect-seven pair has identical
  labelled residues for q=2,...,13, identical full coefficient-height-29
  relation boxes, and the same maximum, but its mod-17 margins are 0 and 2.
  Thus a fixed CRT bank plus bounded relations plus scalar magnitude still
  mixes the next modular certificate.
source: codex-2026-07-24-relation-carry-spectrum
depends_on: []
related:
  - THM-2051
  - THM-2052
  - THM-2072
  - THM-2144
  - THM-2145
  - THM-2161
script: 04-computation/lrc14_radix_relation_carry_fiber_mixing_referee_codex_20260724.py
output: 05-knowledge/results/lrc14_radix_relation_carry_fiber_mixing_referee_codex_20260724.out
script_sha256: aa106b288f44729b5ee6a835801904b83bc2e12437e3695d874e49c7e19f704e
output_sha256: f2adaec56ccf3fa29379fb9803ac1038ef1eb36582708a0c22da8046daecf8b6
hash_basis: working-tree bytes (LF)
---

# THM-2163 -- radix relation-carry descent

The relation coefficients in THM-2144 and THM-2145 are bounded, but the
speeds need not be. This theorem supplies the exact scale coordinate that a
carry dynamic program must retain.

## 1. Exact carry theorem

Let

```text
V=(V_1,...,V_d) in Z_(>=0)^d,
0!=m=(m_1,...,m_d) in Z^d,
q in Z, q>=2.                                           (1)
```

For `j>=0`, take coordinatewise Euclidean divisions

```text
V=q^j Z_j+R_j,          0<=R_(j,i)<q^j,                (2)
```

and let

```text
D_j=Z_j mod q in {0,...,q-1}^d.                        (3)
```

Thus `D_j` is the `j`-th base-`q` digit vector of `V`, and

```text
Z_j=q Z_(j+1)+D_j.                                     (4)
```

Suppose `m.V=0`. Define

```text
kappa_j=(m.R_j)/q^j=-m.Z_j.                            (5)
```

Then every `kappa_j` is an integer and

```text
kappa_0=0,                                             (6)
q divides kappa_j+m.D_j,                               (7)
kappa_(j+1)=(kappa_j+m.D_j)/q,                         (8)
|kappa_j|<||m||_1.                                     (9)
```

If `q^J>max_i V_i`, then

```text
Z_J=0,               kappa_J=0.                        (10)
```

Thus a relation of coefficient `l1` norm `C` travels through at most

```text
2C-1                                                       (11)
```

integer carry states, independently of the speed magnitudes.

### Proof

Taking the dot product of (2) with `m` and using `m.V=0` gives

```text
m.R_j=-q^j m.Z_j,                                      (12)
```

which proves integrality and the second equality in (5). Equations (4)--(5)
give

```text
kappa_j=-q m.Z_(j+1)-m.D_j
       =q kappa_(j+1)-m.D_j,                           (13)
```

which is (7)--(8). Since every coordinate of `R_j/q^j` lies in `[0,1)`,

```text
|kappa_j|
 <=sum_i |m_i| R_(j,i)/q^j
 <sum_i |m_i|=||m||_1.                                (14)
```

Equations (6) and (10) follow from `R_0=0` and `Z_J=0`. QED.

For the positive rows relevant to LRC there is a sign-sharp state count.
Write

```text
P=sum_(m_i>0)m_i,              N=sum_(m_i<0)(-m_i).
```

Both are positive because `V_i>0`, `m.V=0`, and `m!=0`. The first expression
in (5), together with `0<=R_(j,i)/q^j<1`, gives

```text
-N<kappa_j<P.                                         (14a)
```

The integer interval in (14a) contains exactly

```text
P+N-1=||m||_1-1                                       (14b)
```

values. This refinement uses positivity; the symmetric bound (11) remains
the correct statement for the nonnegative generality of (1).

## 2. Converse reconstruction

The carry rule loses no relation information. Let digit vectors

```text
D_0,...,D_(J-1) in {0,...,q-1}^d                       (15)
```

and integer carries `kappa_0,...,kappa_J` obey

```text
q kappa_(j+1)=kappa_j+m.D_j,
kappa_0=kappa_J=0.                                     (16)
```

Reconstruct

```text
V=sum_(j=0)^(J-1) q^j D_j.                             (17)
```

Multiplying (16) by `q^j` and summing telescopes:

```text
m.V=sum_j q^j m.D_j
   =sum_j(q^(j+1)kappa_(j+1)-q^j kappa_j)
   =q^J kappa_J-kappa_0=0.                             (18)
```

Hence (6)--(8), with terminal zero, are equivalent to the relation on the
finite digit word.

At the first scale,

```text
D_0=V mod q.                                           (19)
```

This is exactly the labelled residue word used by a modulus certificate.
The relation adds the affine congruence

```text
m.D_j=-kappa_j mod q                                   (20)
```

and (8) transports it to the next digit.

## 3. The termination sidecar

Define the quotient-owner mask

```text
O_j={i:Z_(j,i)>0}={i:V_i>=q^j}.                        (21)
```

The masks are nested and eventually empty. They record the information that
the bounded carry alone destroys:

```text
source:       the integer row and a bounded relation;
quotient:     the finite carry state kappa_j;
preserved:    the affine digit congruence (20);
destroyed:    which coordinates still have higher digits;
sidecar:      O_j, or the full quotient Z_j;
termination:  O_J=empty and kappa_J=0.                 (22)
```

A zero digit does not imply that its coordinate has terminated: higher
digits can remain. Therefore a carry-only automaton can cycle and cannot
bound the search depth. The owner mask is not optional bookkeeping.

If `|m_i|<=H`, then `||m||_1<=dH`. For the positive thirteen-coordinate
relations of THM-2144 at height `29`, (14b) gives at most `376` carry values
before primitive normalization; for THM-2145's low-core height `298`, it
gives at most `3873`. These are real finite state bounds, but neither supplies
a bound on the number of digit levels.

## 4. A relation-preserving CRT hostile pair

For a row `V`, define its complete height-29 relation box

```text
R_29(V)={c in Z^13: ||c||_infinity<=29 and c.V=0}.     (23)
```

Put `L=lcm(2,...,13)=360360` and consider

```text
A=(1,2,3,4,5,6,
   1441447,44324288,1331890569,39958158250,
   1198739341811,35962180254012,1078865413025413),

B=(1,2,3,4,5,6,
   360367,13693688,412972569,12384492130,
   371535484331,11146064529612,1078865413025413).      (24)
```

Both rows are primitive and have thirteen distinct coordinates. They agree
with `(1,...,13)` label by label modulo `L`, replace exactly the seven AP
entries `7,...,13`, and share the maximum

```text
1078865413025413.                                      (25)
```

For every tail index `j=7,...,13` in either row,

```text
V_j>29 sum_(i<j)V_i.                                   (26)
```

The seven positive differences in (26) are

```text
A: 1440838,2521716,4683645,6124825,719136,718818,6123871,
B: 359758,3242436,5404365,719425,1439856,1439538,
   744483478576591.                                    (27)
```

Suppose `c.V=0`, `||c||_infinity<=29`, and take the largest tail index with
`c_j!=0`. Then

```text
|c_j|V_j>=V_j
 >29 sum_(i<j)V_i
 >=|sum_(i<j)c_i V_i|,                                 (28)
```

a contradiction. Every tail coefficient is zero. Since the first six
coordinates of the two rows coincide,

```text
R_29(A)=R_29(B)
 ={(c_1,...,c_6,0,...,0):
   |c_i|<=29 and sum_(i=1)^6 i c_i=0}.                 (29)
```

The exact companion counts `98,935,667` vectors in this common box, including
zero.

## 5. The adaptive-prime split

Because the rows agree coordinatewise modulo `L`, they have identical
labelled residues at every modulus `2<=q<=13`. In fact their `q`-th labelled
coordinate is zero modulo `q`, so

```text
m_q(A)=m_q(B)=0,                 2<=q<=13.             (30)
```

At modulus `17`, the tail residues are

```text
A: 0,1,1,1,1,1,1,
B: 1,1,1,1,1,1,1.                                    (31)
```

The zero gives `m_17(A)=0`. For `B`, multiplier `2` leaves every core and
tail residue at centered distance at least `2`, while the seven points
`0,a,2a,...,6a` show that some core residue is within `2` of zero for every
nonzero `a mod 17`. Hence

```text
m_17(B)=2,
Gap(B)>=2/17.                                          (32)
```

No assertion is made that `Gap(A)<2/17`; the exact loneliness values are not
being compared.

Equations (25), (29), and (30)--(32) prove that the quotient

```text
V |->
((V mod q)_(2<=q<=13), R_29(V), max(V))                (33)
```

is target-mixed for the next modular certificate. This strengthens
THM-2161's `Q=13` fixed-bank hostile control along an orthogonal
sidecar-preservation axis: even the **complete** bounded relation box and
scalar maximum do not restore the lost adaptive digit. Neither theorem
implies the other in full, since THM-2161 is universal in the bank cutoff
and separates actual gaps.

## 6. Consequences and limits

The theorem gives the exact bridge requested by the relation-and-carry
spectrum:

```text
fixed residue prefix + relation
  -> bounded affine carry path
  -> adaptive next digit + owner/termination sidecar.  (34)
```

It does not make LRC(14) finite. A usable defect-six or defect-seven search
must still control the number of digit levels, connect the digit constraints
to a safe phase, or prove that a repeated carry/owner state can be pumped
without changing the target. The hostile pair proves that replacing the
owner profile by the scalar maximum is insufficient.

The companion script verifies all integers in (24)--(32), counts the common
core relation box by exact dynamic programming, and exhaustively checks both
directions of the carry theorem on `1,412` accepted small digit paths in
bases `2,3,4`. Normal and optimized Python execute the same raising checks
and produce identical output.

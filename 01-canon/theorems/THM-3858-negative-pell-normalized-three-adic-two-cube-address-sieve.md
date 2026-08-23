---
id: THM-3858
title: "Negative-Pell normalized three-adic two-cube address sieve"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the
  negative-Pell successor orbit, dividing the two nontrivial three-adic
  branches by their exact powers of three recovers the unit part of the Pell
  depth modulo nine, with opposite signs on the two adjacent indices.  A
  scaled primitive two-cube address therefore excludes unit parts 4 and 5
  modulo nine, removing two of the six raw unit classes.  This is sharp for
  unconstrained local existence over Z_3: every surviving unit is a sum of
  two Z_3 cubes.  Prescribed summand data, positivity, global integrality,
  other places, LRC owner/phase and LRC(14) remain open.
source: root + pell_depth_extension + pell_mod9_hostile_audit, 2026-08-23
audit: >
  TWO INDEPENDENT HOSTILE AUDITS PASS.  The first uses an order-three
  binomial truncation, opposite residue charts, 4,096 Pell depths, 540
  normalized unit controls, and complete cube images through 3^7.  The
  second avoids that derivation: it stabilizes
  E_s=(B^(3^s)-I)/3^(s+1) modulo nine, checks both direct Pell offsets,
  exhausts the valuation law through depth 1,200, and constructs every local
  target through 3^9.  The 19,218- and 33,234-gate companions agree under
  normal, optimized, and frozen replay.  The audit restricts the stopping
  claim to the unconstrained two-summand Z_3 projection after exact scale
  removal.
related:
  - THM-3742-square-triangular-pell-mod13-central-sign-projective-cycle
  - THM-3793-inert-prime-sum-all-scale-two-cube-singleton
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-3825-prime-colour-valuation-two-cube-decoder
valuation_script: 04-computation/lrc14_negative_pell_three_adic_depth_sieve_20260823.py
valuation_output: 05-knowledge/results/lrc14_negative_pell_three_adic_depth_sieve_20260823.out
valuation_script_sha256: f2b4afc9e6c8491f7105d2de2bd8f2d5259fb674e1bc39a9970537eea436ae05
valuation_output_sha256: f60042e9b76cf5d8b60888b2a696db71fb4c22ebccb279666ce40f112ae0a34e
extension_script: 04-computation/lrc14_negative_pell_normalized_unit_extension_thm3858.py
extension_output: 05-knowledge/results/lrc14_negative_pell_normalized_unit_extension_thm3858.out
extension_script_sha256: aa37c15564920f455e709a40e2b3327467172d81ec2f35c115fa493844105e44
extension_output_sha256: 8c951f7ef5428a9441f39cd728d8b136b3edc16bbe69ad4d19fc9545428d0144
extension_semantic_sha256: de96643f7b05fc26e9d804ff4abe21a37aa6d8d07f446e923560c32da0104d64
independent_script: 04-computation/lrc14_negative_pell_normalized_unit_independent_audit_thm3858.py
independent_output: 05-knowledge/results/lrc14_negative_pell_normalized_unit_independent_audit_thm3858.out
independent_script_sha256: 26d2d4b1c881873037363d8859909f19d0501d5886ba3aaae693f9c8724533bb
independent_output_sha256: 185ea67f317bdd15aacfbbb492d94990f9d9e9fa4273e0144227952f4f8d5d1f
independent_semantic_sha256: fa0d2ddf38d16ac88722ee94bfa52aa8490d291550f4420b87e9db844cb03fe7
hash_basis: raw LF bytes
---

# THM-3858 -- normalized Pell depth retains one final three-adic digit

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Put

```text
A=((2,3),(1,2)),             (X_k,Y_k)^t=A^k(1,1)^t,
Q_k=(Y_k^2-1)/8.                                               (1)
```

For every `j>=1`,

```text
v_3(Q_(3j-1))=v_3(Q_(3j))=1+v_3(j),
v_3(Q_(3j+1))=0.                                               (2)
```

Let `s=v_3(j)>=2` and `h=j/3^s`.  Then the integral normalized values obey

```text
Q_(3j)   /3^(s+1) ==  h  (mod 9),
Q_(3j-1) /3^(s+1) == -h  (mod 9).                              (3)
```

Consequently suppose, on either nontrivial branch, that

```text
Q_k=g^3(a^3+b^3),       v_3(g)=r>=1,
3 does not divide a^3+b^3.                                    (4)
```

Then

```text
s=3r-1,                    h mod 9 in {1,2,7,8}.                (5)
```

The unit classes `4,5` are impossible.  At `r=1`, the first excluded depths
and their normalized residues are

```text
k=107:5,       k=108:4,       k=134:4,       k=135:5.          (6)
```

Thus (5) removes exactly two of the six raw unit residue classes, not
necessarily one third of an already filtered or nonequidistributed subfamily.

## 1. Exact valuation law

Set

```text
C=((9,15),(5,9)),              B=I-3C=-A^3,
w_+=(1,1),                     w_-=A^(-1)w_+=(-1,1).            (7)
```

For either seed `w`, let `ell` select the second coordinate and write

```text
d=ell(B^j w)-1
 =sum_(m>=1)(-3)^m binom(j,m) ell(C^m w).                       (8)
```

The leading moments are

```text
ell(Cw_+)=14,                  ell(Cw_-)=4,                     (9)
```

both three-adic units.  For `m>=2`,

```text
v_3(binomial(j,m))>=v_3(j)-v_3(m),
m-v_3(m)>=2.                                                     (10)
```

The first term of (8) therefore has valuation `s+1`, while every later term
has valuation at least `s+2`.  Since

```text
Q=d(2+d)/8,                                                      (11)
```

and `(2+d)/8` is a three-adic unit, this proves the first line of (2).
For `k=3j+1`, the second coordinate of `Aw_+` is zero modulo three, so
`Y_k^2-1` is a unit.  This proves the remaining line.

The indexing is exact:

```text
A^(3j)w_+   =(-1)^j B^j w_+,
A^(3j-1)w_+=(-1)^j B^j w_-.                                   (12)
```

The global sign disappears after squaring, so the plus and minus branches in
(3) cannot be interchanged by an indexing convention.

## 2. The normalized mod-nine digit

After division by `3^(s+1)`, every order `m>=4` in (8) vanishes modulo nine.
Writing `c_m=ell(C^m w)`, the first three orders reduce to

```text
d/3^(s+1)==h[-c_1+3(c_2-c_3)]  (mod 9).                       (13)
```

The exact moment triples are

```text
(c_1,c_2,c_3)(w_+)=(14,246,4344),
(c_1,c_2,c_3)(w_-)=(4,66,1164).                                (14)
```

Hence

```text
d_+/3^(s+1)==4h,             d_-/3^(s+1)==5h  (mod 9).         (15)
```

Because `s+1>=3`, equation (11) contributes `(2+d)/8==7 mod 9`.
Multiplication by seven changes `4` to `1` and `5` to `-1`, proving (3).

An independent certificate avoids this truncation.  Define

```text
E_s=(B^(3^s)-I)/3^(s+1).                                       (16)
```

Directly,

```text
E_2==((0,3),(4,0))  (mod 9).                                   (17)
```

If `B^(3^s)=I+3^(s+1)E_s`, cubing gives

```text
E_(s+1)=E_s+3^(s+1)E_s^2+3^(2s+1)E_s^3.                       (18)
```

Thus `E_(s+1)==E_s mod 9` for all `s>=2`, and raising to the unit power `h`
recovers independently the coefficients `4h,5h` in (15).

## 3. Scaled two-cube addresses

In (4), equation (2) gives

```text
s+1=v_3(Q_k)=3v_3(g)=3r,
```

which is the first part of (5).  After division by `3^(3r)`, the target is a
unit sum of two cubes, including the cube of the unit part of `g`.  Cubes
modulo nine are `0,+1,-1`, so a unit sum of two cubes lies in

```text
{1,2,7,8}.                                                       (19)
```

This set is stable under negation and multiplication by a unit cube.  Formula
(3) therefore gives the second part of (5).  Taking `s=2` and the first
excluded unit parts `h=4,5` gives (6).

The application is typed.  THM-3818's scaled support-two packet retains the
common cube scale, so `r` is available and (5) genuinely removes an infinite
address subfamily.  It does not restore the other eleven runners, owner,
phase, arrival, or a common lonely time.

## 4. The exact three-adic stopping boundary

For a unit `z in Z_3`,

```text
z=x^3+y^3 for some x,y in Z_3
iff z mod 9 is in {1,2,7,8}.                                   (20)
```

Necessity is (19).  For sufficiency, the three-adic logarithm gives

```text
(Z_3^*)^3={z in Z_3^*:z==+1 or -1 (mod 9)}.                    (21)
```

If `z==+/-1 mod 9`, then `z-27` is a unit cube and
`z=(z-27)+3^3`.  If `z==2 mod 9`, then `z-1` is a unit cube; if
`z==-2 mod 9`, then `z+1` is a unit cube.  This proves (20).

Therefore no modulus `3^n` can further prune the **unconstrained local
existence problem (20)** after exact common-cube-scale removal.  This does
not close prescribed valuations or residues of the individual summands,
labelled pair data, positivity, distinctness as global integers, global
integrality, another prime, height constraints, or any LRC sidecar.  Higher
three-adic information may still matter when one of those coordinates is
retained.

## 5. Exact replay and scope

The first extension companion also independently checks the repaired
valuation-one Poincare residues of THM-3842; that cross-frontier control
changes no Pell or JC statement.  Its Pell half verifies 4,096 full depths,
540 normalized cases and unit cube images through `3^7`.  The independent
block-lift companion verifies 264 lift cases, 240 direct offset cases, the
valuation law through depth 1,200 and 19,680 local targets through `3^9`.
Both are assertion-free and agree with their frozen outputs under normal and
optimized Python.

Reproduction:

```text
python -B 04-computation/lrc14_negative_pell_normalized_unit_extension_thm3858.py
python -B -O 04-computation/lrc14_negative_pell_normalized_unit_extension_thm3858.py
python -B 04-computation/lrc14_negative_pell_normalized_unit_independent_audit_thm3858.py
python -B -O 04-computation/lrc14_negative_pell_normalized_unit_independent_audit_thm3858.py
```

The result is an address sieve and its exact local stopping theorem.  It does
not prove global two-cube uniqueness, LRC(14), JC(2), or a bridge between the
Pell recurrence and the rational cubic Keller tower.  **QED.**

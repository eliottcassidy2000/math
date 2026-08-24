---
id: THM-4024
title: "LRC(14) complete divisor-incidence envelope"
status: >
  PROVED RELATIVE TO THM-4004'S CITED LRCUpTo13 + VERIFIED-EXACT +
  INDEPENDENTLY HOSTILE-AUDITED. The divisor-incidence sequence of every
  unresolved rank-eleven `11+2` row is bounded by `(11,9,8,8,7,7,...)`.
  This strengthens the prime-only profile to every divisor of t and gives
  exact subset-gcd compression. It is a necessary reduction, not LRC(14).
source: root + oeis_transfer + residue_sequence_scout / sequence continuation, 2026-08-24
audit: >
  PASS. A separate gcd/order-profile engine found only the d=3 and d=4
  equality states through d=500 and constructed a primitive distinct d=4
  selector hostile. The tempting d=6 equality was independently realized
  and shown nonprimitive. Primary and audit outputs match in normal and
  optimized Python; all hashes below are raw LF bytes.
depends_on:
  - THM-4004-lrc14-three-detuned-divisor-comb-profile
related:
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse
  - THM-668-detuned-harmonic-dispatch
script: 04-computation/lrc14_complete_divisor_incidence_envelope_thm4024.py
output: 05-knowledge/results/lrc14_complete_divisor_incidence_envelope_thm4024.out
independent_audit_script: 04-computation/lrc14_complete_divisor_incidence_independent_audit_thm4024.py
independent_audit_output: 05-knowledge/results/lrc14_complete_divisor_incidence_independent_audit_thm4024.out
script_sha256: 9e3ffa63b8f8a48cdcfe78b222bf422542cf20ffa518523922700bb2a4d70d50
output_sha256: eab07bed5b66c0751010df0f6ce83ffd60451ce868171a50ea033c7df241843a
independent_audit_script_sha256: 6c8c87ca5d36181ccc4fd9254e636aef50938401de9674236adadb41ed4ac587
independent_audit_output_sha256: 8e4391dfedf27499c7f5f702889c5a1f8a007f4dcabe0775ff912a0242284d61
hash_basis: raw LF bytes
---

# THM-4024 -- the complete divisor-incidence envelope

**PROVED RELATIVE TO THM-4004'S CITED LRCUpTo13 + VERIFIED-EXACT +
INDEPENDENTLY HOSTILE-AUDITED.** The argument below is an all-divisor
consequence of THM-4004, not a finite extrapolation. The exact companion
checks the strict open-comb constants and boundary states; an independent
gcd/order engine verifies the proof by a different representation. LRC(14)
remains open.

## 1. Inheritance pass and the sequence object

Retain THM-4004's exact rank-eleven two-component row

```text
n=s(u_1,...,u_11) direct-sum t(p,q),
s in {1,2}, gcd(s,t)=gcd(p,q)=1.                      (1)
```

For a positive divisor `d|t`, define the divisor-incidence sequence

```text
c_d(u)=#{i in {1,...,11}: d divides u_i}.              (2)
```

The closest proved mechanism is THM-4004's common-branch comb: if `c_d=8`,
the eight divided body owners and the divided pair form a ten-speed harmonic
pack, leaving three labelled detuned owners. Its canonical hostile is the
order-three equality row. The corrected near miss is to inspect only prime
divisors. The least-used sidecar is the lower-divisor incidence `c_g` for
`g=gcd(d,u_j)`.

The sequence `(2)` is antitone in the divisibility poset:

```text
e|d  ==>  c_d<=c_e.                                   (3)
```

An ordinary list of its values forgets this operation. It also forgets the
owner subsets, cited pack phase and branch overlaps; those remain necessary
for any sufficient certificate.

## 2. The inherited two-exception floor

THM-4004 Section 2 proves, for every divisor `d|t`, that `c_d=10` closes and
that `c_d=9` can remain unresolved only when `d=2`. Hence

```text
c_2<=9,                    c_d<=8 for every d>=3.       (4)
```

This all-divisor statement, rather than only its prime specialization, is the
input that makes the complete envelope possible.

## 3. Three exceptions at an arbitrary divisor

Assume `d>=3` and `c_d=8`. Let `delta_1,delta_2,delta_3` be the three actual
detuned body speeds. Because `gcd(s,d)=1`, put equivalently

```text
g_j=gcd(d,delta_j)=gcd(d,u_j),       m_j=d/g_j.         (5)
```

THM-4004's exact number of bad labelled branches is

```text
b_d(delta_j)=g_j m_j/7                   if 7|m_j,
               g_j ceil(m_j/7)           otherwise.   (6)
```

If some `g_j>=3`, then the eight coordinates divisible by `d` together with
`u_j` give `c_(g_j)>=9`. Since `g_j|d|t`, this contradicts `(4)`. Therefore

```text
g_j in {1,2} for j=1,2,3.                            (7)
```

Moreover, `c_2<=9` says that at most one of the three `g_j` can equal two.
Write

```text
f(m)=ceil(m/7)/m.                                      (8)
```

Formula `(8)` is an upper bound for the normalized count `b_d/d`; it is exact
unless `7|m`, when openness makes the true value the smaller `1/7`.

If `d>=5` is odd, all three `g_j` equal one, so

```text
sum_j b_d(delta_j)/d <= 3 f(d) <= 3/4 < 1.             (9)
```

Here `f(m)<=1/4` for every `m>=5`. If `d>=6` is even, at most one owner has
`g_j=2`; hence

```text
sum_j b_d(delta_j)/d
 <= f(d/2)+2f(d) <= 1/3+2/4=5/6 < 1.                 (10)
```

The elementary bounds used in `(9)--(10)` follow directly by writing
`m=7q+r`, `0<=r<7`: for `m>=3`, `f(m)<=1/3`, and for `m>=5`,
`f(m)<=1/4`. THM-4004's strict common-label union certificate now closes
every assumed `c_d=8` row with `d>=5`. Together with `(4)`, this proves

```text
c_d<=7 for every divisor d|t with d>=5.                (11)
```

This is the new gain. It includes composite divisors such as `6` and `12`;
no primality or prime-power assumption is present.

## 4. Sharp small-modulus boundary and prime-power readout

At `d=3`, the possible order profile is `(3,3,3)` and the normalized tariff
is exactly one. THM-4004's typed order-three hostile shows why a phase-free
strict union bound cannot improve this boundary.

At `d=4`, `(4)` and `c_2<=9` allow at most one even exception. If all three
exceptions are odd, their counts sum to `3<4` and the row closes. Therefore
an unresolved equality row must have exactly one exception of 2-adic
valuation one and two odd exceptions, giving order profile `(2,4,4)` and
bad-count sum `2+1+1=4`. Equality is deliberately left unresolved.

This boundary is selector-sharp even in a primitive distinct typed control:

```text
s=1, t=4, (p,q)=(1,4),
u=(8,12,20,24,28,32,36,40,2,9,11).                   (11a)
```

The divided pack is `{1,...,10}`. At its phase `y=1/11`, the bad label sets
of exceptions `(2,9,11)` are respectively

```text
{0,2}, {3}, {1}.                                       (11b)
```

They partition all four lifts. The row itself is not a counterexample: at
`x=21/22` its full clearance is `1/11`. Thus `(11a)--(11b)` is a hostile to
improving the phase-free strict count, not to LRC(14).

Thus the complete natural-number envelope is

```text
d:       1  2  3  4  5  6  7  8 ...
c_d:    11 <=9 <=8 <=8 <=7 <=7 <=7 <=7 ... .          (12)
```

Equivalently, its prime-power valuation profiles are

```text
c_(2^a) <= 9,8,7,7,...       for a=1,2,3,4,...,
c_(3^a) <= 8,7,7,...         for a=1,2,3,...,
c_(ell^a)<=7                 for ell>=5, a>=1.         (13)
```

These are valuation/counting sequences in the useful sense highlighted by
the sequence-work signal: prime powers turn divisibility into nested level
sets. No closed form is needed; the operation `(3)` is the structure.

## 5. Equivalent subset-gcd compression

For `I subset {1,...,11}`, put

```text
g_I=gcd(t,{u_i:i in I}).                               (14)
```

The envelope is equivalently the three statements

```text
|I|=10  ==> g_I=1,
|I|=9   ==> g_I in {1,2},
|I|=8   ==> g_I in {1,2,3,4}.                         (15)
```

Indeed, any divisor of `g_I` has incidence at least `|I|`, so `(12)` proves
`(15)`. Conversely, apply `(15)` to a subset of owners counted by `c_d`.
This is a compact exact filter for future searches: every common divisor at
least five is forbidden from every eight-owner body packet.

## 6. What is preserved and what is lost

The map from a body to `(c_d)_(d|t)` preserves divisor incidence, nested
prime-power levels, and the subset-gcd consequences `(15)`. It destroys
which owners realize different levels, the order of owners, the cited
ten-speed phase, the three bad branch sets and all geometric arrival data.
The missing sidecar is therefore the labelled family

```text
(j,a) |-> 1[||delta_j(y+a)/d||<1/14].                 (16)
```

The cheapest decisive test after `(12)` is not another scalar sequence. It
is a branch-overlap or phase-selection theorem for the equality moduli
`d=2,3,4`, followed by a CRT-compatible synthesis across distinct divisors.

## 7. Scope

This theorem is relative to the same cited `LRCUpTo13` input and exact
`11+2` decomposition as THM-4004. It applies only to divisors of the relative
scale `t`; it asserts nothing about `c_d` when `d` does not divide `t`.
It gives necessary conditions on a hypothetical counterexample, not a
construction of a lonely phase. The three small equality moduli, owner/phase
compatibility, both relative-scale lanes, and LRC(14) remain open.

Reproduce the current exact companion from the repository root with

```text
python3 -B 04-computation/lrc14_complete_divisor_incidence_envelope_thm4024.py
python3 -B -O 04-computation/lrc14_complete_divisor_incidence_envelope_thm4024.py
python3 -B 04-computation/lrc14_complete_divisor_incidence_independent_audit_thm4024.py
python3 -B -O 04-computation/lrc14_complete_divisor_incidence_independent_audit_thm4024.py
```

Both normal/optimized pairs reproduce their frozen outputs. **QED.**

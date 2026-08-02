---
id: THM-3055
title: "LRC grid witnesses form the unhit-divisor upset"
status: >
  PROVED + VERIFIED-EXACT.  For every modulus q>=2 and every nonempty finite
  integer speed set S, the threshold-1/q witnesses on the q-grid are exactly
  the reduced fractions b/d whose denominators d|q divide no speed in S.
  These surviving denominators form an upward-closed divisor upset, generated
  by an antichain of minimal unhit divisors, and the witness count is the sum
  of Euler phi(d) over that upset.  Adjoining a speed v deletes exactly the
  surviving divisors of gcd(v,q).  The q-grid witness set equals the nonempty
  reduced-denominator-q unit layer iff q divides no speed and every proper
  divisor d>1 of q divides at least one speed.  This classifies the grid
  carrier, not the full safe set, its measure, or tightness.
source: codex-2026-08-01-grid-witness-divisor-upset
depends_on: []
related:
  - THM-3043-lrc-tight-instances-are-not-only-APs-and-what-FC2-does-not-transfer
script: 04-computation/lrc_grid_witness_divisor_upset_thm3055.py
output: 05-knowledge/results/lrc_grid_witness_divisor_upset_thm3055.out
script_sha256: fca129b3546af8acb32b74a23babc0bd7100b3ddb642e2ceebdd98c794aae5a7
output_sha256: ffeefa32710d2e746fb59504274c1be063ff888da025447c3fe609f21f72cbe3
hash_basis: LF-normalized bytes
---

# THM-3055 -- LRC grid witnesses form the unhit-divisor upset

**PROVED + VERIFIED-EXACT.**

The small tight instances in THM-3043 suggest that the important object is the
witness set rather than the shape of the speed set.  On the distinguished
`q=n+1` grid, that object has an exact divisor-lattice description.  It also
explains why prime and composite witness denominators behave differently.

## 1. Grid-witness theorem

Let `q>=2`, let `S` be a nonempty finite set of positive integers, and use the
threshold `1/q`.  Define

```text
W_q(S)={a in {0,...,q-1}: ||av/q||>=1/q for every v in S},             (1)

U_q(S)={d:d|q and d does not divide v for every v in S}.              (2)
```

For `a in {0,...,q-1}`, put

```text
d(a)=q/gcd(a,q),                                                       (3)
```

the reduced denominator of `a/q` (with `d(0)=1`).  Then

```text
a in W_q(S)  iff  d(a) in U_q(S).                                     (4)
```

Equivalently, as a disjoint union of reduced fractions,

```text
{a/q:a in W_q(S)}
 = disjoint union over d in U_q(S) of
   {b/d:1<=b<d and gcd(b,d)=1}.                                      (5)
```

Consequently

```text
|W_q(S)|=sum_(d in U_q(S)) phi(d).                                   (6)
```

### Proof

For one speed `v`, the residue `av mod q` is either zero or lies in
`{1,...,q-1}`.  In the second case its circular distance from zero is at least
`1/q`; in the first case it is zero.  Therefore

```text
a in W_q(S)  iff  q does not divide av for each v in S.               (7)
```

Write `g=gcd(a,q)`, `a=gb`, and `q=gd`, so `gcd(b,d)=1`.  Then

```text
q|av  iff  gd|gbv  iff  d|bv  iff  d|v.                              (8)
```

Equations `(7)--(8)` prove `(4)`.  Numerators with `d(a)=d` are exactly the
`phi(d)` unit numerators modulo `d`, proving `(5)--(6)`.  Since `S` is
nonempty, `1` is hit by every speed, so the apparent `a=0` layer never enters
`U_q(S)`.  **QED.**

## 2. Divisor upset and its minimal antichain

Order the divisors of `q` by divisibility.  If `d in U_q(S)` and
`d|e|q`, then `e|v` would imply `d|v`; hence

```text
d in U_q(S), d|e|q  =>  e in U_q(S).                                (9)
```

Thus `U_q(S)` is an upward-closed **upset**.  It is uniquely determined by
the antichain `Min(U_q(S))` of its minimal elements.  This word is deliberate:
`U_q(S)` need not be closed under gcd and therefore need not be a lattice
filter.  For example, with `q=12` and `S={2}`, both `3` and `4` are unhit but
their gcd `1` is hit.

Adjoining a new speed has an equally exact operation law.  For any positive
integer `v`,

```text
U_q(S union {v})
 =U_q(S) setminus {d:d|gcd(v,q)}.                                    (10)
```

Indeed a divisor of `q` divides `v` exactly when it divides `gcd(v,q)`.
So the correct finite state for successive speed insertion is the minimal
unhit antichain, while `gcd(v,q)` is the complete update payload.  Raw speed
magnitudes and residues finer than their divisor incidence are forgotten.

## 3. Exact top-denominator boundary

The `q`-grid witness set is exactly the nonempty top unit layer

```text
{b/q:1<=b<q and gcd(b,q)=1}
```

precisely when

```text
U_q(S)={q}.                                                          (11)
```

By definition, `(11)` is equivalent to the two conditions

```text
q divides no v in S,                                                  (12)
every proper d|q with d>1 divides at least one v in S.                (13)
```

Condition `(12)` leaves the top divisor unhit; `(13)` hits every possible
lower denominator layer.  When `q` is prime, `(13)` is vacuous and every
nonzero numerator is a unit.  When `q` is composite, `(13)` explains the
missing nonunit numerators.

The six corrected THM-3043 rows all satisfy `(12)--(13)`:

```text
q=4, S={1,2,3}:       denominator 2 is hit;
q=5, two S rows:      q is prime;
q=7, S={1,...,6}:     q is prime;
q=6, S={1,3,4,5,9}:  denominators 2 and 3 are hit;
q=8, S={1,...,7}:     denominators 2 and 4 are hit.                    (14)
```

Hence their grid witnesses have reduced denominator exactly `q` for a
structural reason; this is no longer only a six-row observation.  In
particular, the `q=6` numerators are the units `{1,5}` and the `q=8`
numerators are the units `{1,3,5,7}`.

## 4. Sharp hostile and LRC consequence

It is false that every safe numerator on a composite grid must be a unit.
For `q=8` and `S={1}`,

```text
U_8(S)={2,4,8},
W_8(S)={1,2,3,4,5,6,7};                                              (15)
```

the fractions `1/2`, `1/4`, and `3/4` supply lower reduced denominators.
The first failed implication in the unit-only guess is the omission of
unhit proper divisors.

At the LRC threshold `q=n+1`, `(4)` immediately recovers the standard grid
certificate: if no speed is divisible by `n+1`, then `1/(n+1)` is safe.  More
generally, *any* unhit divisor `d|(n+1)` supplies all `phi(d)` unit witnesses
at denominator `d`.  Therefore a counterexample must hit every divisor of
`n+1` greater than one.  For `n=13`, a multiple of `14` already hits
`2,7,14`, so this necessary condition does not close LRC(14).

## 5. Connection contract and scope

The map is

```text
speed set S -> hit divisors {d|q:some d|v} -> minimal unhit antichain. (16)
```

It preserves exactly the safe points on the fixed `q`-grid and their reduced
denominators.  It destroys the arrangement between grid points, arc owners,
safe-set measure, and the maximum of the loneliness function.  The full
rational interval arrangement is the required sidecar for any claim about
tightness or covering.  Thus THM-3043's six *full* safe sets motivate `(16)`,
but the present theorem alone does not prove that an instance is tight, that
its safe set has measure zero, or that no off-grid witness exists.

## 6. Exact companion

The dependency-free companion independently computes the geometric grid
predicate and the divisor upset.  It exhausts all `8,177` nonempty residue
supports modulo every `q=2..12` (`90,035` grid cells), checks the upset and
totient count, checks `18,378` adjoining-speed updates, rebuilds the six
THM-3043 rows, and retains `(15)` as a hostile control.

Reproduce with

```text
python 04-computation/lrc_grid_witness_divisor_upset_thm3055.py
python -O 04-computation/lrc_grid_witness_divisor_upset_thm3055.py
```

Both modes byte-match the stored seven-line transcript.

---
id: THM-3804
title: "Rule 30 all-period amplitude-lattice Smith law"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  spatial period n and iterate r, the integral amplitude-lift image, kernel,
  cokernel, and Smith form are explicit: dyadic period factors create
  saturated free defects, while only the odd core carries one cyclic
  2-power obstruction.  The projective update halves every even declared
  period.  Consequently no physical translation-ray profile has dyadic
  spatial period, and THM-3778's odd-period no-go excludes physical exact
  finite scale-cycles at every finite declared period.  This proves no Rule
  30 center-column prize; all three prizes remain open.
source: root + rule30_lattice_audit / all-frontier session, 2026-08-23
audit: >
  PASS.  An independent derivation proved the odd determinant/index law, the
  even pair-sum/repetition factorization, primitivity of the periodic
  sublattice, and the physical period-halving boundary.  Exact controls
  include 216 Smith forms, 432 power-factorization gates, 812 kernel gates,
  3,168 constructive-lift gates, 898 image hostiles, 21,290 finite-field
  transition gates, and 4,098 physical-gap gates.  Normal and optimized
  executions agree with the frozen transcript and the companion contains no
  Python assert gate.
depends_on:
  - THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary
  - THM-3778-rule30-odd-period-finite-scale-cycle-projective-profile-no-go
related:
  - THM-3511-rule30-orbit-signalizer-gap-renormalization-and-shallow-portrait-hostile
  - THM-3516-rule30-marked-van-der-put-carry-and-power-section-bridge
script: 04-computation/rule30_all_period_amplitude_lattice_thm3804.py
output: 05-knowledge/results/rule30_all_period_amplitude_lattice_thm3804.out
script_sha256: f8cada48d0b12b7d18a51b0621c23fd5af3508527dd89f938f18dc5bf298977d
output_sha256: b4e5d373680fc0cf93b9e430e7800f189428c86d7a62f0b2ddef4a45e59fc051
semantic_sha256: cab8468cf3d6970a0cd5ce0bb17d36597ebb94e9339aecc41e33a3231fcce248
hash_basis: raw LF bytes
---

# THM-3804 -- the amplitude lift has one odd-core carry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The main
lattice theorem below is elementary and unconditional.  THM-3512 is used only
for the physical spatial-period corollary, and THM-3778 only for the final
scale-cycle corollary.  Every Rule 30 center-column prize remains **OPEN**.

For `n>=1`, label coordinates by `Z/nZ` and define the integral
amplitude-lift operator

```text
(T_n A)_j=A_(2j)+A_(2j+1).                            (1)
```

Fix `r>=0` and write

```text
n=2^a m,                 m odd,
b=min(a,r),              d=n/2^b,
s=r-b=max(0,r-a).                                    (2)
```

Let `Per_d(Z^n)` be the primitive sublattice of vectors whose period divides
`d`.

## 1. Exact image, kernel, cokernel, and Smith form

The image is exactly

```text
T_n^r Z^n
 = {B in Per_d(Z^n):
       sum_(j=0)^(d-1) B_j = 0 mod 2^s}.             (3)
```

Put `q=2^b=n/d`.  The kernel is

```text
ker(T_n^r)
 = {A: sum_(t=0)^(q-1) A_(qj+t)=0
       for every 0<=j<d},                            (4)
```

with saturated basis

```text
e_(qj+t)-e_(qj),      0<=j<d, 1<=t<q.               (5)
```

Consequently

```text
rank(T_n^r)=d,
coker(T_n^r) = Z^(n-d) direct-sum Z/(2^s),
SNF(T_n^r) = diag(1^(d-1),2^s,0^(n-d)),             (6)
```

where `Z/(1)=0`.  Equivalently,

```text
0<=r<=a:
  rank=n/2^r,
  SNF=diag(1^(n/2^r),0^(n-n/2^r));

r>a:
  rank=m,
  SNF=diag(1^(m-1),2^(r-a),0^(n-m)).                (7)
```

Thus period peeling creates only primitive free defects until the odd core is
reached.  Thereafter there is exactly one cyclic carry.

### Odd core

Suppose first that `n` is odd.  Let `S(A)_j=A_(j+1)` and
`P_2(A)_j=A_(2j)`.  Doubling is then a coordinate permutation, so

```text
T_n=P_2(I+S),                 |det(T_n)|=2.           (8)
```

Indeed `S` is one odd cycle and
`det(I+S)=1-(-1)^n=2`.  Also

```text
sum(T_n A)=2 sum(A).                                  (9)
```

The image of `T_n^r` therefore lies in the sum-congruence lattice
`sum B_j=0 mod 2^r`.  Both lattices have index `2^r`: the first by
`|det(T_n^r)|), the second because total sum maps onto `Z/(2^r)`.
Inclusion and equal index prove (3) on the odd core and give the cyclic Smith
factor.

The same proof shows a useful covariance statement.  For odd `n`, replace
each `P_2` at successive scales by an arbitrary coordinate permutation
`P_i`.  Every product

```text
P_r(I+S) ... P_1(I+S)                                (10)
```

still has image `{B:sum B_j=0 mod 2^r}`: permutations preserve total sum,
and the determinant index remains `2^r`.

### Peel the two-part

For even `k`, define

```text
Q_k:Z^k -> Z^(k/2),    (Q_k A)_j=A_(2j)+A_(2j+1),
E_k:Z^(k/2) -> Z^k,    (E_k C)_j=C_(j mod k/2).      (11)
```

Here `Q_k` is onto and `E_k` is injective.  Direct substitution gives

```text
T_k=E_k Q_k,
T_k E_k=E_k T_(k/2).                                 (12)
```

Let `Q_(n,b)` take the sums of the `d` consecutive blocks of length
`2^b`, and let `E_(n,d)` repeat a `d)-vector.  Iterating (12) yields the
load-bearing factorization

```text
T_n^r=E_(n,d) T_d^s Q_(n,b).                         (13)
```

The middle map is the identity when `r<=a` and is nonsingular on the odd
core when `r>a`.  Its image is the odd-core sum-congruence lattice already
proved.  Surjectivity of `Q_(n,b)` and injectivity of `E_(n,d)` now prove
(3)--(5).

Finally, `Per_d(Z^n)=E_(n,d)(Z^d)` is primitive: restriction to the first
`d` coordinates is a left inverse to `E_(n,d)`.  Its quotient in
`Z^n` is therefore free of rank `n-d`, while the quotient inside that
periodic summand is `Z/(2^s)`.  This proves (6)--(7), including every
zero-iterate and one-point boundary.

## 2. Even projective periods are nonminimal

For a labeled periodic profile on which the following rational map is
defined, put

```text
R(g)_j=-g_(2j)g_(2j+1)
        (1-g_(2j+2))/(1-g_(2j)).                     (14)
```

If `n=2d`, indices alone give

```text
R(g)_(j+d)=R(g)_j.                                   (15)
```

Thus one defined update removes one factor two from the declared spatial
period.  If `R^h(g)=g` for some `h>=1`, then an even declared period was
not minimal: `R(g)`, hence every later iterate and `g) itself, has period
dividing `n/2`.  Repetition reduces every defined finite projective
scale-cycle to its odd minimal spatial core.  No saturation beyond the
nonzero denominators in (14) is needed for this reduction.

### Physical consequences

Retain THM-3512's physical translation-ray profiles.  Every entry is an odd
2-adic unit, the update is (14), and the exact gap
`d_l=nu_2(1-G_l)` is phase-independent.  THM-3512 proves that no scale has a
phase-constant projective profile: a constant odd value `c` would update to
`-c^2`, and

```text
nu_2(1+c^2)=1                                        (16)
```

would force two consecutive gap-one scales, contradicting the inherited
no-`111` law.

If a physical profile had spatial period dividing `2^a`, repeated use of
(15) would make a later scale constant.  Hence

```text
no physical translation-ray profile has spatial period
dividing a power of two.                              (17)
```

More generally, a finite spatial period `2^a m`, with `m` odd, reduces
after `a` scales to a nonconstant odd period dividing `m`; in particular
`m>1`.

If such a physical profile also returned exactly after finitely many scales,
the preceding minimal-period reduction would place it in the odd-period
fully saturated setting of THM-3778.  That theorem excludes the return.
Therefore

```text
no physical Rule 30 translation-ray projective profile can have
both finite spatial period and an exact finite scale-cycle.             (18)
```

This is a spatial/projective obstruction, not a center-column
nonperiodicity theorem.

## 3. What the cyclic carry forgets

The quotient in (6) preserves integral lift existence exactly, but discards
the projective profile, phase owner, signed gauge, gap locations, scale, and
center-cell chronology.  Two small controls show why the carry cannot replace
those sidecars.

At `n=3,r=2`,

```text
A=(1,1,-2),                 T_3^2 A=A.                (19)
```

Its total sum is zero, so it passes every cyclic carry test, but its even
coordinate violates the all-odd physical type.  Conversely,

```text
A=(1,5,9),                  T_3 A=2(3,5,7).           (20)
```

The normalized parent remains all odd and has the same total sum, yet it is
not proportional to `A`.  Even a correctly typed one-step carry therefore
does not create projective closure.

The exact connection contract is:

```text
source:      labeled integral amplitudes before projectivization
target:      free period-defect coordinates plus one cyclic 2-carry
map:         dyadic block sums, then odd-core total sum modulo 2^s
preserved:   existence of an integral T_n^r-lift
destroyed:   projective ratios, owner, signed gauge, gaps, and chronology
sidecars:    lifted amplitude; THM-3511 owner; THM-3516 gauge/carry
hostile:     require closure under the next actual Rule 30 operation
```

## 4. Verification and scope

The independent companion verifies:

- `216` exact integer Smith forms for `1<=n<=24,0<=r<=8`;
- `432` factorization/rank and `812` kernel-basis gates;
- `3,168` constructive lifts and `898` free/torsion hostiles for
  `1<=n<=32,0<=r<=10`;
- `184` even-boundary and `108` permutation-covariance gates;
- `12` rational projective-halving gates;
- `21,290` saturated finite-field transition hostiles; and
- `4,098` physical odd-square and stopping-control gates.

Run

```bash
python3 -B 04-computation/rule30_all_period_amplitude_lattice_thm3804.py
python3 -B -O 04-computation/rule30_all_period_amplitude_lattice_thm3804.py
```

Both streams equal
`05-knowledge/results/rule30_all_period_amplitude_lattice_thm3804.out`.
The companion imports none of the scout code and contains no Python
`assert` gate.

This theorem proves no Rule 30 center-column nonperiodicity, balance or
density, computational lower bound, bounded innovation gap, finite
signalizer, or automatic quotient.  It makes no literature novelty claim.
Every Rule 30 prize remains **OPEN**.

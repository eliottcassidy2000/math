---
id: THM-2352
title: "q-adic prefix-residue collision spectrum and infinite-depth conditioning boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Canonical
  q-adic prefix residues have exponentially thin distinct support: every
  infinite support Dirichlet profile has abscissa zero. Indexed
  multiplicity equals support plus an exact plateau collision tax, and
  every indexed abscissa in [0,infinity], including infinity, occurs
  densely inside every finite digit cylinder. Ordinary stabilization is
  equivalent to finite support, but neither indexed convergence nor
  divergence characterizes it. Under Haar digits the indexed abscissa is
  zero almost surely, while eventual termination and every positive
  abscissa level are dense null tail sets; finite ANOVA/interaction
  quotients therefore erase them. The whole spectrum lifts with zero
  carry into every THM-2231 relation kernel admitting a binary digit ray.
  This is a formal carry/sequence theorem, not a Mahler Z-number, LRC(14),
  or finite-height result.
source: codex-2026-07-25-q-adic-prefix-collision-spectrum
depends_on:
  - THM-2000-support-harmonic-abel-dini-figurate-surface
  - THM-2005-support-dirichlet-automatic-tournament-atlas
  - THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization
  - THM-2231-relation-carry-completion-and-exact-radix-clock
related:
  - THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum
  - THM-2346-global-allocation-anova-normal-form
  - THM-2348-prime-type-rectangularity-and-target-token-conditioning
  - MISTAKE-209
  - MISTAKE-219
script: 04-computation/q_adic_prefix_collision_spectrum_thm2352.py
output: 05-knowledge/results/q_adic_prefix_collision_spectrum_thm2352.out
script_sha256: 7170c830f6d39c2eb12d0d09e4d94abae1f40f8842deb5f4de09ef1fdb842ce9
output_sha256: 1b3ed8e8ed8336b8e46743c0e62f11f06aeba339d035b779fc05076613a07105
hash_basis: working-tree bytes (LF)
---

# THM-2352 -- q-adic prefix residues have a complete collision spectrum

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The independent audit checked every index convention, the last infinite
plateau, both sides of the abscissa construction, arbitrary finite
cylinders, the Haar/Borel--Cantelli estimate, and the zero-carry relation
lift. It also separated the finite THM-2248/2346 interaction carriers from
the genuinely infinite-depth termination predicate. The exact companion
passes under ordinary and optimized Python and matches the stored
transcript after LF normalization.

THM-2000/2005 separates value support from indexing multiplicity.
THM-2228/2231 separates completed radix compatibility from ordinary
termination. On canonical prefix residues these are the same lost
coordinate: long zero-digit plateaux are precisely the collision tax.

## 1. Canonical prefix residues and their plateaux

Fix `q>=2` and a digit path

```text
x=sum_(n>=0)e_n q^n in Zhat_q,

e_n in {0,...,q-1}.                                  (1)
```

Use the canonical integer representatives

```text
r_m=sum_(0<=n<m)e_n q^n,          0<=r_m<q^m.        (2)
```

For real `s>0`, define the distinct-support and indexed profiles in the
extended nonnegative reals:

```text
S_x={r_m:r_m>0},

D_x(s)=sum_(r in S_x)r^(-s),

I_x(s)=sum_(m:r_m>0)r_m^(-s).                       (3)
```

Let

```text
n_0<n_1<...
```

be the positions of the nonzero digits. Put

```text
R_j=r_(n_j+1),

L_j=n_(j+1)-n_j                                     (4)
```

when a next active position exists. If `n_j` is the final active
position, set `L_j=infinity`.

The prefix sequence is nondecreasing and changes exactly at the active
positions. Hence its positive distinct values are exactly the `R_j`, and
`R_j` occurs with multiplicity `L_j`. Therefore

> **Exact plateau collision law.**
>
> ```text
> D_x(s)=sum_j R_j^(-s),
>
> I_x(s)=sum_j L_j R_j^(-s)
>
>       =D_x(s)+sum_j (L_j-1)R_j^(-s).              (5)
> ```

All identities hold in `[0,infinity]`; in particular an infinite final
plateau makes `I_x(s)=infinity`. Equation (5) is exactly THM-2000's
collision tax, now with the radix plateau length exposed.

## 2. Distinct support is always exponentially thin

At an active position,

```text
R_j>=q^(n_j).
```

Since `n_j>=n_0+j`,

```text
D_x(s)
 <=sum_(j>=0)q^(-s(n_0+j))
 =q^(-s n_0)/(1-q^(-s))<infinity                  (6)
```

for every `s>0`. Thus every **infinite** prefix-residue support has
Dirichlet abscissa exactly zero. Finite supports are Dirichlet
polynomials.

The following conditions are equivalent:

```text
x belongs to Z_(>=0) under the canonical q-adic embedding;

the digit path is eventually zero;

the prefix residues r_m are eventually constant;

S_x is finite;

D_x is a Dirichlet polynomial.                      (7)
```

The equivalence is literal: nonnegative integers are exactly the finite
base-`q` words, and every nonzero digit creates a new larger prefix value.
Consequently convergence of the support profile at any or all positive
arguments cannot detect ordinary stabilization.

## 3. Every indexed abscissa occurs in every cylinder

Define

```text
sigma_I(x)=inf{s>=0:I_x(s)<infinity},
```

with the infimum of the empty set equal to `infinity`.

Fix any finite digit prefix. Append one digit `1`, at or after its end, to
obtain an active position `n_0` and a positive residue `R_0`. For a
prescribed `alpha in [0,infinity)`, continue using only digits zero and
one, recursively setting

```text
L_j=ceil(R_j^alpha),

n_(j+1)=n_j+L_j,

e_(n_(j+1))=1,                                      (8)
```

with zero digits between successive active positions.

If `0<s<=alpha`, then every plateau block satisfies

```text
L_j R_j^(-s)
 >=R_j^(alpha-s)>=1,                                (9)
```

so `I_x(s)` diverges. If `s>alpha`, then

```text
L_j R_j^(-s)
 <=2R_j^(-(s-alpha))
 <=2q^(-(n_0+j)(s-alpha)),                          (10)
```

and the geometric majorant converges. Therefore

```text
sigma_I(x)=alpha.                                   (11)
```

For `alpha=infinity`, use

```text
L_j=R_j^j.                                         (12)
```

For every fixed `s`, the terms `R_j^(j-s)` fail to tend to zero once
`j>s`, so `I_x(s)` diverges. Thus all thresholds in
`[0,infinity]` occur inside **every** finite cylinder. The all-zero
continuation simultaneously supplies an ordinary terminating extension
inside that cylinder.

For `q=2`, THM-2228's cylinder bijection transports this conclusion to
every finite Mahler carry cylinder. It concerns the compatible canonical
residues only; it does not supply THM-2228's strict Archimedean tail
inequalities and therefore proves no `Z`-number.

## 4. Sharp reciprocal controls

For `q=2`, `alpha=1`, beginning with `e_0=1`, recursion (8) gives

```text
(n_j,R_j,L_j)
 =(0,1,1),
  (1,3,3),
  (4,19,19),
  (23,8388627,8388627),...                          (13)
```

Every harmonic plateau contributes

```text
L_j/R_j=1,                                         (14)
```

so `I_x(1)` diverges although the digit path never stabilizes.

At the opposite boundary,

```text
e_n=1 for every n
```

gives

```text
r_m=(q^m-1)/(q-1),             L_j=1,              (15)
```

and `I_x(s)=D_x(s)<infinity` for every `s>0`.

Finally, the ordinary integer `x=1` has

```text
r_m=1 for every m>=1.                               (16)
```

Its support is the singleton `{1}`, while `I_x(s)=infinity` for every
`s>0`. The nonstabilizing construction (12) also diverges for every
positive `s`. Hence:

```text
indexed convergence rules out positive stabilization
but is not equivalent to nonstabilization;

indexed divergence occurs on both terminating and
nonterminating paths.                               (17)
```

No fixed indexed convergence/divergence test characterizes termination.

## 5. Haar-generic paths and the infinite-depth boundary

Give the digits their uniform product, equivalently Haar, measure. Let
`E_n` be the event that the `n` digits beginning at position `n` are all
zero. Then

```text
P(E_n)=q^(-n),

sum_n P(E_n)<infinity.                              (18)
```

By Borel--Cantelli, almost surely only finitely many `E_n` occur. For all
large active positions this gives

```text
L_j<=n_j+1.                                        (19)
```

Together with `R_j>=q^(n_j)`,

```text
sum_j L_j R_j^(-s)
 <=sum_j (n_j+1)q^(-s n_j)<infinity                (20)
```

for every `s>0`. Intersecting the probability-one events for positive
rational `s` and using monotonicity gives

```text
sigma_I(x)=0                  almost surely.        (21)
```

Eventual-zero paths are countable, dense, and Haar-null. For each fixed
`alpha>0`, (8) puts an `alpha`-path in every cylinder, so the
`sigma_I=alpha` level is dense; (21) makes it Haar-null. The same is true
of the `infinity` level, which contains both terminating and
nonterminating paths.

This has two exact consequences.

1. The termination indicator is the zero element of Haar `L^2`; every
   finite-coordinate conditional expectation is identically zero. A naive
   infinite Hoeffding/ANOVA quotient erases the dense pointwise terminal
   set. This does not contradict finite THM-2346/2348.
2. If digit lift-supports are viewed in `P(N)`, termination is the family
   of finite subsets. Its restriction to every finite ground set
   `[0,N)` is the full simplex. There is no finite minimal nonface, so
   THM-2248's finite higher-interaction complex cannot detect this tail
   predicate. The required sidecar is an end marker, eventual-empty owner
   mask, or THM-2231's unbounded radix clock.

This is an infinite-depth/compactness obstruction, not merely an
arbitrarily large finite arity.

## 6. Exact zero-carry relation lift

Let

```text
0!=m in Z^d,

0!=u in {0,...,q-1}^d,

m.u=0.                                               (22)
```

For any **binary-valued** scalar digit path from Section 3, put

```text
D_j=e_j u.                                          (23)
```

These are valid coordinatewise base-`q` digits. In THM-2231's carry
recurrence,

```text
m.D_j=0,

kappa_j=0                 for every j,              (24)
```

and the prefix vector is

```text
V_<(M)=r_M u.                                       (25)
```

Every fixed ray norm merely rescales the profile; for example,

```text
||V_<(M)||_infinity=||u||_infinity r_M.             (26)
```

Thus the whole collision spectrum transfers into the completed relation
kernel while preserving zero carry and the active owner mask.

For THM-2231's radix-`13` support-three relation

```text
m=(1,0,...,0,1,-1),
```

take

```text
u=(1,2,3)
```

on coordinates `(1,12,13)` and zero elsewhere. Then `m.u=0`, so every
threshold above appears on one fixed formal relation ray. This has zero
coordinates and is not asserted to be a distinct-positive LRC row. It is
a carrier hostile proving that bounded carry and owner data do not retain
plateau multiplicity or termination.

## 7. Loss ledger and scope

The exact chain is

```text
digit path
 -> canonical prefix multiset
 -> distinct residue support
 -> Dirichlet evaluation or abscissa.               (27)
```

The first arrow retains the plateau clock as multiplicity. The support
quotient destroys it but still detects termination through **finiteness**.
Any positive reciprocal evaluation destroys even that distinction because
all infinite supports already converge. Indexed multiplicity recovers the
plateaux but still does not distinguish termination from sufficiently long
nonterminating plateaux.

This theorem proves no finite height bound, no Mahler `Z`-number, and no
LRC(14) row. Its reusable warning is sharper: neither finite interaction
data nor almost-everywhere completion may replace a pointwise terminal
sidecar.

## 8. Exact companion

The companion exhausts finite binary and ternary prefix banks, checks the
support/collision identity and exponential support bound with exact
fractions, reconstructs (13), audits finite controls at thresholds
`0,1,2,infinity`, verifies terminating/all-one hostile boundaries, checks
the exact Haar zero-block probabilities, enumerates the full finite
interaction simplex, and replays the radix-`13` zero-carry lift. No
load-bearing check uses `assert`.

Reproduce with

```bash
python3 04-computation/q_adic_prefix_collision_spectrum_thm2352.py
python3 -O 04-computation/q_adic_prefix_collision_spectrum_thm2352.py
```

Both transcripts must match

```text
05-knowledge/results/q_adic_prefix_collision_spectrum_thm2352.out
```

byte-for-byte after LF normalization.

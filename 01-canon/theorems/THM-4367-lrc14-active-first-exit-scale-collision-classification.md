---
id: THM-4367
title: "LRC14 active first-exit scale collision classification"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4365 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED; LRC(14) OPEN. Equal strict-active metric exits in THM-4365's common
  828-quotient tail are exactly primitive-ray classes with an arithmetic
  progression of scale addresses. Their multiplicity is at most 241, with
  all equality classes explicit. The metric quotient is not a state for the
  shift P->P+2 and does not determine the physical binder. No safety
  decrement or LRC(14) conclusion is asserted.
source: root + exit-collision scout + hostile referee + clean-room referee / next-sharp session, 2026-09-03
depends_on:
  - THM-4365-lrc14-cofinite-828-quotient-fibre-and-centered-residue-exit-law
related:
  - THM-4363-lrc14-828-body-completed-chain-first-exit-nonfactorization
mistake_firewall:
  - MISTAKE-300
primary_script: 04-computation/lrc14_active_first_exit_scale_collision_thm4367.py
primary_output: 05-knowledge/results/lrc14_active_first_exit_scale_collision_thm4367.out
primary_script_sha256: 03f8b6f1d9f7930b14b8c30a58d41fb532deae20f42f560174bda1f0659a64d6
primary_output_sha256: 97a32cb62b46f4f3bd1fa25078aedcfa078bbe504259140e88745050d1d1c33c
independent_referee_script: 04-computation/lrc14_active_first_exit_scale_collision_independent_referee_thm4367.py
independent_referee_output: 05-knowledge/results/lrc14_active_first_exit_scale_collision_independent_referee_thm4367.out
independent_referee_script_sha256: 7651b566b8749326ed5c56e64e8ad544b94c6daa357115d2c26fc7ae4e2a27dc
independent_referee_output_sha256: 87e53af1f98245c8ff4a7d48be6e77358b649205e1a600fea93559f60885c956
hash_basis: raw LF bytes
audit: >
  PASS. The 279,091-check primary and 478,870-check import-free referee
  independently verify the normal form, converse, reduced metric quotient,
  exact scale fibres, both leastness notions, all 241 equality cases,
  boundary controls, shift hostile, and physical sidecar. Normal, optimized,
  isolated, hash-seeded, and frozen LF streams agree.
---

# THM-4367 -- LRC14 active first-exit scale collision classification

**PROVED ELEMENTARY RELATIVE TO THM-4365 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THIS CLASSIFIES A METRIC CONSUMER ON ONE EXPLICIT COFINITE FIBRE.
IT DOES NOT DECREMENT THE OPEN LRC(14) OBLIGATION.**

## 1. Statement

Use THM-4365's fixed `h=420` family and restrict to its common 828-status and
completed-chain fibre

```text
P odd,                     P>=11019.                     (1)
```

Put

```text
A=3371,        S=1303,        M=14A=47194,
SP=Mn+rho,                  -M/2<rho<=M/2.               (2)
```

A cell is strict active when `|rho|<A`. THM-4365 then gives its first exit

```text
E_x=S/M+(A-rho)/(MP).                                   (3)
```

The theorem gives the following exact classification.

For every strict-active cell, set

```text
c=A-rho,        g=gcd(c,P),        c=ag,        P=bg.    (4)
```

Then `a` is positive even, `b,g` are positive odd, `gcd(a,b)=1`, and there is
a positive integer `kappa` satisfying

```text
a+Sb=A kappa,                 g kappa=1 (mod 14).        (5)
```

Conversely, every positive integer quadruple `(a,b,kappa,g)` satisfying

```text
a even, b odd, gcd(a,b)=1,      a+1303b=3371kappa,
g kappa=1 (mod 14),             bg>=11019,
2<=ag<=6740                                                   (6)
```

gives exactly one strict-active cell through

```text
P=bg,             rho=3371-ag,             n=(gkappa-1)/14. (7)
```

Its metric exit is the already reduced fraction

```text
E_x=kappa/(14b).                                         (8)
```

Therefore two strict-active cells have the same metric exit if and only if
their pairs `(kappa,b)` agree, equivalently if and only if their reduced
fractions `c/P=a/b` agree.

For a fixed admissible `(a,b,kappa)`, let `g_0` be the unique member of
`{1,3,5,9,11,13}` with `g_0 kappa=1 (mod 14)`. All and only its scales are

```text
g=g_0+14j,                         j in Z,
ceil(11019/b)<=g<=floor(6740/a).                         (9)
```

Every metric collision class has at most 241 cells, and this bound is sharp.
All equality classes are

```text
a=2,       b=401+6742t,       kappa=155+2606t,
t>=1,                         t mod 7 in {0,1,2,5}.      (10)
```

The infinitely many integers `t>=0` with `t!=3 (mod 7)` in `(14)` give
pairwise distinct nontrivial collision classes; the subfamilies in `(10)`
are exactly the 241-cell classes.

## 2. Derivation of the normal form

The integer `3371` is prime. From `(2)` and `c=A-rho`,

```text
SP+c=A(14n+1).                                          (11)
```

Substitute `(4)` into `(11)`. The prime `A` cannot divide `g`: if it did,
then the positive even integer `a` would give `c=ag>=2A`, whereas strict
activity gives `c<=2A-2`. Hence `A` divides `a+Sb`; defining the quotient to
be `kappa` gives `(5)` and then

```text
gkappa=14n+1.                                           (12)
```

Conversely, `(6)`--`(7)` reconstruct `(11)`. Its bounds put
`rho=A-ag` strictly in `(-A,A)`, while the congruence makes `n` integral.
This proves both directions and also proves uniqueness, because `(4)` is the
canonical gcd reduction of `c/P`.

Substitution in `(3)` now gives

```text
E_x=S/(14A)+a/(14Ab)=(Sb+a)/(14Ab)=kappa/(14b).          (13)
```

This fraction is reduced. Indeed `A` cannot divide `b`, since then `(5)` and
`gcd(a,b)=1` would force `A|a`. Also

```text
gcd(Akappa,b)=gcd(a+Sb,b)=1,
```

so `gcd(kappa,b)=1`; and `(5)` gives `gcd(kappa,14)=1` through the scale
congruence. Equality of the reduced positive fractions in `(8)` is therefore
exactly equality of `(kappa,b)`. Equation `(9)` is just `(6)` plus the unique
unit inverse of `kappa` modulo 14.

## 3. The sharp multiplicity and all equality cases

Because `a>=2`, every scale in `(9)` has `g<=3370`. A single unit residue
class modulo 14 contains at most 241 positive integers in that interval.
Thus every metric fibre has size at most 241.

If equality holds, then `a` must be two: for `a>=4`, the upper scale bound is
at most 1685 and no residue class can have 241 elements. Solving `(5)` at
`a=2`, then imposing odd `b`, gives exactly

```text
b=401+6742t,             kappa=155+2606t,       t>=0.   (14)
```

Modulo 14, `kappa=1+2t`. The class is inadmissible for `t=3 (mod 7)`. At
`t=0`, the tail lower bound removes two scales and leaves 239. For `t>=1`,
the six remaining residue cases have respectively

```text
t mod 7:       0   1   2   4   5   6
number:      241 241 241 240 241 240.                   (15)
```

This proves `(10)`. For example, `t=1` gives

```text
(a,b,kappa)=(2,7143,2761),
g=5,19,...,3365                  (241 values),
P=35715,...,24036195,             E_x=2761/100002.       (16)
```

The corrections `2/b` in `(14)` are pairwise distinct, so the admissible
values of `t` also prove the infinitude claim.

## 4. Exact first collisions and hostile controls

The least parameter belonging to any nontrivial metric fibre lies in

```text
(a,b,kappa)=(196,2217,857),
g=(5,19,33),
P=(11085,42123,73161),
rho=(2391,-353,-3097),
n=(306,1163,2020),
E_x=857/31038.                                           (17)
```

The first repeat encountered while scanning `P` upward is `P=12675`, which
repeats the exit at `P=11625`. Its complete fibre is

```text
(a,b,kappa)=(34,75,29),
g=(155,169,183,197),
P=(11625,12675,13725,14775),
rho=(-1899,-2375,-2851,-3327),
n=(321,350,379,408),
E_x=29/1050.                                             (18)
```

This metric quotient is static, not a process state for the natural odd-tail
translation `T(P)=P+2`. The four points in `(18)` have one common exit, while
their successors

```text
P=(11627,12677,13727,14777),
rho=(707,231,-245,-721)
```

have the four distinct exits

```text
4495/162778, 4901/177478, 5307/192178, 5713/206878.       (19)
```

Thus equality of `(kappa,b)` is not an operation congruence under `T`.

These two minimality statements are FINITE-EXACT. The primary enumerates all
odd `P` from 11019 through 174375. This is decisive: if `P_0=bg` is the least
member of a nontrivial class, its next member is `b(g+14)`, at most `15P_0`.
Thus every class whose least member could improve `(17)` appears far enough
to expose its second member. The same enumeration directly proves
the earliest-repeat claim in `(18)`.

Residue phase alone is not the metric quotient. The pair

```text
P=(48679,95873),                rho=(1,1)                (20)
```

has exits `18817/681506` and `37059/1342222`, respectively.

## 5. Metric quotient versus physical binder

Inside one metric fibre, `(kappa,b)` deliberately forgets `g`. The physical
binding tooth is

```text
P@n = bg @ ((gkappa-1)/14).                              (21)
```

Consequently `(kappa,b)` is the exact quotient for the metric consumer, but
not for THM-4363's full physical exit record. Within a fixed metric fibre,
retaining the physical binder requires distinguishing every scale. One
injective scale coordinate is necessary and sufficient, canonically `g`
(equivalently `P` or `n`).
For example, the four cells in `(18)` bind at

```text
11625@321, 12675@350, 13725@379, 14775@408.              (22)
```

The nonactive endpoint needs a separate flag. For `|rho|>3371`, the exit is
`1303/47194` with binder `{3371}`. For `|rho|=3371`, the same metric exit has
binder set `{3371,P}`. The least representatives of the two boundary phases
in the cofinite tail `P>=11019` are

```text
(P,rho,n)=(43823,-3371,1210),       (50565,3371,1396).   (23)
```

The teeth are open, so equality remains safe.

This supplies the precise natural-number lesson. An injective encoding of
`(kappa,b)` is a lossless ordinal name for the metric exit. Giving every
scale in `(9)` the same name is a valid quotient only for that consumer;
the physical binder still requires the scale address.

## 6. Audit and scope

The primary performs 279,091 exact checks. It enumerates all 11,665 active
cells in the minimality universe, matches exact rational fibres with `(9)`,
checks `(10)`--`(23)`, and records 114 bounded collision classes. The
478,870-check import-free referee reconstructs the arithmetic in both forward
and reverse directions over three periods and a separate primitive-pair box.
Normal, optimized, isolated, hash-seeded, and frozen LF streams agree.

The theorem uses THM-4365's already proved geometry rather than rerunning
it. It classifies one consumer on one explicit cofinite quotient fibre. All
rows remain safe controls. It does not improve a safe time, force a missing
component, enter the unresolved seam, decrement any LRC ledger, or prove
LRC(14).

Reproduce from the repository root:

```text
python3 -B 04-computation/lrc14_active_first_exit_scale_collision_thm4367.py
python3 -B -O 04-computation/lrc14_active_first_exit_scale_collision_thm4367.py
python3 -B 04-computation/lrc14_active_first_exit_scale_collision_independent_referee_thm4367.py
python3 -B -O 04-computation/lrc14_active_first_exit_scale_collision_independent_referee_thm4367.py
```

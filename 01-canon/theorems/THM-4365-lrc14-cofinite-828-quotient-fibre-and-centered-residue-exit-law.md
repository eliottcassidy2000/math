---
id: THM-4365
title: "LRC14 cofinite 828-quotient fibre and centered-residue exit law"
status: >
  PROVED FINITE-EXACT PARAMETRIC THEOREM + INDEPENDENTLY AUDITED; LRC(14)
  OPEN. In one explicit primitive h=420 family, every odd P>=11019 has the
  same 828-component status/completed-chain quotient, and P=11017 proves
  this eventual odd threshold sharp. From the earlier sharp threshold
  P>=10141, the first exit is governed by one centered residue modulo 47194;
  P=10139 is the hostile predecessor. The common cofinite quotient fibre has
  3,370 active odd residue cells and infinitely many distinct metric exits.
  Every row is an explicitly safe control. No counterexample family, ledger
  decrement, or arbitrary-bank quotient theorem is asserted.
source: root + cofinite-fibre scout + clean-room referee / next-sharp session, 2026-09-02
depends_on:
  - THM-4363-lrc14-828-body-completed-chain-first-exit-nonfactorization
related:
  - THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal
  - THM-4345-lrc14-halfperiodic-anchor-strip-euclidean-remainder-and-current-envelope
  - THM-4348-lrc14-prefix-envelope-third-tooth-and-nested-wall-shadow
  - THM-2050-period14-top-germs-do-not-determine-global-loneliness
mistake_firewall:
  - MISTAKE-300
primary_script: 04-computation/lrc14_cofinite_828_quotient_fibre_centered_residue_exit_thm4365.py
primary_output: 05-knowledge/results/lrc14_cofinite_828_quotient_fibre_centered_residue_exit_thm4365.out
primary_script_sha256: a7ea9bb5828ba0e3169a066f4d7832443c86af2db3476e735b9eacc28c754bd5
primary_output_sha256: d404f5240a38d3b366a249f51e1942521f62a9a0c4faab117a93fddd769367cf
independent_referee_script: 04-computation/lrc14_cofinite_828_quotient_fibre_centered_residue_exit_independent_referee_thm4365.py
independent_referee_output: 05-knowledge/results/lrc14_cofinite_828_quotient_fibre_centered_residue_exit_independent_referee_thm4365.out
independent_referee_script_sha256: 33d9bb2976ad3f6ee2900c3b80c6c626894941bf0b69a206bbe82d177a2e67ed
independent_referee_output_sha256: 1115c92f2b15a3d49854c8581de5240fb140ca2b9d7fde1bb7446e72b03da2d2
hash_basis: raw LF bytes
audit: >
  PASS WITH THRESHOLD AND ENDPOINT REPAIRS. The 2,772,267-check primary and
  an import-free interval referee collectively verify the fixed quotient,
  extremal reach/gap constants, finite bridges, hostile predecessors,
  centered-residue law, boundary ties, and infinite family. The referee also
  matches THM-4363's serialization hashes and audits the robust 307-row exit
  strip and its minimum margin. Normal/-O/hash-seeded/frozen LF streams match.
---

# THM-4365 -- LRC14 cofinite 828-quotient fibre and centered-residue exit law

**PROVED FINITE-EXACT PARAMETRIC THEOREM + INDEPENDENTLY AUDITED. LRC(14)
REMAINS OPEN. EVERY ROW BELOW HAS AN EXPLICIT SAFE POINT; THIS IS NOT A
COUNTEREXAMPLE FAMILY AND LOWERS NO LRC LEDGER.**

## 1. Family, quotient, and exact statement

Keep THM-4363's open danger teeth, farthest-right selection rule, anchor
components, and collar. Put

```text
h=420,                         anchor=840,
B=(3,39,11,1691,3371,5051,6731,8411,10091,525,945),
W_P=B union {P},
C={19,20,259,260,299,300,539,540,579,580,819,820}.       (1)
```

For odd `P>10091`, the thirteen relative speeds `{840} union W_P` are
distinct and primitive. Let `Q(P)` retain `C`, all 828 labelled residual
statuses, and every complete physical `(speed,tooth-address)` chain on a
`span` or `renew` component, while erasing partial chains on `missing`
components. Let `E(P)` be THM-4363's first missing-component exit record.

The theorem has four parts:

```text
(i)   Q(P)=Q(241) for every odd P>=11019;                (2)
(ii)  P=11017 is outside that fibre, so 11019 is the
      least eventual odd threshold;
(iii) the centered-residue formula (17) below holds for
      every odd P>=10141, and fails at P=10139;
(iv)  on the single fibre (2), E has infinite range.      (3)
```

The inheritance pass was:

- closest proved mechanism: THM-4335/4348's physical-address renewal walk
  and THM-4363's four-point quotient hostile;
- canonical hostile: a status word plus all completed chains can still erase
  the partial trace read by the next metric consumer;
- corrected near miss: the phase-free width threshold is only a sufficient
  tail, not the actual onset;
- least-used sidecar: the largest connected gap in the union of all fixed
  teeth on each missing component.

The live board was

```text
selected reach | union gap | added-tooth diameter | centered residue
Euclidean address | open endpoint | first exit | quotient fibre.           (4)
```

## 2. Fixed-bank protection of completed chains

Remove `P` and run the exact renewal on the 828 residual components. The
fixed-bank census is

```text
missing=546,                 span=276,
renew=6,                     completed=282.              (5)
```

Its status and completed-chain serializations are exactly THM-4363's common
ones: 4,860 bytes with SHA-256
`4d34447b9eca8c8a9302a0f799a56300b4b96135cd0c3a245a97960975f9a347`,
and 10,940 bytes with SHA-256
`64704c712a0e3e70ca0ecb3264834b3f61c4d473ddc20ea8c44ccc5c2c616d11`.

At every selected frontier in every completed fixed trace, the chosen tooth
advances by at least

```text
Delta=1/28665.                                            (6)
```

Equality occurs only at the second selected step (zero-based step index one)
on components `216,496,776`, for the teeth `945@244,945@559,945@874`. An
active `P` tooth advances strictly less than its full diameter `1/(7P)`.
Hence for `P>=4095` it cannot tie or beat the fixed winner. Induction over
selected frontiers proves that every completed physical chain remains
literally unchanged.

This is stronger than preservation of union coverage: it protects the exact
chain retained by `Q`.

## 3. Missing gaps and the sharp eventual odd onset

For each fixed-missing component, merge all clipped fixed open teeth and take
the largest connected complementary gap. The minimum of those 546 maxima is

```text
G=7003/594127807,                                        (7)
```

attained exactly at the reflected components

```text
k=410: [57597/117754,69103/141274],
k=429: [72171/141274,60157/117754].                       (8)
```

Distinct teeth of speed `P` are disjoint. To cross one connected fixed-free
gap, a single open `P` tooth must strictly contain it, so its diameter must
exceed the gap length. But

```text
1/(7G)=84875401/7003,              12119<1/(7G)<12120.   (9)
```

Thus `P>=12120`, and in particular every odd `P>=12121`, preserves all 546
missing statuses. Together with `(6)`, this is a phase-independent sufficient
tail for `(2)`. It is sharp for that diameter argument, not for the actual
parameter family.

The exact finite bridge checks all 551 odd values `11019<=P<=12119`. For each
parameter and fixed-missing component it exhibits a positive fixed gap not
strictly contained by any `P` tooth. Hence all these statuses also remain
missing, while `(6)` protects the completed chains. This proves `(2)`.

At the immediate odd predecessor `P=11017`, exactly two components change
from `missing` to `renew`:

```text
k=147:
1691@296 -> 11017@1929 -> 10091@1767 -> 525@92 -> 11@2;

k=692:
11@9 -> 525@433 -> 10091@8324 -> 11017@9088 -> 1691@1395. (10)
```

The row is still distinct and primitive. Therefore 11019 is the sharp
consecutive-odd onset, without claiming a closed-form classification of all
smaller fibre points.

## 4. The first exit as a centered-residue law

The fixed trace on the least missing component is

```text
k=23,       945@26 -> 3371@93 -> EXIT,
x1=1303/47194.                                           (11)
```

The first two selected reaches are

```text
13/105840,                         92/4459833.            (12)
```

Both exceed `1/(7P)` for `P>=10141`, so the varying tooth cannot alter this
prefix before `x1`. For odd `P`, define the unique centered remainder and
Euclidean address by

```text
1303P=47194n+rho,                  -23597<rho<=23597.     (13)
```

At `x1`, the address-`n` tooth is active exactly when

```text
|P*x1-n|<1/14  iff  |rho|<3371.                          (14)
```

If active, its right wall is

```text
b(P,n)=x1+(3371-rho)/(47194P).                           (15)
```

The next fixed left wall is `left(6731@186)=2603/94234`, at distance

```text
g=2110/158831407                                         (16)
```

from `x1`. Since `1/(7g)=22690201/2110` lies strictly between 10753 and
10754, `(15)` cannot reach the next fixed tooth for odd `P>=10755`. An exact
audit of the remaining odd strip `10141<=P<=10753` gives the same conclusion;
its least positive margin is `19/493079405`, at
`(P,rho,n)=(10465,-3171,289)`.

Consequently the sharp eventual formula is

```text
E_x(P)=x1,                                      if |rho|>=3371,
E_x(P)=x1+(3371-rho)/(47194P),                  if |rho|<3371, (17)
```

for every odd `P>=10141`. At `P=10139`, `(rho,n)=(-3203,280)` and the
predicted `P` wall crosses the next fixed start by `31/68245609`. The actual
chain continues

```text
945@26 -> 3371@93 -> 10139@280 -> 6731@186 -> 10091@279
```

and exits at `3907/141274`. Thus the threshold 10141 is sharp in the same
consecutive-odd sense.

## 5. Exactly 3,370 active cells and infinitely many exits

The factorization and inverse

```text
47194=14*3371,                  gcd(1303,47194)=1,
1303^(-1)=1485 mod 47194                              (18)
```

show that multiplication by 1303 permutes the 23,597 odd residue classes.
Their centered remainders are odd, and `(14)` selects exactly

```text
rho in {+/-1,+/-3,...,+/-3369}:       2*1685=3370       (19)
```

active cells. The other 20,227 are inactive. At the boundary cells
`rho=+/-3371`, the exit is `x1`, both speeds `3371` and `P` bind at distance
exactly `1/14`, and both open teeth are inactive. Thus equality is safe but
the binder is not unique. In strict active cells `P` is the unique binder;
in strict inactive cells `3371` is.

One explicit active class is

```text
P_j=1485+47194j,             n_j=41+1303j,             j>=1. (20)
```

It has `rho=1` and

```text
E_x(P_j)=x1+3370/(47194P_j).                             (21)
```

These values decrease strictly to `x1`, so they are pairwise distinct. Since
every `P_j>=48679` lies in `(2)`, one fixed quotient fibre has infinitely
many metric exits.

This is the useful natural-number structure: finite state is periodic modulo
47194, the tooth address advances by 1303 along one cell, and metric position
is affine in `1/P`. Replacing the physical address by an arbitrary ordinal
would retain the count but destroy `(15)`.

## 6. Safety, audit, and scope

Every odd row in `(2)` is safe at its component-23 exit: the binding distance
is exactly `1/14`, and the binding open tooth or teeth are inactive at their
endpoint. Thus the theorem supplies controls, not LRC counterexamples.

The primary checks the fixed census, extremal constants, 551-row bridge, both
hostile predecessors, a 305-row exit strip followed by a sharper residue-
aware tail, the 3,370-cell count, explicit infinite family, and the earlier
bounded 51-signature scan, in 2,772,267 exact checks. The import-free referee
rebuilds the intervals, matches THM-4363's two serialization hashes, and uses
the robust 307-row strip plus phase-free tail stated in Section 4; it also
checks the minimum strip margin and exact boundary binder sets. Normal,
optimized, hash-seeded, and frozen LF streams agree.

The proved universe is only the family `(1)`. Nothing here classifies an
arbitrary bank or height, restricts a hypothetical minimal counterexample,
proves `Q` optimal, lowers a residual rank, or advances the logical LRC(14)
bound.

Reproduce from the repository root:

```text
python3 -B 04-computation/lrc14_cofinite_828_quotient_fibre_centered_residue_exit_thm4365.py
python3 -B -O 04-computation/lrc14_cofinite_828_quotient_fibre_centered_residue_exit_thm4365.py
python3 -B 04-computation/lrc14_cofinite_828_quotient_fibre_centered_residue_exit_independent_referee_thm4365.py
python3 -B -O 04-computation/lrc14_cofinite_828_quotient_fibre_centered_residue_exit_independent_referee_thm4365.py
```

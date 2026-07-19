---
id: THM-1215
title: THE ALL-RANGE q=14 BEAT-PUNCTURE STALK CLOSURE — a difference beat q=b6-b5=14d with five common gcd-d fast combs leaves a lonely beat point in every slow gap
status: PROVED (elementary all-range residue argument; THM-1192 phase-free contradiction; exact referee replay; Lean integer core)
source: codex-2026-07-18-S80
depends_on: [THM-1192]
related: [THM-1176, THM-1178, THM-1179, THM-1182, THM-1193]
script: 04-computation/lrc14_q14_beat_puncture_stalk_referee_codex_S80.py
output: 05-knowledge/results/lrc14_q14_beat_puncture_stalk_referee_codex_S80.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCBeatPunctureQ14.lean
script_sha256: f712a4781bcb951609f4d6857d935f1c8c250444d4c99bd9eb893db529542ee8
output_sha256: bf1a2525bef60c44c9b86b8f8221c58ebf8b8b915661df1980489dab4689768c
formalization_sha256: 49c2273b1876aa8f964709af5f32663e293dac1cbfe935a5831f5da68c277cfa
---

# THM-1215 — the all-range q=14 beat-puncture stalk closes

## Statement

Let

```text
a < b1 < b2 < b3 < b4 < b5 < b6
```

be positive integers.  Put

```text
q = b6-b5 = 14d,
```

and suppose

```text
q >= 7a,                         (1)
gcd(c,q) = d  for c=b1,...,b5.  (2)
```

For every closed slow-carrier safe gap

```text
G_k(a)=[(14k+1)/(14a),(14k+13)/(14a)],   0<=k<a,       (3)
```

there is an integer `p` such that

```text
t=p/q lies in G_k(a)
and
||a t||, ||b1 t||, ..., ||b6 t|| >= 1/14.             (4)
```

Thus none of these slow gaps can be covered by the six faster dangerous
combs.  This closes the stated `q=14` mixed-gcd stalk at **every** scale; it
does not close the general slow-gap branch or the Lonely Runner Conjecture
for fourteen runners.

## Direct residue proof

The interval of eligible beat numerators is the consecutive block

```text
P_q(k)={p in Z:
  ceil(q(14k+1)/(14a)) <= p <= floor(q(14k+13)/(14a))}. (5)
```

Its underlying real interval has length

```text
q(14k+13)/(14a)-q(14k+1)/(14a)
  = 6q/(7a) >= 6                                                   (6)
```

by (1).  If `p0` is the ceiling of its left endpoint, then
`p0,p0+1,...,p0+5` all belong to (5).  At most one of these six consecutive
integers is divisible by `14`, so choose one with

```text
p != 0 (mod 14).                                           (7)
```

For `c=b1,...,b5`, write `c=d r_c`.  Equation (2) and `q=14d` give

```text
gcd(r_c,14)=1.                                             (8)
```

At `t=p/q`, the strict-danger condition for `c` is

```text
||c p/q|| < 1/14
iff ||r_c p/14|| < 1/14
iff r_c p = 0 (mod 14)
iff p = 0 (mod 14).                                       (9)
```

The middle equivalence is exact because the least circular residue is an
integer and must be strictly below `1`.  Condition (7) therefore makes all
of `b1,...,b5` safe.  Also

```text
gcd(b6,q)=gcd(b5+q,q)=gcd(b5,q)=d,                        (10)
```

so the same argument makes `b6` safe.  Finally (5) puts `t` in the closed
carrier-safe gap (3), hence `a` is safe as well.  This proves (4).

The direct proof exposes more than a density obstruction: after the gcd
quotient, all six fast combs have the identical strict dangerous mask
`{0}` on `Z/14Z`.  The witness is any beat numerator in (5) outside that
single residue.

## THM-1192 phase-free proof

The same closure is already forced by THM-1192's weaker phase-free necessary
law.  Apply that theorem to the defining fast pair `(b5,b6)` and its
difference beat `q`.  Write

```text
N=#P_q(k),
Q_c=q/gcd(c,q).
```

Equations (2) and (10) give

```text
Q_c=14  for c=b1,...,b6.                                  (11)
```

THM-1192 has

```text
A(Q)=2 ceil(Q/14)-1,
U(N,Q)=floor(N/Q) A(Q)+min(N mod Q,A(Q)).                  (12)
```

Consequently

```text
A(14)=1,
U(N,14)=floor(N/14)+min(N mod 14,1)=ceil(N/14).            (13)
```

Deleting the defining pair leaves the four complementary speeds
`b1,...,b4`.  A hypothetical cover of (3) would therefore have to satisfy

```text
N-ceil(N/14) <= 4 ceil(N/14),
or equivalently
N <= 5 ceil(N/14).                                        (14)
```

But (6) gives `N>=6`, and

```text
N > 5 ceil(N/14)  for every N>=6.                         (15)
```

For completeness, write `N=14m+r`, `0<=r<14`.  If `r=0`, then `m>=1`
and the difference in (15) is `9m>0`.  If `r>0,m=0`, then `r=N>=6` and
the difference is `r-5>0`.  If `r>0,m>=1`, the difference is
`9m+r-5>=5`.  This contradicts (14).

This second proof is intentionally retained because it identifies a reusable
coarse certificate: five q=14 period caps cannot absorb even the shortest
slow-gap beat block.  The direct residue proof shows why the certificate has
so much slack—all five masks coincide rather than merely obeying a union
bound.

## Exact example

For

```text
(a;b1,...,b6)=(3;6,10,18,22,26,54),
q=28=14*2,
```

the five gcds are all `2` and `q>=7a`.  The three gap blocks and first
nonzero-mod-14 witnesses are

| `k` | `P_q(k)` | `N` | `U=ceil(N/14)` | `N-U > 4U` | witness `p` |
|---:|:---:|---:|---:|:---:|---:|
| 0 | `1..8` | 8 | 1 | `7>4` | 1 |
| 1 | `10..18` | 9 | 1 | `8>4` | 10 |
| 2 | `20..27` | 8 | 1 | `7>4` | 20 |

Thus `p/28` is a literal lonely point in each slow gap, not merely a failed
capacity inequality.

## Structural and tournament audit

The useful quotient is not a tournament on runners.  It is the reduced
fourteen-section clock with the consecutive block (5) retained as a phase
sidecar.  This quotient preserves the exact strict-danger predicate at every
beat numerator.  If the block sidecar is discarded, it loses which residues
the actual slow gap visits; if the unit multipliers are discarded in a
general period, it would also lose the comb phases.  The exceptional feature
of this stalk is precisely that every unit multiplier has the same singleton
danger mask at radius `1/14`.

For the required Tournament Analysis audit, take equality of reduced danger
masks as the pairwise observable on the six fast speeds.  Every pair ties.
Breaking ties by increasing speed is a gauge, with tie Hamiltonian path

```text
b1 -> b2 -> b3 -> b4 -> b5 -> b6.
```

The resulting tournament is transitive: score histogram
`{0,1,2,3,4,5}`, zero directed cycles, six singleton SCCs, one Hamiltonian
path, and all `15` edges flip under the reverse gauge.  None of these
fingerprints proves (4); they only record that a binary orientation has no
remaining information once the masks coincide.

Vertices were challenged as runners, slow gaps, fixed circle sections,
section boundaries, wall-crossing events, beat numerators, residue masks, and
proof obligations.  Residue sections plus the block-phase sidecar are the
faithful choice.  A runner tournament preserves speed order but destroys the
common singleton mask and the location of the numerator block.

## Formalization boundary

`LRCBeatPunctureQ14.lean` kernel-checks

```text
U(N,14)=ceil(N/14),
N>=6 -> 5 ceil(N/14)<N,
not (N-U(N,14)<=4U(N,14)),
```

and a typed consumer whose explicit hypotheses are the block-size supplier,
the defining/complementary q=14 cap suppliers, and THM-1192's necessary law.
The module deliberately leaves the real-gap-to-integer-block translation and
the speed-gcd reduction as named external bridge assumptions; no geometric
claim is hidden in the arithmetic core.

## Exact replay

```text
python3 04-computation/lrc14_q14_beat_puncture_stalk_referee_codex_S80.py
python3 -O 04-computation/lrc14_q14_beat_puncture_stalk_referee_codex_S80.py
```

The normal and optimized outputs are byte-identical.  Artifact hashes are
recorded in the front matter after the frozen replay.

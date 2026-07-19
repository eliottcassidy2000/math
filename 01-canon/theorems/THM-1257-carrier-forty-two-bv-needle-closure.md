---
id: THM-1257
title: The bounded-variation needle certificate closes carrier 42 exactly
status: PROVED FINITE-EXACT for carrier c=42, every slow-gap phase, and every faster integer speed.  Twenty-one reflection-representative 200-bin probability densities have exact closed-comb load below 1/6 on a finite vertical ladder and an exact BV-tail cap below 1/6 thereafter.  Consequently no six faster combs cover a complete 42-gap, even with repeated speeds.  Together with THM-1182/1197/1255 this closes every carrier 1<=c<=42; carriers c>=43 remain open
source: codex-2026-07-19 carrier-42 exact census
depends_on: [THM-1182, THM-1197, THM-1255]
related: [THM-1176, THM-1198, THM-1232, THM-1250, THM-1254]
script: 04-computation/lrc14_slow_gap_bv_needle_carrier42_thm1257.py
independent_verifier: 04-computation/lrc14_slow_gap_bv_needle_carrier42_thm1257_verify.py
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCCarrier42BVNeedle.lean
certificate: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier42_thm1257_certificate.json
output: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier42_thm1257.out
independent_output: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier42_thm1257_verify.out
---

# THM-1257 -- carrier 42 has no six-comb slow-gap cover

For an integer speed `d`, let

```text
Dbar_d={t in R: ||dt||<=1/14}
```

and put

```text
G_k=[(14k+1)/(14*42),(14k+13)/(14*42)].              (1)
```

> **Theorem.**  For every integer `k` there is an absolutely continuous
> probability measure `mu_k`, supported on `G_k`, such that
>
> ```text
> mu_k(Dbar_d)<1/6                                    (2)
> ```
>
> for every integer `d>42`.  Consequently no six faster danger combs cover
> any complete `42`-gap.

Repeated speeds are allowed, so this also excludes every distinct or
primitive six-comb packet at carrier `42`.

## 1. Exact finite and tail certificate

For each `0<=k<=20`, divide `G_k` into `200` equal bins.  The certificate
stores a nonnegative integer mass on every bin, normalized to a probability
density.  At a row-specific cutoff `B_k=42M_k`, it verifies exactly

```text
max_(43<=d<=B_k) mu_k(Dbar_d)<1/6,                    (3)

1/7+3 TV_R(mu_k)/(49(B_k+1))<1/6.                    (4)
```

The multiplier census is

```text
M          12  16  20  28  44
rows        7   1   1   9   3.                      (5)
```

No phase needs the longer `60`, `76`, or `108` discovery cutoffs.  The exact
worst finite and tail loads are

```text
28499599523095/170999999979651
 =1/6-800947027/341999999959302,                     (6)

65421208326679/392583333284653
 =1/6-56083324579/2355499999707918.                  (7)
```

Both displayed margins are strictly positive.  The LP is used only to
discover integer bin masses; ordinary and optimized replays consume the
stored integers with exact rational arithmetic.

THM-1182 supplies the bounded-variation inequality

```text
|mu_k(Dbar_d)-1/7|<=3 TV_R(mu_k)/(49d).               (8)
```

For `d>B_k`, equation (4), monotonicity of `1/d`, and (8) prove (2).  The
finite range is (3), so every faster speed is covered without using a
projective ratio truncation.

The independent verifier does not import the LP generator or its clipped
tooth integrator.  It recomputes every finite load from

```text
F(y)=floor(y)/7+min(frac(y),1/7),

|[u,v] intersect Dbar_d|
 =[F(dv+1/14)-F(du+1/14)]/d,                         (9)
```

and independently reproduces (6)--(7).

## 2. The exact reflection count is 21

Translation reduces `k` modulo `42`.  Reflection `rho(t)=1-t` gives

```text
rho(G_k)=G_(41-k),                                   (10)
```

because

```text
1-(14k+13)/(14*42)=(14(41-k)+1)/(14*42),
1-(14k+1 )/(14*42)=(14(41-k)+13)/(14*42).            (11)
```

Also `||d(1-t)||=||dt||`, so every closed danger comb is reflection
invariant.  The involution `k -> 41-k` has no integral fixed point: a fixed
point would satisfy `2k=41`.  Thus the `42` phase residues form exactly `21`
two-element orbits, represented by

```text
k=0,1,...,20.                                        (12)
```

There are not `22` reflection representatives; `k=21` is the reflected mate
of `k=20`.  The endpoint identities, absence of a fixed phase, and the fact
that every `0<=k<=41` is represented by (12) are all checked in the sorry-free
Lean consumer.

Finally, a six-comb cover would give

```text
1=mu_k(G_k)
 <=mu_k(union_(i=1)^6 Dbar_(d_i))
 <=sum_(i=1)^6 mu_k(Dbar_(d_i))<1,                   (13)
```

a contradiction.  Closed teeth only strengthen the strict-open lonely
runner formulation.

## 3. Structural fingerprint

Write `d=42q+s`, with `0<=s<42`.  The largest finite load occurs at

```text
(k,d,q,s,M)=(18,171,4,3,28).                         (14)
```

Carrier `41` also had its hardest finite obstruction in residue `s=3`, but
at the first vertical rung.  At carrier `42` the same residue persists while
the sharp rung moves to `q=4`.  This is useful evidence that the residue and
vertical-frequency coordinates must both remain in the phase stalk.

The worst tail is the `k=16,M=16` row.  Its unit-gap zero-extension
variation is

```text
444660714300/83333333323.                            (15)
```

The three longest rows are precisely `k=5,6,7` at `M=44`; there is no long
isolated `M=60+` phase.  Thus the carrier-42 difficulty is a narrow block of
phase-orbit obligations, not a remote-speed tail.

In normalized coordinates `t=(k+u)/42` and `d=42q+s`, the exact orbit law is

```text
dt=(q+s/42)u+sk/42          modulo 1.                (16)
```

THM-1255's failed literal `41 -> 42` density transport and the successful
new certificate show what toothpick self-similarity means here: the
selection operation recurs on a changing orbit, while the literal density
vector does not.

## 4. Tournament and alternate-carrier audit

We challenged runners, teeth, fixed sections, tooth walls, residues, phase
gaps, bins, and proof obligations as vertices.  A runner quotient is too
small: one density row controls infinitely many faster runners.  The faithful
objects are the `21` phase-obligation stalks, each decorated by its full
`200`-bin density and finite/tail split.

For the pairwise loss audit, orient two stalks by their exact bottleneck
`max(finite load,tail load)`, breaking equality by `k`.  This gives a
transitive tournament with score histogram `0,1,...,20`, no directed cycles,
`21` singleton SCCs, and one Hamiltonian path.  It preserves row difficulty
only and destroys the density weights; hence it cannot prove noncoverage.

The assumption challenged by this result is that tournament vertices should
be runners or arcs.  Here they are proof obligations over cyclic phase
orbits, while the proof-bearing sidecar is a Kakeya-needle density stalk.

## 5. Reproducibility, formalization boundary, and frontier

Ordinary and optimized runs are byte-identical for both exact replays.  The
independent cumulative verifier passes all `21` rows.  Frozen hashes are

```text
generator_sha256   = 7f613cfb09514b222b3522528d5a4e97718e25afa1bfa05b15e50756e8a69daa
independent_sha256 = 181c9552848338e2425966cfd8bdc55b6fdd87094190bf350c0c3b40f3c0a942
certificate_sha256 = f24ef5d3e4db9163e1c3b36f5f948fba4893bf57a6d37986a3fe77c8daf641ee
output_sha256      = 9966d58a99eded54d3ccdac087b4a8a4f68049c3b477c77d8c72863fe73e07ce
ind_output_sha256  = 0d8819b4f1159536856e1b8394cf9f87d5a9b3c6cc0c9addcc8e5ccf573a7fdc
lean_sha256        = 9cbdb2bfe1ac413a647a14464196b0962fee1f7fe6a78bc4ae2b429aa092616a
certificate_bytes  = 51265
```

The Lean module checks both rational margins, both endpoint-reflection laws,
the fixed-point-free 21-orbit decomposition, the finite/tail common cap, and
the six-load contradiction.  JSON completeness, exact interval integration,
and the analytic BV estimate (8) remain explicit computation/paper providers.

With THM-1182, THM-1197, and THM-1255, every carrier through `42` is now
closed.  This does not prove uniform six-comb noncoverage or LRC(14).
Carrier `43` is the next bounded exact row, while the higher-leverage target
is a parametric phase-orbit density-selection theorem capable of consuming
THM-1254's bounded `c<=1171` branch without hundreds of thousands of
independent LP rows.

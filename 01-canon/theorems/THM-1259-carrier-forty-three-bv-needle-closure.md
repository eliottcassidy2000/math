---
id: THM-1259
title: The bounded-variation needle certificate closes carrier 43 exactly
status: PROVED FINITE-EXACT for carrier c=43, every slow-gap phase, and every faster integer speed.  Twenty-two reflection-representative 200-bin probability densities have exact closed-comb load below 1/6 on a finite vertical ladder and an exact BV-tail cap below 1/6 thereafter.  Consequently no six faster combs cover a complete 43-gap, even with repeated speeds.  Together with THM-1182/1197/1255/1257 this closes every carrier 1<=c<=43; carriers c>=44 remain open
source: codex-2026-07-19 carrier-43 exact census
depends_on: [THM-1182, THM-1197, THM-1255, THM-1257]
related: [THM-1176, THM-1198, THM-1232, THM-1250, THM-1254]
script: 04-computation/lrc14_slow_gap_bv_needle_carrier43_thm1259.py
independent_verifier: 04-computation/lrc14_slow_gap_bv_needle_carrier43_thm1259_verify.py
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCCarrier43BVNeedle.lean
certificate: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier43_thm1259_certificate.json
output: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier43_thm1259.out
independent_output: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier43_thm1259_verify.out
---

# THM-1259 -- carrier 43 has no six-comb slow-gap cover

Let

```text
Dbar_d={t in R: ||dt||<=1/14},
G_k=[(14k+1)/(14*43),(14k+13)/(14*43)].              (1)
```

> **Theorem.**  For every integer `k` there is an absolutely continuous
> probability measure `mu_k`, supported on `G_k`, such that
>
> ```text
> mu_k(Dbar_d)<1/6                                    (2)
> ```
>
> for every integer `d>43`.  Hence no six faster danger combs cover a
> complete `43`-gap.

The statement permits repeated speeds and therefore excludes all primitive
or distinct six-comb packets as special cases.

## 1. The 22 exact density rows

Each representative gap is divided into `200` equal bins and assigned a
stored nonnegative integer probability density.  For a row-specific
`B_k=43M_k`, exact arithmetic verifies

```text
max_(44<=d<=B_k) mu_k(Dbar_d)<1/6,                    (3)

1/7+3 TV_R(mu_k)/(49(B_k+1))<1/6.                    (4)
```

The cutoff census is

```text
M          12  16  20  28  44
rows        6   2   2   8   4.                      (5)
```

No phase resists the original multiplier bank or needs `M>=60`.  The largest
finite and tail loads are respectively

```text
2331708819707/13999999998460
 =1/6-4873540109/41999999995380,                     (6)

75385062494229/452374999952953
 =1/6-64624987579/2714249999717718.                  (7)
```

The floating LP only discovers the integer weights.  Ordinary and optimized
replays check all stored loads with `Fraction` arithmetic and are
byte-identical.

THM-1182's BV estimate

```text
|mu_k(Dbar_d)-1/7|<=3 TV_R(mu_k)/(49d)                (8)
```

combines with (4) to control every `d>B_k`; equation (3) controls the finite
range.  An independent standard-library verifier recomputes the interval
loads using

```text
F(y)=floor(y)/7+min(frac(y),1/7),
|[u,v] intersect Dbar_d|
 =[F(dv+1/14)-F(du+1/14)]/d.                         (9)
```

It reproduces (6)--(7) exactly.

## 2. Reflection has one fixed phase

After reducing `k` modulo `43`, reflection `t -> 1-t` sends

```text
G_k -> G_(42-k),                                     (10)
```

with endpoint identities

```text
1-(14k+13)/(14*43)=(14(42-k)+1)/(14*43),
1-(14k+1 )/(14*43)=(14(42-k)+13)/(14*43).            (11)
```

Every integer danger comb is reflection invariant.  The phase involution has
the unique fixed point

```text
k=42-k  iff  k=21.                                   (12)
```

Thus the `43` residues split into `21` pairs and the fixed phase `21`, giving
exactly the `22` representatives

```text
k=0,1,...,21.                                        (13)
```

The Lean consumer proves (11)--(13), including representative coverage of
every `0<=k<=42`.

If six combs covered `G_k`, the exact measure would give

```text
1<=mu_k(union_(i=1)^6 Dbar_(d_i))
 <=sum_(i=1)^6 mu_k(Dbar_(d_i))<1,                   (14)
```

a contradiction.

## 3. Structural and tournament fingerprints

The hardest finite row is

```text
(k,d,q,s,M)=(15,56,1,13,28),       d=43q+s.          (15)
```

The sharp residue has moved from `s=3` at carriers `41` and `42` to `s=13`.
The fixed reflection phase `k=21` is not extremal.  The hardest tail is the
short `k=1,M=12` row, with normalized zero-extension variation

```text
500468023300/124999999987.                           (16)
```

Hence neither the central fixed phase nor the deepest cutoff controls this
carrier.  The obstruction is again a low vertical rung on the cyclic orbit

```text
dt=(q+s/43)u+sk/43       modulo 1.                   (17)
```

We challenged runners, teeth, walls, residues, bins, gap phases, and proof
obligations as vertices.  The faithful objects are density-decorated phase
stalks.  Orienting their bottleneck loads, with `k` as gauge, gives a
transitive 22-vertex tournament: score histogram `0,...,21`, no directed
cycles, singleton SCCs, and one Hamiltonian path.  It records difficulty but
forgets the weights, so it is a loss audit rather than the proof carrier.

## 4. Formalization, hashes, and frontier

The independent replay passes all `22` rows, and the Lean module checks the
two rational margins, reflection endpoints, the unique fixed phase,
representative coverage, the common load cap, and the six-load contradiction.
JSON completeness, interval integration, and the analytic BV estimate remain
explicit external providers.

```text
generator_sha256   = de2fc18021419614bf37071d0f3f85177b3585e2b17e449be0ec56ef9c0b44bf
independent_sha256 = f27a813f8f2b7d0e53642cf1d052ccab122d1f6a7545ef9584c7202a42ace90d
certificate_sha256 = def4e2893ac89c8b03affb155fa0da09db7f4d5efa37c391147df03cfd65fe65
output_sha256      = 02f8593615337a240fd9b60d9a1b78d778c61a2d714c4306bd21f02e069eccde
ind_output_sha256  = a4309a4886c2818b9a292e2ddfbdfcf4fe3a611721e5fc5f6a36069e3c86c660
lean_sha256        = 740435b38fe1f8cc9de33a93233aa79f73e08de621f17875ce8a83fc96d32281
certificate_bytes  = 53513
```

Together with THM-1182/1197/1255/1257, this closes all carriers through
`43`.  It does not prove LRC(14) or uniform six-comb noncoverage.  Carrier
`44` is the next exact rung; the structural goal remains a parametric
phase-orbit density-selection theorem that can absorb THM-1254's entire
bounded `c<=1171` branch.

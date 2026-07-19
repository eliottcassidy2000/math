---
id: THM-1265
title: The bounded-variation needle certificate closes carrier 45 exactly
status: PROVED FINITE-EXACT for carrier c=45, every slow-gap phase, and every faster integer speed. Twenty-three reflection-representative 200-bin probability densities have exact closed-comb load below 1/6 on a finite vertical ladder and an exact BV-tail cap below 1/6 thereafter. Consequently no six faster combs cover a complete 45-gap, even with repeated speeds. Together with THM-1182/1197/1255/1257/1259/1261 this closes every carrier 1<=c<=45; carriers c>=46 remain open
source: codex-2026-07-19 carrier-45 exact census
depends_on: [THM-1182, THM-1197, THM-1255, THM-1257, THM-1259, THM-1261]
related: [THM-1176, THM-1198, THM-1232, THM-1250, THM-1254, THM-1264]
script: 04-computation/lrc14_slow_gap_bv_needle_carrier45_thm1265.py
independent_verifier: 04-computation/lrc14_slow_gap_bv_needle_carrier45_thm1265_verify.py
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCCarrier45BVNeedle.lean
certificate: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier45_thm1265_certificate.json
output: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier45_thm1265.out
independent_output: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier45_thm1265_verify.out
generator_sha256: c1a63571fe59ef64e7890afc203b2e43430009a0d1d0177196cf87eede4488ca
independent_sha256: 867d370e114f9d4abdcccb889d9803e4acf30dcac8eb0771b3745560538f4e3a
certificate_sha256: 6c25cbae1eef02f18bfa5d4b3d573519d5c3e1a68b9d38dd0605b5310b738ebc
output_sha256: ce8206bc4bc449d79a44211ccccf7487f884b6a717e6d9a21481b82bfeaa01bd
independent_output_sha256: e14d86e1f6328041e9974a90bf08ecb991e4abfb7738dcf91bfaa72c8e65c17f
lean_sha256: ded8ad6d91eb0294a69db78f484573d2f1b48fd8dcf39ea2f8efc5a28d77898b
---

# THM-1265 -- carrier 45 has no six-comb slow-gap cover

For every positive integer speed `d`, put

```text
Dbar_d={t in R: ||dt||<=1/14},                        (1)
```

and let a complete carrier-`45` safe gap be

```text
G_k=[(14k+1)/(14*45),(14k+13)/(14*45)].              (2)
```

> **Theorem.** For every integer `k` there is an absolutely continuous
> probability measure `mu_k`, supported on `G_k`, such that
>
> ```text
> mu_k(Dbar_d)<1/6                                    (3)
> ```
>
> for every integer `d>45`.  Therefore no six faster danger combs cover any
> complete `45`-gap.

The conclusion permits repeated faster speeds.  It consequently excludes
distinct and primitive six-comb packets as special cases.

## 1. The 23 exact density rows

For each representative phase, divide `G_k` into `200` equal bins and assign
the stored nonnegative integer masses

```text
w_(k,0),...,w_(k,199),               W_k=sum_j w_(k,j)>0.  (4)
```

Give bin `j` constant density `w_(k,j)/(W_k |I_j|)`.  The resulting
probability measure satisfies, for its row-specific multiplier `M_k`,

```text
max_(46<=d<=45M_k) mu_k(Dbar_d)<1/6,                 (5)

1/7+3 TV_R(mu_k)/(49(45M_k+1))<1/6.                 (6)
```

Here `TV_R` is the variation of the zero-extended density, including both
support-endpoint jumps.  The exact cutoff histogram is

```text
M          12  16  28  44  60  108
rows        7   1   4   9   1    1.                 (7)
```

The rows are

```text
M=12:   k=0,1,2,11,14,15,22
M=16:   k=21
M=28:   k=3,4,10,20
M=44:   k=5,8,9,12,13,16,17,18,19
M=60:   k=6
M=108:  k=7.                                         (8)
```

Every phase closes at its predeclared preferred multiplier; no fallback row
is substituted after discovery.  The largest finite and tail loads are
respectively

```text
12991498270615/77999999992356
 =1/6-8501728111/77999999992356,                      (9)

294191299975319/1765399999814633
 =1/6-252199962719/10592399998887798.                (10)
```

Both margins are positive exact rationals.  The floating linear program is
discovery-only: ordinary replay checks all stored rows with `Fraction`
arithmetic.  An optimized replay is byte-identical.

THM-1182's bounded-variation estimate

```text
|mu_k(Dbar_d)-1/7|<=3 TV_R(mu_k)/(49d)               (11)
```

turns (6) into (3) for every `d>45M_k`; equation (5) handles the finite
ladder.  Thus the finite and analytic layers meet with no unverified speed.

## 2. Structurally independent cumulative replay

The independent verifier does not import the generator or its clipped-tooth
integration routine.  It uses the cumulative function

```text
F(y)=floor(y)/7+min(frac(y),1/7)                     (12)
```

and reconstructs every bin load from

```text
|[u,v] intersect Dbar_d|
 =[F(dv+1/14)-F(du+1/14)]/d.                        (13)
```

It reproduces all stored finite maxima, variations, and tails exactly.  The
finite bottleneck is

```text
(k,M,d)=(6,60,156),              156=45*3+21.        (14)
```

The tail bottleneck is the `k=3,M=28` row, whose normalized unit-gap
zero-extension variation is

```text
1866280000080/199999999979.                          (15)
```

Both verifier modes are optimization-safe: all proof-critical checks raise
explicit exceptions, each script parses its own source and rejects any
Python `assert` node, and normal and optimized outputs are byte-identical.

## 3. Reflection and the fixed phase

Reduce `k` modulo `45`.  Reflection `rho(t)=1-t` sends

```text
G_k -> G_(44-k),                                     (16)
```

because

```text
1-(14k+13)/(14*45)=(14(44-k)+1)/(14*45),
1-(14k+1 )/(14*45)=(14(44-k)+13)/(14*45).            (17)
```

Every integer danger comb is reflection invariant.  The `45` phases split
into `22` two-element orbits and the single fixed phase

```text
k=22=44-22.                                          (18)
```

Therefore `k=0,...,22` is a complete representative bank.  The fixed row is
stored explicitly at multiplier `12`; it is not inferred from another phase.
The Lean consumer proves both endpoint identities, the involution, uniqueness
of (18), and representative coverage of every `0<=k<=44`.

Finally, if six combs covered `G_k`, then

```text
1=mu_k(G_k)
 <=mu_k(union_(i=1)^6 Dbar_(d_i))
 <=sum_(i=1)^6 mu_k(Dbar_(d_i))<1,                  (19)
```

a contradiction.  Passing from strict danger teeth to their closed versions
only strengthens the hypothetical cover.

## 4. What the cutoff clock says

Carrier `45` has the first multiplier-`108` row in the exact carrier ladder:

```text
(k,M)=(7,108).                                       (20)
```

It is not the finite bottleneck.  Its finite gap below `1/6` is

```text
10876245421/13999999998642 > 7.7e-4,                 (21)
```

whereas the `k=6,M=60` bottleneck gap in (9) is about `1.09e-4`.  The deep
cutoff primarily finances variation for the infinite tail; it does not
identify a remote speed that nearly covers the selected density.  Likewise,
the fixed phase `k=22` closes already at `M=12` with a finite gap exceeding
`3.8e-3`.

Write

```text
d=45q+s,                         0<=s<45,
t=(k+u)/45.                                           (22)
```

Then the normalized orbit law is

```text
dt=(q+s/45)u+sk/45                  modulo 1.         (23)
```

The multiplier histogram is therefore a clock for the difficulty of
selecting a density on this coupled slope/phase orbit.  The jump to `108`
does not by itself shrink the unbounded carrier frontier.  The structural
target remains a parametric density-selection law for (23) with controlled
variation, capable of replacing one certificate row per phase and carrier.

## 5. Tournament and alternate-carrier audit

Runner vertices are especially lossy: one density row controls infinitely
many faster runners at once.  We challenged runners, individual teeth,
fixed gap sections, section boundaries, wall-crossing events, residues,
Fourier modes, bins, phase gaps, and proof obligations as vertices.  The
proof-bearing carrier is the `23`-vertex family of **density-decorated phase
stalks**, each retaining its `200` weights, finite cutoff, exact worst load,
variation, and tail cap.

For the required pairwise loss audit, orient stalk `k` before stalk `ell`
when its exact bottleneck `max(finite load,tail load)` is smaller, breaking a
tie by the phase label.  This is a transitive tournament with score histogram
`0,1,...,22`, no directed triangles, `23` singleton SCCs, and one Hamiltonian
path:

```text
11,8,9,6,14,0,2,15,5,18,22,10,7,20,16,17,13,1,19,4,21,12,3.  (24)
```

The tournament preserves relative row difficulty and the cutoff clock.  It
destroys every density weight and therefore cannot prove (3).  Section
boundaries and wall events reconstruct individual finite loads, but without
the variation sidecar they cannot control the infinite tail.  This is a loss
audit, not a disguised runner tournament.

## 6. Reproducibility, formalization, and frontier

The exact certificate contains `23` rows and `56,099` bytes.  The independent
cumulative replay passes all rows.  The sorry-free Lean module checks both
rational margins, reflection endpoints, the fixed-phase equivalence,
representative coverage, the common load cap, and the six-load contradiction.
Certificate completeness, exact bin integration, and the analytic BV estimate
(11) remain explicit external providers; there are no proof placeholders or
`native_decide` calls.

```text
generator_sha256   = c1a63571fe59ef64e7890afc203b2e43430009a0d1d0177196cf87eede4488ca
independent_sha256 = 867d370e114f9d4abdcccb889d9803e4acf30dcac8eb0771b3745560538f4e3a
certificate_sha256 = 6c25cbae1eef02f18bfa5d4b3d573519d5c3e1a68b9d38dd0605b5310b738ebc
output_sha256      = ce8206bc4bc449d79a44211ccccf7487f884b6a717e6d9a21481b82bfeaa01bd
ind_output_sha256  = e14d86e1f6328041e9974a90bf08ecb991e4abfb7738dcf91bfaa72c8e65c17f
lean_sha256        = ded8ad6d91eb0294a69db78f484573d2f1b48fd8dcf39ea2f8efc5a28d77898b
certificate_bytes  = 56099
```

Together with THM-1182/1197/1255/1257/1259/1261, this closes every slow
carrier through `45`.  It does not prove uniform six-comb noncoverage or
LRC(14).  Carrier `46` is the next exact finite rung; the higher-leverage
target is still a uniform phase-orbit density-selection theorem that absorbs
the growing cutoff clock rather than extending it one carrier at a time.

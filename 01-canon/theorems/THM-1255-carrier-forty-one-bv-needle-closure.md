---
id: THM-1255
title: The bounded-variation needle certificate closes carrier 41 exactly
status: PROVED FINITE-EXACT for carrier c=41, every slow-gap phase, and every faster integer speed.  Twenty-one reflection-representative 200-bin probability densities have exact closed-comb load below 1/6 on a finite vertical ladder and an exact BV-tail cap below 1/6 thereafter.  Consequently no six faster combs cover a complete 41-gap, even with repeated speeds.  Together with THM-1182/1197 this closes every carrier 1<=c<=41; carriers c>=42 remain open
source: codex-2026-07-19 carrier-41 exact census
depends_on: [THM-1182, THM-1197]
related: [THM-1176, THM-1198, THM-1232, THM-1250, THM-1253, THM-1254]
script: 04-computation/lrc14_slow_gap_bv_needle_carrier41_codex_20260719.py
independent_verifier: 04-computation/lrc14_slow_gap_bv_needle_carrier41_codex_20260719_verify.py
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCCarrier41BVNeedle.lean
certificate: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier41_codex_20260719_certificate.json
output: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier41_codex_20260719.out
independent_output: 05-knowledge/results/lrc14_slow_gap_bv_needle_carrier41_codex_20260719_verify.out
---

# THM-1255 -- carrier 41 has no six-comb slow-gap cover

For an integer speed `d`, let

```text
Dbar_d={t in R: ||dt||<=1/14}
```

be the closed danger comb, and put

```text
G_k=[(14k+1)/(14*41),(14k+13)/(14*41)].              (1)
```

> **Theorem.**  For every integer `k` there is an absolutely continuous
> probability measure `mu_k`, supported on `G_k`, such that
>
> ```text
> mu_k(Dbar_d)<1/6                                    (2)
> ```
>
> for every integer `d>41`.  Consequently no six faster danger combs cover
> any complete `41`-gap.

The statement allows repeated faster speeds.  It therefore excludes, a
fortiori, every distinct or primitive six-comb packet at the first carrier
beyond THM-1197's range.

## 1. The exact 21-row certificate

For each `0<=k<=20`, divide `G_k` into `200` equal bins.  The certificate
stores nonnegative integer masses

```text
w_(k,0),...,w_(k,199),             W_k=sum_j w_(k,j)>0. (3)
```

Give bin `j` constant density `w_(k,j)/(W_k |I_j|)`.  This defines a
probability measure `mu_k`.  For a row-specific multiplier `M_k`, the two
exact checks are

```text
max_(42<=d<=41 M_k) mu_k(Dbar_d)<1/6,                 (4)

1/7+3 TV_R(mu_k)/(49(41 M_k+1))<1/6.                 (5)
```

Here `TV_R(mu_k)` means the total variation of the zero-extended density,
including both support-endpoint jumps.  The multiplier census is

```text
M          12  16  20  28  44  60
rows        7   1   4   5   3   1.                  (6)
```

The worst finite and tail loads are respectively

```text
5498776467337/32999999996568
 =1/6-1223532091/32999999996568,                     (7)

957696499918301/5746999999316107
 =1/6-820999806301/34481999995896642.                (8)
```

Both margins are positive exact rationals.  The floating linear program is
discovery-only: the stored integer masses and every load in (4)--(5) are
replayed with `Fraction` arithmetic.

## 2. Why the finite ledger controls every faster speed

THM-1182 proves the elementary bounded-variation estimate

```text
|mu_k(Dbar_d)-1/7|<=3 TV_R(mu_k)/(49d).               (9)
```

For `d>41M_k`, monotonicity in `d` and (5) give

```text
mu_k(Dbar_d)
 <=1/7+3 TV_R(mu_k)/(49d)
 <=1/7+3 TV_R(mu_k)/(49(41M_k+1))
 <1/6.                                                (10)
```

The finite range is exactly (4), so (2) holds for every `d>41`; no ratio
truncation from THM-1232 is being used.

The finite integrations have a second, independent implementation.  If

```text
F(y)=floor(y)/7+min(frac(y),1/7),                     (11)
```

then on a bin `[u,v]`

```text
|[u,v] intersect Dbar_d|
 =[F(dv+1/14)-F(du+1/14)]/d.                         (12)
```

The standard-library verifier recomputes (4) using (11)--(12), rather than
importing the generator's clipped-tooth integrator, and independently obtains
the same maxima and margins.

## 3. Why 21 phases are all phases

First reduce `k` modulo `41`.  Directly from (1),

```text
G_(k+41)=G_k+1,                                      (13)
```

and every integer danger comb is one-periodic.  Reflection `rho(t)=1-t`
sends

```text
G_k -> G_(40-k),                                     (14)
```

because

```text
1-(14k+13)/(14*41)=(14(40-k)+1)/(14*41),
1-(14k+1 )/(14*41)=(14(40-k)+13)/(14*41).            (15)
```

Moreover `||d(1-t)||=||dt||`, so `rho` preserves every `Dbar_d`.  The
residues `0,...,40` therefore split into the `20` pairs `{k,40-k}` plus the
fixed phase `k=20`.  The rows `0<=k<=20` prove all phases.  Both identities
in (15) and the downstream six-load contradiction are checked in the
sorry-free Lean consumer.

Finally, if six combs covered `G_k`, the union bound against `mu_k` would give

```text
1=mu_k(G_k)
 <=mu_k(union_(i=1)^6 Dbar_(d_i))
 <=sum_(i=1)^6 mu_k(Dbar_(d_i))<1,                   (16)
```

a contradiction.  Closing the teeth only strengthens the strict-open lonely
runner formulation.

## 4. Structural fingerprints beyond the yes/no census

Write every faster speed as

```text
d=41q+s,                     0<=s<41.                (17)
```

The largest finite load occurs at

```text
(k,d,q,s,M)=(8,44,1,3,16).                           (18)
```

Thus the first post-carrier ladder rung, not a remote tooth, is the sharp
finite obligation.  The largest tail cap is the `k=18, M=20` row; its
unit-gap zero-extension variation is

```text
6668121952000/999999999881.                          (19)
```

The longest certificate is `k=15, M=60`, but its worst finite speed is only
`d=81=41+40`.  Its large cutoff finances the normalized variation

```text
4997030487700/249999999973,                          (20)
```

rather than defending against a genuinely remote dangerous speed.  This is
the same functional `H`-drift already visible in THM-1197: a phase row buys
low-frequency separation by spending variation, and the tail charges that
variation at order `1/M`.

The normalized orbit is the underlying object.  With `t=(k+u)/41` and
`d=41q+s`, the phase is

```text
dt=(q+s/41)u+sk/41       modulo 1.                   (21)
```

The certificate is therefore a density-selection result on the cyclic orbit
`{(s/41,sk/41)}` with a vertical `q`-ladder, not a brute-force list of
six-tuples.

## 5. Tournament and alternate-carrier audit

We challenged runners, teeth, fixed circle sections, section boundaries,
wall events, residues, bins, phase gaps, and proof obligations as tournament
vertices.  Runner vertices are particularly lossy here: one row controls
infinitely many possible runners at once.  The faithful carrier is the
`21`-vertex family of **phase-obligation stalks**, each decorated by its full
`200`-bin density and finite/tail split.

For the pairwise loss audit, orient stalk `k` before stalk `ell` when its exact
bottleneck `max(finite load, tail load)` is smaller, breaking equality by the
phase label.  This is a transitive tournament with score histogram
`0,1,...,20`, no directed cycles, `21` singleton SCCs, and one Hamiltonian
path.  It preserves only row difficulty and destroys the density weights;
therefore it cannot prove (2).  The useful tournament conclusion is negative:
the binary relation is merely a clock for the certificate bank, while the
stalk decoration carries the mathematics.

## 6. Reproducibility, formalization boundary, and frontier

Ordinary and optimized generator replays are byte-identical.  The independent
cumulative-formula verifier also passes all `21` rows.  The frozen artifact
hashes are

```text
generator_sha256   = 186c1c4ec9a75c6ccada7a0b98616f6ed147b0fc9d8874ffe296dfa8ef459fe3
independent_sha256 = aef831215129fd9a459641997a12543b8055e04f153edb45bf57d0f63e291a74
certificate_sha256 = e31a44e6dc3cc9559262faeee530b272710ab3388f2f449b6502fa05be44def8
output_sha256      = 15b4b44b57ef537136d9e5c8695295d7491f90edd79d254e6cece1f606abdb84
ind_output_sha256  = c59bed4ace4801a9dc92a95ca6824b063f0d504df60c6026f75d3a86dcbb5c61
lean_sha256        = 5e648bc051afc103201e5df72680ad8e5b1cd2f955826034b8555ea588a3ca76
certificate_bytes  = 51092
```

The Lean module checks the two displayed rational margins, both reflection
endpoint identities, the finite/tail load consumer, and the six-load
contradiction.  The `200`-bin certificate completeness, exact interval
integration, and analytic BV estimate (9) remain explicit paper/computation
providers; there is no hidden `sorry` claiming otherwise.

One post-certificate transport guardrail adjusts the next target.  Reusing
the `k`-th carrier-`41` weight vector unchanged on the `k`-th carrier-`42`
gap passes exact `Fraction` replay only for `k=0,1`; all `19` rows
`2<=k<=20` fail the finite-load `<1/6` test.  A broader floating scout that
allowed arbitrary reselection among all `21` stored carrier-`41` shapes also
found no template for those `19` phases.  The latter is telemetry, not a
theorem, while fresh exact LP discovery succeeds on the sampled carrier-`42`
phases `k=2,8,15,20`.  Thus literal stalk transport and a fixed finite atlas
are the wrong self-similarity targets: the operation that selects a new
density from the changing phase orbit is the object that must become uniform.

Together with THM-1182 and THM-1197, this theorem closes every slow carrier
through `41`.  It does not prove uniform six-comb noncoverage or LRC(14).
The next exact carrier is `42`, while the structural target is a uniform
density-selection lemma on the cyclic orbit (21) whose variation budget is
bounded independently of the carrier.  Equivalently, a well-founded carrier
or address descent extracted from THM-1253/1254's fully paid chronological
word would now have an exact terminal base through `41`; this theorem supplies
that base but not the descent.

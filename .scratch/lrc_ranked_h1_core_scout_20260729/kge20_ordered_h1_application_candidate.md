# K>=20 marked-suffix ordered-H1 application (candidate)

Status: `FINITE-EXACT`; scoped computation, not an LRC(14) theorem.

## Universe and inherited proof mechanism

The universe is the exact 62-root `K(E)>=20` slice of the locked THM-2896
adaptive gate atlas.  It contains 1,289 marked gate branches.  For each
scalar-open branch the data retained are

```text
E, K(E), gate rank j, apex a_j,
P_j={a_1,...,a_j},
C=G_E\D_(a_j), h=|C|, and the literal interval carrier C.
```

The allowed labels are the external speeds `w>=15` outside `P_j`.  Write
their exact singleton coverages in decreasing order as

```text
q_1>=...>=q_5,                  R_4=q_1+...+q_4.
```

On the 181 branches satisfying

```text
R_4<6h/7,
```

THM-2897 supplies a valid four-label complement cap and THM-2893, with
`k=5,s=1`, forces every label of a hypothetical five-cover into

```text
H_1={w>=15:w notin P_j and |C intersect D_w|>=h-R_4}.
```

The strict THM-735(ii) discrepancy tail makes this core finite, with

```text
epsilon=6h/7-R_4>0,
N=ceil((99/70) components(C)/(7 epsilon))-1.
```

The verifier scans every allowed `15<=w<=N`; no finite-window proposal is
treated as a global core.

## Ordered pivot certificate

Order the exact core as

```text
H_1=(x_1,...,x_t)
```

by decreasing coverage on `C`, breaking ties by speed.  A five-subset of
`H_1` has a unique earliest member `x_i`.  Its other four members lie in
the exact later suffix `{x_(i+1),...,x_t}`.  Put

```text
C_i=C\D_(x_i),
u_i=sum of the four largest |C_i intersect D_y|
    over y in {x_(i+1),...,x_t}.
```

If fewer than four labels remain, `x_i` cannot be the earliest member of a
five-set.  Otherwise the strict inequality

```text
u_i<|C_i|
```

excludes every five-set with earliest member `x_i` by subadditivity.
Checking this for every `i` excludes all five-covers.

The implementation first uses the cheaper inherited upper bounds

```text
|C_i intersect D_y|<=|C intersect D_y|.
```

Only if their best-four sum fails does it compute all exact child
coverages on the same later `H_1` suffix.  It never enlarges that suffix
to all allowed labels.  This is exactly the `s=1` specialization of
THM-2893's ordered high-core pivot refinement.

## Exact K>=20 result

With cutoff `N<=15000`, 180 of the 181 eligible branches have exact finite
cores.  Their raw five-subset workloads reach

```text
|H_1|=355,             binom(|H_1|,5)=45,674,610,696.
```

The ordered certificate closes 178 of these 180 branches.  Across those
certificates, at most four pivots per branch require exact local
recomputation, and the minimum positive certified margin is

```text
154787/42146180076.
```

Exactly two canonical boundary controls have one unresolved first pivot:

```text
E=(1,2,5,7,8,9,11), K=20, rank=5, apex=65,
P=(26,19,17,23,65), |H_1|=51, binom(|H_1|,5)=2,349,060,
first unresolved pivot=25;

E=(1,2,5,7,9,11,12), K=21, rank=3, apex=23,
P=(26,16,23), |H_1|=29, binom(|H_1|,5)=118,755,
first unresolved pivot=32.
```

The deterministic exact residual branch-and-bound closes these in `102`
and `58` nodes respectively.  More structurally, they are repaired
immediately by THM-2897's relative Hunter maximum-spanning-tree cap.  Each
branch has exactly one star-hostile five-set:

```text
Q=(20,25,30,40,60):
  pivot-star invoice = 110063/600600,
  maximum-tree invoice Psi = 6673/54600,
  h-Psi = 69337/1146600,
  nonstar extra credit = 47/770,
  maximum-tree edges =
    (20,40):3/140,
    (30,60):19/924,
    (20,30):2/105,
    (25,60):1/350;

Q=(17,25,32,34,50):
  pivot-star invoice = 295612237/1517392800,
  maximum-tree invoice Psi = 12978773/89258400,
  h-Psi = 4404907/89258400,
  nonstar extra credit = 94663/1915900,
  maximum-tree edges =
    (25,50):43/2100,
    (17,34):2141/114954,
    (17,25):92/8925,
    (25,32):1/350.
```

Thus no witness occurs, and the clean structural decomposition is

```text
178 ordered-depth-one closures
+ 2 genuine nonstar Hunter-tree closures
= 180 finite-core Hunter closures.
```

The depth-two branch-and-bound remains an independent exact control, not
the preferred theorem mechanism.

The remaining eligible row is deliberately separate:

```text
E=(2,5,6,8,9,12,14), K=20, rank=3, apex=33,
P=(22,26,33),
epsilon=7081/111831720,
N=121253.
```

It is a cutoff deferment, not a failed certificate and not an open witness.

## Reproduction state

Scratch verifier:

```text
.scratch/lrc_ranked_h1_core_scout_20260729/scout.py
```

Locked imported ranked-suffix engine:

```text
04-computation/lrc14_j6_all_root_ranked_suffix_scalar_census_codex_20260729.py
sha256=e0ff69252870f194549bba61289c1c5b15bef451e37a72836d71f9e71b1016e9
```

Exact K>=20 output:

```text
.scratch/lrc_ranked_h1_core_scout_20260729/highk_depth1_fallback.out
sha256=8265158592fd70153358b9b955d72f296894275360d59b4779cd17ddcce24c69
```

Current verifier source hash for that replay:

```text
sha256=6abc06972cda64bb4f53db26cca62bbea5362ca178bdbf5bf7398e8b0f28317a
```

The output ledger hash internal to the verifier is

```text
7ff1b8234c190b61025a91e75d1e8f391c94b3d850cce0df0ecebe071719faa9.
```

Relative Hunter audit:

```text
.scratch/lrc_ranked_h1_core_scout_20260729/hunter_two_star_exceptions.py
sha256=c452781c7cf6a2ada6be8984b9c9cbe7aab8369c5a9333721ef8ecc8e0207393

.scratch/lrc_ranked_h1_core_scout_20260729/hunter_two_star_exceptions.out
sha256=f6c7c9fe60df589a9f6b7daf036531d04cfeccf95b6a337a2975d1fedd007d7e

aggregate ledger:
be88409549d4af1e3276d1556c09f30b8d4781958edd9c674cbf801b9205f31d
```

Ordinary and optimized replay are byte-identical at the displayed output
hash.  A no-assert audit (`rg '\bassert\b'`) is empty in both the H1 scout
and the Hunter audit; all load-bearing checks use explicit `require`.

These hashes must be regenerated if the verifier is subsequently
generalized or instrumented; the mathematical counts above must not be
silently carried across a source change.

# Seven-body adaptive gate and complement-cap flag battery

**Status: FINITE-EXACT SCOPED RESEARCH NOTE.  No uniform theorem.**

This note records the cheapest tested successor to the THM-2888
eight-body/five-slot atlas.  It uses THM-2893's complement-cap finite-core
flag lemma, but it samples only `36/3432` seven-body roots and three rank-one
apex carriers.  It does not prove the seven-body/six-slot rung or LRC(14).

## 1. The fixed top-twelve proposal fails

For a seven-body root `E`, let `q_i(E)` be its globally ranked external
single-comb coverages.  A top-`K` hitting gate for a six-speed covering set
is valid when

```text
q_(K+1)+...+q_(K+6) < |G_E|.
```

The verifier takes twelve deterministic lexicographic quantiles from each
of the three exact strata

```text
E subset {1,...,12}:                         792 roots;
exactly one of {13,14} in E:                1848 roots;
both 13 and 14 in E:                         792 roots.
```

It globally tail-seals the top thirty coverages.  The fixed `K=12` gate
passes only `20/36` sampled roots: `10/12`, `5/12`, and `5/12` by stratum.
The first exact counter-sample is

```text
E=(1,4,7,8,10,11,12),
|G_E|-(q_13+...+q_18)=-238717/34474440.
```

Adaptive gates exist throughout the sample, with maximum least values
`K=16,20,21` in the three strata.  The largest top-thirty discrepancy
threshold is `1554`; scanning through `1600` seals every sample.  Thus a
fixed top-twelve successor is refuted, while an adaptive top-`K` atlas
remains computationally plausible.

THM-885 does **not** remove the `792` low roots: its `j=6` sweep is explicitly
unfinished.  Treating those roots as inherited terminals would be a scope
error.

## 2. Correct five-slot recursion

After a first apex, six slots become five.  A pair cap `B_2<5h/7` alone
does not reproduce THM-2888's four-slot argument.  THM-2893 gives two
correct alternatives.

1. Compute a global triple cap `B_3<5h/7`.  With
   `theta=h-B_3`, every pair in a hypothetical five-cover is
   `theta`-heavy.  At least four cover vertices lie in the finite high core,
   so a heavy high triangle remains; its literal residual must be excluded
   from every pair.
2. Prove the stronger pair cap `B_2<4h/7`.  With
   `theta=h-B_2`, every triple in a hypothetical five-cover is heavy and at
   least three cover vertices lie in the corresponding high core.  Every
   heavy high triple again leaves a literal residual pair obligation.

For the first route, the finite triple head is exact.  If one speed exceeds
`2500`, the strict discrepancy bound gives

```text
U(x,y,z) < B_2 + h/7 + (99/70)r/(7*2501).
```

Every residual pair cap is independently sealed by an exact scan through
`2500`, a global rank-one tail check, and a direct reconstruction of its
first paid finite maximizer.

## 3. Three hostile rank-one apex cases

The battery selects roots whose adaptive gates require `K=19,20,21`:

```text
E=(2,8,9,10,11,13,14), a=19;
E=(1,3,9,10,11,12,14), a=39;
E=(2,5,9,11,12,13,14), a=16.
```

All three satisfy both strict flag hypotheses.  The cheap estimate
`B_3<=B_2+q_1` fails on the first and third, but exact triple maximization
repairs them with positive margins

```text
31389571261/3605249423160,
810055139/33369996660,
2283219/451250800.
```

For the `(k,s)=(5,2)` route the finite cores have sizes `29,6,54`.
There are `6106` heavy triangles, and all `6106` literal residuals have
global pair cap strictly below their residual mass.

For the `(k,s)=(5,3)` route the finite cores have sizes `33,12,24`.
Among `7700` scalar-admissible triples, exactly `2207` are heavy.  All
`2207` residual pair obligations close; `1072` coincide with obligations
from the first route.  Thus the alternative strong-pair flag is not merely
formal on this hostile battery.

## 4. Exact scope and next probe

The positive battery is evidence for an adaptive seven-body atlas, not a
proof of one.  The cheapest next computation is:

1. exact adaptive top-`K` gates on a larger stratified root sample;
2. pair caps on the retained first apices;
3. exact triple caps only where `B_2+q_1<5h/7` fails;
4. choose per carrier between THM-2893's two flag routes by predicted core
   size and residual workload.

Even a uniform success would prove only the seven-body/six-slot near-AP
rung.  It would not settle configurations with at most six speeds in
`{1,...,14}`, hence would not prove unrestricted LRC(14).

## 5. Reproduction

```text
04-computation/lrc14_j6_adaptive_gate_triple_heavy_battery_codex_20260729.py
SHA-256 f0080ea06c805481ec0fd13ff122f903c52cd389a425ebe91ab6713ec0bcdadb

05-knowledge/results/lrc14_j6_adaptive_gate_triple_heavy_battery_codex_20260729.out
SHA-256 6cd26c47bdfca219e59bfd05e9e9bdc2e25faba00e67092a8e1f48c4f806f9a6
```

The script hash-pins the THM-2888 exact kernel and transcript, uses rational
arithmetic for every decision, and locks both aggregate digests.  Its stored
output ends with `all_exact_controls=PASS`.

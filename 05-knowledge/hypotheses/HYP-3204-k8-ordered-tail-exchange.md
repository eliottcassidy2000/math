# HYP-3204: k=8 Ordered-Tail Exchange

**Status:** EVIDENCE / exact bounded-bank scout; not an LRC14 proof.

## Claim Tested

HYP-3200 and HYP-3161 leave the k=8 hard node in a cleaner form: consecutive
speeds are exact extremizers for several bimodality/covariance signals in the
bounded bank, but the odd Worpitzky/central channel cannot be collapsed to a
scalar `1/7` law.  HYP-3204 tests whether the remaining `L_y` target can be
split into two proof-safe pieces:

1. Replace full convex order by a one-sided high-tail/stop-loss barrier for
   the empty-sector count `N`.
2. Replace raw `q3` maximization by an exchange-rate lemma: any gain in
   central mass `q3` must cost at least as much extreme bimodality
   `q0+q6`.

The exact bank is the same anchored k=8 bounded bank used by HYP-3200:

```text
E = {0} union A,  A subset {1,...,14}, |A|=7.
```

It has `3432` rows and `3431` primitive rows.  The nonprimitive even AP
`(0,2,4,6,8,10,12,14)` has the same miss-count distribution as consecutive
speeds, so the theorem-facing statement should remain in primitive normal
form.

## Exact Readout

For consecutive speeds `E=(0,1,2,3,4,5,6,7)`, the miss-count distribution is

```text
q = [481/1470, 359/1470, 25/147, 26/245, 17/210, 5/98, 1/49].
```

The full stop-loss / increasing-convex-order route is false:

```text
primitive rows with some stop-loss above consec = 3429
primitive rows dominating consec in stop-loss order = 0
low stop-loss beaters: t0=3429, t1=2355, t2=8
```

The surviving route is one-sided and upper-tail:

```text
primitive beaters above consec:
q0=0
q5=0
q6=0
tail_ge_4=0
tail_ge_5=0
tail_ge_6=0
stop_ge_3=0
stop_ge_4=0
stop_ge_5=0
bimod=q0+q6=0
bimod_plus_q3=q0+q6+q3=0
Ly=q0+q6+q3/10=0
```

Two guardrails are exact:

```text
tail_ge_3 has 431 primitive beaters
q3 has 2879 primitive beaters
```

Thus `tail_ge_3` and raw `q3` are not proof targets by themselves.

The positive exchange signal is:

```text
(q3 - q3_consec)_+ <= (q0+q6)_consec - (q0+q6)
```

over all primitive rows in the exact bank.  The scout found `0` violations
among the `2879` primitive rows with `q3 > q3_consec`.  The worst observed
exchange ratio is

```text
12882/17161 = 0.750655556...
at E=(0,1,4,5,9,10,13,14)
```

with

```text
bimod_loss = 17161/114660
q3_gain    = 2147/19110
Ly_margin  = 19841/143325
```

This is stronger than the current `L_y` need: it implies

```text
q0 + q6 + q3 <= (q0 + q6 + q3)_consec
```

and hence also

```text
q0 + q6 + q3/10 <= (q0 + q6 + q3/10)_consec.
```

## Proof-Frontier Use

The new proof angle is an insurance-pricing lemma for the central mass.  The
odd/Worpitzky channel is allowed to increase `q3`, but that increase must buy
its way through a larger loss in the two ordered states `N=0` and `N=6`.
This reframes `L_y` as a dominated exchange rather than as a direct moment,
entropy, or root-radius maximization.

The candidate theorem packet is:

```text
primitive_normal_form
+ q0/q6 bimodality atom
+ upper ordered-tail stop-loss barrier
+ central exchange-rate lemma
=> L_y extremality at consecutive speeds.
```

This connects back to:

- HYP-3210: the Joukowski/Hermite-Biehler bridge is the best current home
  for proving the odd Worpitzky/interlacing sidecar rather than treating
  `q3` as a standalone scalar.
- HYP-3211/HYP-3212: the multiplicative octonion/G2 route is a productive
  negative for `kappa_3`, while the additive cyclotomic/Chebyshev trace is
  the cap route.  This supports keeping HYP-3204 as a coefficient-exchange
  lemma rather than a hidden apex-7 symmetry claim.
- HYP-3221 and the de Moivre-denominator thread reinforce the same guardrail:
  config-blind algebra is not the closure mechanism; the live statement is an
  analytic/equidistribution coefficient inequality.
- HYP-3200: the covariance/ferromagnetic target supplies a nearby even
  proof route.
- HYP-3202: cyclic-distance covariance layers and finite exchange traps are
  the closest current route to proving the needed `q0+q6` atom.
- HYP-3201: law-defect entropy is the quotient-legality guardrail; it is not
  the raw row-extremal entropy route refuted here.
- HYP-3153: `L_y=q0+q6+q3/10` is the bimodality functional in the
  Lee-Yang/Worpitzky/quartic packet.
- HYP-3154/HYP-3162: the Lee-Yang/Joukowski/cyclotomic side explains why
  radius and circle tightness are sidecars, while this lemma prices the
  actual coefficient exchange.
- HYP-3147/HYP-3151: the odd `q3` channel should keep its Worpitzky,
  minority-edge, and ordered-function sidecars rather than being scalarized.

The local-compression route remains a guardrail, not a proof: one-coordinate
single-step local maxima counts are `432` for `L_y`, `420` for `bimod`, and
`363` for `Sigma kappa_2`.  Greedy coordinate descent has many traps.

## Tournament Analysis

Vertices are proof angles and coefficient functionals, not runners, arcs, or
sectors.  The pairwise observable is exact bounded-bank survival score.  The
switch/gauge directs `A -> B` when angle `A` has higher survival score, with
lexical tie-breaking.

The resulting priority path is:

```text
central_exchange_rate_lemma
-> upper_stoploss_barrier
-> q0_bimodality_atom
-> primitive_normal_form
-> full_convex_order_route
-> single_step_compression_gradient
-> raw_q3_maximization
-> raw_entropy_route
```

This deliberately challenges the default modeling assumption.  The quotient
preserves the LRC predicate `L_y` and selected coefficient/tail sidecars
`q0,q3,q5,q6,tail_ge_m,stop_ge_m`; it destroys the full miss-count
distribution and all runner/arc geometry not reconstructed by those sidecars.

## Remaining Obligations

1. Prove the central exchange-rate lemma in primitive normal form.
2. Prove or import the `q0+q6` bimodality atom, likely through the
   covariance/ferromagnetic or reflection-fold machinery.
3. Lift the exact bounded-bank packet to the legal LRC14 bounded-core proof
   setting, with the even-AP dilation handled as nonprimitive normalization.
4. Do not use full convex order, `tail_ge_3`, raw `q3`, entropy, or naive
   coordinate compression as terminal proof routes.

## Artifacts

- Script: `04-computation/lrc_k8_ordered_tail_exchange_codex_20260628.py`
- Result: `05-knowledge/results/lrc_k8_ordered_tail_exchange_codex_20260628.out`

## Links

HYP-3204 refines HYP-3212, HYP-3211, HYP-3210, HYP-3203, HYP-3202, and
HYP-3200 and should be read with HYP-3201, HYP-3163, HYP-3162, HYP-3161,
HYP-3160, HYP-3154, HYP-3153, HYP-3152, HYP-3151, HYP-3150, HYP-3147,
HYP-3144, HYP-3142, HYP-3139, HYP-3138, HYP-3132, HYP-3122, THM-577, T1304,
LTI-304, LTT-204, and OPEN-Q-108.

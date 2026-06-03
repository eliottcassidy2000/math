---
id: HYP-2076
status: PROGRESS - case certificates found; global descent lemma open
source: codex-2026-06-02-S563
related:
  - HYP-1844
  - HYP-1868
  - HYP-2061
  - HYP-2062
  - HYP-2075
  - HYP-2072
  - HYP-2073
  - THM-369
---

# HYP-2076: Recursive multi-sieves should shrink the Tao hard regime to a residual core

## Statement

The coarse THM-369 denominator sieve already proves the full Lonely Runner
bound `1/(k+1)` whenever some denominator `q <= k+1` divides no speed.  Hence
Tao/Bedert-type global lower bounds are operative only on the sieve-covered
core.

HYP-2076 proposes a next proof norm for that core:

```text
each residual proof obligation either
  (1) has a local recursive denominator witness with margin >= 1/(k+1), or
  (2) exports positive weighted frontier mass to child obligations.
```

The S563 ledger tests the first branch using recursive denominator tiers:

```text
T0 small_sieve       : q <= n = k+1
T1 first_fine_window : n < q <= 2n
T2 prime_power_lift  : q is a prime power with base prime <= n
T3 crt_smooth_lift   : every prime factor of q is <= n
T4 external_fine     : remaining q <= Q
```

For a denominator `q`, `t=a/q` certifies a lower bound `theta` when

```text
min_i ||a v_i / q|| >= theta.
```

This is not a new universal asymptotic improvement over Tao or Bedert.  It is a
case/recursive certificate program: where the certificate applies, it reaches
the conjectural `1/(k+1)` scale, far above the general first-moment-plus-surplus
regime.

## Evidence

`lrc_recursive_multisieve_lower_bound_s563.py` runs with denominator cutoff
`Q=220`.

Coarse non-covered examples are settled immediately:

```text
AP_k13:                first witness q=14, margin=1/14
noncovered_random_k13: first witness q=3,  margin=1/3
```

The useful evidence is inside the sieve-covered core, where no `q <= 14`
witness exists:

```text
blind_lcm_2_14:       first full-bound witness t=2/27, margin=2/27
blind_lcm_2_24:       first full-bound witness t=2/27, margin=2/27
S562_packet_n14:      first full-bound witness t=6/23, margin=2/23
S562_packet_n14_lift: first full-bound witness t=3/23, margin=2/23
S562_packet_n17:      first full-bound witness t=9/25, margin=2/25
S562_packet_n18:      first full-bound witness t=2/29, margin=2/29
```

The best witnesses within `Q=220` can move into smooth or external fine
denominators:

```text
blind_lcm_2_14: best t=16/209, margin=16/209
blind_lcm_2_24: best t=15/196, margin=15/196
S562_packet_n18: best t=101/208, margin=15/208
```

Thus the recursive tiers do real work after the coarse sieve goes blind.  The
blind LCM rows are deliberately hostile to all small denominators, yet the
first fine window already restores a full-bound witness.  The HYP-2073
residual packets behave similarly: dyadic or square-payload packets that were
viewed as exported endpoint debt also have local rational witnesses above the
conjectural threshold.

## Relation To Tao And Bedert

Tao proves a global improvement over the first-moment scale, with lower bound
of the form

```text
1/(2k) + c log(k)/(k^2 (log log k)^2)
```

for an absolute constant `c > 0`.

Bedert's later Riesz-product paper gives a polynomial-size global surplus of
the form

```text
1/(2n) + 1/n^(5/3+o(1))
```

in that paper's parameter convention.

S563 does not compete with those as a universal theorem.  Its point is
different: if a recursive sieve tier reaches `1/(k+1)` on a residual class,
then that class is removed from the hard regime entirely.  The remaining target
for a true improvement-beyond-Tao proof is not "all speed sets" but:

```text
sieve-covered residuals with no local tier witness and no positive exported
frontier mass.
```

That is a much narrower object, and it matches the HYP-2073 residual-packet
frontier picture.  Incoming HYP-2075 sharpens the denominator side: the natural
complete pinch moduli are pair-sums, so a future version of this ledger should
extract local witnesses from pair-sum residues before exporting frontier mass.

## Tournament Analysis

Vertices were speed-set proof obligations, not runners.

Pairwise observable:

```text
(best certified margin up to Q, first conjecture-reaching tier, denominator)
```

Switch/gauge:

```text
larger margin wins; ties prefer earlier recursive tier and smaller denominator.
The displayed sample order orients any exact remaining tie.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles=0
sccs=[1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

The transitivity is informative.  A margin-led lower-bound ledger is an ordered
certificate stack.  To see the hidden no-return/cycle-exclusion structure, the
next quotient should keep endpoint owners, wall-crossing events, or residual
frontier signs rather than only best margins.

## Assumption Challenge

Possible vertices considered:

```text
runners, denominators, residues, gaps, fixed circle sections, section
boundaries, wall-crossing events, CRT channels, endpoint owners, cover arcs,
Fourier modes, matroid circuits, residual packets, and proof obligations.
```

Chosen quotient:

```text
vertices = speed-set proof obligations.
```

Predicate preserved:

```text
existence and tier location of a rational witness with margin >= 1/(k+1).
```

Information destroyed:

```text
which runner owns each endpoint, how cover arcs overlap, wall-crossing order,
local cascade factors, and cross-prime cancellation signs.
```

Challenged assumption:

```text
once a set is sieve-covered for q <= n, the only remaining comparison is
global/asymptotic.
```

S563 refutes that working assumption in examples: the first fine window and
smooth denominator lifts can still certify the full conjectural lower bound.

S600 follow-up: HYP-2146 abstracts the Tao/Bedert comparison as a scale-currency
ledger.  The recursive tiers here are not merely denominators; they are places
where the proof earns or pays scale mass.  Ordinary scale savings give `log`,
compressed denominator tiers give `log log`, and meta-tiers of residual packets
or determinant Helly witnesses can introduce `log log log` as either a tax or a
dividend.  This complements HYP-2145's inverse-hyperoperation framing and
reframes the open descent lemma as: identify the unpaid scale currency of any
residual that survives the local tiers.

## Open Work

1. Replace the raw finite `Q` scan by a proof-producing recursive tier rule.
2. Combine this with HYP-2073's frontier product: local witness or positive
   exported product-frontier mass.
3. Keep endpoint ownership so the Tournament Analysis can see no-return cycles
   instead of only a transitive margin order.
4. Test generated residuals from the actual Rosenfeld/Trakulthongchai lift
   workload, not only deterministic stress rows.

## Sources

- Tao, "Some remarks on the lonely runner conjecture": https://arxiv.org/abs/1701.02048
- Bedert, "Riesz products and the Lonely Runner Conjecture: A wider gap of loneliness": https://arxiv.org/abs/2511.16636

## Files

- `04-computation/lrc_recursive_multisieve_lower_bound_s563.py`
- `05-knowledge/results/lrc_recursive_multisieve_lower_bound_s563.out`
- `07-reflections/lrc-recursive-multisieve-lower-bound-ledger-s563.md`

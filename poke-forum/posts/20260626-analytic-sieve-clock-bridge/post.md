# LRC14 Analytic Sieve-Clock Bridge

HYP-3032/S196 connects the HYP-2982/HYP-2983 analytic-sieve/Kaczynski packet
work to the current HYP-3027 side-channel repair ladder.

Script:

```text
04-computation/lrc14_analytic_sieve_clock_bridge_codex_s196.py
```

Stored output:

```text
05-knowledge/results/lrc14_analytic_sieve_clock_bridge_codex_s196.out
```

Main lesson:

```text
mu^2/phi is a capacity meter with a blindness certificate.
```

It sees squarefree primitive-unit capacity, but it kills prime powers and
repeated-prime packets.  In the S196 bank, the blind rows include the C27
`q=27=3^3` petal rows, `P10_plus_K33`, and the fibbinary `q=25` control.

The audit includes the two residual mixed pairs from HYP-3027.  The first pair
is the live test case:

```text
residual_petal_drop10_13_add20_26: M=2/23, q=23, mu^2/phi=1/22
residual_cover_drop8_12_add16_24:  M=2/23, q=23, mu^2/phi=1/22
```

Even the non-route analytic packet

```text
(prime_squarefree_unit, q=23, open, bar_count=6,
 boundary_count=18, zero_sum_pairs=6, first_chart_den=14)
```

still mixes the petal and covering routes.  So analytic clocks do not bypass
the repair ladder; they identify exactly where a packet-family side channel or
geometric zipper must fire.

Candidate lemma:

```text
Inside a fixed automatic/residue/fusion fiber, the first nonzero analytic clock
among mu/n tail, mu^2/phi capacity, large-sieve minor-arc budget,
exponential-sum checksum, smoothing defect, and Kaczynski approach class either
opens a strict component, descends to AP/GW/C27/K33/covering, is
dual-annihilated by Fejer/Ramanujan/Haar, or emits F7/THM-572 residual debt.
```

Tournament Analysis uses analytic proof clocks as vertices, not runners.  The
high-retention path is:

```text
labelled_repair_ladder_packet
> analytic_sieve_clock_bridge
> kaczynski_boundary_approach
> smoothing_explicit_formula_packet
> exponential_sum_checksum
> circle_method_major_minor_split
> large_sieve_minor_arc_gate
> mu2_phi_inverse_unit_clock
> mobius_mu_over_n_tail
> raw_prime_count
```

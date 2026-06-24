---
id: HYP-2982
title: LRC14 analytic sieve packet weights and Goldbach smoothing guardrails
status: EVIDENCE / finite arithmetic atlas and proof-carrier guardrail; not a proof
source: codex-2026-06-24-S163
related:
  - HYP-2981
  - HYP-2979
  - HYP-2978
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2901
  - HYP-2899
  - HYP-2679
  - HYP-1953
  - HYP-1963
  - THM-572
  - OPEN-Q-108
artifacts:
  - 04-computation/lrc14_analytic_sieve_packet_weights_codex_s163.py
  - 05-knowledge/results/lrc14_analytic_sieve_packet_weights_codex_s163.out
  - 07-reflections/lrc14-analytic-sieve-packet-weights-codex-s163.md
---

# HYP-2982: LRC14 Analytic Sieve Packet Weights

This hypothesis reserves and starts the prompt's analytic-number-theory lane:
sums over primes, `sum mu(n)`, `sum mu(n)/n`, `sum mu(n)^2/phi(n)`,
large-sieve and circle-method weights, upper-bound quadratic/Selberg sieve
packets, exponential sums, Goldbach smoothing choices, saddle-point /
explicit-formula terms, and the repo's Kaczynski/Kaczorowski
boundary/exceptional-set threads.

Working thesis:

```text
The LRC14 proof should treat analytic-sieve weights as labelled packet
certificates, not scalar density estimates.  A smoothing or sieve quotient is
allowed only when it declares its kernel, its exceptional-set boundary, and the
packet labels it forgets.
```

Initial anchors:

- HYP-2899/HYP-2900 already keep denominator Mobius/totient data on the product
  ledger `Div(D) x B_r`.
- HYP-2901 shows no fixed finite denominator basis closes LRC14; prime powers
  and exact-period packets remain load-bearing.
- HYP-2978/HYP-2979 show divisor/Ramanujan packets are useful only when
  endpoint, Haar, route, or state-lift labels are retained.
- HYP-2981 supplies the Fejer interval certificate target, where smoothing
  choices become exact interval proof obligations.
- HYP-2679 and the Kaczynski/Bagemihl boundary-function reflections warn that
  approach regions and boundary values must be carried separately.

## Computation

Script:

```text
04-computation/lrc14_analytic_sieve_packet_weights_codex_s163.py
```

Stored output:

```text
05-knowledge/results/lrc14_analytic_sieve_packet_weights_codex_s163.out
```

The script computes, through `N=200000`,

```text
M(x)      = sum_{n<=x} mu(n)
A(x)      = sum_{n<=x} mu(n)/n
G(x)      = sum_{n<=x} mu(n)^2/phi(n)
Phi(x)    = sum_{q<=x} phi(q)
theta(x)  = sum_{p<=x} log p
```

and checks the LRC packet denominators `q in {14,27,41,84,168}`.

Selected readout:

```text
Q        M(Q)     sum_mu/n      G=sum_mu2/phi   G-logQ        Phi(Q)
14       -2      -0.00592741      4.01666667    1.377609          64
27       -1       0.03734102      4.57184343    1.276007         230
41       -1       0.02609263      5.07005772    1.356486         530
8192     22       0.00287450     10.34382878    1.332915    20401610
65536    14       0.00017888     12.42304021    1.332685  1305514926
200000   -1       0.00002507     13.53859018    1.332518 12158598918
```

The important split is:

```text
Phi(z)=sum_{q<=z}phi(q)       grows quadratically as primitive packet capacity.
G(z)=sum_{d<=z}mu(d)^2/phi(d) grows logarithmically as inverse-unit normalizer.
```

For the LRC denominators:

```text
q=14:  phi=6,  mu= 1, mu^2/phi=1/6,  Phi=  64, G=4.0167, 1/G=0.24896
q=27:  phi=18, mu= 0, mu^2/phi=0,    Phi= 230, G=4.5718, 1/G=0.21873
q=41:  phi=40, mu=-1, mu^2/phi=1/40, Phi= 530, G=5.0701, 1/G=0.19724
q=84:  phi=24, mu= 0, mu^2/phi=0,    Phi=2166, G=5.7653, 1/G=0.17345
q=168: phi=48, mu= 0, mu^2/phi=0,    Phi=8610, G=6.4423, 1/G=0.15522
```

This is exactly the analytic-sieve warning HYP-2899/HYP-2901 need.  Killing
primitive packets is a `Phi`-sized local-capacity problem, but upper-bound
quadratic/Selberg and large-sieve normalizations see the inverse-unit weight
`G`.  Those are different quotients.  A proof step that replaces one with the
other must carry the smoothing transform, exceptional-set boundary, and lost
denominator labels.

## Goldbach Transfer

External anchors checked during this session:

```text
Helfgott, "The ternary Goldbach conjecture is true" (arXiv:1312.7748)
Helfgott, "Minor arcs for Goldbach's problem" (arXiv:1205.5252)
Helfgott, "Major arcs for Goldbach's problem" (arXiv:1305.2897)
Kaczorowski-Perelli-Pintz, "A note on the exceptional set for Goldbach's
problem in short intervals" (Monatshefte Math. 116, 1993)
```

The useful facts are structural:

- Helfgott's theorem paper states the proof route as circle method, large
  sieve, and exponential sums, with all estimates explicit.
- The minor-arcs paper highlights reducing the cost of Vaughan's identity and
  exploiting minor-arc tails in the large sieve.
- The major-arcs paper uses explicit formulas and parabolic-cylinder estimates
  to make Gaussian smoothing usable.
- Kaczorowski-Perelli-Pintz is the relevant "Kaczorowski" Goldbach
  exceptional-set thread; the repo's older "Kaczynski" thread is the boundary
  function / curvilinear convergence warning.  They meet in LRC as:

```text
exceptional set = boundary approach class that cannot be averaged away.
```

Thus the prompt's smoothing-function changes are not cosmetic.  In this repo's
language, a smoothing function is part of the proof packet: changing it changes
which transforms, tails, zero terms, and endpoint labels are retained.

## LRC14 Use

The theorem-facing transfer is a middle certificate, not an endpoint proof.

```text
small denominator killing / lcm wall
  -> inverse-unit normalizer G(z)
  -> declared smoothing transform
  -> large-sieve or Selberg family bound
  -> twist/Fejer/exact-period witness or HYP-2908 state-lift handoff
```

This should be used at HYP-2901's Node-3 wall after fixed finite denominator
bases fail.  It does not replace:

- HYP-2978 controlled-kernel/fiber-homogeneity;
- HYP-2979 exact-period Ramanujan packets;
- HYP-2981 Fejer interval packet certificates;
- endpoint owners, safe-measure labels, C27/K33 state labels, or AP/GW
  boundary atoms.

It says where analytic number theory belongs in the proof stack: as a labelled
family estimate for late denominators and exceptional phases, with exact packet
labels still attached.

## Tournament Analysis

Vertices are analytic proof carriers, not runners:

```text
raw_prime_count
mertens_mobius_cancellation
mu_over_n_explicit_tail
squarefree_unit_normalizer_G
selberg_quadratic_upper_packet
large_sieve_minor_arc_packet
circle_method_major_minor_split
explicit_formula_smoothing_packet
labelled_lrc_interval_packet
```

Pairwise observable:

```text
retention vector =
  local data,
  Mobius sign,
  scale tail,
  unit capacity,
  quadratic kernel,
  family phase control,
  major/minor split,
  smoothing/zero packet.
```

Tie Hamiltonian path is the listed vertex order.  Stored fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
hamiltonian_paths=1
Hamiltonian path:
labelled_lrc_interval_packet >
explicit_formula_smoothing_packet >
circle_method_major_minor_split >
large_sieve_minor_arc_packet >
selberg_quadratic_upper_packet >
squarefree_unit_normalizer_G >
mu_over_n_explicit_tail >
mertens_mobius_cancellation >
raw_prime_count
```

Assumption challenge: considered vertices included primes, residues, moduli,
Mobius signs, squarefree divisors, primitive denominator packets, smoothing
functions, zero terms, exceptional sets, LRC endpoint labels, and proof
obligations.  This hypothesis chooses proof-carrier vertices.  The quotient
preserves the role an analytic estimate can play in an LRC implication; it
destroys raw integer/primality identity and must reattach exact-period,
endpoint, and packet-family labels before a theorem step.

## New Inquiry Bucket

- Compute interval-enclosed versions of the `G(z)` normalizer for the exact
  finite denominator ranges used in HYP-2972/HYP-2981.
- Test whether late q packets `{25,27,41}` from HYP-2979 have stable
  squarefree-unit normalizer profiles after endpoint labels are attached.
- Build a toy Selberg upper-bound sieve over HYP-2963 packet families with
  variables `lambda_d`, then audit quotient fibers exactly as HYP-2978 did.
- Compare Gaussian, Fejer, and compactly supported smoothing transforms by
  retained LRC labels, not just numerical decay.
- Revisit the Kaczynski boundary-function lane as the qualitative counterpart
  of Kaczorowski exceptional sets: admissible approach regions versus wild
  exceptional phases.

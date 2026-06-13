---
id: HYP-2152
status: CONCEPTUAL PROOF-PROGRAM plus S601 deterministic sample evidence
source: codex-2026-06-03-S601
related:
  - HYP-2144
  - HYP-2146
  - HYP-2145
  - HYP-2151
  - HYP-2150
  - HYP-2149
  - HYP-2142
  - HYP-2107
  - THM-398
---

# HYP-2152: determinant Helly logarithms measure certificate entropy

## Claim

In a bounded component-language proof, the logarithm attached to a Helly
argument should be read as certificate-search entropy, not necessarily as
denominator height or CRT modulus size.

For `M` live component languages and a Helly rank cutoff `H`, define

```text
B_H(M)      = sum_{1 <= h <= H} binom(M,h)
Lambda_H(M) = log B_H(M).
```

Here `B_H(M)` is the number of candidate empty subfamilies up to rank `H`, and
`Lambda_H(M)` is the proof currency paid by a naive local search or union bound.
For fixed `H`,

```text
Lambda_H(M) = H log M + O_H(1).
```

Thus the same bounded-Helly proof has different iterated-log profiles depending
only on the growth of the live component count:

```text
M(N) ~ N^a        => Lambda_H(M(N)) ~ log N
M(N) ~ log N      => Lambda_H(M(N)) ~ log log N
M(N) ~ log log N  => Lambda_H(M(N)) ~ log log log N.
```

This is the determinant-component specialization of HYP-2151's broader Helly
entropy accounting and HYP-2146's scale-currency ledger.

## S601 Evidence

`04-computation/helly_scale_log_laws_s601.py` reuses the S599 two-block
determinant classifier and computes the certificate entropy
`Lambda_h(M)` for each extracted minimal empty subfamily.

The deterministic sample used `1800` row attempts for each `n` and regime.  In
aggregate:

```text
certified_empty = 1113
high_order_empty = 0
bounded_live = 0
preempted = 16873
h=1 certificates = 1084 (0.974)
h=2 certificates =   29 (0.026)
h=3 certificates =    0
h=4 certificates =    0
```

For the `BC_only` regime after Bprime and Lemma C:

```text
n=6   h1=419 h2=8
n=8   h1=275 h2=2
n=10  h1=143 h2=0
n=12  h1=93  h2=0
n=14  h1=25  h2=0
```

For the full S598+owner stack:

```text
n=6   h1=103 h2=18
n=8   h1=22  h2=1
n=10  h1=4   h2=0
n=12  all preempted
n=14  all preempted
```

So the present n=14 determinant branch is not paying a large Helly logarithm.
It is earning a bounded local certificate dividend before the global CRT scale
appears.

## Tax Versus Dividend

A Helly logarithm can enter in two opposite ways.

**Helly tax:** if a proof has to union-bound over all candidate subfamilies, it
pays

```text
Lambda_H(M) = log sum_{h <= H} binom(M,h).
```

**Helly dividend:** if the residual almost always exposes a singleton wall or
pair incompatibility, the proof avoids the global CRT modulus.  Then the Helly
scale is a compression device: local certificates beat global residue
intersection.

S601 puts the S599 two-block determinant residual in the dividend regime.

## Harmonic Helly Ladder

There is a second, more Tao-like way for Helly scale to generate logarithms.  If
rank-`h` certificates remove `alpha/h` of the remaining residual mass, then

```text
R_H <= R_1 product_{h<=H} (1 - alpha/h)
    <= R_1 exp(-alpha sum_{h<=H} 1/h)
    <= R_1 H^{-alpha}.
```

If `H` itself is an ambient logarithmic cutoff, this produces a `log log`
factor.  If `H` is a logarithm of a logarithmic certificate count, it produces
`log log log`.

This gives a precise form of the phrase "Helly-scale logarithm": the harmonic
variable is certificate rank or certificate count, not denominator size.

## LRC Payoff

For the n=14 Lonely Runner program, the useful proof obligations are:

1. Prove the singleton determinant wall criterion from HYP-2144.
2. Prove the pair determinant incompatibility criterion.
3. Bound the live component count `M` after cap, owner, and bisection gates.
4. Use `Lambda_H(M)` as the accounting variable for any remaining finite search.

If steps 1 and 2 hold uniformly, the determinant branch is bounded-Helly and
does not need Tao-style iterated-log losses.  If a family forces `H` or `M` to
grow, HYP-2152 says exactly which logarithm is being paid.

## Relation To HYP-2151 And The Collatz Log Tower

Opus S598/HYP-2151 gives the broad Helly entropy accounting: clearance entropy,
H-entropy, and Helly number are co-monotone order parameters, and the AP
worry-set is a full-Helly isostatic trap.  HYP-2152 is narrower: it measures the
certificate entropy of the post-owner two-block determinant component languages.

HYP-2145 and the Collatz excursion-tower work (HYP-2149/HYP-2150) locate
logarithmic floors by iterated linearization of an orbit.  HYP-2152 locates a
different source of the same symbols: finite proof certificates.  In Collatz,
the log variable is orbit height or excursion length.  In the two-block LRC
residual, the log variable can be the number of local determinant rows one must
search before finding an empty subfamily.

Both are scale ladders.  They differ in what the ladder counts.

## Tournament Analysis

Vertices are proof currencies:

```text
bounded Helly certificate,
component-count logarithm,
certificate-entropy logarithm,
bounded-w window logarithm,
prime-power CRT modulus,
ambient denominator height.
```

Pair observable:

```text
(certificate locality, compression depth, search entropy).
```

Switch/gauge: earlier bounded certificates beat larger global log currencies.
The S601 fingerprint is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1}
directed_3_cycles=0
SCCs=6 singleton SCCs
Hamiltonian_paths=1
tie_HP=bounded_helly_certificate > component_count_log
       > certificate_entropy_log > bounded_w_window_log
       > prime_power_CRT_modulus > ambient_denominator_height
```

## Assumption Challenge

Candidate vertices considered: runners, safe arcs, component rows, owner pairs,
`w` residues, prime powers, determinant rows, and proof obligations.  The chosen
vertices are Helly certificate budgets `B_H(M)`.

Preserved predicate:

```text
whether a bounded common w can survive the chosen component-language subfamily.
```

Destroyed information: exact phase geometry, row labels outside the bounded
`w` window, and the full CRT modulus.  Challenged assumption: the logarithm need
not come from denominator height; it can come from the number of possible local
certificates.

## Honest Status

This is not a theorem.  The entropy identity for `Lambda_H(M)` is elementary,
and the S601 numbers are deterministic sample evidence from the S599 extractor.
The open mathematical work is to prove the singleton/pair determinant criteria
and to bound `M` or `H` in the true residual families.

**See:** `04-computation/helly_scale_log_laws_s601.py`
(+ `05-knowledge/results/helly_scale_log_laws_s601.out`),
HYP-2144, HYP-2146, HYP-2145, HYP-2151, HYP-2149, HYP-2150, HYP-2142, HYP-2107.

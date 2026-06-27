---
id: HYP-3000
title: Fibonacci path-rank additive-basis Farey bridge
status: SYNTHESIS / proof-interface carrier; not a proof
source: codex-2026-06-24-S169
script: 04-computation/fibonacci_additive_basis_farey_bridge_codex_s169.py
result: 05-knowledge/results/fibonacci_additive_basis_farey_bridge_codex_s169.out
related:
  - HYP-2998
  - HYP-2999
  - HYP-1902
  - HYP-2218
  - HYP-2219
  - HYP-2934
  - HYP-2940
  - HYP-2982
  - HYP-2984
  - LTI-150
---

# HYP-3000: Fibonacci Path-Rank Additive-Basis Farey Bridge

## Claim

The Fibonacci arrangement

```text
1, 1, 1+1, 1+2, 1+3+1, 1+4+3, 1+5+6+1, ...
```

is the rank-vector identity

```text
F_n = sum_k binom(n-k-1,k).
```

For `n >= 2`, the row entries are the numbers of independent sets of size `k`
in the path graph `P_{n-2}`.  Thus the arrangement is not just a Fibonacci
mnemonic.  It is the graded object between two old repo poles:

```text
Fibonacci path-rank polynomial      I(P_{n-2}; x)
Zeckendorf normal form              independent-set support with no adjacent carries
OCF/tournament conflict ledgers     independence polynomials on more complicated graphs
```

This gives a bridge from the additive-basis work to the Farey mutation work.
Goldbach, ternary Goldbach, Fermat polygonal numbers, Zeckendorf, and Farey
payloads are all representation-fiber systems, but they pay for coverage with
different currencies.

Rebase signal: incoming HYP-2998 covers the golden Stern-Brocot/Fibonacci spine
`F_i/F_{i+1}` and the broad representation-economy split; incoming HYP-2999
keeps the Pascal-slope row fields as labelled packet data.  This HYP-3000 lane
is deliberately narrower: it keeps the path-rank identity in the user's exact
row indexing and stresses the LRC14 unit-excess chain `p/(14p-1)` as a
proof-currency classifier.

## The Four Currencies

```text
Goldbach / Hardy-Littlewood:
  high entropy over prime representation fibers, controlled by sieve/singular
  series smoothing.

Ternary Goldbach / Helfgott:
  one extra summand gives a thicker hypergraph, so smoothing has enough room.

Fermat polygonal:
  bounded arity; a fixed number of polynomial atoms absorbs local residues.

Zeckendorf / Fibonacci:
  path-width-one carry graph; local adjacent conflicts collapse to a unique
  normal form.
```

The Fibonacci row identity is the graded count of all legal path supports
before Zeckendorf chooses a single canonical support for an integer.

## Farey Translation

For an LRC14 binding value `M=p/q`, HYP-2984 says to retain

```text
e = 14p - q
```

before applying the prompt's Farey mutations.  On the unit-excess chain
`p/q = p/(14p-1)`, the payloads become:

```text
q        = 14p - 1       binding denominator
p+q      = 15p - 1       additive / Stern-Brocot / affine-safe lane
p*q      = 14p^2 - p     product / incidence / K_{p,q} lane
q^p,p^q                  magnitude stress tests
```

The script verifies again that `p+q` is the least destructive additive payload:
it is just `q(1+M)`.  The product payload is not a scalar proof order; it is an
incidence side channel, most useful after the unit-excess gate.  The power
payloads are deliberately violent and should be used as stress tests for
magnitude-blind quotients.

## Computation

Script:

```text
04-computation/fibonacci_additive_basis_farey_bridge_codex_s169.py
```

Result:

```text
05-knowledge/results/fibonacci_additive_basis_farey_bridge_codex_s169.out
```

Exact row check:

```text
n= 5 F_n=  5 row=1+3+1
n= 6 F_n=  8 row=1+4+3
n= 7 F_n= 13 row=1+5+6+1
n= 8 F_n= 21 row=1+6+10+4
```

Finite bounded-arity check through `250` recovers the Fermat pattern for
`3`- through `8`-gonal atoms: the maximum minimum number of summands observed
equals the number of polygon sides in each case.

Goldbach sample counts show the representation-entropy contrast:

```text
binary unordered prime-pair counts: {20: 2, 50: 4, 100: 6, 200: 8}
ternary unordered prime-triple counts: {101: 38, 501: 283, 999: 769}
```

## Tournament Analysis

Vertices are proof currencies, not runners:

```text
zeckendorf_normal_form
fibonacci_path_rank
farey_sum_affine_check
fermat_polygonal_bounded_arity
ternary_goldbach_smoothing
farey_product_incidence
binary_goldbach_sieve
farey_power_stress_test
```

Pairwise observable:

```text
(quotient_safety, local_to_global, entropy_control, LRC_transfer, exact_normal_form)
```

Gauge: coordinate majority.  Tie Hamiltonian path is the listed order.

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

The transitive order is not claiming that Zeckendorf is "better" than
Goldbach.  It says that for LRC-style quotient safety, normal forms and
path-rank carriers preserve more proof data than high-entropy smoothing
carriers.  Goldbach-style smoothing remains valuable when the problem supplies
enough independent branches.

## LRC Use

This bridge suggests a compact decision rule:

```text
1. If a packet has many independent branches, use Goldbach/sieve smoothing.
2. If a fixed number of atoms can absorb every residue, use Fermat bounded arity.
3. If local repairs conflict along a path, search for a Zeckendorf/Ostrowski
   normal form.
4. If a Farey payload is used, state whether it preserves additive scale
   (`p+q`), incidence/product data (`p*q`), or only stress-tests magnitude
   (`q^p`, `p^q`).
```

For LRC14, the likely value is not a new scalar invariant.  It is a guardrail:
endpoint-debt and packet-debt should be classified by whether they are
high-entropy, bounded-arity, or path-normal-form phenomena before choosing a
certificate.

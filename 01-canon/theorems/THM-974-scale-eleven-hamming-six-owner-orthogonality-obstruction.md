---
id: THM-974
title: Scale-eleven Hamming-six owner-orthogonality obstruction
status: PROVED FINITE-EXACT — the complete 1,636,866,000-context primitive proper AP-centred common-scale-eleven Hamming-six bank reduces to 66 all-order-eleven supports; literal replay of all 66,000,000 remaining unit words gives zero survivors, independently refereed through an exact seven-orbit covariance proof
source: codex-2026-07-17-S66 scale-eleven exact C++ certificate and independent Python referee
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-962, THM-963, THM-969, THM-970, THM-972, HYP-6820]
verification:
  - 04-computation/lrc13_scale_eleven_hamming_six_owner_orthogonality_codex_c11.cpp
  - 05-knowledge/results/lrc13_scale_eleven_hamming_six_owner_orthogonality_codex_c11.out
  - 04-computation/lrc13_scale_eleven_hamming_six_referee_codex_c11.py
  - 05-knowledge/results/lrc13_scale_eleven_hamming_six_referee_codex_c11.out
---

# THM-974 — scale eleven is killed by owner orthogonality

The primitive proper AP-centred common-scale-eleven Hamming-six sheet bank is
empty.

For `c=11`, the effective states are `(1,0)` and `(11,e)` with
`e in {1,...,10}`.  Hereditary leave-one-out lcm gives 57 divisor words and

```text
sum_(k=2)^6 binom(6,k) 10^k = 1,771,500
```

state words per support, hence `924*1,771,500=1,636,866,000` labelled raw
contexts.  Scalar capacity leaves exactly 66 contexts on 66 supports, all
with six order-eleven coordinates; owner-local feasibility leaves the same
66 contexts.  They split into seven multiplication orbits of sizes

```text
12,12,12,4,12,12,2.
```

On the first 64 supports every owner obligation has size 560 and distinct
owner obligations are pairwise disjoint; the unit-word satisfaction profile
is `0:996640, 1:3360`.  The final two supports are the quadratic residues and
nonresidues.  Their six obligations pair into three equal antipodal pairs of
size 600, all cross-pairs are disjoint, and the profile is
`0:998200, 2:1800`.  Thus no `10^6` unit word satisfies all six owners.

## Completeness and independent replay

The C++ certificate literally visits all `66*10^6=66,000,000` remaining unit
words, so its zero verdict uses no symmetry quotient.  O3, O0, and
AddressSanitizer/UndefinedBehaviorSanitizer builds reproduce the frozen output
byte-for-byte.

The independent standard-library Python referee reconstructs the CRT masks by
an algebraic CRT formula, enumerates the actual divisor words, uses Python
sets for owner-local mask reachability, and proves the exact multiplication
covariance before replaying one fibre in each orbit.  The covariance is not
raw same-unit equality: multiplication permutes provider and owner
coordinates, sends the order-eleven unit `e` to `M*e mod 11` for the chosen
lift `M`, and permutes each owner's eleven sheet bits by
`z -> M^(-1)z mod 143`.  The referee checks every generator-level mask before
using this quotient.  Normal and `python -O` runs are byte-identical.

Frozen SHA-256 values:

```text
C++ source     aaa1f1c845f9203fa7d8e26fb3e1363d3a16cad1a3a298a5d0af797e9105cbb7
C++ output     f60af512554df71635d98e23871501c5c3bd2f1f06a9dc073f5021d85c61b56d
Python source  d09b1c405556efeb36c75134c73fcf693abb4188aab92b3886b537ef6e4c3ce5
Python output  3f4950d266b587df678a479f55a1151c32068ef8e1e0d17ddabecb2e75821368
```

## Tournament and vertex audit

Take owner obligations as vertices and nonempty pair intersection as the
binary observable.  Gauge a support to the lexicographically least member of
its multiplication orbit and complete ties along the sorted Hamiltonian path.
Both the generic empty intersection graph and the exceptional `3K2` graph can
then yield the same transitive tournament fingerprint: scores
`0,1,2,3,4,5`, no directed cycles, six singleton SCCs, and one directed
Hamiltonian path.  The tournament therefore destroys even the distinction
between the two terminal species.  Owner obligations—or dually their exact
unit-word sets—are the faithful vertices.

This closes only the primitive proper AP-centred common-scale-eleven
Hamming-six face under THM-860's common-sheet and hereditary-lcm conclusions.
It does not close `c>=12`, H5 ramification, non-AP-centred/deep-sheet packets,
or global sporadic emptiness.  The next common-scale frontier is `c=12`.  ∎

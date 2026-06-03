---
status: comparative hypothesis / proof-program
source: codex-2026-06-03-S580
tags:
  - lonely-runner
  - finite-sieve
  - certificate-compression
  - antipodal-witnesses
  - summand-graph
  - algorithmic-math
---

# HYP-2106: For total n=11,12,13 the paper sieve can be compressed by changing the proof object from prime fibers to certificate clocks

## Claim

The Sungkawichai-Trakulthongchai exact verifier proves total
`n=11,12,13` by certifying `J(k,p)=empty` for many primes, where
`k=n-1`.  The repo's newer LRC methods suggest a more efficient proof object:

```text
primary clocks = n-clock, (2n-1)-antipodal shells, pair-sum pinches, D/U/N endpoint owners
secondary exactness = prime-fiber J(k,p) only after those clocks fail
```

This is not yet a replacement proof.  It is a compression hypothesis: a complete
certificate quotient would avoid enumerating most or all exact prime fibers.

## Quantitative Evidence

`04-computation/lrc_paper_vs_certificate_efficiency_s580.py` compares the paper
prime schedules with the certificate-clock quotient.

```text
total n  paper k  paper raw sum p^k  full exact ansatz  cert scan ops
11       10       10^27.63           10^38.05           10^4.86
12       11       10^31.27           10^55.01           10^5.26
13       12       10^35.31           10^48.68           10^5.61
```

The full exact ansatz column is the scale of direct exact tuple checking at the
paper's terminal lift denominator.  The certificate column is a deliberately
small model: scan all `2^(n-1)` antipodal shell transversals and evaluate the
pair clocks plus D/U/N obligations.

The direct prime-lift penalties also isolate why the paper needs the polynomial
shortcut:

```text
n=11: c=11 lift is 2.53e7 times a c=2 lift for fixed input set
n=13: c=13 lift is 5.69e9 times a c=2 lift for fixed input set
```

THM-398 adds a sharper pre-filter than the S580 scan model originally counted.
For any chosen runner `v`, if the safe set of the other runners `G(S\{v})` has
a connected component longer than `2/(n v)`, then the component cannot be
contained inside one danger arc of `v`, so the whole set is loose.  Therefore a
multiple `v=nw` only reaches the hard residual when every component of
`G(S\{v})` is all-short:

```text
n=11, v=11: component cap 2/121 = 0.016529; n-clock safe gap 9/121
n=12, v=12: component cap 2/144 = 0.013889; n-clock safe gap 10/144
n=13, v=13: component cap 2/169 = 0.011834; n-clock safe gap 11/169
```

For `v=nw` the cap shrinks by `1/w`.  Thus the heavy residual is not "contains a
multiple of `n`"; it is "contains a multiple of `n` and the remaining safe set
is chopped into AP-coverable micro-components."

S573/HYP-2104 names this same split the Vitali handoff: `n !| v` lives on the
constructive `t=1/n` side, while `n | v` creates the periodic Vitali-cover
family that measure can attack.  Its sampled mult-of-`n` computation reports
that Criterion B' already proves `72.4%, 78.7%, 88.9%, 91.5%, 96.8%` at
`n=6,8,10,12,14` respectively, with zero all-short tight cases in the sample.
For the paper comparison, this says the `n=12` multiple branch is already
mostly removed by interval geometry before shell/pair/D/U/N certificates fire.

## Rotated Arithmetic

The paper's decisive denominator is `n=k+1`.  Its polynomial method likes `n`
odd prime; the composite case `n=12` is handled by small-factor lifts
`2,2,2,2,3,3`.

The repo's antipodal/summand denominator is instead

```text
C = 2n - 1.
```

This rotates which cases are arithmetically clean:

```text
n=11: C=21 = 3*7, unit/nonunit shells = 6/4
n=12: C=23 prime, unit/nonunit shells = 11/0
n=13: C=25 = 5^2, unit/nonunit shells = 10/2
```

Thus total `n=12` is the most attractive certificate case even though it is the
paper's composite-lift case: every `C=23` shell is a unit shell, so the
unit-spine residual is zero.

## Clocks That Matter First

1. `n`-clock.  HYP-2102 and THM-398 show the global split: no multiple of
   `n` gives the immediate witness `t=1/n`; the whole residual is the C'
   branch `n | v -> positive measure`.  Incoming HYP-2103 sharpens the large
   branch from "large multiple" to general dominance dodge.
2. `2n-1` antipodal shell clock.  The strict sub-edge witness threshold is
   organized by shells `{a, C-a}`.  Unit shells are visible under multiplication;
   nonunit shells are the residual.
3. Pair-sum/pinch clocks.  HYP-2095 says the hard measure-zero rows are expected
   to have an unblocked small reduced-sum pair, so a pinch witness fires before
   exact grid lifting.
4. D/U/N endpoint-owner obligations.  These preserve the finite-cover content
   that raw tuple symmetries destroy.

## Clocks To Defer

Generic `lp` ansatz denominators, arbitrary time grids, unmarked runner
tournaments, and full prime fibers should be treated as fallback clocks.  They
are exact, but they do not preserve the proof predicate as cheaply as the
endpoint-labelled clocks above.

## Algorithmic Directions

- Insert a THM-398 interval oracle before the certificate scan: build
  `G(S\{v})` by endpoint sweep, discharge every long component, and pass only
  the all-short AP-cover residual to pair/D/U/N machinery.
- Generate endpoint-obligation cores first, then solve realizability over
  speeds; do not generate `p^k` residue tuples first.
- Build a ZDD/BDD over D/U/N obligations and pair-pinch blockers; use prime
  lifts only as clause generators for the obligation automaton.
- Turn every failed paper seed family into a labelled certificate germ
  `(shells, pair blockers, endpoint owners)` and learn cross-prime clauses.
- Prove the HYP-2102 small-multiple residual by an arc-equidistribution
  non-covering lemma.  This would bypass the whole product-sieve proof, not only
  speed it up.
- Use the S573 Vitali framing to make the equidistribution lemma quantitative:
  a hypothetical counterexample is an arc-alignment of short rational intervals
  against the `1/(nw)` lattice, not an arbitrary cover.  Erdos-Turan or
  three-distance bounds should attack this alignment directly.
- Rephrase the all-short residual as an endpoint-owner matching problem: each
  surviving safe component would need to fit inside a single arc of the AP
  `{m/(nw)}`; its two endpoints are owned by other runners, so a cover imposes
  rigid pair-sum/divisor/unit congruences rather than a free time-grid search.
- Reuse HYP-2111's gate vocabulary for the fallback order: the paper-sieve
  residual is `G3_all_short_Cprime_residual`, and the desired replacement for
  exact prime fibres is `G8_endpoint_cover_circuit_positivity`.
- Add HYP-2107's bounded multiplier-residue automaton as the next fallback
  before exact fibres: large-owner cover congruences should be intersected in
  `w` and pruned by dominance.
- For `n=12`, exploit `C=23` prime directly: the unit action is transitive on
  all antipodal shells, so any counterexample-like object must come from
  endpoint owners rather than nonunit shell holes.

## Tournament Analysis

For this comparison, vertices should not be runners.  The useful tournament
vertices are proof clocks or certificate quotients:

```text
n-clock, C-shells, pair-pinch, D/U/N owners, endpoint core, exact prime fallback
```

Observable: residual log-work after applying the clock while preserving the
predicate "does a witness/certificate exist?"  Switch: orient toward the larger
residual burden.  The fingerprint is transitive in the S580 model, so cycles
should be sought only after splitting endpoint-owner fibres or cross-prime
learned clauses.

## Status

Supported as a methodology comparison by S580 scale computations and by
HYP-2095/HYP-2100/HYP-2102 plus THM-398/HYP-2103's dominance-dodge
formalization.  Open as a proof: we still need completeness of the certificate
quotient, especially the all-short small-multiple non-covering lemma and a
realization theorem for endpoint-owner cores.

## Sources

- Sungkawichai and Trakulthongchai, "Eleven, twelve, and thirteen lonely
  runners", arXiv:2604.23906, https://arxiv.org/abs/2604.23906
- Public verifier repository, https://github.com/vzsky/13-lonely-runners
- `04-computation/lrc_paper_vs_certificate_efficiency_s580.py`
- `07-reflections/lrc-paper-sieve-vs-certificate-compression-s580.md`
- HYP-2111, n=14 certificate calculus / local certificate primitives
- HYP-2108, endpoint-cover circuit positivity
- HYP-2107, large-owner residue automaton
- THM-398 / HYP-2103, dominance-dodge and all-short residual formalization
- HYP-2104, Vitali handoff / arc-alignment residual

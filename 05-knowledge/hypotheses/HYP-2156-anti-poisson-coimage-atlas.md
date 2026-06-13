---
id: HYP-2156
status: OPEN synthesis program; anti-Poisson/coimage atlas for LRC collapse, strong tournament residuals, and adjacent arithmetic systems
source:
  - codex-2026-06-03-S604
  - codex-2026-06-03-S605
related:
  - HYP-2153
  - HYP-2154
  - HYP-2155
  - HYP-2151
  - HYP-2152
  - HYP-2144
  - HYP-2146
  - HYP-2145
  - HYP-2143
  - HYP-2114
  - HYP-2104
  - HYP-2135
  - HYP-2157
  - THM-401
  - THM-406
---

# HYP-2156: anti-Poisson coimage atlas

## Core Claim

The anti-Poisson frame is the recurring residual pattern where a free or
independent pushforward baseline predicts a positive ground cell, but arithmetic
or structural correlations force exact all-orders cancellation in the coimage
while preserving a canonical witness floor.

For LRC this is now precise:

```text
free baseline:       independent sparse danger arcs, depth ~ Poisson(lambda)
true coimage:        p_k = meas{t : depth_delta(t)=k}
ground cell:         p_0 = lonely measure
anti-Poisson locus:  p_0 = 0 by full inclusion-exclusion
witness floor:       boundary unit skeleton at t = a/n
```

THM-406 is the key upgrade.  It proves

```text
p_0 = sum_j (-1)^j S_j,
S_j = total j-fold overlap volume,
E[binom(N,j)] = S_j.
```

So anti-Poisson collapse is not "low variance" or "large second moment."  It is
an all-orders overlap cancellation in the spectral/coimage measure.  S599c also
shows that `S_2` does not separate additive chains from generic controls.  The
right object is therefore the whole coimage distribution, not any finite moment
truncation.

## Category / Number-Theory Addendum

The categorical reading and the number-theoretic reading are complementary, not
parallel metaphors.

The categorical side:

```text
N_delta : clock circle -> depth values,
coimage(N_delta) = the minimal quotient retaining the depth predicate,
Yoneda = the coimage is recognized by all natural probes into it.
```

Those probes include the depth distribution, factorial moments, spectral
measure, Helly certificates, partition-function evaluations, CRT witnesses, and
unit-shell observers.  If all these probes keep recovering the same obstruction,
the object is canonical in the Yoneda sense: it represents the problem's
observable functor.

The number-theoretic side says which probes matter at the LRC floor.  At
`delta=1/n`, the sharp resonance modulus is

```text
C = 2n - 1.
```

A missed unit antipodal shell `{a,-a} mod C` gives the inverse clock
`k=a^{-1}` and therefore a witness at `k/C` with gap `2/C > 1/n`.  Thus any
`p_0=0` collapse must pass through the `C`-shell quotient:

```text
V mod C -> antipodal shells -> unit action by (Z/CZ)^x -> gcd strata.
```

This is the THM-401/S571 witness functor.  The `2n-1` resonances are the
arithmetic part of the all-orders cancellation.  They do not replace
`sum_j (-1)^j S_j`; they compress the conditions under which that alternating
sum can vanish.  Prime `C` makes every shell unit-visible.  Composite `C`
creates nonunit lanes where sporadic collapse can hide, such as `C=27` for the
`n=14` route.

## Dictionary

| thread | free or Poisson-like baseline | anti-Poisson residual | coimage / master object | surviving floor |
|---|---|---|---|---|
| LRC covering depth | independent danger arcs, `Poisson(2m delta)` | AP and sporadic additive chains force `p_0=0` | depth law `{p_k}` | unit boundary witnesses |
| Vitali wall | measure/Bonferroni truncations | no finite overlap order decides `p_0` | all factorial moments / spectral measure | exact boundary shell |
| Helly entropy | small subfamily certificates | full-Helly or isostatic rows | minimal covering obstruction | one-stage singleton wall or full cascade |
| two-block determinant | independent component languages | CRT intersections empty by aligned determinant rows | allowed-`w` language intersection | singleton/pair determinant certificates |
| Collatz cycles | drift/density model | cycle equation `2^E - 3^k` is an exact correlated tail | rapidity/linear-form coimage | trivial `3+1=4` floor |
| tournament OCF | random conflict graph / hard-core gas | forbidden `H`, real-root failures, SCC obstructions | partition function / path-homology signature | Hamiltonian path tie order |
| additive bases | high representation entropy | canonical support, carry constraints, thin exact covers | representation-count law | unique/canonical expansion |
| covering systems | random residue coverage | exact cover or no-cover at residue level | residue-depth distribution | uncovered residue class |

This dictionary is not claiming these systems are identical.  It claims that
their difficult residuals share the same shape: the low-order or average
observable sees a free baseline, while the exact coimage sees a structured
ground-cell cancellation.

## Strong Tournaments

The strong-tournament connection is best read after the coimage choice, not
before it.  A strong tournament is the tournament residue that remains when no
single scalar ranking, source, sink, or cheap Hamiltonian elimination order
explains the quotient.  In that sense, strong tournaments are the tournament
analogue of the anti-Poisson residual:

```text
transitive quotient      -> sorted certificate order, finite elimination works;
strong SCC quotient      -> cyclic proof obligations, no one-dimensional rank works.
```

For LRC this matters only if the tournament vertices are chosen correctly.  The
vertices should not automatically be runners or raw arcs.  Useful vertex sets
include:

```text
coverage-depth cells,
overlap orders,
unit/nonunit shell gaps,
section boundaries,
wall-crossing events,
residue classes,
two-block determinant component languages,
cover arcs,
Fourier modes,
matroid circuits,
proof obligations.
```

The preserved predicate is the LRC residual decision: can the quotient certify
positive lonely measure, or can it force `p_0=0` while preserving the boundary
floor?  The destroyed information is endpoint phase order and sometimes the
exact geometry of intervals.  That loss is acceptable when the quotient is used
as a proof-obligation tournament; it is not acceptable if the goal is to
recover the full covering arrangement.

Strong SCCs in such a quotient should be treated as scarce proof targets.  A
transitive proof-lens tournament says "there is a priority order of exits."
A strong proof-lens tournament says "the exits mutually depend on each other,"
which is exactly where anti-Poisson all-orders cancellation can hide.

## Working Formal Subproblem

Attach to a candidate row or residual branch its anti-Poisson signature:

```text
APSig(V) =
  baseline law b_k,
  exact coimage p_k,
  overlap sequence S_j,
  defects D_j = S_j - S_j(free),
  ground-cell value p_0,
  witness floor W,
  Helly rank h,
  tournament fingerprint Tau.
```

Here `Tau` should include score histogram, SCCs, directed 3-cycles, edge flips,
and at least one tie Hamiltonian path for the proof-obligation tournament.  If
the tournament quotient is transitive, the branch likely has a cheap ordered
certificate.  If it has a strong SCC, the branch is an anti-Poisson candidate.

For the current `p_0=0` additive-chain family, the open classifier becomes:

```text
p_0 collapse at delta=1/n
  iff
two-seed additive generation
  + unit-boundary witness floor
  + THM-401 shell compatibility
  + all-orders overlap cancellation.
```

The "iff" is intentionally open.  S602 finds that every collapsed row in its
targeted boxes is a two-seed addition chain, but many chains are false
positives before the shell and floor filters.

## N=14 Payoff

For the `n=14` LRC route, the anti-Poisson frame says to stop asking only
whether a row is AP-like.  The row should be projected through:

```text
coimage depth law,
Res_27(V): unit/nonunit C=27 shell floor,
additive-chain labels,
two-block determinant component languages,
Helly certificate rank,
proof-obligation tournament SCCs.
```

The difficult branch should be exactly the one where these quotients retain a
strong cyclic core after cheap exits.  If every quotient becomes transitive or
small-Helly, then the anti-Poisson threat has been localized.  If a strong SCC
survives, it becomes the next proof object rather than a search artifact.

## Concrete Tests

1. Extend `lrc_p0_collapse_additive_chains_s602.py` with an `APSig` table:
   `p_k`, `S_j`, alternating partial sums, shell floor, and proof-lens SCCs.
   Add `Res_C(V)` columns: antipodal shell coverage, unit-character orbits, gcd
   strata, and additive-chain relations.
2. Compare AP, `(1,3,4,7)`, `(1,3,4,5,9)`, the two `n=8` sporadics, and
   non-chain controls by their full inclusion-exclusion profiles, not only
   `S_2`.
3. Build a determinant analogue: component-language depth distribution,
   singleton/pair Helly rank, and a proof-obligation tournament whose strong
   SCCs are the live CRT intersections.
4. Build a tournament OCF analogue: treat the hard-core/independence
   partition function as the coimage and look for forbidden `H` values as
   ground-cell cancellations against the free conflict-graph baseline.
5. Keep a challenged-assumption block in every exploratory script: list the
   vertex sets considered, the LRC predicate preserved, and the information
   destroyed by the quotient.

## Status

This is a synthesis hypothesis, not a theorem.  Its solid anchors are:

```text
THM-406: p_0 is exact all-orders inclusion-exclusion and {p_k} is spectral.
HYP-2153: the p_0 collapse family is larger than AP and contains sporadic chains.
HYP-2154: the free/depth baseline is Poisson-like.
HYP-2155: the master object is a coimage.
HYP-2151/HYP-2152: Helly rank is certificate entropy.
THM-401/S571: `2n-1` unit-shell resonances give witness exits.
```

The new contribution is the named proof program: anti-Poisson means
coimage-level ground-cell collapse by structured correlation, and strong
tournament residuals are the cyclic proof-obligation subset where that
collapse can hide.  S605 sharpens the number-theory layer: `2n-1` resonances
are the shell/character probes whose coverage decides whether the all-orders
cancellation is even possible.

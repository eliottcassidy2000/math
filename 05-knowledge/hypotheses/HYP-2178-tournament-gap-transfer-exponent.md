---
id: HYP-2178
status: SUPPORTED synthesis and proof-program; no new forbidden H value claimed
source: user-2026-06-03; codex-2026-06-03-S615
related:
  - HYP-2176
  - HYP-2175
  - HYP-2171
  - HYP-2157
  - HYP-2155
  - HYP-2154
  - HYP-2144
  - THM-002
  - THM-029
  - THM-079
tags: [tournaments, forbidden-H, coimage, unit-distance, LRC, Collatz, exponent, tournament-analysis]
---

# HYP-2178: tournament forbidden-H gaps are transferable coimage certificates

## Claim

The hard shared object behind tournament Hamiltonian paths, LRC floor rows,
unit-distance grid disproofs, and Collatz/two-block residues is an
evaluation-spectrum obstruction:

```text
visible quotient asks for an evaluation state
but the retained carrier forces extra compatible structure
so the requested state is unavailable.
```

For tournaments this is literal.  By OCF,

```text
H(T) = I(Omega(T), 2),
```

and the proved permanent gaps `H=7` and `H=21` are unavailable evaluations of
the odd-cycle conflict graph.  The transfer principle is that these gaps should
be read as upper-cap certificates: if a coimage route in another problem
demands the analogue of an unavailable tournament state, then the quotient has
forgotten a forced side-channel.

This hypothesis does **not** claim a new permanent forbidden `H` value.  The
repo guardrail remains: only `7` and `21` are proved permanent gaps; `63` is
achieved at `n=8` and is not a Mersenne-pattern continuation.

## Why This Propagates

The same proof shape appears in four domains:

- Tournament H-spectrum: `H=7` and `H=21` are blocked because compatible
  odd-cycle packets cannot remain exactly as sparse as the scalar evaluation
  demands.
- LRC: `p_0=0` floor rows are all-orders depth cancellations; finite moments
  and low-order density cannot decide the wall, so `Res_(2n-1)` must retain
  owner labels, carry cocycles, pinch witnesses, and endpoint certificates.
- Unit distance: graph-only or visible-grid quotients forget geometry.  At
  small `n=22`, the graph coimage reaches `62` but geometry/unfaithfulness kills
  the candidates; asymptotically, Sawin's explicit `n^1.014` construction uses
  a high-degree arithmetic carrier whose plane projection is only the coimage.
- Collatz/two-block residues: residue density is blind when carries are
  correlated, so the carrier must retain valuation/carry history and determinant
  compatibility.

The shared abstraction is therefore not "everything is a tournament on the
original vertices."  Tournament Analysis should challenge the vertex set.  In
S615, the vertices are transfer mechanisms and proof obligations; alternatives
such as runners, residues, unit-distance points, cover arcs, Fourier modes, and
wall-crossing events are demoted unless they preserve the predicate.

## The 1.014 Exponent

The unit-distance side is normalized: Sawin proves sets of arbitrarily large
`n` planar points with more than `n^1.014` unit-distance pairs, by changing the
carrier from a visible grid to arithmetic number fields with many small split
primes.

The tournament side is not yet normalized.  Raw unlabelled tournament growth
`A000568` and raw scalar `H` comparisons are the wrong coordinates.  If the
user's shared `1.014` is real, it should appear as an entropy dividend or
deficit of the feasible OCF/proof-obligation carrier after forbidden packets and
side-channels are retained:

```text
delta = normalized growth of certificate-bearing carrier
        minus visible/coarse quotient growth.
```

This is the fixed-exponent cousin of the Tao-style log/loglog/logloglog
abstraction already in HYP-2146 and HYP-2152: one must choose the scale ladder
whose product of savings is being measured.  Logs appear when savings thin over
scales; a positive exponent such as `0.014` appears when the retained carrier
keeps a positive entropy dividend across a tower.

## S615 Evidence

`04-computation/tournament_gap_transfer_exponent_s615.py` records:

- the canonical `H=7`/`H=21` guardrail and OCF arithmetic decompositions;
- the cross-domain carrier ledger for tournaments, LRC, unit distances, and
  Collatz/two-block residues;
- a Tournament Analysis over transfer mechanisms, with vertices chosen as proof
  routes rather than runners or points;
- a normalization audit separating Sawin's unit-distance theorem from the
  still-unmeasured tournament-side exponent.

The mechanism tournament is transitive.  Its ranking is:

```text
coimage side-channel retention
> H-gap completion forcing
> two-block determinant walls
> exponent entropy normalization
> Collatz carry-residue transfer
> raw scalar H matching.
```

The transitivity is useful: it says the program should keep side-channels and
completion laws first, then measure exponents, and only then look for scalar
coincidences.

## Next Proof Probes

1. Complete-core `H=21`: prove the `r=10` single-core signature gap or explain
   why THM-079 requires non-core decomposition cases.
2. Unit-distance transfer: for each killed graph-only state, identify the
   missing side-channel, not just the missing edge count.
3. LRC transfer: treat `Res_(2n-1)` floor atoms as forbidden-evaluation tests;
   a lift must preserve owner/carry/pinch probes Yoneda-style.
4. Exponent test: define a common carrier-size variable and measure the
   tournament entropy dividend against Sawin's `n^1.014` unit-distance side.

## Artifacts

- `04-computation/tournament_gap_transfer_exponent_s615.py`
- `05-knowledge/results/tournament_gap_transfer_exponent_s615.out`
- `07-reflections/tournament-gap-transfer-exponent-s615.md`

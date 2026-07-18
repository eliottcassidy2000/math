---
id: THM-988
title: Scale-twenty-one Hamming-six owner-feasibility deficit
status: PROVED STRUCTURAL + FINITE-EXACT — independent literal-search C++ and algebraic-CRT Python certificates exhaust the 77,810,217,408-context bank; an exact cubic-character argument explains the one-sheet deficit on the two scalar-tight rows
source: codex-2026-07-17-S66 scale-twenty-one continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-982, THM-983, THM-986]
related: [THM-978, THM-980, THM-981, THM-989, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twenty_one_hamming_six_frontier_scout_codex_c21.cpp
  - 05-knowledge/results/lrc13_scale_twenty_one_hamming_six_frontier_scout_codex_c21.out
  - 04-computation/lrc13_scale_twenty_one_hamming_six_referee_codex_c21.py
  - 05-knowledge/results/lrc13_scale_twenty_one_hamming_six_referee_codex_c21.out
  - 04-computation/lrc13_scale_twenty_one_cubic_character_overlap_codex_c21.py
  - 05-knowledge/results/lrc13_scale_twenty_one_cubic_character_overlap_codex_c21.out
---

# THM-988 — scale twenty-one has at least two impossible owners

The primitive proper AP-centred common-scale-twenty-one Hamming-six sheet bank
is empty before the global covariant unit-fibre or metric-continuation stages.

For `c=21`, the effective orders are `1,3,7,21`, with twenty-one literal
`(D,e)` states.  Hereditary leave-one-out lcm enumeration gives 3,249 divisor
words and 84,210,192 literal state words per support, hence

```text
924*84,210,192 = 77,810,217,408
```

unquotiented labelled state contexts.  Both certificates leave 350 scalar
contexts on 80 supports across nine multiplicity profiles.  Exact owner-local
reachability gives the census

```text
0 feasible owners: 182 contexts,
1 feasible owner :  96 contexts,
4 feasible owners:  72 contexts.
```

Thus every scalar row has at least two empty owner projections.  The 2,100
owner rows have maximum-union histogram

```text
17:1056, 18:648, 20:12, 21:384.
```

The two all-order-twenty-one scalar rows are especially revealing.  All six
owners are exactly scalar-capacity tight, yet every owner-local union has
maximum twenty.  This is a pure overlap/projection obstruction, invisible to
scalar slack: every owner owes exactly one sheet, but the owed sheet is visible
only before cardinality compresses the mask family.

## Symbolic overlap obstruction on the two tight rows

The two all-order supports are the quadratic-residue cosets

```text
Q  = {1,3,4,9,10,12},
2Q = {2,5,6,7,8,11}
```

in `F_13^*`.  For an owner `o`, the affine sheet coordinate
`t=o^-1+13s (mod 21)` is a bijection.  Both supports then have ratio set `Q`.
Writing `u=e^-1` for a provider's order-twenty-one unit, its mask is `u A_r`,
where the strict-left/closed-right interval gives

```text
A_1  = {0,3,8,16}       A_3  = {6,11,19}
A_4  = {1,6,14,19}      A_9  = {2,7,15,20}
A_10 = {2,10,15}        A_12 = {5,13,18}.
```

Their total cardinality is exactly twenty-one, split into one zero, two
gcd-seven residues, six gcd-three residues and twelve units.  Consequently a
full union would have to be a disjoint partition preserving each of these
multiplicative strata.

Let `chi : F_7^* -> {1,omega,omega^2}` be the cubic character with
`chi(3)=omega`.  Every `u A_r` contains one gcd-three point `3d_r` and two
units.  Exact Eisenstein arithmetic gives the unit-pair character sum

```text
-omega chi(d_r)    for r in {+/-1,+/-3},
-omega^2 chi(d_r)  for r in {+/-4}.
```

If the six masks partitioned `Z/21Z`, let `S_1,S_3,S_4` be the sums of the two
selected `chi(d_r)` marks in the three opposite-ratio pairs.  Partitioning the
gcd-three stratum and the unit stratum would respectively force

```text
S_1 + S_3 + S_4 = 0,
-omega(S_1 + S_3) - omega^2 S_4 = 0.
```

Hence `(omega-omega^2)S_4=0`, so `S_4=0`.  This is impossible: `S_4` is a sum
of two cube roots of unity, and no two cube roots are negatives.  Thus every
local union has size at most twenty.  The explicit multipliers
`(2,17,2,17,2,2)` cover every residue except one, so the maximum is exactly
twenty.  This proves the tight-row deficit symbolically; the exhaustive
certificate below checks the normalized-mask derivation and every exact
character identity without floating-point arithmetic.

The same verifier classifies all `12^6=2,985,984` multiplier tuples.  Exactly
672 attain size twenty, producing the eighteen complements
`Z/21Z - {t}` with `7` not dividing `t`; residues `0,7,14` are never missed.
Of the maximizers, 288 miss a unit `x` and duplicate `4x`, while 384 miss a
gcd-three residue `x` and duplicate one of `3x,4x`.  In the original owner
coordinate the missed sheet therefore satisfies
`7` not dividing `o^-1+13s`.

## Independent replay

The primary certificate constructs CRT representatives by literal search,
checks each sheet cardinality against an independent period formula, and uses
a packed exact union-mask table.  The Python referee was derived before reading
the primary.  It instead solves CRT algebraically, enumerates Cartesian grammar
fibres, stores each reachable owner bank as an immutable set, and hashes every
mask in increasing order.  After freezing its answer, it replayed the C++ FNV
serialization and matched the primary's mask, order, multiplicity, scalar-bank,
capacity-vector, row-capacity, and owner-profile digests exactly.

The implementations agree on all decisive counts: the nine multiplicity
profiles, support-context histogram `1:14,2:48,8:12,24:6`, multiplication
orbits `2:1,6:1,12:6`, 131 capacity vectors, 34 owner maximum vectors, all
owner histograms above, and zero all-six contexts.

| artifact | SHA-256 |
|---|---|
| `04-computation/lrc13_scale_twenty_one_hamming_six_frontier_scout_codex_c21.cpp` | `08772ed1a3e5acf59eb1f33ac0420de6394fa6a7648b8c6efcdb79a3d36268bc` |
| `05-knowledge/results/lrc13_scale_twenty_one_hamming_six_frontier_scout_codex_c21.out` | `65a55c93aaea9d53ee98ff2d3b36fa61afcc4511eecf90904d621a415e4fb284` |
| `04-computation/lrc13_scale_twenty_one_hamming_six_referee_codex_c21.py` | `14a5d181a93befc5516711a3340b1ceaa4899d3e478d649cb4504a7e3b027ae7` |
| `05-knowledge/results/lrc13_scale_twenty_one_hamming_six_referee_codex_c21.out` | `fbf1680761bf9d42e2b95e4ca4f37e239a7e9840b40fad46816870c641b1987b` |
| `04-computation/lrc13_scale_twenty_one_cubic_character_overlap_codex_c21.py` | `159f34f851909f0545da6ba45f7e87147adea781b90373f21f253cef55fd709f` |
| `05-knowledge/results/lrc13_scale_twenty_one_cubic_character_overlap_codex_c21.out` | `3b6f5ffbcb1764a83684766db4335c50921a282029ff55ab5a6650cee89dfdcc` |

Optimized, unoptimized, and ASan/UBSan C++ outputs agree byte-for-byte.  Normal
and `python -O` referee outputs also agree byte-for-byte.

## Faithful carrier and terminal implication

For the complete bank, the faithful terminal carrier is the labelled owner
capacity/feasibility/max-union vector.  On the two tight rows its explanatory
refinement is a colored exact-cover/group-ring system: six gcd-three marks
coupled to six two-point unit edges.  The three vertices obtained by pairing
opposite ratios have coefficients `(-omega,-omega,-omega^2)`; their completed
tournament is acyclic and forgets the additive character cancellation that
does the work.  Owner, provider, individual-residue, gap, boundary, wall-event
and cover-arc vertices likewise lose either the multiplicative stratum or the
shared multiplier incidence without this sidecar.

For every context a global unit word would induce a feasible choice in all six
owner projections, so the exact maximum of four is terminal.  The converse is
not used: separately feasible owner witnesses need not glue to one global
word.  This theorem concerns only the primitive proper AP-centred
common-scale-twenty-one Hamming-six face.  THM-989 separately closes scale 22
and THM-990 claims scale 24; the H5 bank, non-AP/deep-sheet continuations, and
global sporadic emptiness remain open.

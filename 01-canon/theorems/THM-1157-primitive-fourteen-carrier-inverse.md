---
id: THM-1157
title: Primitive q=14 carrier inverse and the three-exit LRC14 trichotomy
status: PROVED FORMAL REDUCTION.  Under no Lonely13, Covering through 14 is equivalent to the single q=14 carrier bit; the primitive Covering supplier, primitive carrier supplier, and three-exit trichotomy are equivalent.  Given the carrier supplier and the cited LRC through 13, gcd normalization and dilation prove the global LRC14 statement.  The structural carrier supplier itself, crown collapse, and n=12 equality rigidity remain OPEN
source: codex-2026-07-18-S67 continuation
depends_on: [THM-1158, THM-1149]
related: [THM-1131, HYP-4382, MISTAKE-170]
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCPrimitiveCarrierINV.lean
---

# THM-1157 -- the primitive inverse has one carrier bit

For a positive thirteen-speed family `v`, write

```text
L13(v)  := exists t, Lonely 13 v t,
C14(v)  := exists i, 14 divides v_i,
D13(v)  := exists j, for every i != j, 13 v_i <= v_j.
```

The literal `INVcov` was refuted by the doubled AP because Covering is not
preserved when the common gcd is divided out.  The correct normalized supplier
does not transport Covering.  It retains the only divisor obligation that is
not already forced at the `1/13` threshold.

## 1. The q=14 atom

> **Theorem A.**  If `not L13(v)`, then
>
> ```text
> Covering(2..14,v)  iff  C14(v).                       (1)
> ```

For every `2<=q<=13`, absence of a `Lonely 13` time and the denominator sieve
force some speed to be divisible by `q`: otherwise `t=1/q` is already a
`Lonely 13` time.  Hence only `q=14` remains.  The reverse implication in (1)
combines those twelve forced obligations with the supplied 14-carrier.

This exactly identifies what dilation destroyed in THM-1158.  The doubled AP
has a 14-carrier, but its primitive core `{1,...,13}` does not.  The obligations
`q=2,...,13` survive because they are forced again by `not L13`; the q=14 bit
must be tested after normalization.

## 2. The honest primitive supplier

Define

```text
PrimitiveCarrierINV:
  positive(u) and gcd(u)=1 and C14(u) and not L13(u)
  => D13(u).                                            (2)
```

The operational form is the three-exit trichotomy

```text
positive(u) and gcd(u)=1
=> L13(u) or not C14(u) or D13(u).                      (3)
```

> **Theorem B.**  `PrimitiveCarrierINV` is equivalent both to (3) and to the
> primitive version of `INVcov`, in which `Covering(2..14)` is assumed.

The equivalence with (3) is classical case analysis.  The equivalence with the
primitive Covering form is precisely Theorem A.  Unlike `ResidualINV`, (2) is
not obtained by replacing `not L13` with the actual no-`Lonely14` branch.  It is
a genuinely stronger structural supplier and therefore a noncircular target.
It remains open.

## 3. Global normalization and the three exits

Let `v` be any positive family, put

```text
g = gcd(v),             u_i = v_i/g.
```

Then `g>0`, `u` is positive and primitive, and `v=g*u`.  Apply (3).

1. If `u` is `Lonely 13` at `s`, band monotonicity makes it `Lonely 14` at
   `s`; dilation makes `v` lonely at `s/g`.
2. If `u` has no 14-carrier, the exact sieve time is `s=1/14`; the transported
   time for `v` is `1/(14g)`.
3. If `u` has a 13-dominant speed, `ap_core_bridge` and the cited LRC through
   13 give a `Lonely 14` time `s`; dilation again transports it to `s/g`.

> **Theorem C.**  The cited LRC through 13 together with
> `PrimitiveCarrierINV` implies both the positive working LRC(14) statement and
> the signed/nonzero `LRC14Statement`.

Lean proves the normalization through `tupleGcd`, `primPart`, and
`exists_lonely_of_dvd`.  It also records the calibration

```text
not exists Lonely14(v) => Covering(primPart(v)),         (4)
```

because the no-`Lonely14` premise transports through dilation and the
q=2,...,14 sieve may then be rerun.  Replacing `14` by `13` in (4) is invalid:
only q through 13 is recovered.  Thus the safe order is **normalize first,
then split on q=14**, never “transport Covering through division.”

## 4. Exact connection to the essential crown

A failure of (3) is exactly a positive primitive family satisfying

```text
not L13(u),       C14(u),       not D13(u).              (5)
```

By Theorem A it is a Cover14 strict cover.  For pairwise-distinct speeds,
`not D13` is the compact top-ratio condition.  THM-1149 therefore places (5)
in the tight-deletion/all-loose essential-crown dichotomy.

The missing bridge can be stated without another global reformulation:

```text
CrownCollapse14:
  every family in (5) has some deletion W with M(W)=1/13.
```

Assume `CrownCollapse14` and the open n=12 equality classification

```text
M(W)=1/13  =>  W=d*{1,...,12}.                          (R12)
```

THM-1149's Farey regeneration gives `13d | w` for the deleted speed `w`.
Primitivity of the full family forces `d=1`.  Since `{1,...,12}` has no
14-carrier, (5) then forces `14|w`; hence `182|w`.  Consequently

```text
w >= 182 > 13*12,
```

which is `D13`, a contradiction.  Therefore

```text
CrownCollapse14 + R12 + Farey regeneration
  => PrimitiveCarrierINV => LRC14.                      (6)
```

The genuinely missing crown statement is exclusion of the **all-loose
Cover14 crown**.  The exact non-Cover14 all-loose family `V0` from THM-1149 is
not a counterexample to (2): it lies in the `not C14` exit and is dispatched at
`t=1/14`.  This is why the carrier bit and the integer lift sidecar are
load-bearing.

## 5. Tournament and alternate-vertex audit

The runner-order tournament remains unfaithful.  `{1,...,13}` and its doubled
dilation both produce the same transitive tournament: score histogram
`0,1,...,12`, no directed cycle, thirteen singleton SCCs, and one Hamiltonian
path.  It cannot see the q=14 carrier that changes the truth of Covering.

The smallest proof-preserving normalized object is instead

```text
primitive scale orbit
  + q=14 carrier mask
  + the three labelled exits {L13, no-C14, D13}.
```

If these three proof obligations are forced into a tournament, orient them by
discharge order

```text
L13 -> no-C14 -> D13,
```

with the displayed labels as the tie switch.  Its fingerprint is the
transitive three-vertex tournament: scores `0,1,2`, no cycle, singleton SCCs,
and one Hamiltonian path.  This preserves certificate flow, not the full
runner geometry.

On the crown residual, even that object is too coarse.  Faithful vertices are
deletion/private-stalk obligations, with hyperedges for multiplicity-at-least-
three chambers and sidecars carrying the owner tooth, q=14 carrier, deletion
margin, and integer lift.  Pairwise midpoint orientation again becomes
transitive and destroys the triple incidence that pays for private mass.  The
challenged assumption is therefore explicit: neither runners nor pairwise
stalk comparisons are the underlying vertices of the hard object; the carrier-
labelled chamber hypergraph is.

## 6. Formal audit

`LRCPrimitiveCarrierINV.lean` proves Theorems A--C, both supplier
equivalences, gcd/loneliness invariance specialized to `primPart`, and (4).
It contains no `sorry`, `admit`, or `native_decide`; public theorem axioms are
audited in the file.

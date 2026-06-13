---
id: HYP-2235
status: OPEN finite-field carrier transfer hypothesis with exact S659 scout
source: codex-2026-06-05-S659
related:
  - HYP-2234
  - HYP-2233
  - HYP-2231
  - HYP-2230
  - HYP-2222
  - HYP-2197
  - HYP-2176
  - HYP-1804
---

# HYP-2235: Finite-Field Kakeya/Falconer Gives an LRC No-Leak Carrier

## Claim

Finite-field Kakeya/Falconer should be imported into the repo as a
predicate-preserving carrier atlas, not as a famous-problem analogy.

The useful dictionary is:

```text
Kakeya direction cover       -> every clock/direction obligation is hit
line concurrency packet      -> owner/multiplicity side channel
Falconer distance support    -> scalar safe/distance support
pinned distance fibers       -> observer-coupled side channel
unit-distance direction spine -> unit-spine/direction-support side channel
Paley/Fourier character data  -> finite-field phase orientation
```

The LRC `n=14` translation is that after a scalar wall is hit, the proof cannot
forget the owner/concurrency label.  In finite fields that label is line
multiplicity and pinned distance data; in LRC `n=14` it is the `C=27` gcd shell,
carry, and owner route.

## External Anchors

The finite-field Kakeya anchor is Dvir's polynomial-method theorem
([arXiv:0803.2336](https://arxiv.org/abs/0803.2336)): a Kakeya set in
`F_q^n` contains a line in every direction and has size at least `C_n q^n`.

The finite-field Falconer anchor is the Iosevich-Rudnev distance-set program
([arXiv:math/0509005](https://arxiv.org/abs/math/0509005)): for
`E subset F_q^d`, estimate

```text
Delta(E) = { (x_1-y_1)^2 + ... + (x_d-y_d)^2 : x,y in E }
```

from `|E|` using finite-field Fourier analysis and exponential sums.  The
pinned-distance and distance-graph variants, such as Hart-Iosevich-Koh-Senger-
Uriarte-Tuero ([arXiv:0804.3036](https://arxiv.org/abs/0804.3036)), are
important because they retain the observer fiber that raw distance support can
erase.

## S659 Evidence

`04-computation/finite_field_kakeya_falconer_s659.py` builds affine line
families in `F_p^2`, one line in every projective direction, and records union
size, point multiplicity, distance support, pinned distance counts,
unit-distance direction support, Paley distance signs, and proof-carrier
Tournament Analysis.

Exact branch-and-bound over all one-line-per-direction families gives:

```text
p=3: exact_min=7,  min_family_count=72,    branch_nodes=121
p=5: exact_min=17, min_family_count=3000,  branch_nodes=16781
p=7: exact_min=31, min_family_count=16464, branch_nodes=2221668
```

These match the odd-prime parabola carrier size

```text
(p^2 + 2p - 1)/2.
```

For parabola carriers at `p=3,5,7,11,13`, the distance support is already full:

```text
p=3:  distances=3/3
p=5:  distances=5/5
p=7:  distances=7/7
p=11: distances=11/11
p=13: distances=13/13.
```

This is the key scalar-collapse warning.  Full Falconer distance support is
too coarse to identify the load-bearing Kakeya/Falconer packet.

The pinned and unit channels still vary:

```text
p=5:  pinned_min/max/full = 4/5/16, unit_edges=24
p=7:  pinned_min/max/full = 6/7/30, unit_edges=79
p=11: pinned_min/max/full = 10/11/69, unit_edges=245
p=13: pinned_min/max/full = 13/13/97, unit_edges=344.
```

The exhaustive `p=5` scalar-leak audit is especially sharp.  Across all
`5^6=15625` one-line-per-direction families:

```text
union_size_hist={17:3000, 18:4000, 19:8000, 21:600, 25:25}
distance_count_hist={5:15625}
pinned_min_hist={4:3600, 5:12025}.
```

Every row has full distance support, but union/concurrency and pinned fibers
still distinguish rows.  Thus distance support is not the sufficient statistic.

Tournament Analysis ranks the proof carriers, not raw point sets:

```text
score_hist={0:1,2:3,4:1,6:2,7:2}
directed_3cycles=3
scc_sizes=[4,3,1,1]
hamiltonian_paths=15.
```

The top carriers are pinned distance fibers and concurrency multiplicity
packets, followed by unit-distance direction spine and direction coverage.
Raw cardinality density is last.

## LRC14 Transfer

HYP-2231 says the active `n=14` wall is off-diagonal odd complement pairs
`(1,13),(3,11),(5,9)` plus `C=27` gcd-shell labels.  HYP-2230 says the carry
coordinate simultaneously controls parity and the `mod 14` obstruction.
HYP-2222 says the `C=27` gcd shell is a fixed clock.

HYP-2235 adds a finite-field incidence model for the same proof shape:

```text
direction cover is necessary but too scalar;
distance support can saturate early;
the proof object is pinned/concurrency/owner data.
```

The proposed LRC move is therefore a Kakeya-style no-leak lemma:

```text
If the n=14 scalar wall/pair-sum cover is fixed, then preserving
odd wall pairs + C=27 gcd shells + carry/owner labels should force
AP, Vstar, nonprimitive 2AP, or strict looseness.
```

The finite-field audit suggests how to search for failures: hold the scalar
support fixed, jackknife the hidden owner/concurrency labels, and watch which
rows lose the target predicate.  This is the S656 forcing-jackknife protocol in
incidence form.

## Assumption Challenge

The vertices in the S659 tournament are not points, lines, runners, or
distances.  They are proof carriers:

```text
direction cover, concurrency, Falconer support, pinned fibers, unit spine,
Paley character orientation, polynomial certificate, raw cardinality,
LRC shell/owner lift.
```

The quotient preserves direction coverage, distance support, pinned support,
unit-distance support, and owner/concurrency salience.  It destroys exact
embedding and exact line placement, which must be reattached before any
external Kakeya or Falconer theorem claim.

**See also:** HYP-2234, HYP-2233, HYP-2231, HYP-2230, HYP-2222, HYP-2197,
HYP-2176, HYP-1804, `04-computation/finite_field_kakeya_falconer_s659.py`,
`05-knowledge/results/finite_field_kakeya_falconer_s659.out`,
`07-reflections/finite-field-kakeya-falconer-lrc-carrier-s659.md`.

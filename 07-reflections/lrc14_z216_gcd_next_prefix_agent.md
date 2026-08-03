# LRC(14) `z1=216`: the next gcd-family prefix and a two-threshold overlap cut

**Status: FINITE-EXACT INDEPENDENT AUDIT OF
[THM-3313](../01-canon/theorems/THM-3313-projected-k3-z216-third-ruler-cost-prefix-multicover-closure.md),
PROJECTED-ATLAS SCOPE ONLY.**  Derived concurrently by a separate script and
certificate taxonomy, the next complete intrinsic-family prefix after the two
audited cost-prefix probes contains four rows, and its full exact necessary
upper screen is empty.  Thus the maintained projected ledger changes by
exactly

```text
373161 -> 373157,
z1=216 wall rows: 357 -> 353,
complete live families: 33 -> 31.                         (1)
```

This is not a physical-cover classification and makes no LRC(14) claim.  It
does not restore endpoint, owner, phase, current, chronology, arbitrary
`k<=1`, or the final rung.

THM-3313 is the current canonical proof source.  Its THM-3308 multicover
compiler and this note's coordinate-union/two-fan/two-threshold taxonomy both
produce the same unique gap-`510` two-layer obstruction, providing an
independent implementation and presentation audit rather than a second
ledger decrement.

Companions:

- [`lrc14_j7_k3_z216_gcd_next_prefix_agent.py`](../04-computation/lrc14_j7_k3_z216_gcd_next_prefix_agent.py);
- its [frozen deterministic output](../05-knowledge/results/lrc14_j7_k3_z216_gcd_next_prefix_agent.out).

## Inheritance and exact selection

The closest proved screen is [THM-3139](../01-canon/theorems/THM-3139-projected-k3-z225-terminal-and-z224-screen-double-layer-descent.md).
[THM-3281](../01-canon/theorems/THM-3281-projected-k3-z216-three-natural-wall-family-screen-descent.md)
licenses complete `(gcd(216,L),L)` fibres as bounded work units.  The two
[independently audited intrinsic prefixes](lrc14-z216-intrinsic-ruler-cost-prefix-independent-audit-codex-20260803.md)
leave exactly `357` wall rows in `33` such families.

The new program does not start from four pinned atlas addresses.  It rebuilds
the two old atlas parsers, reconstructs the disjoint gcd-eight, gcd-eighteen,
order, and THM-3281 closures, derives the first two cost prefixes again, and
then reranks all `357` live rows by complete-family invoice `sum_E L(E)r(E)`.
The next prefix through the first nonsingleton family is thereby forced to be

| family | row | body `E` | components `r` | invoice `Lr` |
|---|---:|---|---:|---:|
| `gcd72/L72072` | 157 | `(1,4,9,11,12,13)` | 24 | 1,729,728 |
| `gcd24/L25872` | 53 | `(1,2,6,8,11,14)` | 26 | 672,672 |
| `gcd24/L25872` | 293 | `(2,4,6,8,11,14)` | 26 | 672,672 |
| `gcd24/L25872` | 357 | `(2,6,8,11,12,14)` | 30 | 776,160 |

The total invoice is `3,851,232`; the boundary before the next family is
strict.  Cost remains only a scheduling order, never a safety invariant.

## Complete exact screen

Every selected row closes before a terminal step:

| row | states | crude | exact status | residual |
|---:|---:|---:|---:|---:|
| 157 | 111 | 93 | 18 | 0 |
| 53 | 27 | 17 | 10 | 0 |
| 293 | 14 | 9 | 5 | 0 |
| 357 | 271 | 170 | 101 | 0 |
| **total** | **423** | **289** | **134** | **0** |

All `134` solver-returned Farkas certificates were separately rebuilt from
the divisor tuple and verified over exact rationals.  Their vectors and
normalization-dependent contradiction magnitudes are not persisted, following
MISTAKE-331 and MISTAKE-333.  One negative-`alpha` mutation per row was
deliberately rejected.

## A new two-threshold overlap lemma

The earlier prefix taxonomy left a small weighted exact-Farkas core.  In the
present prefix, ordinary coordinate unions, zero-marginal reductions, and
two-fans explain `133/134` statuses.  The remaining row-157 state has a short
two-threshold incidence obstruction.

Let `T_t` denote the histogram mass with load at least `t`, let `m_j` be the
four status marginals, and let `M=sum_j m_j` be total bit incidence.  Suppose:

1. every pattern capable of supporting the lower threshold `u` is nonempty;
2. outside the coordinate union `S`, every pattern capable of supporting the
   higher threshold `v>u` contains at least `1+e` bits.

At most `sum_(j in S)m_j` high-tail points can lie in that coordinate union.
Every remaining high-tail point contributes at least `e` incidences beyond
the one incidence already forced at the lower threshold.  Therefore

```text
M >= T_u + e max(0, T_v - sum_(j in S) m_j).              (2)
```

For the sole formerly weighted state,

```text
divisors = (99,616,1001,8008),       q=9009,
marginals = (1365,1287,1287,1287),  M=5226,
capacities by mask = (2,8,3,8,8,8,8,8,3,8,4,8,8,8,8,8),
T_3=5190, T_4=3198, S={0,2}, e=1.                       (3)
```

At threshold four, the only capable mask outside `{bit 0 or bit 2}` is mask
`10`, containing bits one and three.  Hence at least

```text
3198 - (1365+1287) = 546
```

high-tail points have a second bit.  Formula `(2)` would require at least
`5190+546=5736` incidences, but only `5226` exist: an exact contradiction of
gap `510`.

The full solver-free taxonomy is consequently

```text
128 coordinate-union statuses,
  2 zero-reduced coordinate unions,
  3 two-fans,
  1 two-threshold overlap,
  0 weighted-only states.                                  (4)
```

This suggests a useful connection with Boolean rank and submodular cover
inequalities: instead of classifying each load threshold separately, couple
nested tail events through total incidence.  The cheapest next theoretical
test is to enumerate all four-bit two-threshold rank inequalities and compare
them with the eleven weighted states from the preceding sixteen-row family.

## Hostile controls and direction

The same run reconstructs three known nonempty screens, derived from the
pinned low-cost and costly gcd-eight evidence:

| row | states | crude | status | residual | residual quotient |
|---:|---:|---:|---:|---:|---|
| 64 | 115 | 50 | 57 | 8 | `8^8` |
| 238 | 1,035 | 237 | 683 | 115 | `8^115` |
| 370 | 354 | 213 | 122 | 19 | `8^19` |

Every control status certificate is exactly audited by the inherited screen.
These controls rule out a vacuous convention that labels every screen empty.
Their known terminal closures are not used in the new four-row argument.

The only logical direction used is

```text
empty exact necessary upper screen
    => corresponding projected wall row is empty.          (5)
```

No converse from screen feasibility is asserted.  Since the four selected
screens have no residual, no carrier, high-slot, or terminal claim is needed.

## Consequence and next queue head

Removing the two complete families gives `(1)`.  The next family queue begins

```text
gcd24/L76440: one row, invoice 2,293,200;
gcd24/L30576: three rows, invoice 2,629,536.                (6)
```

The first is a cheap sentinel before another coherent three-row fibre.  A
parallel theory lane is now available: test whether the two-threshold lemma
eliminates the earlier eleven weighted cores before expanding the atlas queue.

## Reproduction and hashes

```text
python 04-computation/lrc14_j7_k3_z216_gcd_next_prefix_agent.py --processes 3
python -O 04-computation/lrc14_j7_k3_z216_gcd_next_prefix_agent.py --processes 3
```

The normal and optimized computations have empty standard error and identical
LF-normalized output.  The source has no `assert` truth gates and no float
literals.

```text
semantic sha256  45a88091b7a1d6963807c81cde810ded1b3d0dc14f2fdc468dba6669295a5ef2
script sha256    89b0f1ff38ffa2ed254d05e824bf166f709ca603d620b1b74628883bf30ef7ca
output sha256    0a7884faa1b661833a4aed0156b04cbcfbdff25923967b1cb3f25be586ab0407
```

---
id: HYP-2044
status: OPEN; synthesis from existing Goldbach/Lemoine, n=4 parity, and LRC gauge-stack notes
source: codex-2026-06-01-S551
related:
  - HYP-1963
  - HYP-1965
  - HYP-1966
  - HYP-1976
  - HYP-1983
  - HYP-1984
  - HYP-2039
  - HYP-2043
  - THM-393
---

# HYP-2044: doubled primes are first-even parity gates in the LRC add/multiply stack

## Statement

For additive prime fibers inside the LRC denominator stack, the distinction

```text
even N:  N = p + q
odd N:   N = p + 2q
```

is not only the usual Goldbach/Lemoine parity fix.  The doubled prime `2q` is
the minimal dyadic bridge that carries an odd prime core through the first-even
row of the natural-number address

```text
N = 2^h * odd_core.
```

Thus a Lemoine representation should be read as a two-layer object:

```text
odd prime atom p
+ first-even prime-core bridge 2q.
```

The bridge spends exactly one factor of `2` to hit an odd target while
preserving prime-core arithmetic in the even leg.  In LRC terms, doubled
primes are candidate parity gates: they break all-odd boundary symmetries but
do not behave like generic composite even speeds.

## Why this belongs to LRC

HYP-1984 already puts denominators on the `x+2` / `x*2` grid: addition moves
horizontally through odd cores, multiplication moves vertically through dyadic
height, and product-sum equations label collisions.  Goldbach pairs occupy
even fixed-sum fibers.  Lemoine pairs occupy odd fixed-sum fibers, but only by
placing one summand on the `h=1` branch.

At `n=4`, HYP-2043/THM-393 make the parity role sharper.  The safe indicator
has Fourier support only on odd harmonics,

```text
g_k = -chi4(k)/(pi k),
```

so a pairwise correction vanishes whenever one reduced cofactor is even.  A
doubled-prime leg is therefore pairwise-neutral in the mod-4 character
ledger, while still being able to destroy the all-odd witness `t=1/4`.  This
suggests the doubled-prime gate moves obstruction out of the pairwise `R2`
layer and into boundary, full-support, or higher-order resonance layers.

## Relation to odd and even cycles

Tournament parity has the same shape at a different level.  In the OCF,
directed odd cycles are the surviving atoms:

```text
H(T) = I(Omega(T), 2) = 1 + 2 alpha_1 + 4 alpha_2 + ...
```

Even-cycle structure is mostly cancellation, pairing, or quotient data; odd
cycles survive as the obstruction-carrying pieces.  In additive number theory,
all primes except `2` are odd atoms.  Goldbach's even target uses two odd
atoms.  Lemoine's odd target cannot use two odd atoms; it needs one dyadic
operator.  The doubled prime is exactly that operator with the odd atom still
visible underneath.

This is the proposed analogy:

```text
odd cycle survives quotienting        <-> odd prime core survives doubling
even cycle cancels/pairs              <-> factor 2 acts as parity transport
Goldbach p+q                          <-> pairwise odd-atom smoothing
Lemoine p+2q                          <-> odd-atom smoothing plus dyadic gate
```

## Tournament Analysis declaration

The right vertices should not automatically be runners.  Candidate vertex
sets:

```text
1. denominators N,
2. additive representation pair-cells (p,q) or (p,2q),
3. doubled-prime cores q,
4. first-even bridge events N-p=2q,
5. endpoint owners touched by the bridge,
6. residue channels, Fourier modes, or proof obligations.
```

For the first implementation, use vertices `b=(N,p,q)` representing Lemoine
bridges `N=p+2q` and, for comparison, Goldbach bridges `N=p+q`.

Pairwise observable:

```text
debt(b) = (
  target parity,
  v2(even_leg),
  odd_core(even_leg),
  local residue obstruction count,
  endpoint-debt proxy phi(N)/N,
  pairwise chi4-neutrality flag,
  representation abundance in the same N-fiber
)
```

Switch/gauge:

```text
b1 -> b2 iff b1 has lower boundary/higher-repair debt lexicographically;
ties follow increasing (N, q, p).
```

Tie Hamiltonian path:

```text
sort by N, then q, then p.
```

Required fingerprints:

```text
score histogram, directed 3-cycles, SCCs, Hamiltonian-path count,
edge flips against the S513 add/multiply majority tournament,
and whether bridge-rich odd N have lower endpoint-pressure debt.
```

## Assumption challenge

Challenged assumption: doubled primes are merely a number-theory way to repair
odd parity.

Alternative reading: doubled primes are LRC first-even bridge events.  They
preserve the LRC-relevant predicate "can this additive fiber provide a parity
gate without losing prime-core smoothing?"  They destroy exact time geometry,
endpoint ownership, and pressure-core information unless those labels are
reattached.  Therefore a runner-only or denominator-only tournament is too
coarse; a useful quotient needs either bridge pair-cells or endpoint-owner
decorations.

## Predictions

1. In `n=4` style character decompositions, doubled-prime legs should be
   pairwise quiet but boundary-active: they kill the all-odd `t=1/4` witness
   while contributing no `chi4` pair correction when the reduced cofactor is
   even.
2. Odd denominators with many Lemoine bridges should have more first-even
   repair channels; Lemoine-scarce odd denominators with high endpoint debt
   are better counterexample-shaped candidates.
3. A useful odd-denominator LRC stack should separate "arbitrary even leg"
   from "doubled prime leg"; the latter is the minimal dyadic bridge with a
   prime odd core.
4. The next computation should compare bridge-cell tournaments against the
   existing S513 denominator tournament and report whether doubled-prime
   bridges explain any edge flips between additive abundance, endpoint debt,
   and A000568 odd survival.

## See also

- `05-knowledge/hypotheses/HYP-1983-lrc-operation-arc-zoo.md`
- `05-knowledge/hypotheses/HYP-1984-lrc-add-mult-gauge-stack.md`
- `05-knowledge/hypotheses/HYP-2043-lrc-n4-measure-gap-unique-tight.md`
- `07-reflections/lrc-doubled-primes-as-first-even-parity-gates-s551.md`
- `04-computation/lrc_add_mult_gauge_stack_s513.py`
- `04-computation/lrc_operation_arc_zoo_s511.py`

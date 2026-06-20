# LRC14 Glaisher/Witt/Even-Graph Bridge - S37 Reflection

Date: 2026-06-19
Session: codex-2026-06-19-S37
Claimed: HYP-2660 / T905

## Why This Session

The prompt asked to keep aiming at LRC14 while thinking through three
equinumerosities:

- odd partitions versus distinct partitions,
- tournaments versus even graphs,
- Witt algebras and Laurent polynomials.

The existing repo state made the first link unusually live.  KPS HYP-2656 had
already found the CRT/halving structure in the LRC14 one-hole AP-window layer:
odd speeds are a rigid base, and even speeds are dyadic refinements.  The
Glaisher bijection is the exact combinatorial language for that sentence.

## Main Finding

The common object is not equality of counts.  It is closure after gauge choice.

Odd partitions close to distinct partitions by expanding multiplicities in base
2.  Tournaments close to even graphs by adding a root whose incident edges force
even degrees.  Euler products close to Witt ghost coordinates after applying
`D=x d/dx` to the logarithm.  LRC14 should close raw speed data into a tower/CRT
/endpoint-owner address before asking scalar questions.

That gives the current proof hierarchy:

```text
glaisher_tower_word
> crt_residue_cell
> endpoint_owner_ledger
> even_graph_parity_closure
> witt_ghost_divisor_sum
> laurent_wall_polynomial
> raw_sumset_excess
> raw_speed
```

## Useful Correction

The new single-deletion ledger caught a small but important reading issue in KPS
HYP-2656.  The low-even block is `{2,4,6,10,12}`; drop `8=2^3` is the high-level
even outlier:

```text
drop  6: 7/858
drop 12: 426/35035
drop 10: 1520/63063
drop  4: 97/4004
drop  2: 11/364
drop 13: 6617/194040
drop  5: 103/2860
drop  3: 11/273
drop  8: 950/21021
```

So the Glaisher tower is an address, not a one-number ranking.  That is exactly
the kind of nuance the LRC14 proof route keeps needing: the dyadic mechanism
spots the right layer, but the CRT cell and endpoint-owner ledger decide the
actual inequality.

## Witt/Laurent Reading

The script verifies

```text
m [x^m] log prod(1+x^n) = sigma_odd(m)
```

through `m<=40`.  This is the practical Witt-vector/ghost-coordinate face of
the partition identity.  For this proof project, the useful Laurent-polynomial
operator is not the whole Witt algebra at once; it is the logarithmic derivation
`x d/dx`, because it changes product walls into additive divisor ledgers.

That is a good analogy for LRC wall functions.  Product descriptions and raw
wall arrangements can look too big.  The ghost map asks for divisor/address
coordinates, and those are closer to the CRT/tower packets already appearing in
the LRC14 work.

## Tournament/Even-Graph Reading

The tournament/even-graph bijection is exact:

```text
tournaments on n vertices = even graphs on n+1 vertices.
```

The constructive map is more useful than the count.  Free tournament bits become
internal graph edges; root edges are then forced by parity.  That mirrors the
LRC situation: one can choose local wall defects, but the safe-set predicate
forces CRT and endpoint closures.  A proof that treats all pair bits as
independent will miss the root closure.

## What Failed

Raw equinumerosity failed as a proof invariant.  Knowing that two sets have the
same size does not preserve `G_C`, component endpoints, or endpoint owners.

Raw dyadic level also failed as a total order.  Drop 8 is the warning.  The
right lesson from Glaisher is not "even beats odd"; it is "record the dyadic
tower word before projecting to CRT cells and owners."

## Next Target

For OPEN-Q-108, the next concrete computation should extend the Glaisher tower
word from one-hole AP-window rows to the THM-543/544 AP-tail layers.  The first
question is whether all sub-`426/35035` tails have a small finite list of tower
defect words plus endpoint-owner bubble types.  If yes, this gives a finite
state-template lemma compatible with HYP-2658's far-survival recursion and
HYP-2659's odd-shell carry carrier.

---
id: HYP-2110
status: NEW local lemmas E/F PROVED; sampled n=14 multiple branch closed; global large-owner residual OPEN
source: codex-2026-06-03-S581
related:
  - THM-398
  - HYP-2108
  - HYP-2107
  - HYP-2105
  - HYP-2104
  - HYP-2103
  - HYP-2102
---

# HYP-2110: one small endpoint owner forces a pin-or-peel descent

## Claim

HYP-2105's small-owner theorem does not stop at the both-small-owner case.
A single small endpoint owner already creates a rigid pin.  For a component of
`G(S')` and a multiple runner `v=nw`, a small endpoint owner `u<n` must satisfy

```text
w(k n +/- 1) = j u
```

to sit in a `v`-danger arc.  Therefore the component peels in either of two
new ways:

```text
Lemma E: if u does not divide w(k n +/- 1), the endpoint is off the v-centre
         lattice, so the component is uncoverable.

Lemma F: if the endpoint is pinned to a v-centre but the component length is
         > 1/(n^2 w), the component is too long for a single half-arc, so it
         is uncoverable.
```

Together with Bprime and Lemma C, this is a much sharper endpoint-cover
circuit positivity sieve:

```text
long interval -> loose
both small owners -> loose
one small owner off lattice -> loose
one small owner pinned but longer than half-radius -> loose
remaining residual -> large-owner or mixed tiny half-radius window
```

## S581 Evidence

`04-computation/lrc_small_owner_descent_s581.py` samples multiple-of-`n` rows and
routes each row by the first component-level proof criterion.

```text
n= 6: rows=4984; proved=4607 (92.44%); residual=377
n= 8: rows=5000; proved=4863 (97.26%); residual=137
n=10: rows=5000; proved=5000 (100.00%); residual=0
n=12: rows=5000; proved=5000 (100.00%); residual=0
n=14: rows=5000; proved=5000 (100.00%); residual=0
```

For `n=14`, the route split is:

```text
Bprime_long_interval           3407
LemmaC_both_small               843
LemmaE_small_pin_off_lattice    748
LemmaF_small_pin_half_radius      2
large_owner_residual              0
```

This does not prove Cprime at `n=14`, but it turns the previous S574 `~1%`
large-binding-owner sample residual into zero in the S581 sample.  The remaining
general obstruction is no longer "large owner somewhere"; it is:

```text
every component is short,
no component has a both-small owner pair,
every one-small endpoint is either pinned and half-radius tiny or absent,
all remaining binding data are large-owner congruence windows.
```

That is a strictly smaller CRT/Diophantine target.

## Endpoint-Cover Circuit Positivity

The endpoint-cover positivity route is now:

```text
cover assignment
  -> endpoint owner congruences
  -> small owners become exact pins
  -> off-lattice or half-radius failure peels
  -> only large-owner tiny windows remain
```

This is why endpoint-cover circuits are a better proof object than raw speed
rows.  The circuit carries the exact obstruction type.  A failed circuit is not
just "not covered"; it has a named local reason that produces a positive
measure interval.

## Tournament Analysis

Vertices are proof criteria:

```text
Bprime_long_interval,
LemmaC_both_small,
LemmaE_small_pin_off_lattice,
LemmaF_small_pin_half_radius,
large_owner_residual.
```

Pair observable:

```text
(proof_rank, coverage)
```

Switch:

```text
stronger proved criterion beats weaker criterion; residual loses.
```

Fingerprint:

```text
score_hist={0:1, 1:1, 2:1, 3:1, 4:1}
directed_3_cycles=0
hamiltonian_path =
  Bprime_long_interval -> LemmaC_both_small ->
  LemmaE_small_pin_off_lattice -> LemmaF_small_pin_half_radius ->
  large_owner_residual
```

This remains a transitive proof ledger.  Nontrivial cycles should be sought
inside the remaining large-owner congruence-window graph, not at this route
level.

## Assumption Challenge

Possible vertices considered:

```text
runners, endpoints, safe components, owner pairs, congruence windows,
cover assignments, and proof criteria.
```

Chosen quotient: proof criteria.  It preserves the predicate:

```text
some component is certified uncoverable, hence S is loose.
```

It destroys exact endpoint locations and speed identities after they have been
summarized by route.  That loss is deliberate: the current question is no
longer which runner is large, but which endpoint-cover circuit can still avoid
pinning or peeling.

## Honest Status

Lemmas E and F are local proofs and have been folded into THM-398's status
ledger.  The S581 `n=14` sample closes, but this is not a proof of LRC(14).
The next theorem target is the large-owner tiny-window lemma: prove that the
remaining large-owner congruence windows cannot be simultaneously assigned
around the whole endpoint-cover circuit.

**See:** `04-computation/lrc_small_owner_descent_s581.py`
(+ `05-knowledge/results/lrc_small_owner_descent_s581.out`),
`07-reflections/lrc-small-owner-descent-past-lemma-c-s581.md`,
THM-398, HYP-2108, HYP-2107, HYP-2105, HYP-2104, HYP-2103, HYP-2102.

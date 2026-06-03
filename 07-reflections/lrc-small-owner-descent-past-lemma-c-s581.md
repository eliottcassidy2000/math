# Small-owner descent past Lemma C (S581)

The S574 cover-to-congruence translator gave the right door:

```text
endpoint cover geometry -> endpoint owner congruences.
```

Lemma C walked through the easiest part of that door: if both endpoint owners of
a safe component are small (`<n`), then both endpoints have to be the same
`v=nw` arc centre.  That forces `a=b`, impossible for a component.

Today's extension is that one small owner already matters.

## The pin

For an endpoint owner `u`, the congruence window is

```text
|w(k n +/- 1) - j u| < u/n.
```

If `u<n`, the right side is `<1`, and the left side is an integer.  So it must
be `0`.  The small endpoint owner pins the endpoint exactly to a `v=nw` arc
centre:

```text
w(k n +/- 1) = j u.
```

That gives two exits.

1. If `u` does not divide `w(k n +/- 1)`, the endpoint is not even on the
   centre lattice.  The component cannot be covered.
2. If the endpoint is pinned, the component can use only one side of the
   danger arc.  So it must fit within a half-radius `1/(n^2 w)`, not the full
   Bprime radius `2/(n^2 w)`.

Those are Lemma E and Lemma F in the updated THM-398 ledger.

## What changed

Before S581, the S574 residual was:

```text
every component is short,
and every component has at least one large binding owner.
```

After the one-small-owner pin, the residual is smaller:

```text
every component is short,
no component has two small owners,
every one-small component is centre-pinned and half-radius tiny,
and the rest is large-owner congruence slack.
```

In the deterministic S581 sample this closes the whole multiple-of-14 sample:

```text
n=14, rows=5000:
  Bprime_long_interval           3407
  LemmaC_both_small               843
  LemmaE_small_pin_off_lattice    748
  LemmaF_small_pin_half_radius      2
  residual                          0
```

It also closes every sampled row at `n=10` and `n=12`; residuals remain only in
the smaller `n=6,8` samples, where they are large-owner or mixed tiny windows.

## Why this is endpoint-cover circuit positivity

The proof route is no longer "find a witness time" directly.  It is:

```text
cover assignment fails locally
  -> a component survives
  -> positive measure interval exists.
```

Small owners are exact pins.  Exact pins are brittle.  They either miss the
centre lattice or consume only half an arc.  This is exactly the kind of local
peeling law an endpoint-cover circuit positivity theorem needs.

## Next target

The remaining theorem should be stated over large-owner congruence windows:

```text
No closed endpoint-cover circuit can be made entirely of
large-owner or centre-pinned half-radius windows.
```

That is not a speed-set enumeration.  It is a finite CRT/Diophantine
compatibility statement on endpoint owners, arc centres, and multiplier `w`.
The current script says the n=14 sampled version has no survivors; the next
script should build the large-owner window graph directly and test cycle
incompatibility.

## Tournament Analysis

Vertices used here are proof criteria rather than runners:

```text
Bprime, Lemma C, Lemma E, Lemma F, residual.
```

The tournament is transitive.  That is fine.  It says the route ledger is
ordered.  Any future cycle worth caring about should appear inside the
large-owner window graph itself, where cover choices can compete.

## Handoff

Do not go back to broad row enumeration yet.  The next useful computation is:

```text
for every residual component, build all feasible large-owner arc-centre windows;
then form the compatibility graph across components and search for a global
cover assignment.
```

If that graph has no full transversal in the n=14 owner families, the Cprime
residual becomes an owner-window Hall obstruction.
